/*
 * Copyright 2024 The Android Open Source Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifdef UHDR_ENABLE_HEIF

#include "ultrahdr/heifr.h"
#include "ultrahdr/gainmapmath.h"
#include "ultrahdr/gainmapmetadata.h"

namespace ultrahdr {

#define HEIF_ERR_CHECK(x)                                               \
  {                                                                     \
    heif_error err = (x);                                               \
    if (err.code != heif_error_Ok) {                                    \
      status.error_code = UHDR_CODEC_ERROR;                             \
      status.has_detail = 1;                                            \
      snprintf(status.detail, sizeof status.detail, "%s", err.message); \
      goto CleanUp;                                                     \
    }                                                                   \
  }

class MemoryWriter {
 public:
  MemoryWriter() : data_(nullptr), size_(0), capacity_(0) {}

  ~MemoryWriter() { free(data_); }

  const uint8_t* data() const { return data_; }

  size_t size() const { return size_; }

  void write(const void* data, size_t size) {
    if (capacity_ - size_ < size) {
      size_t new_capacity = capacity_ + size;
      uint8_t* new_data = static_cast<uint8_t*>(malloc(new_capacity));
      if (data_) {
        memcpy(new_data, data_, size_);
        free(data_);
      }
      data_ = new_data;
      capacity_ = new_capacity;
    }
    memcpy(&data_[size_], data, size);
    size_ += size;
  }

 public:
  uint8_t* data_;
  size_t size_;
  size_t capacity_;
};

static struct heif_error writer_write([[maybe_unused]] struct heif_context* ctx, const void* data,
                                      size_t size, void* userdata) {
  MemoryWriter* writer = static_cast<MemoryWriter*>(userdata);
  writer->write(data, size);
  struct heif_error err {
    heif_error_Ok, heif_suberror_Unspecified, nullptr
  };
  return err;
}

static struct heif_error fill_img_plane(heif_image* img, heif_channel channel, void* srcBuffer,
                                        int width, int height, int srcRowStride) {
  heif_error err = heif_image_add_plane(img, channel, width, height, 8 /*bitDepth */);
  if (err.code != heif_error_Ok) return err;

  int dstStride;
  uint8_t* dstBuffer = heif_image_get_plane(img, channel, &dstStride);
  uint8_t* srcRow = static_cast<uint8_t*>(srcBuffer);

  for (int y = 0; y < height; y++) {
    memcpy(dstBuffer, srcRow, width);
    srcRow += srcRowStride;
    dstBuffer += dstStride;
  }
  return err;
}

static heif_color_primaries map_cg_cicp_primaries(uhdr_color_gamut_t gamut) {
  heif_color_primaries color_primaries = heif_color_primaries_unspecified;
  if (gamut == UHDR_CG_BT_709) {
    color_primaries = heif_color_primaries_ITU_R_BT_709_5;
  } else if (gamut == UHDR_CG_DISPLAY_P3) {
    color_primaries = heif_color_primaries_SMPTE_EG_432_1;
  } else if (gamut == UHDR_CG_BT_2100) {
    color_primaries = heif_color_primaries_ITU_R_BT_2020_2_and_2100_0;
  }
  return color_primaries;
}

static heif_transfer_characteristics map_ct_cicp_transfer_characteristics(
    uhdr_color_transfer_t tf) {
  heif_transfer_characteristics transfer_characteristics = heif_transfer_characteristic_unspecified;
  if (tf == UHDR_CT_SRGB) {
    transfer_characteristics = heif_transfer_characteristic_ITU_R_BT_709_5;
  } else if (tf == UHDR_CT_LINEAR) {
    transfer_characteristics = heif_transfer_characteristic_linear;
  } else if (tf == UHDR_CT_PQ) {
    transfer_characteristics = heif_transfer_characteristic_ITU_R_BT_2100_0_PQ;
  } else if (tf == UHDR_CT_HLG) {
    transfer_characteristics = heif_transfer_characteristic_ITU_R_BT_2100_0_HLG;
  }
  return transfer_characteristics;
}

static heif_matrix_coefficients map_cg_cicp_matrix_coefficients(uhdr_color_gamut_t gamut) {
  heif_matrix_coefficients matrix_coefficients = heif_matrix_coefficients_unspecified;
  if (gamut == UHDR_CG_BT_709) {
    matrix_coefficients = heif_matrix_coefficients_ITU_R_BT_709_5;
  } else if (gamut == UHDR_CG_DISPLAY_P3) {
    matrix_coefficients = heif_matrix_coefficients_chromaticity_derived_non_constant_luminance;
  } else if (gamut == UHDR_CG_BT_2100) {
    matrix_coefficients = heif_matrix_coefficients_ITU_R_BT_2020_2_non_constant_luminance;
  } else if (gamut == UHDR_CG_BT_601) {
    matrix_coefficients = heif_matrix_coefficients_ITU_R_BT_601_6;
  }
  return matrix_coefficients;
}

static heif_error map_fmt_to_heif_chroma_vars(uhdr_img_fmt_t fmt, unsigned int w, unsigned h,
                                              enum heif_chroma& heif_img_fmt, int& chromaWd,
                                              int& chromaHt) {
  if (fmt == UHDR_IMG_FMT_12bppYCbCr420) {
    heif_img_fmt = heif_chroma_420;
    chromaWd = w / 2;
    chromaHt = h / 2;
  } else if (fmt == UHDR_IMG_FMT_16bppYCbCr422) {
    heif_img_fmt = heif_chroma_422;
    chromaWd = w;
    chromaHt = h / 2;
  } else if (fmt == UHDR_IMG_FMT_24bppYCbCr444) {
    heif_img_fmt = heif_chroma_444;
    chromaWd = w;
    chromaHt = h;
  } else {
    heif_error status;
    status.code = heif_error_Unsupported_feature;
    status.subcode = heif_suberror_Unsupported_color_conversion;
    status.message = "mapping uhdr image format to heif chroma format failed";
    return status;
  }
  return {heif_error_Ok, heif_suberror_Unspecified, nullptr};
}

static heif_error set_internal_color_format(heif_encoder* encoder, enum heif_chroma heif_img_fmt) {
  struct heif_error err {
    heif_error_Ok, heif_suberror_Unspecified, nullptr
  };
  switch (heif_img_fmt) {
    case heif_chroma_420:
    case heif_chroma_monochrome:
      err = heif_encoder_set_parameter(encoder, "chroma", "420");
      break;
    case heif_chroma_422:
      err = heif_encoder_set_parameter(encoder, "chroma", "422");
      break;
    case heif_chroma_444:
      err = heif_encoder_set_parameter(encoder, "chroma", "444");
      break;
    default:
      err.code = heif_error_Invalid_input;
      err.subcode = heif_suberror_Unsupported_parameter;
      err.message = "unsupported heif image format setting";
      break;
  }
  return err;
}

HeifR::HeifR(void* uhdrGLESCtxt, int mapDimensionScaleFactor, int mapCompressQuality,
             bool useMultiChannelGainMap, float gamma, uhdr_enc_preset_t preset,
             float minContentBoost, float maxContentBoost, float targetDispPeakBrightness,
             uhdr_codec_t codec)
    : UltraHdr(uhdrGLESCtxt, mapDimensionScaleFactor, mapCompressQuality, useMultiChannelGainMap,
               gamma, preset, minContentBoost, maxContentBoost, targetDispPeakBrightness),
      mCodec(codec) {}

/* Encode API-0 */
uhdr_error_info_t HeifR::encodeHEIFR(uhdr_raw_image_t* hdr_intent, uhdr_compressed_image_t* dest,
                                     int quality, uhdr_mem_block_t* exif) {
  uhdr_img_fmt_t sdr_intent_fmt;
  if (hdr_intent->fmt == UHDR_IMG_FMT_24bppYCbCrP010) {
    sdr_intent_fmt = UHDR_IMG_FMT_12bppYCbCr420;
  } else if (hdr_intent->fmt == UHDR_IMG_FMT_30bppYCbCr444) {
    sdr_intent_fmt = UHDR_IMG_FMT_24bppYCbCr444;
  } else if (hdr_intent->fmt == UHDR_IMG_FMT_32bppRGBA1010102 ||
             hdr_intent->fmt == UHDR_IMG_FMT_64bppRGBAHalfFloat) {
    sdr_intent_fmt = UHDR_IMG_FMT_32bppRGBA8888;
  } else {
    uhdr_error_info_t status;
    status.error_code = UHDR_CODEC_INVALID_PARAM;
    status.has_detail = 1;
    snprintf(status.detail, sizeof status.detail, "unsupported hdr intent color format %d",
             hdr_intent->fmt);
    return status;
  }
  std::unique_ptr<uhdr_raw_image_ext_t> sdr_intent = std::make_unique<uhdr_raw_image_ext_t>(
      sdr_intent_fmt, UHDR_CG_UNSPECIFIED, UHDR_CT_UNSPECIFIED, UHDR_CR_UNSPECIFIED, hdr_intent->w,
      hdr_intent->h, 64);

  // tone map
  UHDR_ERR_CHECK(toneMap(hdr_intent, sdr_intent.get()));

  // If hdr intent is tonemapped internally, it is observed from quality pov,
  // generateGainMapOnePass() is sufficient
  mEncPreset = UHDR_USAGE_REALTIME;  // overriding the config option

  // generate gain map
  uhdr_gainmap_metadata_ext_t metadata(kJpegrVersion);
  std::unique_ptr<uhdr_raw_image_ext_t> gainmap;
  UHDR_ERR_CHECK(generateGainMap(sdr_intent.get(), hdr_intent, &metadata, gainmap,
                                 sdr_intent->cg == UHDR_CG_DISPLAY_P3 /* sdr_is_601 */,
                                 false /* use_luminance */));

  // compress sdr image
  std::unique_ptr<uhdr_raw_image_ext_t> sdr_intent_yuv_ext;
  uhdr_raw_image_t* sdr_intent_yuv = sdr_intent.get();
  if (isPixelFormatRgb(sdr_intent->fmt)) {
#if (defined(UHDR_ENABLE_INTRINSICS) && (defined(__ARM_NEON__) || defined(__ARM_NEON)))
    sdr_intent_yuv_ext = convert_raw_input_to_ycbcr_neon(
        sdr_intent.get(), sdr_intent->cg == UHDR_CG_DISPLAY_P3 /* use bt601 */);
#else
    sdr_intent_yuv_ext = convert_raw_input_to_ycbcr(
        sdr_intent.get(), sdr_intent->cg == UHDR_CG_DISPLAY_P3 /* use bt601 */);
#endif
    sdr_intent_yuv = sdr_intent_yuv_ext.get();
  }

  std::shared_ptr<DataStruct> baseIcc = IccHelper::writeIccProfile(UHDR_CT_SRGB, sdr_intent->cg);
  std::shared_ptr<DataStruct> alternateIcc =
      IccHelper::writeIccProfile(gainmap->ct, gainmap->cg, true);

  return encodeHEIFR(sdr_intent_yuv, gainmap.get(), &metadata, dest, quality, exif, baseIcc.get(),
                     alternateIcc.get());
}

/* Encode API-1 */
uhdr_error_info_t HeifR::encodeHEIFR(uhdr_raw_image_t* hdr_intent, uhdr_raw_image_t* sdr_intent,
                                     uhdr_compressed_image_t* dest, int quality,
                                     uhdr_mem_block_t* exif) {
  // generate gain map
  uhdr_gainmap_metadata_ext_t metadata(kJpegrVersion);
  std::unique_ptr<uhdr_raw_image_ext_t> gainmap;
  UHDR_ERR_CHECK(generateGainMap(sdr_intent, hdr_intent, &metadata, gainmap,
                                 sdr_intent->cg == UHDR_CG_DISPLAY_P3 /* sdr_is_601 */));

  std::unique_ptr<uhdr_raw_image_ext_t> sdr_intent_yuv_ext;
  uhdr_raw_image_t* sdr_intent_yuv = sdr_intent;
  if (isPixelFormatRgb(sdr_intent->fmt)) {
#if (defined(UHDR_ENABLE_INTRINSICS) && (defined(__ARM_NEON__) || defined(__ARM_NEON)))
    sdr_intent_yuv_ext = convert_raw_input_to_ycbcr_neon(
        sdr_intent, sdr_intent->cg == UHDR_CG_DISPLAY_P3 /* use bt601 */);
#else
    sdr_intent_yuv_ext = convert_raw_input_to_ycbcr(
        sdr_intent, sdr_intent->cg == UHDR_CG_DISPLAY_P3 /* use bt601 */);
#endif
    sdr_intent_yuv = sdr_intent_yuv_ext.get();
  }

  std::shared_ptr<DataStruct> baseIcc = IccHelper::writeIccProfile(UHDR_CT_SRGB, sdr_intent->cg);
  std::shared_ptr<DataStruct> alternateIcc =
      IccHelper::writeIccProfile(gainmap->ct, gainmap->cg, false);

  return encodeHEIFR(sdr_intent_yuv, gainmap.get(), &metadata, dest, quality, exif, baseIcc.get(),
                     alternateIcc.get());
}

uhdr_error_info_t HeifR::encodeHEIFR(uhdr_raw_image_t* sdr_intent, uhdr_raw_image_t* gainmap_img,
                                     uhdr_gainmap_metadata_ext_t* metadata,
                                     uhdr_compressed_image_t* dest, int quality,
                                     uhdr_mem_block_t* exif, DataStruct* baseIcc,
                                     DataStruct* alternateIcc) {
  uhdr_error_info_t status = g_no_error;
  heif_encoder* encoder = nullptr;
  heif_encoding_options* options = nullptr;
  heif_color_profile_nclx* sdrNclx = nullptr;
  heif_color_profile_nclx* gainmapNclx = nullptr;
  heif_color_profile_nclx* hdrNclx = nullptr;
  heif_image* baseImage = nullptr;
  heif_image_handle* baseHandle = nullptr;
  heif_image* secondaryImage = nullptr;
  heif_image_handle* secondaryHandle = nullptr;

  uhdr_raw_image_t* gainmap_img_yuv = gainmap_img;
  std::unique_ptr<uhdr_raw_image_ext_t> gainmap_yuv_ext;

  MemoryWriter writer;
  struct heif_writer w;
  w.writer_api_version = 1;
  w.write = writer_write;

  heif_context* ctx = heif_context_alloc();
  if (!ctx) {
    status.error_code = UHDR_CODEC_MEM_ERROR;
    status.has_detail = 1;
    snprintf(status.detail, sizeof status.detail, "failed to allocate heif context");
    return status;
  }

  HEIF_ERR_CHECK(heif_context_get_encoder_for_format(
      ctx, mCodec == UHDR_CODEC_AVIF ? heif_compression_AV1 : heif_compression_HEVC, &encoder));

  // set the encoder parameters
  HEIF_ERR_CHECK(heif_encoder_set_lossy_quality(encoder, quality));

  // encode the primary image
  enum heif_chroma heif_img_fmt;
  int chromaWd, chromaHt;
  HEIF_ERR_CHECK(map_fmt_to_heif_chroma_vars(sdr_intent->fmt, sdr_intent->w, sdr_intent->h,
                                             heif_img_fmt, chromaWd, chromaHt));
  HEIF_ERR_CHECK(heif_image_create(sdr_intent->w, sdr_intent->h, heif_colorspace_YCbCr,
                                   heif_img_fmt, &baseImage));
  HEIF_ERR_CHECK(fill_img_plane(baseImage, heif_channel_Y, sdr_intent->planes[UHDR_PLANE_Y],
                                sdr_intent->w, sdr_intent->h, sdr_intent->stride[UHDR_PLANE_Y]));
  HEIF_ERR_CHECK(fill_img_plane(baseImage, heif_channel_Cb, sdr_intent->planes[UHDR_PLANE_U],
                                chromaWd, chromaHt, sdr_intent->stride[UHDR_PLANE_U]));
  HEIF_ERR_CHECK(fill_img_plane(baseImage, heif_channel_Cr, sdr_intent->planes[UHDR_PLANE_V],
                                chromaWd, chromaHt, sdr_intent->stride[UHDR_PLANE_V]));
  if (baseIcc != nullptr) {
    void* ptr = (uint8_t*)baseIcc->getData() + kICCIdentifierSize;
    HEIF_ERR_CHECK(heif_image_set_raw_color_profile(baseImage, "prof", ptr,
                                                    baseIcc->getLength() - kICCIdentifierSize));
  }
  sdrNclx = heif_nclx_color_profile_alloc();
  if (!sdrNclx) {
    status.error_code = UHDR_CODEC_MEM_ERROR;
    status.has_detail = 1;
    snprintf(status.detail, sizeof status.detail,
             "failed to allocate nclx color profile for base image");
    goto CleanUp;
  }
  HEIF_ERR_CHECK(set_internal_color_format(encoder, heif_img_fmt));
  HEIF_ERR_CHECK(
      heif_nclx_color_profile_set_color_primaries(sdrNclx, map_cg_cicp_primaries(sdr_intent->cg)));
  HEIF_ERR_CHECK(heif_nclx_color_profile_set_transfer_characteristics(
      sdrNclx, map_ct_cicp_transfer_characteristics(sdr_intent->ct)));
  HEIF_ERR_CHECK(heif_nclx_color_profile_set_matrix_coefficients(
      sdrNclx, map_cg_cicp_matrix_coefficients(sdr_intent->cg == UHDR_CG_DISPLAY_P3
                                                   ? (uhdr_color_gamut_t)UHDR_CG_BT_601
                                                   : sdr_intent->cg)));
  sdrNclx->full_range_flag = sdr_intent->range == UHDR_CR_FULL_RANGE ? true : false;
  options = heif_encoding_options_alloc();
  if (!options) {
    status.error_code = UHDR_CODEC_MEM_ERROR;
    status.has_detail = 1;
    snprintf(status.detail, sizeof status.detail,
             "failed to allocate nclx color profile for base image");
    goto CleanUp;
  }
  options->save_two_colr_boxes_when_ICC_and_nclx_available = 1;
  options->output_nclx_profile = sdrNclx;
  HEIF_ERR_CHECK(heif_image_set_nclx_color_profile(baseImage, sdrNclx));
  HEIF_ERR_CHECK(heif_context_encode_image(ctx, baseImage, encoder, options, &baseHandle));
  if (exif != nullptr) {
    HEIF_ERR_CHECK(heif_context_add_exif_metadata(ctx, baseHandle, exif->data, exif->data_sz));
  }

  // encode the gain map image
  if (isUsingMultiChannelGainMap()) {
    gainmap_yuv_ext = convert_raw_input_to_ycbcr(gainmap_img, true /* use bt601 */);
    gainmap_img_yuv = gainmap_yuv_ext.get();
  }
  if (isUsingMultiChannelGainMap()) {
    HEIF_ERR_CHECK(map_fmt_to_heif_chroma_vars(gainmap_img_yuv->fmt, gainmap_img_yuv->w,
                                               gainmap_img_yuv->h, heif_img_fmt, chromaWd,
                                               chromaHt));
    HEIF_ERR_CHECK(heif_image_create(gainmap_img_yuv->w, gainmap_img_yuv->h, heif_colorspace_YCbCr,
                                     heif_img_fmt, &secondaryImage));
    HEIF_ERR_CHECK(fill_img_plane(secondaryImage, heif_channel_Y,
                                  gainmap_img_yuv->planes[UHDR_PLANE_Y], gainmap_img_yuv->w,
                                  gainmap_img_yuv->h, gainmap_img_yuv->stride[UHDR_PLANE_Y]));
    HEIF_ERR_CHECK(fill_img_plane(secondaryImage, heif_channel_Cb,
                                  gainmap_img_yuv->planes[UHDR_PLANE_U], chromaWd, chromaHt,
                                  gainmap_img_yuv->stride[UHDR_PLANE_U]));
    HEIF_ERR_CHECK(fill_img_plane(secondaryImage, heif_channel_Cr,
                                  gainmap_img_yuv->planes[UHDR_PLANE_V], chromaWd, chromaHt,
                                  gainmap_img_yuv->stride[UHDR_PLANE_V]));
    HEIF_ERR_CHECK(set_internal_color_format(encoder, heif_img_fmt));
  } else {
    HEIF_ERR_CHECK(heif_image_create(gainmap_img_yuv->w, gainmap_img_yuv->h,
                                     heif_colorspace_monochrome, heif_chroma_monochrome,
                                     &secondaryImage));
    HEIF_ERR_CHECK(fill_img_plane(secondaryImage, heif_channel_Y,
                                  gainmap_img_yuv->planes[UHDR_PLANE_Y], gainmap_img_yuv->w,
                                  gainmap_img_yuv->h, gainmap_img_yuv->stride[UHDR_PLANE_Y]));
    HEIF_ERR_CHECK(heif_encoder_set_parameter(encoder, "chroma", "420"));
  }
  gainmapNclx = heif_nclx_color_profile_alloc();
  if (!gainmapNclx) {
    status.error_code = UHDR_CODEC_MEM_ERROR;
    status.has_detail = 1;
    snprintf(status.detail, sizeof status.detail,
             "failed to allocate nclx color profile for base image");
    goto CleanUp;
  }
  HEIF_ERR_CHECK(
      heif_nclx_color_profile_set_color_primaries(gainmapNclx, heif_color_primaries_unspecified));
  HEIF_ERR_CHECK(heif_nclx_color_profile_set_transfer_characteristics(
      gainmapNclx, heif_transfer_characteristic_unspecified));
  HEIF_ERR_CHECK(heif_nclx_color_profile_set_matrix_coefficients(
      gainmapNclx, map_cg_cicp_matrix_coefficients((uhdr_color_gamut_t)UHDR_CG_BT_601)));
  gainmapNclx->full_range_flag = gainmap_img_yuv->range == UHDR_CR_FULL_RANGE ? true : false;
  options->output_nclx_profile = gainmapNclx;
  HEIF_ERR_CHECK(heif_image_set_nclx_color_profile(secondaryImage, gainmapNclx));

  heif_gain_map_metadata heif_metadata;
  uhdr_gainmap_metadata_frac iso_secondary_metadata;
  status =
      uhdr_gainmap_metadata_frac::gainmapMetadataFloatToFraction(metadata, &iso_secondary_metadata);
  if (status.error_code != UHDR_CODEC_OK) goto CleanUp;
  for (int c = 0; c < 3; ++c) {
    heif_metadata.gainMapMinN[c] = iso_secondary_metadata.gainMapMinN[c];
    heif_metadata.gainMapMinD[c] = iso_secondary_metadata.gainMapMinD[c];
    heif_metadata.gainMapMaxN[c] = iso_secondary_metadata.gainMapMaxN[c];
    heif_metadata.gainMapMaxD[c] = iso_secondary_metadata.gainMapMaxD[c];
    heif_metadata.gainMapGammaN[c] = iso_secondary_metadata.gainMapGammaN[c];
    heif_metadata.gainMapGammaD[c] = iso_secondary_metadata.gainMapGammaD[c];
    heif_metadata.baseOffsetN[c] = iso_secondary_metadata.baseOffsetN[c];
    heif_metadata.baseOffsetD[c] = iso_secondary_metadata.baseOffsetD[c];
    heif_metadata.alternateOffsetN[c] = iso_secondary_metadata.alternateOffsetN[c];
    heif_metadata.alternateOffsetD[c] = iso_secondary_metadata.alternateOffsetD[c];
  }
  heif_metadata.baseHdrHeadroomN = iso_secondary_metadata.baseHdrHeadroomN;
  heif_metadata.baseHdrHeadroomD = iso_secondary_metadata.baseHdrHeadroomD;
  heif_metadata.alternateHdrHeadroomN = iso_secondary_metadata.alternateHdrHeadroomN;
  heif_metadata.alternateHdrHeadroomD = iso_secondary_metadata.alternateHdrHeadroomD;
  heif_metadata.backwardDirection = iso_secondary_metadata.backwardDirection;
  heif_metadata.useBaseColorSpace = iso_secondary_metadata.useBaseColorSpace;
  // TODO: pass alternateIcc write color properties of tmap item
  // TODO: write clli of tmap item
  hdrNclx = heif_nclx_color_profile_alloc();
  if (!hdrNclx) {
    status.error_code = UHDR_CODEC_MEM_ERROR;
    status.has_detail = 1;
    snprintf(status.detail, sizeof status.detail,
             "failed to allocate nclx color profile for tmap image");
    goto CleanUp;
  }
  HEIF_ERR_CHECK(
      heif_nclx_color_profile_set_color_primaries(hdrNclx, map_cg_cicp_primaries(gainmap_img->cg)));
  HEIF_ERR_CHECK(heif_nclx_color_profile_set_transfer_characteristics(
      hdrNclx, map_ct_cicp_transfer_characteristics(gainmap_img->ct)));
  HEIF_ERR_CHECK(heif_nclx_color_profile_set_matrix_coefficients(
      hdrNclx, map_cg_cicp_matrix_coefficients(gainmap_img->cg)));
  hdrNclx->full_range_flag = gainmap_img->range == UHDR_CR_FULL_RANGE ? true : false;
  (void)alternateIcc;
  HEIF_ERR_CHECK(heif_encoder_set_lossy_quality(encoder, mMapCompressQuality));
  HEIF_ERR_CHECK(heif_context_encode_gain_map_image(ctx, secondaryImage, baseHandle, encoder,
                                                    options, &heif_metadata, hdrNclx,
                                                    &secondaryHandle));

  heif_context_write(ctx, &w, &writer);
  memcpy(dest->data, writer.data(), writer.size());
  dest->data_sz = dest->capacity = writer.size();

CleanUp:
  if (baseImage) heif_image_release(baseImage);
  baseImage = nullptr;
  if (baseHandle) heif_image_handle_release(baseHandle);
  baseHandle = nullptr;
  if (secondaryImage) heif_image_release(secondaryImage);
  secondaryImage = nullptr;
  if (secondaryHandle) heif_image_handle_release(secondaryHandle);
  secondaryHandle = nullptr;
  if (options) heif_encoding_options_free(options);
  options = nullptr;
  if (hdrNclx) heif_nclx_color_profile_free(hdrNclx);
  hdrNclx = nullptr;
  if (gainmapNclx) heif_nclx_color_profile_free(gainmapNclx);
  gainmapNclx = nullptr;
  if (sdrNclx) heif_nclx_color_profile_free(sdrNclx);
  sdrNclx = nullptr;
  if (encoder) heif_encoder_release(encoder);
  encoder = nullptr;
  if (ctx) heif_context_free(ctx);
  ctx = nullptr;

  return status;
}

}  // namespace ultrahdr

#endif
