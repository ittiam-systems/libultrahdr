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

#ifndef ULTRAHDR_HEIFR_H
#define ULTRAHDR_HEIFR_H

#ifdef UHDR_ENABLE_HEIF

#include <cfloat>

#include "ultrahdr_api.h"
#include "ultrahdr/ultrahdrcommon.h"
#include "ultrahdr/icc.h"

#include "libheif/heif.h"
#include "libheif/heif_experimental.h"

namespace ultrahdr {

class HeifR : public UltraHdr {
 public:
  HeifR(void* uhdrGLESCtxt = nullptr,
        int mapDimensionScaleFactor = kMapDimensionScaleFactorAndroidDefault,
        int mapCompressQuality = kMapCompressQualityAndroidDefault,
        bool useMultiChannelGainMap = kUseMultiChannelGainMapAndroidDefault,
        float gamma = kGainMapGammaDefault,
        uhdr_enc_preset_t preset = kEncSpeedPresetAndroidDefault, float minContentBoost = FLT_MIN,
        float maxContentBoost = FLT_MAX, float targetDispPeakBrightness = -1.0f,
        uhdr_codec_t codec = UHDR_CODEC_AVIF);

  /*!\brief Encode API-0.
   *
   * Create ultrahdr heif/avif image from raw hdr intent.
   *
   * Input hdr image is tonemapped to sdr image. A gainmap coefficient is computed between hdr and
   * sdr intent. sdr intent and gain map coefficient are compressed using heif/avif encoding.
   * compressed sdr intent is signalled as primary item and compressed gainmap is signalled as a
   * tone map derived item
   *
   * \param[in]       hdr_intent        hdr intent raw input image descriptor
   * \param[in, out]  dest              output image descriptor to store compressed ultrahdr image
   * \param[in]       quality           quality factor for sdr intent heif/avif compression
   * \param[in]       exif              optional exif metadata that needs to be inserted in
   *                                    compressed output
   *
   * \return uhdr_error_info_t #UHDR_CODEC_OK if operation succeeds, uhdr_codec_err_t otherwise.
   */
  uhdr_error_info_t encodeHEIFR(uhdr_raw_image_t* hdr_intent, uhdr_compressed_image_t* dest,
                                int quality, uhdr_mem_block_t* exif);

  /*!\brief Encode API-1.
   *
   * Create ultrahdr heif/avif image from raw hdr intent and raw sdr intent.
   *
   * A gainmap coefficient is computed between hdr and sdr intent. sdr intent and gain map
   * coefficient are compressed using heif/avif encoding. compressed sdr intent is signalled as
   * primary item and compressed gainmap is signalled as a tone map derived item
   * NOTE: Color transfer of sdr intent is expected to be sRGB.
   *
   * \param[in]       hdr_intent        hdr intent raw input image descriptor
   * \param[in]       sdr_intent        sdr intent raw input image descriptor
   * \param[in, out]  dest              output image descriptor to store compressed ultrahdr image
   * \param[in]       quality           quality factor for sdr intent heif/avif compression
   * \param[in]       exif              optional exif metadata that needs to be inserted in
   *                                    compressed output
   *
   * \return uhdr_error_info_t #UHDR_CODEC_OK if operation succeeds, uhdr_codec_err_t otherwise.
   */
  uhdr_error_info_t encodeHEIFR(uhdr_raw_image_t* hdr_intent, uhdr_raw_image_t* sdr_intent,
                                uhdr_compressed_image_t* dest, int quality, uhdr_mem_block_t* exif);

 private:
  /*!\brief Encode API.
   *
   * Create ultrahdr image from sdr intent, gainmap image and gainmap metadata
   *
   * A gainmap coefficient is computed between hdr and sdr intent. sdr intent and gain map
   * coefficient are compressed using heif/avif encoding. compressed sdr intent is signalled as
   * primary item and compressed gainmap is signalled as a tone map derived item
   *
   * \param[in]       sdr_intent        sdr intent raw input image descriptor
   * \param[in]       gainmap_img       gainmap raw image descriptor
   * \param[in]       metadata          gainmap metadata descriptor
   * \param[in, out]  dest              output image descriptor to store compressed ultrahdr image
   * \param[in]       quality           quality factor for sdr intent heif/avif compression
   * \param[in]       exif              optional exif metadata that needs to be inserted in
   *                                    compressed output
   * \param[in]       baseIcc           base image icc
   * \param[in]       alternateIcc      alternate image icc
   *
   * \return uhdr_error_info_t #UHDR_CODEC_OK if operation succeeds, uhdr_codec_err_t otherwise.
   */
  uhdr_error_info_t encodeHEIFR(uhdr_raw_image_t* sdr_intent, uhdr_raw_image_t* gainmap_img,
                                uhdr_gainmap_metadata_ext_t* metadata,
                                uhdr_compressed_image_t* dest, int quality, uhdr_mem_block_t* exif,
                                DataStruct* baseIcc, DataStruct* alternateIcc);

  uhdr_codec_t mCodec;
};

}  // namespace ultrahdr

#endif

#endif  // ULTRAHDR_HEIFR_H
