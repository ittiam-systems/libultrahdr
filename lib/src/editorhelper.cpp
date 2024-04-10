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

#include <cstring>
#include <cstdint>

#include "ultrahdr/editorhelper.h"
#include "ultrahdr/gainmapmath.h"

namespace ultrahdr {

template <typename T>
void rotate_buffer_clockwise(T* src_buffer, T* dst_buffer, int src_w, int src_h, int src_stride,
                             int dst_stride, int degree) {
  if (degree == 90) {
    int dst_w = src_h;
    int dst_h = src_w;
    for (int i = 0; i < dst_h; i++) {
      for (int j = 0; j < dst_w; j++) {
        dst_buffer[i * dst_stride + j] = src_buffer[(src_h - j - 1) * src_stride + i];
      }
    }
  } else if (degree == 180) {
    int dst_w = src_w;
    int dst_h = src_h;
    for (int i = 0; i < dst_h; i++) {
      for (int j = 0; j < dst_w; j++) {
        dst_buffer[i * dst_stride + j] = src_buffer[(src_h - i - 1) * src_stride + (src_w - j - 1)];
      }
    }
  } else if (degree == 270) {
    int dst_w = src_h;
    int dst_h = src_w;
    for (int i = 0; i < dst_h; i++) {
      for (int j = 0; j < dst_w; j++) {
        dst_buffer[i * dst_stride + j] = src_buffer[j * src_stride + (src_w - i - 1)];
      }
    }
  }
}

template <typename T>
void mirror_buffer(T* src_buffer, T* dst_buffer, int src_w, int src_h, int src_stride,
                   int dst_stride, uhdr_mirror_direction_t direction) {
  if (direction == UHDR_MIRROR_VERTICAL) {
    for (int i = 0; i < src_h; i++) {
      memcpy(&dst_buffer[(src_h - i - 1) * dst_stride], &src_buffer[i * src_stride],
             src_w * sizeof(T));
    }
  } else if (direction == UHDR_MIRROR_HORIZONTAL) {
    for (int i = 0; i < src_h; i++) {
      for (int j = 0; j < src_w; j++) {
        dst_buffer[i * dst_stride + j] = src_buffer[i * src_stride + (src_w - j - 1)];
      }
    }
  }
}

template <typename T>
void resize_buffer(T* src_buffer, T* dst_buffer, int src_w, int src_h, int dst_w, int dst_h,
                   int src_stride, int dst_stride) {
  for (int i = 0; i < dst_h; i++) {
    for (int j = 0; j < dst_w; j++) {
      dst_buffer[i * dst_stride + j] =
          src_buffer[i * (src_h / dst_h) * src_stride + j * (src_w / dst_w)];
    }
  }
}

std::unique_ptr<uhdr_raw_image_ext_t> apply_rotate(uhdr_raw_image_t* src, int degree) {
  std::unique_ptr<uhdr_raw_image_ext_t> dst;

  if (degree == 90 || degree == 270) {
    dst = std::make_unique<uhdr_raw_image_ext_t>(src->fmt, src->cg, src->ct, src->range, src->h,
                                                 src->w, 1);
  } else if (degree == 180) {
    dst = std::make_unique<uhdr_raw_image_ext_t>(src->fmt, src->cg, src->ct, src->range, src->w,
                                                 src->h, 1);
  } else {
    return nullptr;
  }

  if (src->fmt == UHDR_IMG_FMT_24bppYCbCrP010) {
    uint16_t* src_buffer = static_cast<uint16_t*>(src->planes[UHDR_PLANE_Y]);
    uint16_t* dst_buffer = static_cast<uint16_t*>(dst->planes[UHDR_PLANE_Y]);
    rotate_buffer_clockwise(src_buffer, dst_buffer, src->w, src->h, src->stride[UHDR_PLANE_Y],
                            dst->stride[UHDR_PLANE_Y], degree);
    uint32_t* src_uv_buffer = static_cast<uint32_t*>(src->planes[UHDR_PLANE_UV]);
    uint32_t* dst_uv_buffer = static_cast<uint32_t*>(dst->planes[UHDR_PLANE_UV]);
    rotate_buffer_clockwise(src_uv_buffer, dst_uv_buffer, src->w / 2, src->h / 2,
                            src->stride[UHDR_PLANE_UV] / 2, dst->stride[UHDR_PLANE_UV] / 2, degree);
  } else if (src->fmt == UHDR_IMG_FMT_12bppYCbCr420 || src->fmt == UHDR_IMG_FMT_8bppYCbCr400) {
    uint8_t* src_buffer = static_cast<uint8_t*>(src->planes[UHDR_PLANE_Y]);
    uint8_t* dst_buffer = static_cast<uint8_t*>(dst->planes[UHDR_PLANE_Y]);
    rotate_buffer_clockwise(src_buffer, dst_buffer, src->w, src->h, src->stride[UHDR_PLANE_Y],
                            dst->stride[UHDR_PLANE_Y], degree);
    if (src->fmt == UHDR_IMG_FMT_12bppYCbCr420) {
      for (int i = 1; i < 3; i++) {
        src_buffer = static_cast<uint8_t*>(src->planes[i]);
        dst_buffer = static_cast<uint8_t*>(dst->planes[i]);
        rotate_buffer_clockwise(src_buffer, dst_buffer, src->w / 2, src->h / 2, src->stride[i],
                                dst->stride[i], degree);
      }
    }
  } else if (src->fmt == UHDR_IMG_FMT_32bppRGBA1010102 || src->fmt == UHDR_IMG_FMT_32bppRGBA8888) {
    uint32_t* src_buffer = static_cast<uint32_t*>(src->planes[UHDR_PLANE_PACKED]);
    uint32_t* dst_buffer = static_cast<uint32_t*>(dst->planes[UHDR_PLANE_PACKED]);
    rotate_buffer_clockwise(src_buffer, dst_buffer, src->w, src->h, src->stride[UHDR_PLANE_PACKED],
                            dst->stride[UHDR_PLANE_PACKED], degree);
  } else if (src->fmt == UHDR_IMG_FMT_64bppRGBAHalfFloat) {
    uint64_t* src_buffer = static_cast<uint64_t*>(src->planes[UHDR_PLANE_PACKED]);
    uint64_t* dst_buffer = static_cast<uint64_t*>(dst->planes[UHDR_PLANE_PACKED]);
    rotate_buffer_clockwise(src_buffer, dst_buffer, src->w, src->h, src->stride[UHDR_PLANE_PACKED],
                            dst->stride[UHDR_PLANE_PACKED], degree);
  }
  return std::move(dst);
}

std::unique_ptr<uhdr_raw_image_ext_t> apply_mirror(uhdr_raw_image_t* src,
                                                   uhdr_mirror_direction_t direction) {
  std::unique_ptr<uhdr_raw_image_ext_t> dst = std::make_unique<uhdr_raw_image_ext_t>(
      src->fmt, src->cg, src->ct, src->range, src->w, src->h, 1);

  if (src->fmt == UHDR_IMG_FMT_24bppYCbCrP010) {
    uint16_t* src_buffer = static_cast<uint16_t*>(src->planes[UHDR_PLANE_Y]);
    uint16_t* dst_buffer = static_cast<uint16_t*>(dst->planes[UHDR_PLANE_Y]);
    mirror_buffer(src_buffer, dst_buffer, src->w, src->h, src->stride[UHDR_PLANE_Y],
                  dst->stride[UHDR_PLANE_Y], direction);
    uint32_t* src_uv_buffer = static_cast<uint32_t*>(src->planes[UHDR_PLANE_UV]);
    uint32_t* dst_uv_buffer = static_cast<uint32_t*>(dst->planes[UHDR_PLANE_UV]);
    mirror_buffer(src_uv_buffer, dst_uv_buffer, src->w / 2, src->h / 2,
                  src->stride[UHDR_PLANE_UV] / 2, dst->stride[UHDR_PLANE_UV] / 2, direction);
  } else if (src->fmt == UHDR_IMG_FMT_12bppYCbCr420 || src->fmt == UHDR_IMG_FMT_8bppYCbCr400) {
    uint8_t* src_buffer = static_cast<uint8_t*>(src->planes[UHDR_PLANE_Y]);
    uint8_t* dst_buffer = static_cast<uint8_t*>(dst->planes[UHDR_PLANE_Y]);
    mirror_buffer(src_buffer, dst_buffer, src->w, src->h, src->stride[UHDR_PLANE_Y],
                  dst->stride[UHDR_PLANE_Y], direction);
    if (src->fmt == UHDR_IMG_FMT_12bppYCbCr420) {
      for (int i = 1; i < 3; i++) {
        src_buffer = static_cast<uint8_t*>(src->planes[i]);
        dst_buffer = static_cast<uint8_t*>(dst->planes[i]);
        mirror_buffer(src_buffer, dst_buffer, src->w / 2, src->h / 2, src->stride[i],
                      dst->stride[i], direction);
      }
    }
  } else if (src->fmt == UHDR_IMG_FMT_32bppRGBA1010102 || src->fmt == UHDR_IMG_FMT_32bppRGBA8888) {
    uint32_t* src_buffer = static_cast<uint32_t*>(src->planes[UHDR_PLANE_PACKED]);
    uint32_t* dst_buffer = static_cast<uint32_t*>(dst->planes[UHDR_PLANE_PACKED]);
    mirror_buffer(src_buffer, dst_buffer, src->w, src->h, src->stride[UHDR_PLANE_PACKED],
                  dst->stride[UHDR_PLANE_PACKED], direction);
  } else if (src->fmt == UHDR_IMG_FMT_64bppRGBAHalfFloat) {
    uint64_t* src_buffer = static_cast<uint64_t*>(src->planes[UHDR_PLANE_PACKED]);
    uint64_t* dst_buffer = static_cast<uint64_t*>(dst->planes[UHDR_PLANE_PACKED]);
    mirror_buffer(src_buffer, dst_buffer, src->w, src->h, src->stride[UHDR_PLANE_PACKED],
                  dst->stride[UHDR_PLANE_PACKED], direction);
  }
  return std::move(dst);
}

void apply_crop(uhdr_raw_image_t* src, int left, int top, int wd, int ht) {
  if (src->fmt == UHDR_IMG_FMT_24bppYCbCrP010) {
    uint16_t* src_buffer = static_cast<uint16_t*>(src->planes[UHDR_PLANE_Y]);
    src->planes[UHDR_PLANE_Y] = &src_buffer[top * src->stride[UHDR_PLANE_Y] + left];
    uint32_t* src_uv_buffer = static_cast<uint32_t*>(src->planes[UHDR_PLANE_UV]);
    src->planes[UHDR_PLANE_UV] =
        &src_uv_buffer[(top / 2) * (src->stride[UHDR_PLANE_UV] / 2) + (left / 2)];
  } else if (src->fmt == UHDR_IMG_FMT_12bppYCbCr420 || src->fmt == UHDR_IMG_FMT_8bppYCbCr400) {
    uint8_t* src_buffer = static_cast<uint8_t*>(src->planes[UHDR_PLANE_Y]);
    src->planes[UHDR_PLANE_Y] = &src_buffer[top * src->stride[UHDR_PLANE_Y] + left];
    if (src->fmt == UHDR_IMG_FMT_12bppYCbCr420) {
      for (int i = 1; i < 3; i++) {
        src->planes[i] = &src_buffer[(top / 2) * src->stride[i] + (left / 2)];
      }
    }
  } else if (src->fmt == UHDR_IMG_FMT_32bppRGBA1010102 || src->fmt == UHDR_IMG_FMT_32bppRGBA8888) {
    uint32_t* src_buffer = static_cast<uint32_t*>(src->planes[UHDR_PLANE_PACKED]);
    src->planes[UHDR_PLANE_PACKED] = &src_buffer[top * src->stride[UHDR_PLANE_PACKED] + left];
  } else if (src->fmt == UHDR_IMG_FMT_64bppRGBAHalfFloat) {
    uint64_t* src_buffer = static_cast<uint64_t*>(src->planes[UHDR_PLANE_PACKED]);
    src->planes[UHDR_PLANE_PACKED] = &src_buffer[top * src->stride[UHDR_PLANE_PACKED] + left];
  }
  src->w = wd;
  src->h = ht;
}

std::unique_ptr<uhdr_raw_image_ext_t> apply_resize(uhdr_raw_image_t* src, int dst_w, int dst_h) {
  std::unique_ptr<uhdr_raw_image_ext_t> dst = std::make_unique<uhdr_raw_image_ext_t>(
      src->fmt, src->cg, src->ct, src->range, dst_w, dst_h, 1);

  if (src->fmt == UHDR_IMG_FMT_24bppYCbCrP010) {
    uint16_t* src_buffer = static_cast<uint16_t*>(src->planes[UHDR_PLANE_Y]);
    uint16_t* dst_buffer = static_cast<uint16_t*>(dst->planes[UHDR_PLANE_Y]);
    resize_buffer(src_buffer, dst_buffer, src->w, src->h, dst->w, dst->h, src->stride[UHDR_PLANE_Y],
                  dst->stride[UHDR_PLANE_Y]);
    uint32_t* src_uv_buffer = static_cast<uint32_t*>(src->planes[UHDR_PLANE_UV]);
    uint32_t* dst_uv_buffer = static_cast<uint32_t*>(dst->planes[UHDR_PLANE_UV]);
    resize_buffer(src_uv_buffer, dst_uv_buffer, src->w / 4, src->h / 2, dst->w / 4, dst->h / 2,
                  src->stride[UHDR_PLANE_UV] / 2, dst->stride[UHDR_PLANE_UV] / 2);
  } else if (src->fmt == UHDR_IMG_FMT_12bppYCbCr420 || src->fmt == UHDR_IMG_FMT_8bppYCbCr400) {
    uint8_t* src_buffer = static_cast<uint8_t*>(src->planes[UHDR_PLANE_Y]);
    uint8_t* dst_buffer = static_cast<uint8_t*>(dst->planes[UHDR_PLANE_Y]);
    resize_buffer(src_buffer, dst_buffer, src->w, src->h, dst->w, dst->h, src->stride[UHDR_PLANE_Y],
                  dst->stride[UHDR_PLANE_Y]);
    if (src->fmt == UHDR_IMG_FMT_12bppYCbCr420) {
      for (int i = 1; i < 3; i++) {
        src_buffer = static_cast<uint8_t*>(src->planes[i]);
        dst_buffer = static_cast<uint8_t*>(dst->planes[i]);
        resize_buffer(src_buffer, dst_buffer, src->w / 2, src->h / 2, dst->w / 2, dst->h / 2,
                      src->stride[i], dst->stride[i]);
      }
    }
  } else if (src->fmt == UHDR_IMG_FMT_32bppRGBA1010102 || src->fmt == UHDR_IMG_FMT_32bppRGBA8888) {
    uint32_t* src_buffer = static_cast<uint32_t*>(src->planes[UHDR_PLANE_PACKED]);
    uint32_t* dst_buffer = static_cast<uint32_t*>(dst->planes[UHDR_PLANE_PACKED]);
    resize_buffer(src_buffer, dst_buffer, src->w, src->h, dst->w, dst->h,
                  src->stride[UHDR_PLANE_PACKED], dst->stride[UHDR_PLANE_PACKED]);
  } else if (src->fmt == UHDR_IMG_FMT_64bppRGBAHalfFloat) {
    uint64_t* src_buffer = static_cast<uint64_t*>(src->planes[UHDR_PLANE_PACKED]);
    uint64_t* dst_buffer = static_cast<uint64_t*>(dst->planes[UHDR_PLANE_PACKED]);
    resize_buffer(src_buffer, dst_buffer, src->w, src->h, dst->w, dst->h,
                  src->stride[UHDR_PLANE_PACKED], dst->stride[UHDR_PLANE_PACKED]);
  }
  return std::move(dst);
}

const float BT601RGBtoYUVMatrix[9] = {
    0.299,           0.587, 0.114, (-0.299 / 1.772), (-0.587 / 1.772), 0.5, 0.5, (-0.587 / 1.402),
    (-0.114 / 1.402)};

const float BT709RGBtoYUVMatrix[9] = {0.2126,
                                      0.7152,
                                      0.0722,
                                      (-0.2126 / 1.8556),
                                      (-0.7152 / 1.8556),
                                      0.5,
                                      0.5,
                                      (-0.7152 / 1.5748),
                                      (-0.0722 / 1.5748)};
const float BT2020RGBtoYUVMatrix[9] = {0.2627,
                                       0.6780,
                                       0.0593,
                                       (-0.2627 / 1.8814),
                                       (-0.6780 / 1.8814),
                                       0.5,
                                       0.5,
                                       (-0.6780 / 1.4746),
                                       (-0.0593 / 1.4746)};

std::unique_ptr<uhdr_raw_image_ext_t> convert_input_to_YCbCr(uhdr_raw_image_t* src) {
  std::unique_ptr<uhdr_raw_image_ext_t> dst = nullptr;
  const float* coeffs = BT601RGBtoYUVMatrix;

  if (src->fmt == UHDR_IMG_FMT_32bppRGBA1010102 || src->fmt == UHDR_IMG_FMT_32bppRGBA8888) {
    if (src->cg == UHDR_CG_BT_709) {
      coeffs = BT709RGBtoYUVMatrix;
    } else if (src->cg == UHDR_CG_BT_2100) {
      coeffs = BT2020RGBtoYUVMatrix;
    } else if (src->cg == UHDR_CG_DISPLAY_P3) {
      coeffs = BT601RGBtoYUVMatrix;
    } else {
      return dst;
    }
  }

  if (src->fmt == UHDR_IMG_FMT_32bppRGBA1010102) {
    dst = std::make_unique<uhdr_raw_image_ext_t>(UHDR_IMG_FMT_24bppYCbCrP010, src->cg, src->ct,
                                                 UHDR_CR_LIMITED_RANGE, src->w, src->h, 64);

    uint32_t* rgbData = static_cast<uint32_t*>(src->planes[UHDR_PLANE_PACKED]);
    unsigned int srcStride = src->stride[UHDR_PLANE_PACKED];

    uint16_t* yData = static_cast<uint16_t*>(dst->planes[UHDR_PLANE_Y]);
    uint16_t* uData = static_cast<uint16_t*>(dst->planes[UHDR_PLANE_UV]);
    uint16_t* vData = uData + 1;

    for (size_t i = 0; i < dst->h; i += 2) {
      for (size_t j = 0; j < dst->w; j += 2) {
        float r0 = float(rgbData[srcStride * i + j] & 0x3ff);
        float g0 = float((rgbData[srcStride * i + j] >> 10) & 0x3ff);
        float b0 = float((rgbData[srcStride * i + j] >> 20) & 0x3ff);

        float r1 = float(rgbData[srcStride * i + j + 1] & 0x3ff);
        float g1 = float((rgbData[srcStride * i + j + 1] >> 10) & 0x3ff);
        float b1 = float((rgbData[srcStride * i + j + 1] >> 20) & 0x3ff);

        float r2 = float(rgbData[srcStride * (i + 1) + j] & 0x3ff);
        float g2 = float((rgbData[srcStride * (i + 1) + j] >> 10) & 0x3ff);
        float b2 = float((rgbData[srcStride * (i + 1) + j] >> 20) & 0x3ff);

        float r3 = float(rgbData[srcStride * (i + 1) + j + 1] & 0x3ff);
        float g3 = float((rgbData[srcStride * (i + 1) + j + 1] >> 10) & 0x3ff);
        float b3 = float((rgbData[srcStride * (i + 1) + j + 1] >> 20) & 0x3ff);

        r0 /= 1023.0f;
        g0 /= 1023.0f;
        b0 /= 1023.0f;

        r1 /= 1023.0f;
        g1 /= 1023.0f;
        b1 /= 1023.0f;

        r2 /= 1023.0f;
        g2 /= 1023.0f;
        b2 /= 1023.0f;

        r3 /= 1023.0f;
        g3 /= 1023.0f;
        b3 /= 1023.0f;

        float y = coeffs[0] * r0 + coeffs[1] * g0 + coeffs[2] * b0;
        y = (y * 876.0f) + 64.0f + 0.5f;
        y = CLIP3(y, 64.0f, 940.0f);
        yData[dst->stride[UHDR_PLANE_Y] * i + j] = uint16_t(y) << 6;

        y = coeffs[0] * r1 + coeffs[1] * g1 + coeffs[2] * b1;
        y = (y * 876.0f) + 64.0f + 0.5f;
        y = CLIP3(y, 64.0f, 940.0f);
        yData[dst->stride[UHDR_PLANE_Y] * i + j + 1] = uint16_t(y) << 6;

        y = coeffs[0] * r2 + coeffs[1] * g2 + coeffs[2] * b2;
        y = (y * 876.0f) + 64.0f + 0.5f;
        y = CLIP3(y, 64.0f, 940.0f);
        yData[dst->stride[UHDR_PLANE_Y] * (i + 1) + j] = uint16_t(y) << 6;

        y = coeffs[0] * r3 + coeffs[1] * g3 + coeffs[2] * b3;
        y = (y * 876.0f) + 64.0f + 0.5f;
        y = CLIP3(y, 64.0f, 940.0f);
        yData[dst->stride[UHDR_PLANE_Y] * (i + 1) + j + 1] = uint16_t(y) << 6;

        float u0 = coeffs[3] * r0 + coeffs[4] * g0 + coeffs[5] * b0;
        float u1 = coeffs[3] * r1 + coeffs[4] * g1 + coeffs[5] * b1;
        float u2 = coeffs[3] * r2 + coeffs[4] * g2 + coeffs[5] * b2;
        float u3 = coeffs[3] * r3 + coeffs[4] * g3 + coeffs[5] * b3;

        float v0 = coeffs[6] * r0 + coeffs[7] * g0 + coeffs[8] * b0;
        float v1 = coeffs[6] * r1 + coeffs[7] * g1 + coeffs[8] * b1;
        float v2 = coeffs[6] * r2 + coeffs[7] * g2 + coeffs[8] * b2;
        float v3 = coeffs[6] * r3 + coeffs[7] * g3 + coeffs[8] * b3;

        u0 = (u0 + u1 + u2 + u3) / 4;
        v0 = (v0 + v1 + v2 + v3) / 4;

        u0 = (u0 * 896.0f) + 512.0f + 0.5f;
        v0 = (v0 * 896.0f) + 512.0f + 0.5f;

        u0 = CLIP3(u0, 64.0f, 960.0f);
        v0 = CLIP3(v0, 64.0f, 960.0f);

        uData[dst->stride[UHDR_PLANE_UV] * (i / 2) + j] = uint16_t(u0) << 6;
        vData[dst->stride[UHDR_PLANE_UV] * (i / 2) + j] = uint16_t(v0) << 6;
      }
    }
  } else if (src->fmt == UHDR_IMG_FMT_32bppRGBA8888) {
    dst = std::make_unique<uhdr_raw_image_ext_t>(UHDR_IMG_FMT_12bppYCbCr420, src->cg, src->ct,
                                                 UHDR_CR_FULL_RANGE, src->w, src->h, 64);
    uint32_t* rgbData = static_cast<uint32_t*>(src->planes[UHDR_PLANE_PACKED]);
    unsigned int srcStride = src->stride[UHDR_PLANE_PACKED];

    uint8_t* yData = static_cast<uint8_t*>(dst->planes[UHDR_PLANE_Y]);
    uint8_t* uData = static_cast<uint8_t*>(dst->planes[UHDR_PLANE_U]);
    uint8_t* vData = static_cast<uint8_t*>(dst->planes[UHDR_PLANE_V]);
    for (size_t i = 0; i < dst->h; i += 2) {
      for (size_t j = 0; j < dst->w; j += 2) {
        float r0 = float(rgbData[srcStride * i + j] & 0xff);
        float g0 = float((rgbData[srcStride * i + j] >> 8) & 0xff);
        float b0 = float((rgbData[srcStride * i + j] >> 16) & 0xff);

        float r1 = float(rgbData[srcStride * i + (j + 1)] & 0xff);
        float g1 = float((rgbData[srcStride * i + (j + 1)] >> 8) & 0xff);
        float b1 = float((rgbData[srcStride * i + (j + 1)] >> 16) & 0xff);

        float r2 = float(rgbData[srcStride * (i + 1) + j] & 0xff);
        float g2 = float((rgbData[srcStride * (i + 1) + j] >> 8) & 0xff);
        float b2 = float((rgbData[srcStride * (i + 1) + j] >> 16) & 0xff);

        float r3 = float(rgbData[srcStride * (i + 1) + (j + 1)] & 0xff);
        float g3 = float((rgbData[srcStride * (i + 1) + (j + 1)] >> 8) & 0xff);
        float b3 = float((rgbData[srcStride * (i + 1) + (j + 1)] >> 16) & 0xff);

        r0 /= 255.0f;
        g0 /= 255.0f;
        b0 /= 255.0f;

        r1 /= 255.0f;
        g1 /= 255.0f;
        b1 /= 255.0f;

        r2 /= 255.0f;
        g2 /= 255.0f;
        b2 /= 255.0f;

        r3 /= 255.0f;
        g3 /= 255.0f;
        b3 /= 255.0f;

        float y = coeffs[0] * r0 + coeffs[1] * g0 + coeffs[2] * b0;
        y = y * 255.0f + 0.5f;
        y = CLIP3(y, 0.0f, 255.0f);
        yData[dst->stride[UHDR_PLANE_Y] * i + j] = uint8_t(y);

        y = coeffs[0] * r1 + coeffs[1] * g1 + coeffs[2] * b1;
        y = y * 255.0f + 0.5f;
        y = CLIP3(y, 0.0f, 255.0f);
        yData[dst->stride[UHDR_PLANE_Y] * i + j + 1] = uint8_t(y);

        y = coeffs[0] * r2 + coeffs[1] * g2 + coeffs[2] * b2;
        y = y * 255.0f + 0.5f;
        y = CLIP3(y, 0.0f, 255.0f);
        yData[dst->stride[UHDR_PLANE_Y] * (i + 1) + j] = uint8_t(y);

        y = coeffs[0] * r3 + coeffs[1] * g3 + coeffs[2] * b3;
        y = y * 255.0f + 0.5f;
        y = CLIP3(y, 0.0f, 255.0f);
        yData[dst->stride[UHDR_PLANE_Y] * (i + 1) + j + 1] = uint8_t(y);

        float u0 = coeffs[3] * r0 + coeffs[4] * g0 + coeffs[5] * b0;
        float u1 = coeffs[3] * r1 + coeffs[4] * g1 + coeffs[5] * b1;
        float u2 = coeffs[3] * r2 + coeffs[4] * g2 + coeffs[5] * b2;
        float u3 = coeffs[3] * r3 + coeffs[4] * g3 + coeffs[5] * b3;

        float v0 = coeffs[6] * r0 + coeffs[7] * g0 + coeffs[8] * b0;
        float v1 = coeffs[6] * r1 + coeffs[7] * g1 + coeffs[8] * b1;
        float v2 = coeffs[6] * r2 + coeffs[7] * g2 + coeffs[8] * b2;
        float v3 = coeffs[6] * r3 + coeffs[7] * g3 + coeffs[8] * b3;

        u0 = (u0 + u1 + u2 + u3) / 4;
        v0 = (v0 + v1 + v2 + v3) / 4;

        u0 = u0 * 255.0f + 0.5 + 128.0f;
        v0 = v0 * 255.0f + 0.5 + 128.0f;

        u0 = CLIP3(u0, 0.0f, 255.0f);
        v0 = CLIP3(v0, 0.0f, 255.0f);

        uData[dst->stride[UHDR_PLANE_U] * (i / 2) + (j / 2)] = uint8_t(u0);
        vData[dst->stride[UHDR_PLANE_V] * (i / 2) + (j / 2)] = uint8_t(v0);
      }
    }
  } else if (src->fmt == UHDR_IMG_FMT_12bppYCbCr420) {
    dst = std::make_unique<ultrahdr::uhdr_raw_image_ext_t>(src->fmt, src->cg, src->ct, src->range,
                                                           src->w, src->h, 64);

    uint8_t* y_dst = static_cast<uint8_t*>(dst->planes[UHDR_PLANE_Y]);
    uint8_t* y_src = static_cast<uint8_t*>(src->planes[UHDR_PLANE_Y]);
    uint8_t* u_dst = static_cast<uint8_t*>(dst->planes[UHDR_PLANE_U]);
    uint8_t* u_src = static_cast<uint8_t*>(src->planes[UHDR_PLANE_U]);
    uint8_t* v_dst = static_cast<uint8_t*>(dst->planes[UHDR_PLANE_V]);
    uint8_t* v_src = static_cast<uint8_t*>(src->planes[UHDR_PLANE_V]);

    // copy y
    for (size_t i = 0; i < src->h; i++) {
      memcpy(y_dst, y_src, src->w);
      y_dst += dst->stride[UHDR_PLANE_Y];
      y_src += src->stride[UHDR_PLANE_Y];
    }
    // copy cb & cr
    for (size_t i = 0; i < src->h / 2; i++) {
      memcpy(u_dst, u_src, src->w / 2);
      memcpy(v_dst, v_src, src->w / 2);
      u_dst += dst->stride[UHDR_PLANE_U];
      v_dst += dst->stride[UHDR_PLANE_V];
      u_src += src->stride[UHDR_PLANE_U];
      v_src += src->stride[UHDR_PLANE_V];
    }
  } else if (src->fmt == UHDR_IMG_FMT_24bppYCbCrP010) {
    dst = std::make_unique<ultrahdr::uhdr_raw_image_ext_t>(src->fmt, src->cg, src->ct, src->range,
                                                           src->w, src->h, 64);

    int bpp = 2;
    uint8_t* y_dst = static_cast<uint8_t*>(dst->planes[UHDR_PLANE_Y]);
    uint8_t* y_src = static_cast<uint8_t*>(src->planes[UHDR_PLANE_Y]);
    uint8_t* uv_dst = static_cast<uint8_t*>(dst->planes[UHDR_PLANE_UV]);
    uint8_t* uv_src = static_cast<uint8_t*>(src->planes[UHDR_PLANE_UV]);

    // copy y
    for (size_t i = 0; i < src->h; i++) {
      memcpy(y_dst, y_src, src->w * bpp);
      y_dst += (dst->stride[UHDR_PLANE_Y] * bpp);
      y_src += (src->stride[UHDR_PLANE_Y] * bpp);
    }
    // copy cbcr
    for (size_t i = 0; i < src->h / 2; i++) {
      memcpy(uv_dst, uv_src, src->w * bpp);
      uv_dst += (dst->stride[UHDR_PLANE_UV] * bpp);
      uv_src += (src->stride[UHDR_PLANE_UV] * bpp);
    }
  }
  return std::move(dst);
}

}  // namespace ultrahdr
