/*********************************************************************/
/* DataUnit.h                                                        */
/* Author: James Weaver and Galen Reich                              */
/* This version July 2020                                            */
/*                                                                   */
/* Declares stuff related to VC2 Data Units.                         */
/* Copyright (c) BBC 2011-2020 -- For license see the LICENSE file   */
/*********************************************************************/

#ifndef DATAUNIT_17JUN15
#define DATAUNIT_17JUN15

#include <iosfwd>

#include "Utils.h"
#include "Picture.h"
#include "Slices.h"
#include "WaveletTransform.h"
#include "VLC.h"

enum DataUnitType {
  UNKNOWN_DATA_UNIT,
  SEQUENCE_HEADER,
  END_OF_SEQUENCE,
  AUXILIARY_DATA,
  PADDING_DATA,
  HQ_PICTURE,
  LD_PICTURE,
  HQ_FRAGMENT,
  LD_FRAGMENT
};

class DataUnit {
  public:
    DataUnit();

    int length();
    DataUnitType type;
    Bytes next_parse_offset;
    Bytes prev_parse_offset;

    friend std::istream& operator >> (std::istream& stream, DataUnit &d);

};

class PicturePreamble;
class FragmentSlices;

class Fragment {
 public:
  Fragment()
    : mNSlices(0)
    , mSliceOffsetX(0)
    , mSliceOffsetY(0) {}

  int n_slices() { return mNSlices; }
  int slice_offset_x() { return mSliceOffsetX; }
  int slice_offset_y() { return mSliceOffsetY; }

  friend std::istream& operator >> (std::istream& stream, Fragment &f);

 protected:
  int mNSlices;
  int mSliceOffsetX;
  int mSliceOffsetY;
};

class WrappedPicture {
  public:
    WrappedPicture(const unsigned long picture_number,
                   const WaveletKernel wavelet_kernel,
                   const int depth,
                   const int slices_x,
                   const int slices_y,
                   const int slice_prefix,
                   const int slice_size_scalar,
                   const Slices &slices);

    WrappedPicture(const unsigned long picture_number,
                   const WaveletKernel wavelet_kernel,
                   const int depth,
                   const int slices_x,
                   const int slices_y,
                   const utils::Rational slice_bytes,
                   const Slices &slices);
    unsigned long picture_number;
    WaveletKernel wavelet_kernel;
    int depth;
    int slices_x;
    int slices_y;
    int slice_prefix;
    int slice_size_scalar;
    utils::Rational slice_bytes;
    Slices slices;
};

enum FrameRate { FR_UNSET = -1, FR0, FR24000_1001, FR24, FR25, FR30000_1001, FR30, FR50, FR60000_1001, FR60, FR15000_1001, FR25_2, FR48, FR48_1001, FR96, FR100, FR120_1001, FR120 };
//enum ColourFormat {CF_UNSET=-1, CF444, CF422, CF420};
enum PixelAspectRatio { AR_UNSET = -1, AR0, AR1_1, AR10_11, AR12_11, AR40_33, AR16_11, AR4_3 };
enum SignalRange { SR_UNSET = -1, SR0, SR8_FULL, SR8, SR10, SR12, SR10_FULL, SR12_FULL, SR16, SR16_FULL};
enum ColorSpec { CS_UNSET = -1, CS_CUSTOM, CS_SDTV_525, CS_SDTV_625, CS_HDTV, CS_D_CINEMA, CS_UHDTV, CS_HDRTV_PQ, CS_HDRTV_HLG};
const FrameRate MAX_V2_FRAMERATE = FR48;
enum Profile { PROFILE_UNKNOWN, PROFILE_LD, PROFILE_HQ };

class SequenceHeader {
public:
  SequenceHeader();
  SequenceHeader( Profile profile, 
                  int height,
                  int width,
                  ColourFormat chromaFormat,
                  bool interlace,
                  FrameRate frameRate,
                  bool topFieldFirst,
                  int bitdepth,

                  // Optional Args
                  PixelAspectRatio pixelAspectRatio = AR_UNSET,
                  int cleanWidth = -1,
                  int cleanHeight = -1,
                  int leftOffset = -1,
                  int topOffset = -1,
                  
                  ColorSpec colorSpec = CS_UNSET,
                  int colorPrimaries = 0, // HDTV
                  int colorMatrix = 0, // HDTV
                  int transferFunction = 0, //TV Gamma
                  
                  bool use_v3 = false
                );

  int major_version;
  int minor_version;
  Profile profile;
  int width;
  int height;
  ColourFormat chromaFormat;
  bool interlace;
  FrameRate frameRate;
  bool topFieldFirst;  
  unsigned int frameRateNumer;
  unsigned int frameRateDenom;
  int bitdepth;
  unsigned int lumaExcursion;
  unsigned int lumaOffset;
  unsigned int colorDiffExcursion;
  unsigned int colorDiffOffset;

  PixelAspectRatio pixelAspectRatio;
  unsigned int pixelAspectRatioNumer;
  unsigned int pixelAspectRatioDenom;
  int cleanWidth;
  int cleanHeight;
  int leftOffset;
  int topOffset;
  
  ColorSpec colorSpec;
  int colorPrimaries;
  int colorMatrix;
  int transferFunction;
};

SequenceHeader getDefaultSourceParameters(const int base_video_format_index);

struct video_format {

  video_format();
  video_format(const SequenceHeader &fmt);

  int major_version;
  int minor_version;
  int profile;
  int level;
  int base_video_format;
  
  bool custom_dimensions_flag;
  int frame_width;
  int frame_height;
  bool custom_color_diff_format_flag;
  int color_diff_format;
  bool custom_scan_format_flag;
  int source_sampling;
  bool custom_frame_rate_flag;
  FrameRate frame_rate;
  unsigned int frame_rate_numer;
  unsigned int frame_rate_denom;    
  bool custom_pixel_aspect_ratio_flag;
  int pixel_aspect_ratio;      
  unsigned int pixel_aspect_ratio_numer;
  unsigned int pixel_aspect_ratio_denom;    
  bool custom_clean_area_flag;
  int clean_width;
  int clean_height;
  int left_offset;
  int top_offset;
  bool custom_signal_range_flag;
  int bitdepth;
  unsigned int luma_excursion;
  unsigned int luma_offset;
  unsigned int color_diff_excursion;
  unsigned int color_diff_offset;
  bool custom_color_spec_flag;
  int color_spec;
  bool custom_color_primaries_flag;
  int color_primaries;
  bool custom_color_matrix_flag;
  int color_matrix;
  bool custom_transfer_function_flag;
  int transfer_function;

  bool top_field_first;
};

class PictureHeader {
public:
  unsigned long picture_number;
};

std::istream &operator >> (std::istream &stream, PictureHeader &h);

class PicturePreamble {
 public:
  PicturePreamble();

  WaveletKernel wavelet_kernel;
  WaveletKernel wavelet_kernel_ho;
  int depth;
  int depth_ho;
  int slices_x;
  int slices_y;
  int slice_prefix;
  int slice_size_scalar;
  utils::Rational slice_bytes;
};

namespace dataunitio {
  using sliceio::highQualityCBR;
  using sliceio::highQualityVBR;

  int &majorVersionNum(std::ios_base& stream);

  std::istream& lowDelay(std::istream& stream);

  std::istream& synchronise(std::istream& stream);

  std::ostream& start_sequence(std::ostream& stream);
  std::ostream& end_sequence(std::ostream& stream);

  class fragmentedPictures {
  public:
    fragmentedPictures(const int f) : mFL(f) {};
    void operator() (std::ios_base &stream) const;
  private:
    const int mFL;
  };
};

std::ostream& operator << (std::ostream& stream, dataunitio::fragmentedPictures arg);

std::ostream& operator << (std::ostream& stream, const WrappedPicture& d);

std::ostream& operator << (std::ostream& stream, const SequenceHeader& s);

std::ostream& operator << (std::ostream& stream, const DataUnitType& t);

std::ostream& operator << (std::ostream& stream, const FrameRate& t);

std::istream& operator >> (std::istream& stream, DataUnit &d);

std::istream& operator >> (std::istream& stream, Fragment &d);

std::istream& operator >> (std::istream& stream, SequenceHeader &hdr);

std::istream& operator >> (std::istream& stream, PicturePreamble &hdr);

#endif /* DATAUNIT_17JUN15 */
