/*********************************************************************/
/* DataUnit.cpp                                                      */
/* Author: James Weaver                                              */
/* This version 17th June 2015                                       */
/*                                                                   */
/* Defines stuff relating to data units                              */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

#include <iostream> //For cin, cout, cerr
#include <sstream>

#include "DataUnit.h"
#include "Slices.h"
#include "VLC.h"
#include "Utils.h"

DataUnit::DataUnit()
  : type (UNKNOWN_DATA_UNIT)
  , next_parse_offset(4)
  , prev_parse_offset(4)
  {}

int DataUnit::length() { return next_parse_offset-13; }

WrappedPicture::WrappedPicture(const unsigned long p,
                               const WaveletKernel w,
                               const int d,
                               const int x,
                               const int y,
                               const int sp,
                               const int ss,
                               const Slices &s)
  : picture_number (p),
    wavelet_kernel (w),
    depth (d),
    slices_x (x),
    slices_y (y),
    slice_prefix (sp),
    slice_size_scalar (ss),
    slice_bytes (),
    slices (s) {
}

WrappedPicture::WrappedPicture(const unsigned long p,
                               const WaveletKernel w,
                               const int d,
                               const int x,
                               const int y,
                               const utils::Rational sb,
                               const Slices &s)
  : picture_number (p),
    wavelet_kernel (w),
    depth (d),
    slices_x (x),
    slices_y (y),
    slice_prefix (0),
    slice_size_scalar (0),
    slice_bytes(sb),
    slices (s) {
}

namespace {
  long& prev_parse_offset(std::ios_base& stream) {
    static const int i = std::ios_base::xalloc();
    return stream.iword(i);
  }

  long& major_version_number(std::ios_base& stream) {
    static const int i = std::ios_base::xalloc();
    return stream.iword(i);
  }

  long& fragment_length(std::ios_base& stream) {
    static const int i = std::ios_base::xalloc();
    return stream.iword(i);
  }
};

class ParseInfoIO {
public:
  ParseInfoIO(const DataUnitType du)
    : type (du)
    , next_parse_offset (0) {}
  ParseInfoIO(const DataUnitType du, const unsigned int data_size)
    : type (du)
    , next_parse_offset (data_size + 13) {}
  
  DataUnitType type;
  unsigned int next_parse_offset;

  unsigned char parse_code() const {
    switch (type) {
    case SEQUENCE_HEADER:
      return 0x00;
    case END_OF_SEQUENCE:
      return 0x10;
    case LD_PICTURE:
      return 0xC8;
    case HQ_PICTURE:
      return 0xE8;
    case HQ_FRAGMENT:
      return 0xEC;
    case LD_FRAGMENT:
      return 0xCC;
    default:
      return 0x20;
    }
  }
};

std::ostream& operator << (std::ostream& stream, const ParseInfoIO &piio) {
  stream << Bytes(1, 0x42)
         << Bytes(1, 0x42)
         << Bytes(1, 0x43)
         << Bytes(1, 0x44)
         << Bytes(1, piio.parse_code())
         << Bytes(4, piio.next_parse_offset)
         << Bytes(4, prev_parse_offset(stream));

  prev_parse_offset(stream) = piio.next_parse_offset;
  return stream;
}

std::ostream& LDWrappedPictureIO(std::ostream& stream, const WrappedPicture& d) {
  if (fragment_length(stream) == 0) {
    std::ostringstream ss;
    ss.copyfmt(stream);

    // Picture Header
    ss << Bytes(4, d.picture_number);

    // Transform Params
    ss << vlc::unbounded
       << UnsignedVLC(d.wavelet_kernel)
       << UnsignedVLC(d.depth);

    if (major_version_number(stream) >= 3) {
      ss << Boolean(false)  // asym_transform_index_flag
         << Boolean(false); // asym_transform_flag
    }

    ss << UnsignedVLC(d.slices_x)
       << UnsignedVLC(d.slices_y)
       << UnsignedVLC(d.slice_bytes.numerator)
       << UnsignedVLC(d.slice_bytes.denominator)
       << Boolean(false)
       << vlc::align;

    // Transform Data
    ss << d.slices;

    stream << ParseInfoIO(LD_PICTURE, ss.str().size());

    return (stream << ss.str());
  } else {
    {
      std::ostringstream ss;
      ss.copyfmt(stream);

      // Transform Params
      ss << vlc::unbounded
         << UnsignedVLC(d.wavelet_kernel)
         << UnsignedVLC(d.depth);
    
      ss << Boolean(false)  // asym_transform_index_flag
         << Boolean(false); // asym_transform_flag
    
      ss << UnsignedVLC(d.slices_x)
         << UnsignedVLC(d.slices_y)
         << UnsignedVLC(d.slice_bytes.numerator)
         << UnsignedVLC(d.slice_bytes.denominator)
         << Boolean(false)
         << vlc::align;

      stream << ParseInfoIO(LD_FRAGMENT, ss.str().size() + 8)
             << Bytes(4, d.picture_number)
             << Bytes(2, ss.str().size())
             << Bytes(2, 0);
      stream << ss.str();
    }

    int slice_x = 0;
    int slice_y = 0;
    int slice_offset_x = 0;
    int slice_offset_y = 0;
    int nslices = 0;
    std::ostringstream fragstream;
    fragstream.copyfmt(stream);
    const Array2D& bytes = *sliceio::SliceSizes(stream);
    const bool bytes_valid = (sliceio::SliceSizes(stream)!=0);
    const PictureArray& yuvSlices = d.slices.yuvSlices;
    const Array2D& qIndices = d.slices.qIndices;
    const int waveletDepth = d.slices.waveletDepth;
    while (slice_y*d.slices_x + slice_x < d.slices_y*d.slices_x) {
      std::ostringstream slicestream;
      slicestream.copyfmt(stream);
      if (bytes_valid) slicestream << sliceio::setBytes(bytes[slice_y][slice_x]);
      slicestream << Slice(yuvSlices[slice_y][slice_x], waveletDepth, qIndices[slice_y][slice_x]);

      if (nslices > 0 && (int)(fragstream.str().size() + slicestream.str().size()) > fragment_length(stream)) {
        stream << ParseInfoIO(LD_FRAGMENT, fragstream.str().size() + 12)
               << Bytes(4, d.picture_number)
               << Bytes(2, fragstream.str().size())
               << Bytes(2, nslices)
               << Bytes(2, slice_offset_x)
               << Bytes(2, slice_offset_y)
               << fragstream.str();
        slice_offset_x = slice_x;
        slice_offset_y = slice_y;
        nslices = 0;
        fragstream.str("");
        fragstream.clear();
        fragstream.copyfmt(stream);
      }

      fragstream << slicestream.str();
      nslices++;
      slice_x++;
      if (slice_x == d.slices_x) {
        slice_x = 0;
        slice_y++;
      }
    }
    stream << ParseInfoIO(LD_FRAGMENT, fragstream.str().size() + 12)
           << Bytes(4, d.picture_number)
           << Bytes(2, fragstream.str().size())
           << Bytes(2, nslices)
           << Bytes(2, slice_offset_x)
           << Bytes(2, slice_offset_y)
           << fragstream.str();
    return stream;
  }
}

std::ostream& HQWrappedPictureIO(std::ostream& stream, const WrappedPicture& d) {
  if (fragment_length(stream) == 0) {
    std::ostringstream ss;
    ss.copyfmt(stream);

    // Picture Header
    ss << Bytes(4, d.picture_number);

    // Transform Params
    ss << vlc::unbounded
       << UnsignedVLC(d.wavelet_kernel)
       << UnsignedVLC(d.depth);
    
    if (major_version_number(stream) >= 3) {
      ss << Boolean(false)  // asym_transform_index_flag
         << Boolean(false); // asym_transform_flag
    }
    
    ss << UnsignedVLC(d.slices_x)
       << UnsignedVLC(d.slices_y)
       << UnsignedVLC(d.slice_prefix)
       << UnsignedVLC(d.slice_size_scalar)
       << Boolean(false)
       << vlc::align;

    // Transform Data
    ss << d.slices;

    stream << ParseInfoIO(HQ_PICTURE, ss.str().size());

    return (stream << ss.str());
  } else {
    {
      std::ostringstream ss;
      ss.copyfmt(stream);

      // Transform Params
      ss << vlc::unbounded
         << UnsignedVLC(d.wavelet_kernel)
         << UnsignedVLC(d.depth);
    
      ss << Boolean(false)  // asym_transform_index_flag
         << Boolean(false); // asym_transform_flag
    
      ss << UnsignedVLC(d.slices_x)
         << UnsignedVLC(d.slices_y)
         << UnsignedVLC(d.slice_prefix)
         << UnsignedVLC(d.slice_size_scalar)
         << Boolean(false)
         << vlc::align;

      stream << ParseInfoIO(HQ_FRAGMENT, ss.str().size() + 8)
             << Bytes(4, d.picture_number)
             << Bytes(2, ss.str().size())
             << Bytes(2, 0);
      stream << ss.str();
    }

    int slice_x = 0;
    int slice_y = 0;
    int slice_offset_x = 0;
    int slice_offset_y = 0;
    int nslices = 0;
    std::ostringstream fragstream;
    fragstream.copyfmt(stream);
    const Array2D& bytes = *sliceio::SliceSizes(stream);
    const bool bytes_valid = (sliceio::SliceSizes(stream)!=0);
    const PictureArray& yuvSlices = d.slices.yuvSlices;
    const Array2D& qIndices = d.slices.qIndices;
    const int waveletDepth = d.slices.waveletDepth;
    while (slice_y*d.slices_x + slice_x < d.slices_y*d.slices_x) {
      std::ostringstream slicestream;
      slicestream.copyfmt(stream);
      if (bytes_valid) slicestream << sliceio::setBytes(bytes[slice_y][slice_x]);
      slicestream << Slice(yuvSlices[slice_y][slice_x], waveletDepth, qIndices[slice_y][slice_x]);

      if (nslices > 0 && (int)(fragstream.str().size() + slicestream.str().size()) > fragment_length(stream)) {
        stream << ParseInfoIO(HQ_FRAGMENT, fragstream.str().size() + 12)
               << Bytes(4, d.picture_number)
               << Bytes(2, fragstream.str().size())
               << Bytes(2, nslices)
               << Bytes(2, slice_offset_x)
               << Bytes(2, slice_offset_y)
               << fragstream.str();
        slice_offset_x = slice_x;
        slice_offset_y = slice_y;
        nslices = 0;
        fragstream.str("");
      }

      fragstream << slicestream.str();
      nslices++;
      slice_x++;
      if (slice_x == d.slices_x) {
        slice_x = 0;
        slice_y++;
      }
    }
    stream << ParseInfoIO(HQ_FRAGMENT, fragstream.str().size() + 12)
           << Bytes(4, d.picture_number)
           << Bytes(2, fragstream.str().size())
           << Bytes(2, nslices)
           << Bytes(2, slice_offset_x)
           << Bytes(2, slice_offset_y)
           << fragstream.str();
    return stream;
  }
}

std::ostream& operator << (std::ostream& stream, const WrappedPicture& d) {
  switch (sliceio::sliceIOMode(stream)) {
    case sliceio::LD:
      return LDWrappedPictureIO(stream, d);
      break;
    case sliceio::HQVBR:
    case sliceio::HQCBR:
      return HQWrappedPictureIO(stream, d);
      break;
    default:
      throw std::logic_error("DataUnitIO: Unknown Output Format");
  }
}

std::ostream& dataunitio::start_sequence(std::ostream& stream) {
  prev_parse_offset(stream) = 0;
  return stream;
}

std::ostream& dataunitio::end_sequence(std::ostream& stream) {
  stream << ParseInfoIO(END_OF_SEQUENCE);
  prev_parse_offset(stream) = 0;
  return stream;
}

SequenceHeader::SequenceHeader()
  : major_version(1)
  , minor_version(0)
  , profile (PROFILE_UNKNOWN)
  , width(0)
  , height(0)
  , chromaFormat(CF444)
  , interlace (false)
  , frameRate (FR0)
  , topFieldFirst (false)
  , bitdepth (0) {}

SequenceHeader::SequenceHeader( Profile profile, 
                                int height,
                                int width,
                                ColourFormat chromaFormat,
                                bool interlace,
                                FrameRate frameRate,
                                bool topFieldFirst,
                                int bitdepth,

                                // Optional Args
                                PixelAspectRatio pixelAspectRatio,
                                int cleanWidth,
                                int cleanHeight,
                                int leftOffset,
                                int topOffset,
                                
                                ColorSpec colorSpec,
                                int colorPrimaries,
                                int colorMatrix,
                                int transferFunction,

                                bool use_v3
                              )
  : major_version(1)
  , minor_version(0)
  , profile (profile)

  , width(width)
  , height(height)
  , chromaFormat(chromaFormat)
  , interlace (interlace)
  , topFieldFirst (topFieldFirst)
  , frameRate (frameRate)  
  , bitdepth (bitdepth)
  , pixelAspectRatio (pixelAspectRatio)
  , cleanWidth (cleanWidth)
  , cleanHeight (cleanHeight)
  , leftOffset (leftOffset)
  , topOffset (topOffset)
  , colorSpec (colorSpec)
  , colorPrimaries (colorPrimaries)
  , colorMatrix (colorMatrix)
  , transferFunction (transferFunction) {
  if (profile == PROFILE_HQ) {
    major_version = 2;
  }
  if (use_v3 ||
      frameRate > MAX_V2_FRAMERATE ||
      bitdepth > 12) {
    major_version = 3;
  }
}

SequenceHeader getDefaultSourceParameters(const int base_video_format_index){
  switch (base_video_format_index) {
    // To Do: refactor bitdepth to be the index
    case  0: return  SequenceHeader(PROFILE_UNKNOWN, 480,640,  CF420, false, FR24000_1001,  false,  8, AR1_1,   640,480,0,0,   CS_CUSTOM   ); break;
    case  1: return  SequenceHeader(PROFILE_UNKNOWN, 120,176,  CF420, false, FR15000_1001,  false,  8, AR10_11, 176,120,0,0,   CS_SDTV_525 ); break;
    case  2: return  SequenceHeader(PROFILE_UNKNOWN, 144,176,  CF420, false, FR25_2,        true,   8, AR12_11, 176,144,0,0,   CS_SDTV_625 ); break;
    case  3: return  SequenceHeader(PROFILE_UNKNOWN, 240,352,  CF420, false, FR15000_1001,  false,  8, AR10_11, 352,240,0,0,   CS_SDTV_525 ); break;
    case  4: return  SequenceHeader(PROFILE_UNKNOWN, 288,352,  CF420, false, FR25_2,        true,   8, AR12_11, 352,288,0,0,   CS_SDTV_625 ); break;
    case  5: return  SequenceHeader(PROFILE_UNKNOWN, 480,704,  CF420, false, FR15000_1001,  false,  8, AR10_11, 704,480,0,0,   CS_SDTV_525 ); break;
    case  6: return  SequenceHeader(PROFILE_UNKNOWN, 576,704,  CF420, false, FR25_2,        true,   8, AR12_11, 704,576,0,0,   CS_SDTV_625 ); break;
    case  7: return  SequenceHeader(PROFILE_UNKNOWN, 480,720,  CF422, true,  FR30000_1001,  false, 10, AR10_11, 704,480,8,0,   CS_SDTV_525 ); break;
    case  8: return  SequenceHeader(PROFILE_UNKNOWN, 576,720,  CF422, true,  FR25,          true,  10, AR12_11, 704,576,8,0,   CS_SDTV_625 ); break;
    case  9: return  SequenceHeader(PROFILE_UNKNOWN, 720,1280, CF422, false, FR60000_1001,  true,  10, AR1_1,   1280,720,0,0,  CS_HDTV     ); break;
    case 10: return  SequenceHeader(PROFILE_UNKNOWN, 720,1280, CF422, false, FR50,          true,  10, AR1_1,   1280,720,0,0,  CS_HDTV     ); break;
    case 11: return  SequenceHeader(PROFILE_UNKNOWN, 1080,1920, CF422, true,  FR30000_1001, true,  10, AR1_1,   1920,1080,0,0, CS_HDTV     ); break;
    case 12: return  SequenceHeader(PROFILE_UNKNOWN, 1080,1920, CF422, true,  FR25,         true,  10, AR1_1,   1920,1080,0,0, CS_HDTV     ); break;
    case 13: return  SequenceHeader(PROFILE_UNKNOWN, 1080,1920, CF422, false, FR60000_1001, true,  10, AR1_1,   1920,1080,0,0, CS_HDTV     ); break;
    case 14: return  SequenceHeader(PROFILE_UNKNOWN, 1080,1920, CF422, false, FR50,         true,  10, AR1_1,   1920,1080,0,0, CS_HDTV     ); break;
    case 15: return  SequenceHeader(PROFILE_UNKNOWN, 1080,2048, CF444, false, FR24,         true,  12, AR1_1,   2048,1080,0,0, CS_D_CINEMA ); break;
    case 16: return  SequenceHeader(PROFILE_UNKNOWN, 2160,4096, CF444, false, FR24,         true,  12, AR1_1,   4096,2160,0,0, CS_D_CINEMA ); break;
    case 17: return  SequenceHeader(PROFILE_UNKNOWN, 2160,3840, CF422, false, FR60000_1001, true,  10, AR1_1,   3840,2160,0,0, CS_UHDTV    ); break;
    case 18: return  SequenceHeader(PROFILE_UNKNOWN, 2160,3840, CF422, false, FR50,         true,  10, AR1_1,   3840,2160,0,0, CS_UHDTV    ); break;
    case 19: return  SequenceHeader(PROFILE_UNKNOWN, 4320,7680, CF422, false, FR60000_1001, true,  10, AR1_1,   7680,4320,0,0, CS_UHDTV    ); break;
    case 20: return  SequenceHeader(PROFILE_UNKNOWN, 4320,7680, CF422, false, FR50,         true,  10, AR1_1,   7680,4320,0,0, CS_UHDTV    ); break;
    case 21: return  SequenceHeader(PROFILE_UNKNOWN, 1080,1920, CF422, false, FR24000_1001, true,  10, AR1_1,   1920,1080,0,0, CS_HDTV     ); break;
    case 22: return  SequenceHeader(PROFILE_UNKNOWN, 486,720,  CF422, true,  FR30000_1001, false,  10, AR10_11, 720,486,0,0,   CS_HDTV     ); break; 
  default:
    throw std::logic_error("DataUnitIO: unknown base video format");
  }
}


bool PictureFormatMatches(const SequenceHeader &fmt,
                          const int w,
                          const int h,
                          const ColourFormat cf,
                          const FrameRate r,
                          const int bd,
                          const bool topFieldFirst) {
  return ((fmt.width == w) &&
          (fmt.height == h) &&
          (fmt.chromaFormat == cf) &&
          (fmt.frameRate == r) &&
          (fmt.bitdepth == bd) &&
          (fmt.topFieldFirst == topFieldFirst));
}

bool PictureFormatMatches(const SequenceHeader &fmt,
                          const int index) {

  const SequenceHeader base = getDefaultSourceParameters(index);

return ((fmt.width == base.width) &&
        (fmt.height == base.height) &&
        (fmt.chromaFormat == base.chromaFormat) &&
        (fmt.frameRate == base.frameRate)&&
        (fmt.bitdepth == base.bitdepth)&&
        (fmt.interlace == base.interlace) &&
        (fmt.topFieldFirst == base.topFieldFirst) &&
        // Optional
        ((fmt.pixelAspectRatio==-1)||(fmt.pixelAspectRatio == base.pixelAspectRatio)) &&
        ((fmt.cleanWidth==-1)||(fmt.cleanWidth == base.cleanWidth))&&
        ((fmt.cleanHeight==-1)||(fmt.cleanHeight == base.cleanHeight))&&
        ((fmt.leftOffset==-1)||(fmt.leftOffset == base.leftOffset))&&
        ((fmt.topOffset==-1)||(fmt.topOffset == base.topOffset))&&
        ((fmt.colorSpec==-1)||(fmt.colorSpec == base.colorSpec)));
}

int CheckMatch(const SequenceHeader &fmt,
                          const int index) {

  const SequenceHeader base = getDefaultSourceParameters(index);

  int non_matching_fields =
                          (fmt.width != base.width) +
                          (fmt.height != base.height) +
                          (fmt.chromaFormat != base.chromaFormat) +
                          (fmt.frameRate != base.frameRate) +
                          (fmt.bitdepth != base.bitdepth)+
                          (fmt.interlace != base.interlace) +
                          // Optional
                          ((fmt.pixelAspectRatio!=-1)&&(fmt.pixelAspectRatio != base.pixelAspectRatio)) +
                          ((fmt.cleanWidth!=-1)&&(fmt.cleanWidth != base.cleanWidth))+
                          ((fmt.cleanHeight!=-1)&&(fmt.cleanHeight != base.cleanHeight))+
                          ((fmt.leftOffset!=-1)&&(fmt.leftOffset != base.leftOffset))+
                          ((fmt.topOffset!=-1)&&(fmt.topOffset != base.topOffset))+
                          ((fmt.colorSpec!=-1)&&(fmt.colorSpec != base.colorSpec));


  bool valid = (fmt.topFieldFirst == base.topFieldFirst);

  return valid ? non_matching_fields : -1;
}

video_format::video_format()
    : major_version(0)
    , minor_version(0)
    , profile (0)
    , level (0)
    , base_video_format (0)
    , custom_dimensions_flag (false)
    , frame_width (0)
    , frame_height (0)
    , custom_color_diff_format_flag (false)
    , color_diff_format(0)
    , custom_scan_format_flag (false)
    , source_sampling (0)
    , custom_frame_rate_flag (false)
    , frame_rate (FR0)
    , custom_pixel_aspect_ratio_flag (false) 
    , pixel_aspect_ratio (0)            
    , custom_clean_area_flag (false)
    , clean_width (0)
    , clean_height (0)
    , left_offset (0)
    , top_offset (0)
    , custom_signal_range_flag (false)
    , bitdepth (0)
    , custom_color_spec_flag (false)
    , color_spec (0)
    , custom_color_primaries_flag (0)
    , color_primaries (0)
    , custom_color_matrix_flag (false)
    , color_matrix (0)
    , custom_transfer_function_flag (false)
    , transfer_function (0)
    , top_field_first (false) {}

video_format::video_format(const SequenceHeader &fmt)
    : major_version(0)
    , minor_version(0)
    , profile (0)
    , level (0)
    , base_video_format (0)
    , custom_dimensions_flag (false)
    , frame_width (0)
    , frame_height (0)
    , custom_color_diff_format_flag(false)
    , color_diff_format(0)
    , custom_scan_format_flag (false)
    , source_sampling (0)
    , custom_frame_rate_flag (false)
    , frame_rate (FR0)
    , custom_pixel_aspect_ratio_flag (false) 
    , pixel_aspect_ratio (0)            
    , custom_clean_area_flag (false)
    , clean_width (0)
    , clean_height (0)
    , left_offset (0)
    , top_offset (0)
    , custom_signal_range_flag (false)
    , bitdepth (0)
    , custom_color_spec_flag (false)
    , color_spec (0)
    , custom_color_primaries_flag (0)
    , color_primaries (0)
    , custom_color_matrix_flag (false)
    , color_matrix (0)
    , custom_transfer_function_flag (false)
    , transfer_function (0)
    , top_field_first (false)

     {

  major_version = fmt.major_version;
  minor_version = fmt.minor_version;
  switch (fmt.profile) {
  case PROFILE_LD:
    profile = 0;
    break;
  case PROFILE_HQ:
    profile = 3;
    break;
  default:
    profile = 0;
    break;
  }

  if (fmt.interlace) {
    // Level 2
    if      (PictureFormatMatches(fmt, 7))      { base_video_format =  7; level = 2; }
    else if (PictureFormatMatches(fmt, 8)) { base_video_format =  8; level = 2; }
    else if (PictureFormatMatches(fmt, 22)) { base_video_format = 22; level = 2; }
    else if (fmt.chromaFormat == CF422 &&
             fmt.width == 720 &&
             fmt.height >= 480 &&
             fmt.height <= 486 &&
             fmt.frameRate == FR30000_1001 &&
             fmt.bitdepth == 10) {
      base_video_format = 7;
      level = 2;
      custom_dimensions_flag = true;
      frame_width = fmt.width;
      frame_height = fmt.height;
    }

    // Level 3
    else if (PictureFormatMatches(fmt, 11)) { base_video_format = 11; level = 3; }
    else if (PictureFormatMatches(fmt, 12)) { base_video_format = 12; level = 3; }
  } else {
    // Level 1
    if      (PictureFormatMatches(fmt, 1)) { base_video_format = 1; level = 1; }
    else if (PictureFormatMatches(fmt, 2)) { base_video_format = 2; level = 1; }
    else if (PictureFormatMatches(fmt, 3)) { base_video_format = 3; level = 1; }
    else if (PictureFormatMatches(fmt, 4)) { base_video_format = 4; level = 1; }
    else if (PictureFormatMatches(fmt, 5)) { base_video_format = 5; level = 1; }
    else if (PictureFormatMatches(fmt, 6)) { base_video_format = 6; level = 1; }

    // Level 2
    else if (PictureFormatMatches(fmt, 720, 480, CF422, FR30000_1001, 10, false)) { base_video_format =  7; level = 2; custom_scan_format_flag = true; source_sampling = 0; }
    else if (PictureFormatMatches(fmt, 720, 576, CF422, FR25,         10, true)) { base_video_format =  8; level = 2; custom_scan_format_flag = true; source_sampling = 0; }
    else if (PictureFormatMatches(fmt, 720, 486, CF422, FR30000_1001, 10, false)) { base_video_format = 22; level = 2; custom_scan_format_flag = true; source_sampling = 0; }

    // Level 3
    else if (PictureFormatMatches(fmt, 9)) { base_video_format =  9; level = 3; }
    else if (PictureFormatMatches(fmt, 10)) { base_video_format = 10; level = 3; }
    else if (PictureFormatMatches(fmt, 1920, 1080, CF422, FR30000_1001, 10, true)) { base_video_format = 11; level = 3; custom_scan_format_flag = true; source_sampling = 0; }
    else if (PictureFormatMatches(fmt, 1920, 1080, CF422, FR25,         10, true)) { base_video_format = 12; level = 3; custom_scan_format_flag = true; source_sampling = 0; }
    else if (PictureFormatMatches(fmt, 13)) { base_video_format = 13; level = 3; }
    else if (PictureFormatMatches(fmt, 14)) { base_video_format = 14; level = 3; }
    else if (PictureFormatMatches(fmt, 21)) { base_video_format = 21; level = 3; }

    // Level 4
    else if (PictureFormatMatches(fmt, 15)) { base_video_format = 15; level = 4; }
    else if (PictureFormatMatches(fmt, 2048, 1080, CF444, FR48, 12, true)) { base_video_format = 15; level = 4; custom_frame_rate_flag = true; frame_rate = FR48; }

    // Level 5
    else if (PictureFormatMatches(fmt, 16)) { base_video_format = 16; level = 5; }

    // Level 6
    else if (PictureFormatMatches(fmt, 17)) { base_video_format = 17; level = 6; }
    else if (PictureFormatMatches(fmt, 18)) { base_video_format = 18; level = 6; }

    // Level 7
    else if (PictureFormatMatches(fmt, 19)) { base_video_format = 19; level = 7; }
    else if (PictureFormatMatches(fmt, 20)) { base_video_format = 20; level = 7; }
  }
  
  if (base_video_format == 0) {
    level = 0 ;
    int non_matching_fields;
    int non_matching_fields_prev = 999;
    for (int base_format = 1; base_format<=22; base_format++){
      non_matching_fields = CheckMatch(fmt, base_format);
      if (non_matching_fields == -1) continue;
      if (non_matching_fields < non_matching_fields_prev){
        base_video_format = base_format;
        non_matching_fields_prev = non_matching_fields;
      }
    }

    SequenceHeader base = getDefaultSourceParameters(base_video_format);

    if (fmt.interlace != base.interlace) {
      custom_scan_format_flag = true;
      source_sampling = fmt.interlace;
    }
    if (fmt.width != base.width || fmt.height != base.height) {
      custom_dimensions_flag = true;
      frame_width  = fmt.width;
      frame_height = fmt.height;
    }
    if (fmt.chromaFormat != base.chromaFormat) {
      custom_color_diff_format_flag = true;
      color_diff_format = fmt.chromaFormat;
    }
    if (fmt.frameRate != base.frameRate) {
      custom_frame_rate_flag = true;
      frame_rate = fmt.frameRate;
    }
    if (fmt.bitdepth != base.bitdepth) {
      custom_signal_range_flag = true;
      switch (fmt.bitdepth) {
      case  8: bitdepth = 1; break;
      case 10: bitdepth = 3; break;
      case 12: bitdepth = 4; break;
      case 16: bitdepth = 7; break;
      default:
        throw std::logic_error("DataUnitIO: invalid bit depth");
      }
    }
    if ((fmt.pixelAspectRatio != AR_UNSET)&&(fmt.pixelAspectRatio != base.pixelAspectRatio)) {
      custom_pixel_aspect_ratio_flag = true;
      pixel_aspect_ratio = fmt.pixelAspectRatio;
    }
    if (
       (fmt.cleanHeight != -1 ||
        fmt.cleanWidth  != -1 ||
        fmt.leftOffset  != -1 ||
        fmt.topOffset   != -1 )
        &&
       (fmt.cleanHeight != base.cleanHeight ||
        fmt.cleanWidth  != base.cleanWidth  ||
        fmt.leftOffset  != base.leftOffset  ||
        fmt.topOffset   != base.topOffset)) {
      custom_clean_area_flag = true;
      clean_height = fmt.cleanHeight;
      clean_width = fmt.cleanWidth;
      left_offset = fmt.leftOffset;
      top_offset = fmt.topOffset;
    }
    if ((fmt.colorSpec != CS_UNSET)&&(fmt.colorSpec != base.colorSpec)) {
      custom_color_spec_flag = true;      
      color_spec = fmt.colorSpec;
    }

    if (fmt.colorSpec == CS_CUSTOM){
      if (fmt.colorPrimaries != base.colorPrimaries) {
        custom_color_primaries_flag = true;
        color_primaries = fmt.colorPrimaries;
      }
      if (fmt.colorMatrix != base.colorMatrix) {
        custom_color_matrix_flag = true;
        color_matrix = fmt.colorMatrix;  
      }
      if (fmt.transferFunction != base.transferFunction) {
        custom_transfer_function_flag = true;
        transfer_function = fmt.transferFunction;  
      }
    }
    
  }
}

std::ostream& operator << (std::ostream& ss, const video_format& fmt) {

  major_version_number(ss) = fmt.major_version;
  ss << vlc::unbounded;
  
  ss << UnsignedVLC(fmt.major_version)
     << UnsignedVLC(fmt.minor_version)
     << UnsignedVLC(fmt.profile)
     << UnsignedVLC(fmt.level);

  ss << UnsignedVLC(fmt.base_video_format);

  ss << Boolean(fmt.custom_dimensions_flag);
  if (fmt.custom_dimensions_flag) {
    ss << UnsignedVLC(fmt.frame_width)
       << UnsignedVLC(fmt.frame_height);
  }

  ss << Boolean(fmt.custom_color_diff_format_flag);
  if (fmt.custom_color_diff_format_flag) {
      ss << UnsignedVLC((int)fmt.color_diff_format);
  }

  ss << Boolean(fmt.custom_scan_format_flag);
  if (fmt.custom_scan_format_flag) {
    ss << UnsignedVLC(fmt.source_sampling);
  }

  // To do: handle custom frame rate (index 0)
  ss << Boolean(fmt.custom_frame_rate_flag);
  if (fmt.custom_frame_rate_flag) {
      ss << UnsignedVLC((int) fmt.frame_rate);
  }

  // To do: handle custom pixel aspect ratio (index 0)
  ss << Boolean(fmt.custom_pixel_aspect_ratio_flag);
  if (fmt.custom_pixel_aspect_ratio_flag){
      ss << UnsignedVLC((int) fmt.pixel_aspect_ratio);
  }

  // To do: handle restrictions on clean area (elsewhere)
  ss << Boolean(fmt.custom_clean_area_flag);
  if (fmt.custom_clean_area_flag){
    ss << UnsignedVLC(fmt.clean_width);
    ss << UnsignedVLC(fmt.clean_height);
    ss << UnsignedVLC(fmt.left_offset);
    ss << UnsignedVLC(fmt.top_offset);
  }

  // To do: handle custom signal range (index 0)
  // To do: replace bitdepth with signal range num
  ss << Boolean(fmt.custom_signal_range_flag);
  if (fmt.custom_signal_range_flag) {
    // bitdepth here has been replaced by the index 
    // (in copy_video_fmt_to_hdr)
    ss << UnsignedVLC(fmt.bitdepth);  
  }

  ss << Boolean(fmt.custom_color_spec_flag);
  if (fmt.custom_color_spec_flag) {
    ss << UnsignedVLC((int)fmt.color_spec);
    if (fmt.color_spec == CS_CUSTOM){
      ss << Boolean(fmt.custom_color_primaries_flag);
      if (fmt.custom_color_primaries_flag) {
        ss << UnsignedVLC(fmt.color_primaries);
      }
      ss << Boolean(fmt.custom_color_matrix_flag);
      if (fmt.custom_color_matrix_flag) {
        ss << UnsignedVLC(fmt.color_matrix);
      }
      ss << Boolean(fmt.custom_transfer_function_flag);
      if (fmt.custom_transfer_function_flag) {
        ss << UnsignedVLC(fmt.transfer_function);
      }
    }
  }
  // Use the source sampling to determine the picture coding mode
  // This means that progressive video is always encoded as frames
  // and interlaced video is always encoded as fields
  // To do: add a separate parameter for the picture coding mode (11.5)
  ss << UnsignedVLC(fmt.source_sampling);
  ss << vlc::align;

  return ss;
}

std::istream& operator >> (std::istream& stream, video_format& fmt) {
  stream >> vlc::unbounded;

  UnsignedVLC major_version, minor_version, profile, level;
  stream >> major_version >> minor_version >> profile >> level;
  fmt.major_version = major_version;
  major_version_number(stream) = fmt.major_version;
  fmt.minor_version = minor_version;
  fmt.profile = profile;
  fmt.level = level;

  UnsignedVLC base_video_format;
  stream >> base_video_format;
  fmt.base_video_format = base_video_format;

  Boolean custom_dimensions_flag;
  stream >> custom_dimensions_flag;
  fmt.custom_dimensions_flag = custom_dimensions_flag;
  if (custom_dimensions_flag) {
    UnsignedVLC frame_width, frame_height;
    stream >> frame_width >> frame_height;
    fmt.frame_width  = frame_width;
    fmt.frame_height = frame_height;
  }

  Boolean custom_color_diff_format_flag;
  stream >> custom_color_diff_format_flag;
  fmt.custom_color_diff_format_flag = custom_color_diff_format_flag;
  if (custom_color_diff_format_flag) {
    UnsignedVLC color_diff_format;
    stream >> color_diff_format;
    try {
      fmt.color_diff_format = (ColourFormat)(int)color_diff_format;
    } catch(const std::exception& e) {
      std::stringstream ss;
      ss << "DataUnitIO: Invalid Frame Rate on Input: " << (int)color_diff_format;
      ss << e.what();
      throw std::logic_error(ss.str());
    }
  }

  Boolean custom_scan_format_flag;
  stream >> custom_scan_format_flag;
  fmt.custom_scan_format_flag = custom_scan_format_flag;
  if (custom_scan_format_flag) {
    UnsignedVLC source_sampling;
    stream >> source_sampling;
    fmt.source_sampling  = source_sampling;
  }

  // To do: handle custom frame rate (index 0)
  Boolean custom_frame_rate_flag;
  stream >> custom_frame_rate_flag;
  fmt.custom_frame_rate_flag = custom_frame_rate_flag;
  if (custom_frame_rate_flag) {
    UnsignedVLC index;
    stream >> index;

    try {
      fmt.frame_rate = (FrameRate)(int)index;
    } catch(const std::exception& e) {
      std::stringstream ss;
      ss << "DataUnitIO: Invalid Frame Rate on Input: " << (int)index;
      ss << e.what();
      throw std::logic_error(ss.str());
    }
  }

  // To do: handle custom pixel aspect ratio (index 0)
  Boolean custom_pixel_aspect_ratio_flag;
  stream >> custom_pixel_aspect_ratio_flag;
  if (custom_pixel_aspect_ratio_flag) {
    UnsignedVLC index;
    stream >> index;
    try {
      fmt.pixel_aspect_ratio = (PixelAspectRatio)(int)index;
    } catch(const std::exception& e) {
      std::stringstream ss;
      ss << "DataUnitIO: Invalid Pixel Aspect Ratio on Input: " << (int)index;
      ss << e.what();
      throw std::logic_error(ss.str());
    }
  }

  // To do: handle restrictions on clean area (elsewhere?)
  Boolean custom_clean_area_flag;
  stream >> custom_clean_area_flag;
  if (custom_clean_area_flag) {
    UnsignedVLC clean_width, clean_height, left_offset, top_offset;
    stream >> clean_width >> clean_height >> left_offset >> top_offset;
    
    fmt.clean_height = clean_height;
    fmt.clean_width = clean_width;
    fmt.left_offset = left_offset;
    fmt.top_offset = top_offset;
  }

  // To do: handle custom signal range (index 0)
  Boolean custom_signal_range_flag;
  stream >> custom_signal_range_flag;
  fmt.custom_signal_range_flag = custom_signal_range_flag;
  if (custom_signal_range_flag) {
    UnsignedVLC bitdepth;
    stream >> bitdepth;
    fmt.bitdepth = bitdepth;
  }

  Boolean custom_color_spec_flag;
  stream >> custom_color_spec_flag;
  fmt.custom_color_spec_flag = custom_color_spec_flag;
  if (custom_color_spec_flag) {
    UnsignedVLC custom_color_spec_index;
    stream >> custom_color_spec_index;
    fmt.color_spec = custom_color_spec_index;
    if ((ColorSpec)(int)custom_color_spec_index == CS_CUSTOM) {
      Boolean custom_color_primaries_flag, custom_color_matrix_flag, custom_transfer_function_flag;

      stream >> custom_color_primaries_flag;
      fmt.custom_color_primaries_flag = custom_color_primaries_flag;
      if (custom_color_primaries_flag) {
        UnsignedVLC color_primaries;
        stream >> color_primaries;
        fmt.color_primaries = color_primaries;
      }

      stream >> custom_color_matrix_flag;
      fmt.custom_color_matrix_flag = custom_color_matrix_flag;
      if (custom_color_matrix_flag) {
        UnsignedVLC color_matrix;
        stream >> color_matrix;
        fmt.color_matrix = color_matrix;
      }

      stream >> custom_transfer_function_flag;
      fmt.custom_transfer_function_flag = custom_transfer_function_flag;
      if (custom_transfer_function_flag) {
        UnsignedVLC transfer_function;
        stream >> transfer_function;
        fmt.transfer_function = transfer_function;
      }
    }
  }

  // Determine the picture coding mode from the source_sampling
  // This means that progressive video is always encoded as frames
  // and interlaced video is always encoded as fields
  // To do: add a separate parameter for the picture coding mode (11.5)
  UnsignedVLC source_sampling;
  stream >> source_sampling;
  fmt.source_sampling = source_sampling;

  stream >> vlc::align;

  return stream;
}

std::ostream& operator << (std::ostream& stream, const SequenceHeader& s) {
  video_format fmt(s);

  if (fragment_length(stream) > 0)
    if (s.major_version < 3)
      fmt.major_version = 3;

  std::stringstream ss;
  ss.copyfmt(stream);
  
  ss << fmt;

  stream.copyfmt(ss);

  stream << ParseInfoIO(SEQUENCE_HEADER, ss.str().size()) << ss.str();

  return stream;
}

std::istream& dataunitio::lowDelay(std::istream& stream) {
  sliceio::sliceIOMode(stream) = sliceio::LD;
  return stream;
}

std::istream& dataunitio::synchronise(std::istream &stream) {
  while (stream) {
    Bytes b(1);
    stream >> b;
    if (0x42 != (unsigned char)b)
      continue;

    stream >> b;
    if (0x42 != (unsigned char)b)
      continue;

    stream >> b;
    if (0x43 != (unsigned char)b)
      continue;

    stream >> b;
    if (0x44 != (unsigned char)b)
      continue;

    break;
  }
  
  return stream;
}

std::istream& operator >> (std::istream& stream, DataUnit &d) {

  // Read and check the 4 parse info prefix bytes
  Bytes prefix1(1),prefix2(1),prefix3(1),prefix4(1);
  stream >> prefix1 >> prefix2 >> prefix3 >> prefix4 ;
  if (0x42 != (unsigned char)prefix1||
      0x42 != (unsigned char)prefix2||
      0x43 != (unsigned char)prefix3||
      0x44 != (unsigned char)prefix4
      ){
        throw std::logic_error("Read bytes do not match expected parse_info_header.");
      }

  Bytes type(1);
  stream >> type;

  switch ((unsigned char)type) {
  case 0x00: d.type = SEQUENCE_HEADER; break;
  case 0x10: d.type = END_OF_SEQUENCE; break;
  case 0x20: d.type = AUXILIARY_DATA;  break;
  case 0x30: d.type = PADDING_DATA;    break;
  case 0xC8: d.type = LD_PICTURE;      break;
  case 0xE8: d.type = HQ_PICTURE;      break;
  case 0xCC: d.type = LD_FRAGMENT;     break;
  case 0xEC: d.type = HQ_FRAGMENT;     break;
  default:
    d.type = UNKNOWN_DATA_UNIT;
    throw std::logic_error("Stream Error: Unknown data unit type.");
  }

  stream >> d.next_parse_offset >> d.prev_parse_offset;

  return stream;
}

std::istream& operator >> (std::istream& stream, Fragment &d) {
  Bytes picnum(4);
  Bytes fraglen(2);
  Bytes slice_count(2);
  stream >> fraglen >> slice_count;

  d.mNSlices       = slice_count;

  if (d.mNSlices != 0) {
    Bytes slice_offset_x(2);
    Bytes slice_offset_y(2);
    stream >> slice_offset_x >> slice_offset_y;
    d.mSliceOffsetX = (int)slice_offset_x;
    d.mSliceOffsetY = (int)slice_offset_y;
  }

  return stream;
}

std::ostream& operator << (std::ostream& stream, const DataUnitType& t) {
  switch(t) {
  case SEQUENCE_HEADER: stream << "Sequence Header"; break;
  case END_OF_SEQUENCE: stream << "End of Sequence"; break;
  case AUXILIARY_DATA:  stream << "Auxiliary Data";  break;
  case PADDING_DATA:    stream << "Padding Data";    break;
  case LD_PICTURE:      stream << "LD Picture";      break;
  case HQ_PICTURE:      stream << "HQ Picture";      break;
  case LD_FRAGMENT:     stream << "LD Fragment";     break;
  case HQ_FRAGMENT:     stream << "HQ Fragment";     break;
  default:
    stream << "Unknown Data Unit";
  }
  
  return stream;
}

std::ostream& operator << (std::ostream& stream, const FrameRate& r) {
  switch(r) {
  case FR24000_1001: stream << "24/1.001 fps"; break;
  case FR24: stream << "24 fps"; break;
  case FR25: stream << "25 fps"; break;
  case FR30000_1001: stream << "30/1.001 fps"; break;
  case FR30: stream << "30 fps"; break;
  case FR50: stream << "50 fps"; break;
  case FR60000_1001: stream << "60/1.001 fps"; break;
  case FR60: stream << "60 fps"; break;
  case FR15000_1001: stream << "50/1.001 fps"; break;
  case FR25_2: stream << "25/2 fps"; break;
  case FR48: stream << "48 fps"; break;
  default:
    stream << "unknown";
    break;
  }

  return stream;
}

void copy_video_fmt_to_hdr (SequenceHeader *hdr, video_format &fmt) {

  SequenceHeader base = getDefaultSourceParameters(fmt.base_video_format);

  hdr->profile       = base.profile;
  hdr->width         = base.width;
  hdr->height        = base.height;
  hdr->chromaFormat  = base.chromaFormat;
  hdr->interlace     = base.interlace;
  hdr->frameRate     = base.frameRate;
  hdr->topFieldFirst = base.topFieldFirst;
  hdr->bitdepth      = base.bitdepth;
  // Optional Args
  hdr->pixelAspectRatio = base.pixelAspectRatio;
  hdr->cleanWidth       = base.cleanWidth;
  hdr->cleanHeight      = base.cleanHeight;
  hdr->leftOffset       = base.leftOffset;
  hdr->topOffset        = base.topOffset;
  
  hdr->colorSpec        = base.colorSpec;
  hdr->colorPrimaries   = base.colorPrimaries;
  hdr->colorMatrix      = base.colorMatrix;
  hdr->transferFunction = base.transferFunction;

  hdr->major_version = fmt.major_version;
  hdr->minor_version = fmt.minor_version;

  if (fmt.profile == 0)
    hdr->profile = PROFILE_LD;
  else if (fmt.profile == 3)
    hdr->profile = PROFILE_HQ;

  if (fmt.custom_dimensions_flag) {
    hdr->width = fmt.frame_width;
    hdr->height = fmt.frame_height;
  }
  if (fmt.custom_color_diff_format_flag) {
    hdr->chromaFormat = (ColourFormat)fmt.color_diff_format;
  }
  if (fmt.custom_scan_format_flag) {
    if (fmt.source_sampling == 0)
      hdr->interlace = false;
    else
      hdr->interlace = true;
  }
  if (fmt.custom_frame_rate_flag) {
    hdr->frameRate = fmt.frame_rate;

    if (fmt.frame_rate > MAX_V2_FRAMERATE)
      if (hdr->major_version < 3)
        hdr->major_version = 3;
  }
  if (fmt.custom_pixel_aspect_ratio_flag){
    hdr->pixelAspectRatio = (PixelAspectRatio)fmt.pixel_aspect_ratio;
  }
  if (fmt.custom_clean_area_flag){
    hdr->cleanWidth   = fmt.clean_width;
    hdr->cleanHeight  = fmt.clean_height;
    hdr->leftOffset   = fmt.left_offset;
    hdr->topOffset    = fmt.top_offset;
  }
  if (fmt.custom_signal_range_flag) {
    switch (fmt.bitdepth) {
    case 1: hdr->bitdepth =  8; break;
    case 2: hdr->bitdepth =  8; break;
    case 3: hdr->bitdepth = 10; break;
    case 4: hdr->bitdepth = 12; break;
    case 5: hdr->bitdepth = 10; break;
    case 6: hdr->bitdepth = 12; break;
    case 7: hdr->bitdepth = 16; break;
    case 8: hdr->bitdepth = 16; break;
    }

    if (fmt.bitdepth > 4) {
      if (hdr->major_version < 3)
        hdr->major_version = 3;
    }
  }
  if (fmt.custom_color_spec_flag){
    hdr->colorSpec   = (ColorSpec)fmt.color_spec;
    if ((ColorSpec)fmt.color_spec == CS_CUSTOM){
      if (fmt.custom_color_primaries_flag){
        hdr->colorPrimaries = fmt.color_primaries;
      }
      if (fmt.custom_color_matrix_flag){
        hdr->colorMatrix = fmt.color_matrix;
      }
      if (fmt.custom_transfer_function_flag){
        hdr->transferFunction = fmt.transfer_function;
      }
    }
  }

  //  return hdr;
}

std::istream& operator >> (std::istream& stream, SequenceHeader &hdr) {
  video_format fmt;
  stream >> fmt;
  copy_video_fmt_to_hdr(&hdr, fmt);
  major_version_number(stream) = hdr.major_version;
  return stream;
}

PicturePreamble::PicturePreamble()
  : wavelet_kernel(NullKernel)
  , depth (0)
  , slices_x (0)
  , slices_y (0)
  , slice_prefix (0)
  , slice_size_scalar (0)
  , slice_bytes() {
}

std::istream& operator >> (std::istream& stream, PictureHeader &h) {
  Bytes picnum(4);
  stream >> picnum;
  h.picture_number = (int)picnum;
  return stream;
}

std::istream& operator >> (std::istream& stream, PicturePreamble &hdr) {
  UnsignedVLC wavelet_index, depth;
  stream >> wavelet_index >> depth;
  switch(wavelet_index) {
  case 0: hdr.wavelet_kernel = DD97;     break;
  case 1: hdr.wavelet_kernel = LeGall;   break;
  case 2: hdr.wavelet_kernel = DD137;    break;
  case 3: hdr.wavelet_kernel = Haar0;    break;
  case 4: hdr.wavelet_kernel = Haar1;    break;
  case 5: hdr.wavelet_kernel = Fidelity; break;
  case 6: hdr.wavelet_kernel = Daub97;   break;
  }
  hdr.depth = depth;
  hdr.wavelet_kernel_ho = hdr.wavelet_kernel;
  hdr.depth_ho = 0;

  if (major_version_number(stream) >= 3) {
    Boolean asym_transform_index_flag, asym_transform_flag;
    UnsignedVLC wavelet_index_ho;
    UnsignedVLC dwt_depth_ho;
    stream >> asym_transform_index_flag;
    if (asym_transform_index_flag) {
      stream >> wavelet_index_ho;
      switch(wavelet_index_ho) {
      case 0: hdr.wavelet_kernel_ho = DD97;     break;
      case 1: hdr.wavelet_kernel_ho = LeGall;   break;
      case 2: hdr.wavelet_kernel_ho = DD137;    break;
      case 3: hdr.wavelet_kernel_ho = Haar0;    break;
      case 4: hdr.wavelet_kernel_ho = Haar1;    break;
      case 5: hdr.wavelet_kernel_ho = Fidelity; break;
      case 6: hdr.wavelet_kernel_ho = Daub97;   break;
      }
    }
    stream >> asym_transform_flag;
    if (asym_transform_flag) {
      stream >> dwt_depth_ho;
      hdr.depth_ho = dwt_depth_ho;
    }
  }

  if (sliceio::sliceIOMode(stream) == sliceio::HQVBR || sliceio::sliceIOMode(stream) == sliceio::HQCBR) {
    UnsignedVLC slices_x, slices_y, slice_prefix, slice_size_scalar;
    stream >> slices_x >> slices_y >> slice_prefix >> slice_size_scalar;

    hdr.slices_x          = slices_x;
    hdr.slices_y          = slices_y;
    hdr.slice_prefix      = slice_prefix;
    hdr.slice_size_scalar = slice_size_scalar;
    hdr.slice_bytes       = utils::rationalise(0, 1);
  } else if (sliceio::sliceIOMode(stream) == sliceio::LD) {
    UnsignedVLC slices_x, slices_y, slice_bytes_numerator, slice_bytes_denominator;
    stream >> slices_x >> slices_y >> slice_bytes_numerator >> slice_bytes_denominator;

    hdr.slices_x          = slices_x;
    hdr.slices_y          = slices_y;
    hdr.slice_prefix      = 0;
    hdr.slice_size_scalar = 0;
    hdr.slice_bytes       = utils::rationalise(slice_bytes_numerator, slice_bytes_denominator);
  } else {
    throw std::logic_error("DataUnitIO: Not in HQ or LD Format");
  }
  Boolean custom_quant_matrix;
  stream >> custom_quant_matrix;

  if (custom_quant_matrix) {
    throw std::logic_error("DataUnitIO: Custom Quantisation Matrix flag not supported");
  }

  stream >> vlc::align;

  return stream;
}

void dataunitio::fragmentedPictures::operator() (std::ios_base &stream) const {
  fragment_length(stream) = (long)mFL;
  if (major_version_number(stream) < 3)
    major_version_number(stream) = 3;
}

std::ostream& operator << (std::ostream& stream, dataunitio::fragmentedPictures arg) {
  arg(stream);
  return stream;
}
