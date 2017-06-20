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
  , strm () {}

std::istream &DataUnit::stream() { return strm; }

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

SequenceHeader::SequenceHeader(Profile p, int h, int w, ColourFormat c, bool i, FrameRate f, bool tff, int bd, bool use_v3)
  : major_version(1)
  , minor_version(0)
  , profile (p)
  , width(w)
  , height(h)
  , chromaFormat(c)
  , interlace (i)
  , frameRate (f)
  , topFieldFirst (tff)
  , bitdepth (bd) {
  if (p == PROFILE_HQ) {
    major_version = 2;
  }
  if (use_v3 ||
      f > MAX_V2_FRAMERATE ||
      bd > 12) {
    major_version = 3;
  }
}


struct video_format {
  video_format()
    : major_version(0)
    , minor_version(0)
    , profile (0)
    , level (0)
    , base_video_format (0)
    , custom_dimensions_flag (false)
    , frame_width (0)
    , frame_height (0)
    , custom_scan_format_flag (false)
    , source_sampling (0)
    , custom_signal_range_flag (false)
    , bitdepth (0)
    , custom_frame_rate_flag (false)
    , frame_rate (FR0) {}

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
  bool custom_signal_range_flag;
  int bitdepth;
  bool custom_frame_rate_flag;
  FrameRate frame_rate;
  int picture_coding_mode;
};

bool PictureFormatMatches(const SequenceHeader &fmt,
                          const int w,
                          const int h,
                          const ColourFormat cf,
                          const FrameRate r,
                          const int bd) {
  return ((fmt.width == w) &&
          (fmt.height == h) &&
          (fmt.chromaFormat == cf) &&
          (fmt.frameRate == r) &&
          (fmt.bitdepth == bd));
}


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
    , custom_signal_range_flag (false)
    , bitdepth (0)
    , custom_frame_rate_flag (false)
    , frame_rate (FR0) {

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
    if (PictureFormatMatches(fmt, 720, 480, CF422, FR30000_1001, 10))      { base_video_format =  7; level = 2; }
    else if (PictureFormatMatches(fmt, 720, 576, CF422, FR25,         10)) { base_video_format =  8; level = 2; }
    else if (PictureFormatMatches(fmt, 720, 486, CF422, FR30000_1001, 10)) { base_video_format = 22; level = 2; }
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
    else if (PictureFormatMatches(fmt, 1920, 1080, CF422, FR30000_1001, 10)) { base_video_format = 11; level = 3; }
    else if (PictureFormatMatches(fmt, 1920, 1080, CF422, FR25,         10)) { base_video_format = 12; level = 3; }
  } else {
    // Level 1
    if (PictureFormatMatches(fmt, 176, 120, CF420, FR15000_1001,      8)) { base_video_format = 1; level = 1; }
    else if (PictureFormatMatches(fmt, 176, 144, CF420, FR25_2,       8)) { base_video_format = 2; level = 1; }
    else if (PictureFormatMatches(fmt, 352, 240, CF420, FR15000_1001, 8)) { base_video_format = 3; level = 1; }
    else if (PictureFormatMatches(fmt, 352, 288, CF420, FR25_2,       8)) { base_video_format = 4; level = 1; }
    else if (PictureFormatMatches(fmt, 704, 480, CF420, FR15000_1001, 8)) { base_video_format = 5; level = 1; }
    else if (PictureFormatMatches(fmt, 704, 576, CF420, FR25_2,       8)) { base_video_format = 6; level = 1; }

    // Level 2
    else if (PictureFormatMatches(fmt, 720, 480, CF422, FR30000_1001, 10)) { base_video_format =  7; level = 2; custom_scan_format_flag = true; source_sampling = 0; }
    else if (PictureFormatMatches(fmt, 720, 576, CF422, FR25,         10)) { base_video_format =  8; level = 2; custom_scan_format_flag = true; source_sampling = 0; }
    else if (PictureFormatMatches(fmt, 720, 486, CF422, FR30000_1001, 10)) { base_video_format = 22; level = 2; custom_scan_format_flag = true; source_sampling = 0; }

    // Level 3
    else if (PictureFormatMatches(fmt, 1280, 720, CF422, FR60000_1001,  10)) { base_video_format =  9; level = 3; }
    else if (PictureFormatMatches(fmt, 1280, 720, CF422, FR50,          10)) { base_video_format = 10; level = 3; }
    else if (PictureFormatMatches(fmt, 1920, 1080, CF422, FR30000_1001, 10)) { base_video_format = 11; level = 3; custom_scan_format_flag = true; source_sampling = 0; }
    else if (PictureFormatMatches(fmt, 1920, 1080, CF422, FR25,         10)) { base_video_format = 12; level = 3; custom_scan_format_flag = true; source_sampling = 0; }
    else if (PictureFormatMatches(fmt, 1920, 1080, CF422, FR60000_1001, 10)) { base_video_format = 13; level = 3; }
    else if (PictureFormatMatches(fmt, 1920, 1080, CF422, FR50,         10)) { base_video_format = 14; level = 3; }
    else if (PictureFormatMatches(fmt, 1920, 1080, CF422, FR24000_1001, 10)) { base_video_format = 21; level = 3; }

    // Level 4
    else if (PictureFormatMatches(fmt, 2048, 1080, CF444, FR24, 12)) { base_video_format = 15; level = 4; }
    else if (PictureFormatMatches(fmt, 2048, 1080, CF444, FR48, 12)) { base_video_format = 15; level = 4; custom_frame_rate_flag = true; frame_rate = FR48; }

    // Level 5
    else if (PictureFormatMatches(fmt, 4096, 2160, CF444, FR24, 12)) { base_video_format = 16; level = 5; }

    // Level 6
    else if (PictureFormatMatches(fmt, 3840, 2160, CF422, FR60000_1001, 10)) { base_video_format = 17; level = 6; }
    else if (PictureFormatMatches(fmt, 3840, 2160, CF422, FR50,         10))         { base_video_format = 18; level = 6; }

    // Level 7
    else if (PictureFormatMatches(fmt, 7680, 4320, CF422, FR60000_1001, 10)) { base_video_format = 19; level = 7; }
    else if (PictureFormatMatches(fmt, 7680, 4320, CF422, FR50,         10)) { base_video_format = 20; level = 7; }
  }

  if (base_video_format == 0) {
    if (fmt.interlace) {
      custom_scan_format_flag = true;
      source_sampling = 1;
    }
    if (fmt.width != 640 || fmt.height != 480) {
      custom_dimensions_flag = true;
      frame_width  = fmt.width;
      frame_height = fmt.height;
    }
    if (fmt.chromaFormat != CF420) {
      custom_color_diff_format_flag = true;
      color_diff_format = fmt.chromaFormat;
    }
    if (fmt.frameRate != FR24000_1001) {
      custom_frame_rate_flag = true;
      frame_rate = fmt.frameRate;
    }
    if (fmt.bitdepth != 8) {
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
  }

  if (fmt.interlace)
    picture_coding_mode = 1;
  else
    picture_coding_mode = 0;
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
    switch (fmt.color_diff_format) {
    case CF422:
      ss << UnsignedVLC(1);
      break;
    case CF420:
      ss << UnsignedVLC(2);
      break;
    case CF444:
    case RGB:
    default:
      ss << UnsignedVLC(0);
      break;
    }
  }

  ss << Boolean(fmt.custom_scan_format_flag);
  if (fmt.custom_scan_format_flag) {
    ss << UnsignedVLC(fmt.source_sampling);
  }

  ss << Boolean(fmt.custom_frame_rate_flag);
  if (fmt.custom_frame_rate_flag) {
    switch (fmt.frame_rate) {
    case FR24000_1001:
      ss << UnsignedVLC(1);
      break;
    case FR24:
      ss << UnsignedVLC(2);
      break;
    case FR25:
      ss << UnsignedVLC(3);
      break;
    case FR30000_1001:
      ss << UnsignedVLC(4);
      break;
    case FR30:
      ss << UnsignedVLC(5);
      break;
    case FR50:
      ss << UnsignedVLC(6);
      break;
    case FR60000_1001:
      ss << UnsignedVLC(7);
      break;
    case FR60:
      ss << UnsignedVLC(8);
      break;
    case FR15000_1001:
      ss << UnsignedVLC(9);
      break;
    case FR25_2:
      ss << UnsignedVLC(10);
      break;
    case FR48:
      ss << UnsignedVLC(11);
      break;
    case FR48_1001:
      ss << UnsignedVLC(12);
      break;
    case FR96:
      ss << UnsignedVLC(13);
      break;
    case FR100:
      ss << UnsignedVLC(14);
      break;
    case FR120_1001:
      ss << UnsignedVLC(15);
      break;
    case FR120:
      ss << UnsignedVLC(16);
      break;
    default:
      {
        std::stringstream ss;
        ss << "DataUnitIO: Invalid Frame Rate on output: " << fmt.frame_rate;
        throw std::logic_error(ss.str());
      }
    }
  }

  ss << Boolean(false); // custom_pixel_aspect_ratio_flag

  ss << Boolean(false); // custom_clean_area_flag

  ss << Boolean(fmt.custom_signal_range_flag);
  if (fmt.custom_signal_range_flag) {
    ss << UnsignedVLC(fmt.bitdepth);
  }

  if (fmt.custom_color_diff_format_flag && fmt.color_diff_format == RGB) {
    ss << Boolean(true);
    ss << UnsignedVLC(0); // Custom
    ss << Boolean(false); // No custom primaries
    ss << Boolean(true) << UnsignedVLC(3); // RGB Matrix
    ss << Boolean(false); // No custom gamma
  } else 
    ss << Boolean(false); // custom_color_spec_flag

  ss << UnsignedVLC(fmt.picture_coding_mode);

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
    switch ((unsigned int)color_diff_format) {
    case 0:
      fmt.color_diff_format = CF444;
      break;
    case 1:
      fmt.color_diff_format = CF422;
      break;
    case 2:
      fmt.color_diff_format = CF420;
      break;
    default:
      throw std::logic_error("DataUnitIO: Invalid Color Format");
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

  Boolean custom_frame_rate_flag;
  stream >> custom_frame_rate_flag;
  fmt.custom_frame_rate_flag = custom_frame_rate_flag;
  if (custom_frame_rate_flag) {
    UnsignedVLC index;
    stream >> index;
    switch(index) {
    case 1: fmt.frame_rate = FR24000_1001; break;
    case 2: fmt.frame_rate = FR24; break;
    case 3: fmt.frame_rate = FR25; break;
    case 4: fmt.frame_rate = FR30000_1001; break;
    case 5: fmt.frame_rate = FR30; break;
    case 6: fmt.frame_rate = FR50; break;
    case 7: fmt.frame_rate = FR60000_1001; break;
    case 8: fmt.frame_rate = FR60; break;
    case 9: fmt.frame_rate = FR15000_1001; break;
    case 10: fmt.frame_rate = FR25_2; break;
    case 11: fmt.frame_rate = FR48; break;
    case 12: fmt.frame_rate = FR48_1001; break;
    case 13: fmt.frame_rate = FR96; break;
    case 14: fmt.frame_rate = FR100; break;
    case 15: fmt.frame_rate = FR120_1001; break;
    case 16: fmt.frame_rate = FR120; break;
    default:
      {
        std::stringstream ss;
        ss << "DataUnitIO: Invalid Frame Rate on Input: " << (int)index;
        throw std::logic_error(ss.str());
      }
    }
  }

  Boolean custom_pixel_aspect_ratio_flag;
  stream >> custom_pixel_aspect_ratio_flag;
  if (custom_pixel_aspect_ratio_flag) {
    throw std::logic_error("DataUnitIO: custom_pixel_aspect_ratio_flag set, shouldn't be");
  }

  Boolean custom_clean_area_flag;
  stream >> custom_clean_area_flag;
  if (custom_clean_area_flag) {
    throw std::logic_error("DataUnitIO: custom_clean_area_flag set, shouldn't be");
  }

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
  if (custom_color_spec_flag) {
    UnsignedVLC custom_color_spec_index;
    stream >> custom_color_spec_index;
    if (custom_color_spec_index == 0) {
      Boolean custom_color_primaries_flag, custom_color_matrix_flag, custom_transfer_function_flag;
      stream >> custom_color_primaries_flag;
      if (custom_color_primaries_flag) {
        UnsignedVLC color_primaries;
        stream >> color_primaries;
      }

      stream >> custom_color_matrix_flag;
      if (custom_color_matrix_flag) {
        UnsignedVLC color_matrix;
        stream >> color_matrix;

        if (color_matrix == 3) {
          fmt.custom_color_diff_format_flag = true;
          fmt.color_diff_format = RGB;
        }
      }

      stream >> custom_transfer_function_flag;
      if (custom_transfer_function_flag) {
        UnsignedVLC transfer_function;
        stream >> transfer_function;
      }
    }
  }

  
  UnsignedVLC picture_coding_mode;
  stream >> picture_coding_mode;
  fmt.picture_coding_mode = picture_coding_mode;

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
  }

  Bytes next_parse_offset(4);
  Bytes prev_parse_offset(4);

  stream >> next_parse_offset >> prev_parse_offset;

  if ((unsigned long) next_parse_offset == 0) {
    d.strm.str(std::string());
  } else {
    int bufsize = ((unsigned long) next_parse_offset) - 13;
    char *buf = new char[bufsize];
    if (buf == NULL)
      throw std::logic_error("DataUnitIO: Couldn't allocate memory for Data Unit");
    stream.read(buf, bufsize);
    d.strm.str(std::string(buf, bufsize));
    delete[] buf;
  }

  Bytes prefix(4);
  stream >> prefix;

  d.strm.copyfmt(stream);

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
  SequenceHeader *other;

  switch (fmt.base_video_format) {
  case  0: other = new SequenceHeader(PROFILE_UNKNOWN, 480,  640,  CF420, false, FR24000_1001, false,  8); break;
  case  1: other = new SequenceHeader(PROFILE_UNKNOWN, 120,  176,  CF420, false, FR15000_1001, false,  8); break;
  case  2: other = new SequenceHeader(PROFILE_UNKNOWN, 144,  176,  CF420, false, FR25_2,       true,   8); break;
  case  3: other = new SequenceHeader(PROFILE_UNKNOWN, 240,  352,  CF420, false, FR15000_1001, false,  8); break;
  case  4: other = new SequenceHeader(PROFILE_UNKNOWN, 288,  352,  CF420, false, FR25_2,       true,   8); break;
  case  5: other = new SequenceHeader(PROFILE_UNKNOWN, 480,  704,  CF420, false, FR15000_1001, false,  8); break;
  case  6: other = new SequenceHeader(PROFILE_UNKNOWN, 576,  704,  CF420, false, FR25_2,       true,   8); break;
  case  7: other = new SequenceHeader(PROFILE_UNKNOWN, 480,  720,  CF422, true,  FR30000_1001, false, 10); break;
  case  8: other = new SequenceHeader(PROFILE_UNKNOWN, 576,  720,  CF422, true,  FR25,         true,  10); break;
  case  9: other = new SequenceHeader(PROFILE_UNKNOWN, 720,  1280, CF422, false, FR60000_1001, true,  10); break;
  case 10: other = new SequenceHeader(PROFILE_UNKNOWN, 720,  1280, CF422, false, FR50,         true,  10); break;
  case 11: other = new SequenceHeader(PROFILE_UNKNOWN, 1080, 1920, CF422, true,  FR30000_1001, true,  10); break;
  case 12: other = new SequenceHeader(PROFILE_UNKNOWN, 1080, 1920, CF422, true,  FR25,         true,  10); break;
  case 13: other = new SequenceHeader(PROFILE_UNKNOWN, 1080, 1920, CF422, false, FR60000_1001, true,  10); break;
  case 14: other = new SequenceHeader(PROFILE_UNKNOWN, 1080, 1920, CF422, false, FR50,         true,  10); break;
  case 15: other = new SequenceHeader(PROFILE_UNKNOWN, 1080, 2048, CF444, false, FR24,         true,  12); break;
  case 16: other = new SequenceHeader(PROFILE_UNKNOWN, 2160, 4096, CF444, false, FR24,         true,  12); break;
  case 17: other = new SequenceHeader(PROFILE_UNKNOWN, 2160, 3840, CF422, false, FR60000_1001, true,  10); break;
  case 18: other = new SequenceHeader(PROFILE_UNKNOWN, 2160, 3840, CF422, false, FR50,         true,  10); break;
  case 19: other = new SequenceHeader(PROFILE_UNKNOWN, 4320, 7680, CF422, false, FR60000_1001, true,  10); break;
  case 20: other = new SequenceHeader(PROFILE_UNKNOWN, 4320, 7680, CF422, false, FR50,         true,  10); break;
  case 21: other = new SequenceHeader(PROFILE_UNKNOWN, 1080, 1920, CF422, false, FR24000_1001, true,  10); break;
  case 22: other = new SequenceHeader(PROFILE_UNKNOWN, 486,  720,  CF422, true,  FR30000_1001, false, 10); break;
  default:
    throw std::logic_error("DataUnitIO: unknown base video format");
  }

  hdr->profile       = other->profile;
  hdr->width         = other->width;
  hdr->height        = other->height;
  hdr->chromaFormat  = other->chromaFormat;
  hdr->interlace     = other->interlace;
  hdr->frameRate     = other->frameRate;
  hdr->topFieldFirst = other->topFieldFirst;
  hdr->bitdepth      = other->bitdepth;
  delete other;

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
  if (fmt.custom_frame_rate_flag) {
    hdr->frameRate = fmt.frame_rate;

    if (fmt.frame_rate > MAX_V2_FRAMERATE)
      if (hdr->major_version < 3)
        hdr->major_version = 3;
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
