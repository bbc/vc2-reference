/*********************************************************************/
/* DataUnit.h                                                        */
/* Author: James Weaver                                              */
/* This version 17th June 2015                                       */
/*                                                                   */
/* Declares stuff related to VC2 Data Units.                         */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

#ifndef DATAUNIT_17JUN15
#define DATAUNIT_17JUN15

#include <iosfwd>

#include "Utils.h"
#include "Picture.h"
#include "Slices.h"
#include "WaveletTransform.h"

enum DataUnitType {
  UNKNOWN_DATA_UNIT,
  SEQUENCE_HEADER,
  END_OF_SEQUENCE,
  AUXILIARY_DATA,
  PADDING_DATA,
  HQ_PICTURE,
  LD_PICTURE
};

class DataUnit {
  public:
    DataUnit();

    std::istream &stream();
    DataUnitType type;

    friend std::istream& operator >> (std::istream& stream, DataUnit &d);

  protected:
    std::istringstream strm;
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

enum FrameRate { FR0, FR24000_1001, FR24, FR25, FR30000_1001, FR30, FR50, FR60000_1001, FR60, FR15000_1001, FR25_2, FR48 };
enum Profile { PROFILE_UNKNOWN, PROFILE_LD, PROFILE_HQ };

class SequenceHeader {
public:
  SequenceHeader();
  SequenceHeader(Profile profile, int height, int width, ColourFormat chromaFormat, bool interlace, FrameRate frameRate, bool topFieldFirst, int bitdepth);

  int major_version;
  int minor_version;
  Profile profile;
  int width;
  int height;
  ColourFormat chromaFormat;
  bool interlace;
  FrameRate frameRate;
  bool topFieldFirst;
  int bitdepth;
};

class PicturePreamble {
 public:
  PicturePreamble();
  PicturePreamble(const unsigned long picture_number,
                  const WaveletKernel wavelet_kernel,
                  const int depth,
                  const int slices_x,
                  const int slices_y,
                  const int slice_prefix,
                  const int slice_size_scalar);

  unsigned long picture_number;
  WaveletKernel wavelet_kernel;
  int depth;
  int slices_x;
  int slices_y;
  int slice_prefix;
  int slice_size_scalar;
  utils::Rational slice_bytes;
};

namespace dataunitio {
  using sliceio::highQualityCBR;
  using sliceio::highQualityVBR;

  std::istream& lowDelay(std::istream& stream);

  std::istream& synchronise(std::istream& stream);

  std::ostream& start_sequence(std::ostream& stream);
  std::ostream& end_sequence(std::ostream& stream);
};

std::ostream& operator << (std::ostream& stream, const WrappedPicture& d);

std::ostream& operator << (std::ostream& stream, const SequenceHeader& s);

std::ostream& operator << (std::ostream& stream, const DataUnitType& t);

std::ostream& operator << (std::ostream& stream, const FrameRate& t);

std::istream& operator >> (std::istream& stream, DataUnit &d);

std::istream& operator >> (std::istream& stream, SequenceHeader &hdr);

std::istream& operator >> (std::istream& stream, PicturePreamble &hdr);

#endif /* DATAUNIT_17JUN15 */
