/*********************************************************************/
/* Picture.cpp                                                       */
/* Author: Tim Borer                                                 */
/* This version 7th June 2012                                        */
/*                                                                   */
/* Defines stuff related to pictures.                                */
/* Pictures, in this context, are three components (YUV/RGB), plus   */
/* some metadata specifying the size and colour subsampling (if any) */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

#include <iostream>
#include <string>
#include <stdexcept> // For invalid_argument

#include "Picture.h"
#include "FrameResolutions.h" //List of frame resolutions (frameResolutions)
#include "Utils.h"

std::ostream& operator<<(std::ostream& os, ColourFormat format) {
  const char* s;
  switch (format) {
    case CF444:
      s = "4:4:4";
      break;
    case CF422:
      s = "4:2:2";
      break;
    case CF420:
      s = "4:2:0";
      break;
    case RGB:
      s = "RGB";
      break;
    default:
      s = "Unknown colour format!";
      break;
  }
  return os<<s;
}

std::istream& operator>>(std::istream& strm, ColourFormat& format) {
        std::string text;
        strm >> text;
        if (text == "4:4:4") format = CF444;
        else if (text == "4:2:2") format = CF422;
        else if (text == "4:2:0") format = CF420;
        else if (text == "RGB") format = RGB;
        else format = UNKNOWN;
        // Alternatively
        // else strm.setstate(std::ios_base::badbit||std::ios_base::failbit);
        // else throw std::invalid_argument("invalid colour format");
        return strm;
}

void PictureFormat::construct(int height, int width, ColourFormat cFormat) {
  yHeight = height;
  yWidth = width;
  uvFormat = cFormat;
  switch (uvFormat) {
    case RGB:
    case CF444: 
      uvHeight = yHeight;
      uvWidth = yWidth;
      break;
    case CF422:
      uvHeight = yHeight;
      uvWidth = yWidth/2;
      break;
    case CF420:
      uvHeight = yHeight/2;
      uvWidth = yWidth/2;
      break;
    case UNKNOWN:
      uvHeight = 0;
      uvWidth = 0;
      break;
    default:
      throw std::invalid_argument("Invalid colour format");
  }
}

PictureFormat::PictureFormat() {construct(0, 0, UNKNOWN);} // Used for arrays of Pictures

PictureFormat::PictureFormat(int height, int width, ColourFormat cFormat) {
  construct(height, width, cFormat);
}

PictureFormat::PictureFormat(int lumaHeight, int lumaWidth,
                             int chromaHeight, int chromaWidth,
                             ColourFormat format):
  yHeight(lumaHeight), yWidth(lumaWidth),
  uvHeight(chromaHeight), uvWidth(chromaWidth),
  uvFormat(format) {
}


PictureFormat::operator const bool() const {
  return (yHeight*yWidth*uvFormat)!=0;
}

const ColourFormat PictureFormat::chromaFormat() const {return uvFormat;}

const Shape2D PictureFormat::lumaShape() const {
  Shape2D result = {{yHeight, yWidth}};
  return result; }

const Shape2D PictureFormat::chromaShape() const {
  Shape2D result = {{uvHeight, uvWidth}};
  return result; }

const int PictureFormat::lumaHeight() const {return yHeight;}

const int PictureFormat::lumaWidth() const {return yWidth;}

const int PictureFormat::chromaHeight() const {
  return uvHeight;
}

const int PictureFormat::chromaWidth() const {
  return uvWidth;
}

const int PictureFormat::samples() const {
  return (lumaHeight()*lumaWidth()+2*chromaHeight()*chromaWidth());
}

void PictureFormat::guessFormat(int imageSamples, ColourFormat cFormat) {
  const int resolutions = sizeof(frameResolutions)/sizeof(int);
  for (int r=0; r<resolutions; ++r) {
    construct(frameResolutions[r][0], frameResolutions[r][1], cFormat);
    if (samples()==imageSamples) return;
  }
  construct (0,0,UNKNOWN); //No matching format found
}

PictureFormat::PictureFormat(int height, int width, ColourFormat cFormat, int imageSamples) {
  if (height && width && cFormat) { //frame format explicity defined
    construct(height, width, cFormat);
    if (samples()!=imageSamples) construct(0,0,UNKNOWN);
    return;
  }
  if (height && width) { //colour format unknown (guess it)
    construct(height, width, CF444);
    if (samples()!=imageSamples) construct(height, width, CF422);
    if (samples()!=imageSamples) construct(height, width, CF420);
    if (samples()!=imageSamples) construct(0,0,UNKNOWN);
    return;
  }
  if (cFormat) { //only colour format available (no height or no width)
    guessFormat(imageSamples, cFormat);
    return;
  }
  // guess colour format too (priority to 4:4:4)
  guessFormat(imageSamples, CF444);
  if (!(*this)) guessFormat(imageSamples, CF422);
  if (!(*this)) guessFormat(imageSamples, CF420);
}

Picture::Picture() {}

Picture::Picture(const PictureFormat& f):
  picFormat(f),
  luma(f.lumaShape()),
  chroma1(f.chromaShape()),
  chroma2(f.chromaShape()) {
}

Picture::Picture(const int height, const int width, const ColourFormat format):
  picFormat(height, width, format),
  luma(picFormat.lumaShape()),
  chroma1(picFormat.chromaShape()),
  chroma2(picFormat.chromaShape()) {
}


Picture::Picture(const PictureFormat& f,
                 const Array2D& yArray,
                 const Array2D& c1Array,
                 const Array2D& c2Array):
  picFormat(f) {
    y(yArray);
    c1(c1Array);
    c2(c2Array);
}

PictureFormat Picture::format() const {
  return picFormat;
}

const Array2D& Picture::y() const {
  return luma;
}

const Array2D& Picture::c1() const {
  return chroma1;
}

const Array2D& Picture::c2() const {
  return chroma2;
}

void Picture::y(const Array2D& arg) {
  if (shape(arg)[0]!=picFormat.lumaHeight()) {
    throw std::invalid_argument("wrong luma height");
  }
  if (shape(arg)[1]!=picFormat.lumaWidth()) {
    throw std::invalid_argument("wrong luma width");
  }
  luma = arg;
}

void Picture::c1(const Array2D& arg) {
  if (shape(arg)[0]!=picFormat.chromaHeight()) {
    throw std::invalid_argument("wrong chroma height");
  }
  if (shape(arg)[1]!=picFormat.chromaWidth()) {
    throw std::invalid_argument("wrong chroma width");
  }
  chroma1 = arg;
}

void Picture::c2(const Array2D& arg) {
  if (shape(arg)[0]!=picFormat.chromaHeight()) {
    throw std::invalid_argument("wrong chroma height");
  }
  if (shape(arg)[1]!=picFormat.chromaWidth()) {
    throw std::invalid_argument("wrong chroma width");
  }
  chroma2 = arg;
}

// Get the shape of a 2D PictureArray
const Shape2D shape(const PictureArray& arg) {
  const Shape2D result = {{static_cast<Index>(arg.shape()[0]), static_cast<Index>(arg.shape()[1])}};
  return result;
}

const PictureArray split_into_blocks(const Picture& picture, int ySlices, int xSlices) {
  const Shape2D shape = {{ySlices, xSlices}};
  PictureArray slices(shape);
  const BlockArray luma = split_into_blocks(picture.y(), ySlices, xSlices);
  const BlockArray chroma1 = split_into_blocks(picture.c1(), ySlices, xSlices);
  const BlockArray chroma2 = split_into_blocks(picture.c2(), ySlices, xSlices);
  const ColourFormat colourFormat = picture.format().chromaFormat();
  for (int y=0; y<ySlices; ++y) {
    for (int x=0; x<xSlices; ++x) {
      const int sliceHeight = luma[y][x].shape()[0];
      const int sliceWidth = luma[y][x].shape()[1];
      const PictureFormat sliceformat(sliceHeight, sliceWidth, colourFormat);
      slices[y][x] = Picture(sliceformat, luma[y][x], chroma1[y][x], chroma2[y][x]);
    }
  }
  return slices;
}

const Picture merge_blocks(const PictureArray& blocks) {
  const Shape2D blocksShape = shape(blocks);
  BlockArray lumaBlocks(blocksShape);
  BlockArray chroma1Blocks(blocksShape);
  BlockArray chroma2Blocks(blocksShape);
  const int ySlices = blocksShape[0];
  const int xSlices = blocksShape[1];
  for (int y=0; y<ySlices; ++y) {
    for (int x=0; x<xSlices; ++x) {
      lumaBlocks[y][x] = blocks[y][x].y();
      chroma1Blocks[y][x] = blocks[y][x].c1();
      chroma2Blocks[y][x] = blocks[y][x].c2();
    }
  }
  const Array2D luma = merge_blocks(lumaBlocks);
  const Array2D chroma1 = merge_blocks(chroma1Blocks);
  const Array2D chroma2 = merge_blocks(chroma2Blocks);
  const int height = luma.shape()[0];
  const int width = luma.shape()[1];
  const ColourFormat colourFormat = blocks[0][0].format().chromaFormat();
  const PictureFormat pictureFormat(height, width, colourFormat);
  return Picture(pictureFormat, luma, chroma1, chroma2);
}

// Clip a Picture to specified limits
// First function clips all components to the same values (good for RGB)
const Picture clip(const Picture& picture, const int min_value, const int max_value) {
  Picture result(picture.format());
  result.y(clip(picture.y(), min_value, max_value));
  result.c1(clip(picture.c1(), min_value, max_value));
  result.c2(clip(picture.c2(), min_value, max_value));
  return result;
}

// Second function clips luma and chroma values separately (good for YUV)
const Picture clip(const Picture& picture,
                   const int luma_min, const int luma_max,
                   const int chroma_min, const int chroma_max) {
  Picture result(picture.format());
  result.y(clip(picture.y(), luma_min, luma_max));
  result.c1(clip(picture.c1(), chroma_min, chroma_max));
  result.c2(clip(picture.c2(), chroma_min, chroma_max));
  return result;
}

//**************** IO functions ****************//

namespace {

  // Functions to extract stored values from streams
  // Note the default values are zero if they have not been explicitly set

  // Choice of UNSIGNED, SIGNED, OFFSET, TEXT (see ArrayIO.h)
  long& luma_format(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  // Choice of UNSIGNED, SIGNED, OFFSET, TEXT (see ArrayIO.h)
  long& chroma_format(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  // Returns true if luma and chroma format have been set
  long& picture_format(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  long& luma_bit_depth(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  long& chroma_bit_depth(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  // Returns true if luma and chroma bit depth have been set
  long& picture_bit_depth(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  long& luma_offset(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  long& chroma_offset(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  // Returns true if luma and chroma offset have been set
  long& picture_offset(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

} // end unnamed namespace

// io format manipulator
void pictureio::format::operator()(std::ios_base& stream) const {
  luma_format(stream) = static_cast<long>(lumaFormat);
  chroma_format(stream) = static_cast<long>(chromaFormat);
  picture_format(stream) = static_cast<long>(true);
}

// ostream io format manipulator
std::ostream& operator << (std::ostream& stream, pictureio::format f) {
  f(stream);
  return stream;
}

// istream io format manipulator
std::istream& operator >> (std::istream& stream, pictureio::format f) {
  f(stream);
  return stream;
}

// bit depth format manipulator
void pictureio::bitDepth::operator()(std::ios_base& stream) const {
  luma_bit_depth(stream) = static_cast<long>(lumaBitDepth);
  chroma_bit_depth(stream) = static_cast<long>(chromaBitDepth);
  picture_bit_depth(stream) = static_cast<long>(true);
}

// ostream bit depth format manipulator
std::ostream& operator << (std::ostream& stream, pictureio::bitDepth d) {
  d(stream);
  return stream;
}

// istream bit depth format manipulator
std::istream& operator >> (std::istream& stream, pictureio::bitDepth d) {
  d(stream);
  return stream;
}

// offset format manipulator
void pictureio::offset::operator()(std::ios_base& stream) const {
  luma_offset(stream) = static_cast<long>(lumaOffset);
  chroma_offset(stream) = static_cast<long>(chromaOffset);
  picture_offset(stream) = static_cast<long>(true);
}

// ostream offset format manipulator
std::ostream& operator << (std::ostream& stream, pictureio::offset o) {
  o(stream);
  return stream;
}

// istream offset format manipulator
std::istream& operator >> (std::istream& stream, pictureio::offset o) {
  o(stream);
  return stream;
}

std::istream& operator >> (std::istream& stream, Picture& frame) {
  using arrayio::ioFormat;
  using arrayio::format;
  using arrayio::bitDepth;
  using arrayio::offset;
  // Set luma data format, bit depth and offset
  if (picture_format(stream)) stream >> format(static_cast<ioFormat>(luma_format(stream)));
  if (picture_bit_depth(stream)) stream >> bitDepth(luma_bit_depth(stream));
  if (picture_offset(stream)) stream >> offset(luma_offset(stream));
  stream >> frame.luma;
  // Set chroma data format, bit depth and offset
  if (picture_format(stream)) stream >> format(static_cast<ioFormat>(chroma_format(stream)));
  if (picture_bit_depth(stream)) stream >> bitDepth(chroma_bit_depth(stream));
  if (picture_offset(stream)) stream >> offset(chroma_offset(stream));
  stream >> frame.chroma1 >> frame.chroma2;
  return stream;
}

std::ostream& operator << (std::ostream& stream, const Picture& frame) {
  using arrayio::ioFormat;
  using arrayio::format;
  using arrayio::bitDepth;
  using arrayio::offset;
  // Set luma data format, bit depth and offset
  if (picture_format(stream)) stream << format(static_cast<ioFormat>(luma_format(stream)));
  if (picture_bit_depth(stream)) stream << bitDepth(luma_bit_depth(stream));
  if (picture_offset(stream)) stream << offset(luma_offset(stream));
  stream << frame.y();
  // Set chroma data format, bit depth and offset
  if (picture_format(stream)) stream << format(static_cast<ioFormat>(chroma_format(stream)));
  if (picture_bit_depth(stream)) stream << bitDepth(chroma_bit_depth(stream));
  if (picture_offset(stream)) stream << offset(chroma_offset(stream));
  stream << frame.c1() << frame.c2();
  return stream;
}
