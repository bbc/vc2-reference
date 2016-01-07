/***********************************************************************/
/* Arrays.cpp                                                          */
/* Author: Tim Borer                                                   */
/* This version 8th August 2011                                        */
/*                                                                     */
/* Defines stuff related to mulitdimensional arrays.                   */
/* 20th June 2011: File copied from version 7th October 2010, IO added */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/***********************************************************************/

#include <ostream>
#include <istream>
#include <stdexcept>

#include "Arrays.h"
#include "Utils.h"

using utils::pow;

// Get the shape of a 2D array
const Shape2D shape(const Array2D& arg) {
  const Shape2D result = {{static_cast<Index>(arg.shape()[0]), static_cast<Index>(arg.shape()[1])}};
  return result;
}

const Shape2D shape(const View2D& arg) {
  const Shape2D result = {{static_cast<Index>(arg.shape()[0]), static_cast<Index>(arg.shape()[1])}};
  return result;
}

const Shape2D shape(const ConstView2D& arg) {
  const Shape2D result = {{static_cast<Index>(arg.shape()[0]), static_cast<Index>(arg.shape()[1])}};
  return result;
}

const Shape2D shape(const BlockArray& arg) {
  const Shape2D result = {{static_cast<Index>(arg.shape()[0]), static_cast<Index>(arg.shape()[1])}};
  return result;
}

const Array2D clip(const Array2D& values, const int min_value, const int max_value) {
  Array2D result(values.ranges());
  const Index height = values.shape()[0];
  const Index width = values.shape()[1];
  for (int y=0; y<height; ++y) {
    for (int x=0; x<width; ++x) {
      if (values[y][x]<min_value) result[y][x] = min_value;
      else if (values[y][x]>max_value) result[y][x] =  max_value;
      else result[y][x] = values[y][x];
    }
  }
  return result;
}

// Splits a large 2D array into an array of smaller 2D arrays (blocks)
// Note that if the number of blocks is not a sub-multiple of the input array dimensions then
// the blocks will have different sizes!
// yBlocks and xBlocks are the number of blocks in the vertical and horizontal dimension respectively.
// Splits a picture into slices or a subband into codeblocks.
const BlockArray split_into_blocks(const Array2D& picture, int yBlocks, int xBlocks) {
  // Define array of yBlocks by xBlocks
  BlockArray blocks(extents[yBlocks][xBlocks]);
  const int pictureHeight = picture.shape()[0];
  const int pictureWidth = picture.shape()[1];
  // Note Range(left, right) defines the half open range [left, right),
  // i.e. the rightmost element is not included
  for (int y=0, top=0, bottom=pictureHeight/yBlocks;
       y<yBlocks;
       ++y, top=bottom, bottom=((y+1)*pictureHeight/yBlocks) ) {
    for (int x=0, left=0, right=pictureWidth/xBlocks;
         x<xBlocks;
         ++x, left=right, right=((x+1)*pictureWidth/xBlocks) ) {
      // Assign region of picture to block
      blocks[y][x] = picture[indices[Range(top,bottom)][Range(left,right)]];
    }
  }
  return blocks;
}

// Converts an array of blocks back to a single 2D array
// This is the inverse of "split_into_blocks()"
// The array of blocks might represent slices or codeblocks
const Array2D merge_blocks(const BlockArray& blocks) {
  // First find picture dimensions
  int pictureHeight = 0;
  int pictureWidth = 0;
  const int yBlocks = blocks.shape()[0];
  const int xBlocks = blocks.shape()[1];
  for (int y=0; y<yBlocks; ++y) pictureHeight += blocks[y][0].shape()[0];
  for (int x=0; x<xBlocks; ++x) pictureWidth += blocks[0][x].shape()[1];
  // Define Array2D pictureHeight by PictureWidth
  Array2D picture(extents[pictureHeight][pictureWidth]);
  // Now merge the block together
  // Note Range(top, bottom) and Range(left, right) define half open ranges either
  // [top, bottom) or [left, right), i.e. the bottom/rightmost element is not included
  int bottom;
  for (int y=0, top=0;
       y<yBlocks;
       ++y, top=bottom) {
    bottom = top + blocks[y][0].shape()[0];
    int right;
    for (int x=0, left=0;
         x<xBlocks;
         ++x, left=right) {
      right = left + blocks[0][x].shape()[1];
      // Copy block into a view of the picture
      picture[indices[Range(top,bottom)][Range(left,right)]] = blocks[y][x];
    }
  }
  return picture;
}

//**************** IO functions ****************//

namespace {

  // Functions to extract stored values from streams
  // Note the default values are zero if they have not been explicitly set

  long& word_width(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  long& bit_depth(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  long& is_signed(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  long& is_offset(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  long& is_right_justified(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

#if 0
  long& is_text(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }
#endif

  long& zero_level_set(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  long& zero_level(std::ios_base& stream) {
      static const int i = std::ios_base::xalloc();
      return stream.iword(i);
  }

  const int ioBytes(std::ios_base& stream) {
    return (word_width(stream) ? word_width(stream) : sizeof(int));
  }

  const int ioDepth(std::ios_base& stream) {
    return (bit_depth(stream) ? bit_depth(stream) : 8*ioBytes(stream));
  }

  const int ioShift(std::ios_base& stream) {
    return (is_right_justified(stream) ? 0 : 8*ioBytes(stream)-ioDepth(stream));
  }

  const int ioZero(std::ios_base& stream) {
    const int zeroLevel =
      (zero_level_set(stream) ? zero_level(stream) : pow(2, ioDepth(stream)-1));
    return (is_offset(stream) ? zeroLevel : 0);
  }

} // end unnamed namespace

// io format manipulator to set the io data format
void arrayio::format::operator()(std::ios_base& stream) const {
  switch (iof) {
    case UNSIGNED:
      is_signed(stream) = static_cast<long>(false);
      is_offset(stream) = static_cast<long>(false);
      break;
    case SIGNED:
      is_signed(stream) = static_cast<long>(true);
      is_offset(stream) = static_cast<long>(false);
      break;
    case OFFSET:
      is_signed(stream) = static_cast<long>(false);
      is_offset(stream) = static_cast<long>(true);
      break;
    default:
      throw std::invalid_argument("invalid array io format");
      break;
  }
}

// ostream ioformat format manipulator
std::ostream& operator << (std::ostream& stream, arrayio::format f) {
  f(stream);
  return stream;
}

// istream ioformat format manipulator
std::istream& operator >> (std::istream& stream, arrayio::format f) {
  f(stream);
  return stream;
}

// wordWidth manipulator to set number of bytes per word
void arrayio::wordWidth::operator()(std::ios_base& stream) const {
  word_width(stream) = static_cast<long>(no_of_bytes);
}

// ostream wordWidth format manipulator
std::ostream& operator << (std::ostream& stream, arrayio::wordWidth w) {
  w(stream);
  return stream;
}

// istream wordWidth format manipulator
std::istream& operator >> (std::istream& stream, arrayio::wordWidth w) {
  w(stream);
  return stream;
}

// bitDepth manipulator to set number of bits per word
void arrayio::bitDepth::operator()(std::ios_base& stream) const {
  bit_depth(stream) = static_cast<long>(no_of_bits);
}

// ostream bitDepth format manipulator
std::ostream& operator << (std::ostream& stream, arrayio::bitDepth d) {
  d(stream);
  return stream;
}

// istream bitDepth format manipulator
std::istream& operator >> (std::istream& stream, arrayio::bitDepth d) {
  d(stream);
  return stream;
}

// bitDepth manipulator to set number of bits per word
void arrayio::offset::operator()(std::ios_base& stream) const {
  zero_level_set(stream) = static_cast<long>(true);
  zero_level(stream) = static_cast<long>(zeroValue);
}

// ostream bitDepth format manipulator
std::ostream& operator << (std::ostream& stream, arrayio::offset z) {
  z(stream);
  return stream;
}

// istream bitDepth format manipulator
std::istream& operator >> (std::istream& stream, arrayio::offset z) {
  z(stream);
  return stream;
}

// ostream left justified format manipulator
std::ostream& arrayio::left_justified(std::ostream& stream) {
  is_right_justified(stream) = static_cast<long>(false);
  return stream;
}

// istream left justified format manipulator
std::istream& arrayio::left_justified(std::istream& stream) {
  is_right_justified(stream) = static_cast<long>(false);
  return stream;
}

// ostream right justified format manipulator
std::ostream& arrayio::right_justified(std::ostream& stream) {
  is_right_justified(stream) = static_cast<long>(true);
  return stream;
}

// istream right justified format manipulator
std::istream& arrayio::right_justified(std::istream& stream) {
  is_right_justified(stream) = static_cast<long>(true);
  return stream;
}

// Set data format to offset binary
std::ostream& arrayio::offset_binary(std::ostream& stream) {
  is_signed(stream) = static_cast<long>(false);
  is_offset(stream) = static_cast<long>(true);
  return stream;
}

// Set data format to offset binary
std::istream& arrayio::offset_binary(std::istream& stream) {
  is_signed(stream) = static_cast<long>(false);
  is_offset(stream) = static_cast<long>(true);
  return stream;
}

// Set data format to signed binary
std::ostream& arrayio::signed_binary(std::ostream& stream) {
  is_signed(stream) = static_cast<long>(true);
  is_offset(stream) = static_cast<long>(false);
  return stream;
}

// Set data format to signed binary
std::istream& arrayio::signed_binary(std::istream& stream) {
  is_signed(stream) = static_cast<long>(true);
  is_offset(stream) = static_cast<long>(false);
  return stream;
}

// Set data format to unsigned binary
std::ostream& arrayio::unsigned_binary(std::ostream& stream) {
  is_signed(stream) = static_cast<long>(false);
  is_offset(stream) = static_cast<long>(false);
  return stream;
}

// Set data format to unsigned binary
std::istream& arrayio::unsigned_binary(std::istream& stream) {
  is_signed(stream) = static_cast<long>(false);
  is_offset(stream) = static_cast<long>(false);
  return stream;
}

std::istream& operator >> (std::istream& stream, Array2D& array) {
  // Get word width and bit depth from stream
  // Default word width is size of int, default bit depth fills word width
  const int wordBytes = ioBytes(stream);

  //Create reference for input stream buffer (input via stream buffer for efficiency).
  std::streambuf& inbuf = *(stream.rdbuf());
  // Read wordBytes bytes per array element
  const int size = wordBytes*array.num_elements();
  unsigned char *inBuffer = new unsigned char[size];
  std::istream::sentry s(stream, true);
  if (s) {
    if ( inbuf.sgetn(reinterpret_cast<char*>(inBuffer), size) < size )
      stream.setstate(std::ios_base::eofbit|std::ios_base::failbit);
  }

  const int height = array.shape()[0];
  const int width = array.shape()[1];
  const int shift = ioShift(stream);
  const int offset = ioZero(stream);
  unsigned int value, byte = 0;
  for (int y=0; y<height; ++y) {
    for (int x=0; x<width; ++x) {
      value = 0;
      switch (wordBytes) {
        // Only allowed 4 bytes word width in 32 bit systems
        case 4:
          value |= (inBuffer[byte++]<<24);
        case 3:
          value |= (inBuffer[byte++]<<16);
        case 2:
          value |= (inBuffer[byte++]<<8);
        case 1:
          value |= inBuffer[byte++];
          break;
        default:
          throw std::domain_error("Word width of input stream must be in range 1 to 4");
      }
      if (!is_signed(stream)) value >>= shift; //Use logical shift for unsigned data (value is unsigned int)
      array[y][x] = value;
      if (is_signed(stream)) array[y][x] >>= shift; //Use arithmetic shift for signed data (array[y][x] is int)
      if (is_offset(stream)) array[y][x] -= offset;
    }
  }
  delete[] inBuffer;
  return stream;
}

std::ostream& operator << (std::ostream& stream, const Array2D& array) {
  // Get word width and bit depth from stream
  // Default word width is size of int, default bit depth fills word width
  const int wordBytes = ioBytes(stream);
  
  // Write wordBytes bytes per array element
  const int size = wordBytes*array.num_elements();
  unsigned char *outBuffer = new unsigned char[size];
  const int height = array.shape()[0];
  const int width = array.shape()[1];
  const int shift = ioShift(stream);
  const int offset = ioZero(stream);
  unsigned int value, byte = 0;
  for (int y=0; y<height; ++y) {
    for (int x=0; x<width; ++x) {
      value = array[y][x]+offset;
      value <<= shift;
      switch (wordBytes) {
        // Only allowed 4 bytes word width in 32 bit systems
        case 4:
          outBuffer[byte++] = (value>>24);
        case 3:
          outBuffer[byte++] = (value>>16);
        case 2:
          outBuffer[byte++] = (value>>8);
        case 1:
          outBuffer[byte++] = value;
          break;
        default:
          throw std::domain_error("Word width of input stream must be in range 1 to 4");
          ;
      }
    }
  }

  //Create reference for output stream buffer (output via stream buffer for efficiency).
  std::streambuf& outbuf = *(stream.rdbuf());
  std::ostream::sentry s(stream);
  if (s) {
    if ( outbuf.sputn(reinterpret_cast<char*>(outBuffer), size) < size )
      stream.setstate(std::ios_base::eofbit|std::ios_base::failbit);
  }
  delete[] outBuffer;

  return stream;
}
