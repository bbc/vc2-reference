/*********************************************************************/
/* Arrays.h                                                          */
/* Author: Tim Borer                                                 */
/* This version 11th July 2011                                       */
/*                                                                   */
/* Declares stuff related to mulitdimensional arrays.                */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

#ifndef ARRAYS_25FEB10
#define ARRAYS_25FEB10

#include <iosfwd>
#include "boost/multi_array.hpp"
#include "boost/array.hpp"

using boost::extents; // An extents generator object used to define array sizes
using boost::indices; // An index generator object used to define array views

typedef boost::multi_array_types::index_range Range; // Range(bottom, top+1, stride)

typedef boost::multi_array_types::index Index;

// Shape2D is a 1D Array used to store the size of a picture
typedef boost::array<Index, 2> Shape2D;

// Array1D is (surprise) a one dimensional array or vector
typedef boost::multi_array<int, 1> Array1D;

// Array2D is a 2D array for holding picture or coefficient samples
typedef boost::multi_array<int, 2> Array2D;

// View2D is a view of a Array2D,
// that is a synonym for a subset of array values
typedef Array2D::array_view<2>::type View2D;

// ConstView2D is a constant view of a Array2D,
typedef Array2D::const_array_view<2>::type ConstView2D;

// BlockVector is a 1D Array (i.e. vector) of 2D arrays (each element is, itself, a little 2D array)
// A BlockVector may be used for storing an array of wavelet transform subbands (each subband is,
// itself, a 2D array of coefficients).
typedef boost::multi_array<Array2D, 1> BlockVector;

// BlockArray is a 2D array of 2D arrays (each block is a little 2D array)
typedef boost::multi_array<Array2D, 2> BlockArray;

// ArrayIndices2D specifies the subset of array elements which define a view,
// that is it specifies the subsampling factor and subsampling phase.
typedef boost::detail::multi_array::index_gen<2, 2> ArrayIndices2D;

// Get the shape of a 2D array
const Shape2D shape(const Array2D&);
const Shape2D shape(const View2D&);
const Shape2D shape(const ConstView2D&);
const Shape2D shape(const BlockArray&);

// Splits a large 2D array into an array of smaller 2D arrays (blocks)
// Note that if the number of blocks is not a sub-multiple of the input array dimensions then
// the blocks will have different sizes!
// yBlocks and xBlocks are the number of blocks in the vertical and horizontal dimension respectively.
// Splits a picture into slices or a subband into codeblocks.
const BlockArray split_into_blocks(const Array2D& picture, int yBlocks, int xBlocks);

// Converts an array of blocks back to a single 2D array
// This is the inverse of "split_into_blocks()"
// The array of blocks might represent slices or codeblocks
const Array2D merge_blocks(const BlockArray& blocks);

// Clip an array to specified limits
const Array2D clip(const Array2D& values, const int min_value, const int max_value);

//**************** Array IO declarations ****************//

namespace arrayio {

  enum ioFormat {UNSIGNED, SIGNED, OFFSET};

  // Format manipulator - sets io format
  class format {
    public:
      format(ioFormat f): iof(f) {}; 
      void operator () (std::ios_base& stream) const;
    private:
      const ioFormat iof;
  };

  // Format manipulator - sets number of bytes per element
  class wordWidth {
    public:
      wordWidth(int b): no_of_bytes(b) {}; 
      void operator () (std::ios_base& stream) const;
    private:
      const int no_of_bytes;
  };

  // Format manipulator - sets number of bits used for data
  class bitDepth {
    public:
      bitDepth(int b): no_of_bits(b) {}; 
      void operator () (std::ios_base& stream) const;
    private:
      const int no_of_bits;
  };

  // Format manipulator - sets the offset for offset binary
  // (Default offset, if not explicitly set, is half range)
  class offset {
    public:
      offset(int z): zeroValue(z) {}; 
      void operator () (std::ios_base& stream) const;
    private:
      const int zeroValue;
  };

  // Set data to be left justified within the data word
  std::ostream& left_justified(std::ostream& stream);

  // Set data to be left justified within the data word
  std::istream& left_justified(std::istream& stream);

  // Set data to be right justified within the data word
  std::ostream& right_justified(std::ostream& stream);

  // Set data to be right justified within the data word
  std::istream& right_justified(std::istream& stream);

  // Set data format to offset binary
  std::ostream& offset_binary(std::ostream& stream);

  // Set data format to offset binary
  std::istream& offset_binary(std::istream& stream);

  // Set data format to signed binary (two's complement)
  std::ostream& signed_binary(std::ostream& stream);

  // Set data format to signed binary (two's complement)
  std::istream& signed_binary(std::istream& stream);

  // Set data format to unsigned binary
  std::ostream& unsigned_binary(std::ostream& stream);

  // Set data format to unsigned binary
  std::istream& unsigned_binary(std::istream& stream);

  // Set data format to text
  std::ostream& text(std::ostream& stream);

  // Set data format to text
  std::istream& text(std::istream& stream);

} // end namespace arrayio

// ostream io format manipulator
std::ostream& operator << (std::ostream& stream, arrayio::format f);

// istream io format manipulator
std::istream& operator >> (std::istream& stream, arrayio::format f);

// ostream wordWidth format manipulator
std::ostream& operator << (std::ostream& stream, arrayio::wordWidth w);

// istream wordWidth format manipulator
std::istream& operator >> (std::istream& stream, arrayio::wordWidth w);

// ostream bit depth format manipulator
std::ostream& operator << (std::ostream& stream, arrayio::bitDepth d);

// istream bit depth format manipulator
std::istream& operator >> (std::istream& stream, arrayio::bitDepth d);

// ostream offset format manipulator
std::ostream& operator << (std::ostream& stream, arrayio::offset z);

// istream offset format manipulator
std::istream& operator >> (std::istream& stream, arrayio::offset z);

std::istream& operator >> (std::istream& stream, Array2D& array);

std::ostream& operator << (std::ostream& stream, const Array2D& array);
  
#endif //ARRAYS_25FEB10
