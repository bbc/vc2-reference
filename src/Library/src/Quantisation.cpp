/*********************************************************************/
/* Quantisation.cpp                                                  */
/* Author: Tim Borer,  BBC Research                                  */
/* This version 7th July 2011                                        */
/*                                                                   */
/* Defines stuff related to quantisation                             */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

#include "Quantisation.h"
#include "WaveletTransform.h"
#include "Utils.h"

using utils::pow;

const int adjust_quant_index(const int qIndex, const int qMatrix) {
  int aQIndex = qIndex-qMatrix;
  if (aQIndex<0) return 0;
  return aQIndex;
}

const Array1D adjust_quant_indices(const Array1D& qIndices, const int qMatrix) {
  Array1D aQIndices(qIndices.ranges());
  // Adjust all the quantisers in qIndices
  std::transform(qIndices.data(), qIndices.data()+qIndices.num_elements(),
                 aQIndices.data(),
                 std::bind2nd(std::ptr_fun(adjust_quant_index), qMatrix) );
  return aQIndices;
}

const Array2D adjust_quant_indices(const Array2D& qIndices, const int qMatrix) {
  Array2D aQIndices(qIndices.ranges());
  // Adjust all the quantisers in qIndices
  std::transform(qIndices.data(), qIndices.data()+qIndices.num_elements(),
                 aQIndices.data(),
                 std::bind2nd(std::ptr_fun(adjust_quant_index), qMatrix) );
  return aQIndices;
}

const int quant_factor(int q) {
  // Lookup table valid q<120. For q>= 120 quant_factor(q) requires more than 32 bits.
  static const unsigned int lookup[120] = {
    0x000000004, 0x000000005, 0x000000006, 0x000000007, 0x000000008, 0x00000000A, 0x00000000B, 0x00000000D,
    0x000000010, 0x000000013, 0x000000017, 0x00000001B, 0x000000020, 0x000000026, 0x00000002D, 0x000000036,
    0x000000040, 0x00000004C, 0x00000005B, 0x00000006C, 0x000000080, 0x000000098, 0x0000000B5, 0x0000000D7,
    0x000000100, 0x000000130, 0x00000016A, 0x0000001AF, 0x000000200, 0x000000261, 0x0000002D4, 0x00000035D,
    0x000000400, 0x0000004C2, 0x0000005A8, 0x0000006BA, 0x000000800, 0x000000983, 0x000000B50, 0x000000D74,
    0x000001000, 0x000001307, 0x0000016A1, 0x000001AE9, 0x000002000, 0x00000260E, 0x000002D41, 0x0000035D1,
    0x000004000, 0x000004C1C, 0x000005A82, 0x000006BA2, 0x000008000, 0x000009838, 0x00000B505, 0x00000D745,
    0x000010000, 0x000013070, 0x000016A0A, 0x00001AE8A, 0x000020000, 0x0000260E0, 0x00002D414, 0x000035D14,
    0x000040000, 0x00004C1C0, 0x00005A828, 0x00006BA28, 0x000080000, 0x00009837F, 0x0000B504F, 0x0000D7450,
    0x000100000, 0x0001306FE, 0x00016A09E, 0x0001AE8A0, 0x000200000, 0x000260DFC, 0x0002D413D, 0x00035D13F,
    0x000400000, 0x0004C1BF8, 0x0005A827A, 0x0006BA27E, 0x000800000, 0x0009837F0, 0x000B504F3, 0x000D744FD,
    0x001000000, 0x001306FE1, 0x0016A09E6, 0x001AE89FA, 0x002000000, 0x00260DFC1, 0x002D413CD, 0x0035D13F3,
    0x004000000, 0x004C1BF83, 0x005A8279A, 0x006BA27E6, 0x008000000, 0x009837F05, 0x00B504F33, 0x00D744FCD,
    0x010000000, 0x01306FE0A, 0x016A09E66, 0x01AE89F99, 0x020000000, 0x0260DFC14, 0x02D413CCD, 0x035D13F33,
    0x040000000, 0x04C1BF829, 0x05A82799A, 0x06BA27E65, 0x080000000, 0x09837F052, 0x0B504F334, 0x0D744FCCB,
//  0x100000000, 0x1306FE0A3, 0x16A09E668, 0x1AE89F996, 0x200000000, 0x260DFC146, 0x2D413CCD0, 0x35D13F32B
  };
  if (q<0) q=0;
  return static_cast<int>(lookup[q]);
}

// Quantise according to parameter q
const int quant(int value, int q) {
  bool negative = (value<0);
  if (negative) value *= -1;
  value <<= 2;
  value /= quant_factor(q);
  if (negative) value *= -1;
  return value;
}

const int quant_offset(int q) {
  if (q<0) q=0;
  if (q==0) return 1;
  else if (q==1) return 2;
  else return (quant_factor(q)+1)/2;
}

// Inverse quantise a value according to quantisation index q
const int scale(int value, int q) {
  bool negative = (value<0);
  if (negative) value *= -1;
  value *= quant_factor(q);
  if ( value>0 ) value += quant_offset(q);
  value += 2; //This addition could be included in quant_offset
  value /= 4;
  if (negative) value *= -1;
  return value;
}

// Quantise all the coefficients in a block using
// the same quantiser index
const Array2D quantise_block(const ConstView2D& block, int q) {
  // Construct a new array with same size as block
  Array2D quantisedBlock(block.ranges());
  const int blockHeight = block.shape()[0];
  const int blockWidth = block.shape()[1];
  for (int y=0; y<blockHeight; ++y) {
    for (int x=0; x<blockWidth; ++x) {
      quantisedBlock[y][x] = quant(block[y][x], q);
    }
  }
  return quantisedBlock;
}

// Quantise a block of coefficients using an array of quantisers
// The block to be quantised may either be the transform of the whole picture
// or a subband. In the former case this function will quantise slices, in the 
// latter case this function will quantise codeblocks
const Array2D quantise_block(const ConstView2D& block,
                             const Array2D& qIndices) {
  // Construct a new array with same size as block
  Array2D quantisedBlock(block.ranges());
  const int blockHeight = block.shape()[0];
  const int blockWidth = block.shape()[1];
  const int yBlocks = qIndices.shape()[0];
  const int xBlocks = qIndices.shape()[1];
  // Loop through the slices or codeblocks
  // Note Range(left, right) defines the half open range [left, right),
  // i.e. the rightmost element is not included
  for (int y=0, top=0, bottom=blockHeight/yBlocks;
       y<yBlocks;
       ++y, top=bottom, bottom=((y+1)*blockHeight/yBlocks) ) {
    for (int x=0, left=0, right=blockWidth/xBlocks;
         x<xBlocks;
         ++x, left=right, right=((x+1)*blockWidth/xBlocks) ) {
           const ArrayIndices2D sliceIndices = // Define the samples within the current slice/codeblock
             indices[Range(top,bottom)][Range(left,right)];
           quantisedBlock[sliceIndices] =
             quantise_block(block[sliceIndices], qIndices[y][x]);
    }
  }
  return quantisedBlock;
}

// Inverse quantise a block of quantised coefficients using an array of quantisers.
// The block to be inverse quantised may correspond either to the transform of the whole picture
// or to a subband. In the former case this function will inverse quantise slices, in the 
// latter case this function will inverse quantise codeblocks
const Array2D inverse_quantise_block(const ConstView2D& block,
                                     const Array2D& qIndices) {
  // Construct a new array with same size as block
  Array2D invQuantisedBlock(block.ranges());
  const int blockHeight = block.shape()[0];
  const int blockWidth = block.shape()[1];
  const int yBlocks = qIndices.shape()[0];
  const int xBlocks = qIndices.shape()[1];
  // Loop through the slices or codeblocks
  // Note Range(left, right) defines the half open range [left, right),
  // i.e. the rightmost element is not included
  for (int y=0, top=0, bottom=blockHeight/yBlocks;
       y<yBlocks;
       ++y, top=bottom, bottom=((y+1)*blockHeight/yBlocks) ) {
    for (int x=0, left=0, right=blockWidth/xBlocks;
         x<xBlocks;
         ++x, left=right, right=((x+1)*blockWidth/xBlocks) ) {
           const ArrayIndices2D sliceIndices = // Define the samples within the curent slice/codeblock
             indices[Range(top,bottom)][Range(left,right)];
           invQuantisedBlock[sliceIndices] =
             inverse_quantise_block(block[sliceIndices], qIndices[y][x]);
    }
  }
  return invQuantisedBlock;
}

// Inverse quantise all the coefficients in a block using
// the same quantiser index
const Array2D inverse_quantise_block(const ConstView2D& block, int q) {
  // Construct a new array with same size as block
  Array2D invQuantisedBlock(block.ranges());
  const int blockHeight = block.shape()[0];
  const int blockWidth = block.shape()[1];
  for (int y=0; y<blockHeight; ++y) {
    for (int x=0; x<blockWidth; ++x) {
      invQuantisedBlock[y][x] = scale(block[y][x], q);
    }
  }
  return invQuantisedBlock;
}

/***** Predictive Quantisation for Simple, Main and Low Delay Profiles *****/

// Predict LL subband coefficient, at position [y][x],
// from its neighbours above and to the left.
const int predictDC(const Array2D& llSubband, int y, int x) {
  if (y>0 && x>0) {
    int result = llSubband[y-1][x-1];
    result += llSubband[y-1][x];
    result += llSubband[y][x-1];
    if (result>=0) return (result+1)/3;
    else return (result-1)/3;
  }
  else if (y>0) {
    return llSubband[y-1][x];
  }
  else if (x>0) {
    return llSubband[y][x-1];
  }
  else {
    return 0;
  }
}

// Quantise an LL (DC) subband, including prediction
// This version either quantises the LL subband for low delay mode or
// codes the LL subband for core syntax using codeblocks
const Array2D quantise_LLSubband(const ConstView2D& llSubband,
                                 const Array2D& qIndices) {
  const int LLHeight = llSubband.shape()[0]; // Height of the LL subband
  const int LLWidth = llSubband.shape()[1]; // Width of the LL subband
  const int yBlocks = qIndices.shape()[0]; // Number of vertical slices/codeblocks in the LL subband
  const int xBlocks = qIndices.shape()[1]; // Number of horizontal slices/codeblocks in the LL subband
  Array2D quantisedLL(llSubband.ranges());
  Array2D restoredLL(llSubband.ranges());
  for (int y=0; y<LLHeight; ++y) {
    for (int x=0; x<LLWidth; ++x) {
      // TO DO: Implement more efficient calculation of yb/xb.
      // Calculate (y+1)*yBlocks by incrementing previous version by yBlocks
      // Do division by shift (width is always a power of 2)
      const int yb = ((y+1)*yBlocks-1)/LLHeight; // vertical slice/codeblock number
      const int xb = ((x+1)*xBlocks-1)/LLWidth; // horizontal slice/codeblock number
      const int prediction = predictDC(restoredLL, y, x);
      quantisedLL[y][x] = quant(llSubband[y][x]-prediction, qIndices[yb][xb]);
      restoredLL[y][x] = scale(quantisedLL[y][x], qIndices[yb][xb])+prediction;
    }
  }
  return quantisedLL;
}

// Quantise a subband in in-place transform order
// This version of quantise_subbands assumes multiple quantisers per subband.
// It may be used for either quantising slices or for quantising subbands with codeblocks
const Array2D quantise_subbands(const Array2D& coefficients, const BlockVector& qIndices) {
  const Index transformHeight = coefficients.shape()[0];
  const Index transformWidth = coefficients.shape()[1];
  // TO DO: Check numberOfSubbands=3n+1 ?
  const int numberOfSubbands = qIndices.size();
  const int waveletDepth = (numberOfSubbands-1)/3;
  Index stride, offset; // stride is subsampling factor, offset is subsampling phase
  Array2D result(coefficients.ranges());

  // Create a view of the coefficients, representing the LL subband, quantise it,
  // then assign the result a view of the results array. This puts the quantised
  // LL subband into the result array in in-place transform order.
  // ArrayIndices2D objects specify the subset of array elements within a view,
  // that is they specify the subsampling factor and subsampling phase.
  stride = pow(2, waveletDepth);
  const ArrayIndices2D LLindices = // LLindices specifies the samples in the LL subband
    indices[Range(0,transformHeight,stride)][Range(0,transformWidth,stride)];
  result[LLindices] =
    quantise_LLSubband(coefficients[LLindices], qIndices[0]);

  // Next quantise the other subbands
  // Note: Level numbers go from zero for the lowest ("DC") frequencies to depth for
  // the high frequencies. This corresponds to the convention in the VC-2 specification.
  // Subands go from zero ("DC") to numberOfSubbands-1 for HH at the highest level
  for (char level=1, band=1; level<=waveletDepth; ++level) {
    stride = pow(2, waveletDepth+1-level);
    offset = stride/2;
    // Create a view of coefficients corresponding to a subband, then quantise it
    //Quantise HL subband
    const ArrayIndices2D HLindices = // HLindices specifies the samples in the HL subband
      indices[Range(0,transformHeight,stride)][Range(offset,transformWidth,stride)];
    result[HLindices] = quantise_block(coefficients[HLindices], qIndices[band++]);
    //Quantise LH subband
    const ArrayIndices2D LHindices = // LHindices specifies the samples in the LH subband
      indices[Range(offset,transformHeight,stride)][Range(0,transformWidth,stride)];
    result[LHindices] = quantise_block(coefficients[LHindices], qIndices[band++]);
    //Quantise HH subband
    const ArrayIndices2D HHindices = // HHindices specifies the samples in the HH subband
      indices[Range(offset,transformHeight,stride)][Range(offset,transformWidth,stride)];
    result[HHindices] = quantise_block(coefficients[HHindices], qIndices[band++]);
  }

  return result;
}

// Inverse quantise an LL (DC) subband, including prediction
// This version either inverse quantises the LL subband for low delay mode or
// inverse quantises the LL subband for core syntax using codeblocks
const Array2D inverse_quantise_LLSubband(const ConstView2D& llSubband,
                                         const Array2D& qIndices) {
  const int LLHeight = llSubband.shape()[0]; // Height of the LL subband
  const int LLWidth = llSubband.shape()[1]; // Width of the LL subband
  const int yBlocks = qIndices.shape()[0]; // Number of vertical slices/codeblocks in the LL subband
  const int xBlocks = qIndices.shape()[1]; // Number of horizontal slices/codeblocks in the LL subband
  Array2D invQuantisedLL(llSubband.ranges());
  for (int y=0; y<LLHeight; ++y) {
    for (int x=0; x<LLWidth; ++x) {
      // TO DO: Implement more efficient calculation of yb/xb.
      // Calculate (y+1)*yBlocks by incrementing previous version by yBlocks
      // Do division by shift (width is always a power of 2)
      const int yb = ((y+1)*yBlocks-1)/LLHeight; // vertical slice/codeblock number
      const int xb = ((x+1)*xBlocks-1)/LLWidth; // horizontal slice/codeblock number
      const int prediction = predictDC(invQuantisedLL, y, x);
      invQuantisedLL[y][x] = scale(llSubband[y][x], qIndices[yb][xb])+prediction;
    }
  }
  return invQuantisedLL;
}

// Inverse quantise a subband in in-place transform order
// This version of inverse_quantise_subbands assumes mulitple quantisers per subband.
// It may be used for either inverse quantising slices or for inverse quantising subbands with codeblocks
const Array2D inverse_quantise_subbands(const Array2D& coefficients, const BlockVector& qIndices) {
  const Index transformHeight = coefficients.shape()[0];
  const Index transformWidth = coefficients.shape()[1];
  // TO DO: Check numberOfSubbands=3n+1 ?
  const int numberOfSubbands = qIndices.size();
  const int waveletDepth = (numberOfSubbands-1)/3;
  Index stride, offset; // stride is subsampling factor, offset is subsampling phase
  Array2D result(coefficients.ranges());

  // Create a view of the coefficients, representing the LL subband, quantise it,
  // then assign the result a view of the results array. This puts the quantised
  // LL subband into the result array in in-place transform order.
  // ArrayIndices2D objects specify the subset of array elements within a view,
  // that is they specify the subsampling factor and subsampling phase.
  stride = pow(2, waveletDepth);
  const ArrayIndices2D LLindices = // LLindices specifies the samples in the LL subband
    indices[Range(0,transformHeight,stride)][Range(0,transformWidth,stride)];
  result[LLindices] = inverse_quantise_LLSubband(coefficients[LLindices], qIndices[0]);

  // Next quantise the other subbands
  // Note: Level numbers go from zero for the lowest ("DC") frequencies to depth for
  // the high frequencies. This corresponds to the convention in the VC-2 specification.
  // Subands go from zero ("DC") to numberOfSubbands-1 for HH at the highest level
  for (char level=1, band=1; level<=waveletDepth; ++level) {
    stride = pow(2, waveletDepth+1-level);
    offset = stride/2;
    // Create a view of coefficients corresponding to a subband, then quantise it
    //Quantise HL subband
    const ArrayIndices2D HLindices = // HLindices specifies the samples in the HL subband
      indices[Range(0,transformHeight,stride)][Range(offset,transformWidth,stride)];
    result[HLindices] = inverse_quantise_block(coefficients[HLindices], qIndices[band++]);
    //Quantise LH subband
    const ArrayIndices2D LHindices = // LHindices specifies the samples in the LH subband
      indices[Range(offset,transformHeight,stride)][Range(0,transformWidth,stride)];
    result[LHindices] = inverse_quantise_block(coefficients[LHindices], qIndices[band++]);
    //Quantise HH subband
    const ArrayIndices2D HHindices = // HHindices specifies the samples in the HH subband
      indices[Range(offset,transformHeight,stride)][Range(offset,transformWidth,stride)];
    result[HHindices] = inverse_quantise_block(coefficients[HHindices], qIndices[band++]);
  }

  return result;
}

// Quantise in-place transformed coefficients of a whole picture as slices
// Uses a quantisation matrix
const Array2D quantise_transform(const Array2D& coefficients,
                                 const Array2D& qIndices,
                                 const Array1D& qMatrix) {
  // TO DO: Check numberOfSubbands=3n+1 ?
  BlockVector aQIndices(qMatrix.ranges());
  const int numberOfSubbands = qMatrix.size();
  for (char band=0; band<numberOfSubbands; ++band) {
    aQIndices[band] = adjust_quant_indices(qIndices, qMatrix[band]);
  }
  return quantise_subbands(coefficients, aQIndices);
}

const Array2D inverse_quantise_transform(const Array2D& qCoeffs,
                                         const Array2D& qIndices,
                                         const Array1D& qMatrix) {
  // TO DO: Check numberOfSubbands=3n+1 ?
  BlockVector aQIndices(qMatrix.ranges());
  const int numberOfSubbands = qMatrix.size();
  for (char band=0; band<numberOfSubbands; ++band) {
    aQIndices[band] = adjust_quant_indices(qIndices, qMatrix[band]);
  }
  return inverse_quantise_subbands(qCoeffs, aQIndices);
}

/***** Non-predictive Quantisation for High Quality Profile *****/

// Quantise a subband in in-place transform order (without LL subband prediction)
// This version of quantise_subbands assumes multiple quantisers per subband.
// It may be used for either quantising slices or for quantising subbands with codeblocks
const Array2D quantise_subbands_np(const Array2D& coefficients, const BlockVector& qIndices) {
  const Index transformHeight = coefficients.shape()[0];
  const Index transformWidth = coefficients.shape()[1];
  // TO DO: Check numberOfSubbands=3n+1 ?
  const int numberOfSubbands = qIndices.size();
  const int waveletDepth = (numberOfSubbands-1)/3;
  Index stride, offset; // stride is subsampling factor, offset is subsampling phase
  Array2D result(coefficients.ranges());

  // Create a view of the coefficients, representing the LL subband, quantise it,
  // then assign the result a view of the results array. This puts the quantised
  // LL subband into the result array in in-place transform order.
  // ArrayIndices2D objects specify the subset of array elements within a view,
  // that is they specify the subsampling factor and subsampling phase.
  stride = pow(2, waveletDepth);
  const ArrayIndices2D LLindices = // LLindices specifies the samples in the LL subband
    indices[Range(0,transformHeight,stride)][Range(0,transformWidth,stride)];
  result[LLindices] = quantise_block(coefficients[LLindices], qIndices[0]);

  // Next quantise the other subbands
  // Note: Level numbers go from zero for the lowest ("DC") frequencies to depth for
  // the high frequencies. This corresponds to the convention in the VC-2 specification.
  // Subands go from zero ("DC") to numberOfSubbands-1 for HH at the highest level
  for (char level=1, band=1; level<=waveletDepth; ++level) {
    stride = pow(2, waveletDepth+1-level);
    offset = stride/2;
    // Create a view of coefficients corresponding to a subband, then quantise it
    //Quantise HL subband
    const ArrayIndices2D HLindices = // HLindices specifies the samples in the LL subband
      indices[Range(0,transformHeight,stride)][Range(offset,transformWidth,stride)];
    result[HLindices] = quantise_block(coefficients[HLindices], qIndices[band++]);
    //Quantise LH subband
    const ArrayIndices2D LHindices = // LHindices specifies the samples in the LL subband
      indices[Range(offset,transformHeight,stride)][Range(0,transformWidth,stride)];
    result[LHindices] = quantise_block(coefficients[LHindices], qIndices[band++]);
    //Quantise HH subband
    const ArrayIndices2D HHindices = // HHindices specifies the samples in the LL subband
      indices[Range(offset,transformHeight,stride)][Range(offset,transformWidth,stride)];
    result[HHindices] = quantise_block(coefficients[HHindices], qIndices[band++]);
  }

  return result;
}

// Inverse quantise a subband in in-place transform order (without LL subband prediction)
// This version of inverse_quantise_subbands assumes mulitple quantisers per subband.
// It may be used for either inverse quantising slices or for inverse quantising subbands with codeblocks
const Array2D inverse_quantise_subbands_np(const Array2D& coefficients, const BlockVector& qIndices) {
  const Index transformHeight = coefficients.shape()[0];
  const Index transformWidth = coefficients.shape()[1];
  // TO DO: Check numberOfSubbands=3n+1 ?
  const int numberOfSubbands = qIndices.size();
  const int waveletDepth = (numberOfSubbands-1)/3;
  Index stride, offset; // stride is subsampling factor, offset is subsampling phase
  Array2D result(coefficients.ranges());

  // Create a view of the coefficients, representing the LL subband, quantise it,
  // then assign the result a view of the results array. This puts the quantised
  // LL subband into the result array in in-place transform order.
  // ArrayIndices2D objects specify the subset of array elements within a view,
  // that is they specify the subsampling factor and subsampling phase.
  stride = pow(2, waveletDepth);
  const ArrayIndices2D LLindices = // LLindices specifies the samples in the LL subband
    indices[Range(0,transformHeight,stride)][Range(0,transformWidth,stride)];
  result[LLindices] = inverse_quantise_block(coefficients[LLindices], qIndices[0]);

  // Next quantise the other subbands
  // Note: Level numbers go from zero for the lowest ("DC") frequencies to depth for
  // the high frequencies. This corresponds to the convention in the VC-2 specification.
  // Subands go from zero ("DC") to numberOfSubbands-1 for HH at the highest level
  for (char level=1, band=1; level<=waveletDepth; ++level) {
    stride = pow(2, waveletDepth+1-level);
    offset = stride/2;
    // Create a view of coefficients corresponding to a subband, then quantise it
    //Quantise HL subband
    const ArrayIndices2D HLindices = // HLindices specifies the samples in the LL subband
      indices[Range(0,transformHeight,stride)][Range(offset,transformWidth,stride)];
    result[HLindices] = inverse_quantise_block(coefficients[HLindices], qIndices[band++]);
    //Quantise LH subband
    const ArrayIndices2D LHindices = // LHindices specifies the samples in the LL subband
      indices[Range(offset,transformHeight,stride)][Range(0,transformWidth,stride)];
    result[LHindices] = inverse_quantise_block(coefficients[LHindices], qIndices[band++]);
    //Quantise HH subband
    const ArrayIndices2D HHindices = // HHindices specifies the samples in the LL subband
      indices[Range(offset,transformHeight,stride)][Range(offset,transformWidth,stride)];
    result[HHindices] = inverse_quantise_block(coefficients[HHindices], qIndices[band++]);
  }

  return result;
}

// Quantise in-place transformed coefficients of a whole picture as slices
// Uses a quantisation matrix
const Array2D quantise_transform_np(const Array2D& coefficients,
                                    const Array2D& qIndices,
                                    const Array1D& qMatrix) {
  // TO DO: Check numberOfSubbands=3n+1 ?
  BlockVector aQIndices(qMatrix.ranges());
  const int numberOfSubbands = qMatrix.size();
  for (char band=0; band<numberOfSubbands; ++band) {
    aQIndices[band] = adjust_quant_indices(qIndices, qMatrix[band]);
  }
  return quantise_subbands_np(coefficients, aQIndices);
}

// Quantise all the coefficients in a block using
// the same quantiser index
const Array2D quantise_block(const Array2D& block, int q) {
  // Construct a new array with same size as block
  Array2D quantisedBlock(block.ranges());
  const int blockHeight = block.shape()[0];
  const int blockWidth = block.shape()[1];
  for (int y=0; y<blockHeight; ++y) {
    for (int x=0; x<blockWidth; ++x) {
      quantisedBlock[y][x] = quant(block[y][x], q);
    }
  }
  return quantisedBlock;
}

// Quantise in-place transformed coefficients (without LL subband prediction)
const Array2D quantise_transform_np(const Array2D& coefficients,
                                    const int qIndex,
                                    const Array1D& qMatrix) {
  // TO DO: Check numberOfSubbands=3n+1 ?
  const int numberOfSubbands = qMatrix.size();
  const int waveletDepth = (numberOfSubbands-1)/3;
  BlockVector subbands = split_into_subbands(coefficients, waveletDepth);
  for (int band=0; band<numberOfSubbands; ++band) {
    const int aQIndex = adjust_quant_index(qIndex, qMatrix[band]);
    subbands[band] = quantise_block(subbands[band], aQIndex);
  }
  return merge_subbands(subbands);
}

const Array2D inverse_quantise_transform_np(const Array2D& qCoeffs,
                                            const Array2D& qIndices,
                                            const Array1D& qMatrix) {
  // TO DO: Check numberOfSubbands=3n+1 ?
  BlockVector aQIndices(qMatrix.ranges());
  const int numberOfSubbands = qMatrix.size();
  for (char band=0; band<numberOfSubbands; ++band) {
    aQIndices[band] = adjust_quant_indices(qIndices, qMatrix[band]);
  }
  return inverse_quantise_subbands_np(qCoeffs, aQIndices);
}

// Quantise in-place transformed coefficients of a whole picture as slices
// Using LL (DC) subband prediction
// Uses a quantisation matrix
const Picture quantise_transform(const Picture& transform,
                                 const Array2D& qIndices,
                                 const Array1D& qMatrix) {
  Picture result(transform.format());
  result.y(quantise_transform(transform.y(), qIndices, qMatrix));
  result.c1(quantise_transform(transform.c1(), qIndices, qMatrix));
  result.c2(quantise_transform(transform.c2(), qIndices, qMatrix));
  return result;
}

const Picture inverse_quantise_transform(const Picture& qCoeffs,
                                         const Array2D& qIndices,
                                         const Array1D& qMatrix) {
  Picture result(qCoeffs.format());
  result.y(inverse_quantise_transform(qCoeffs.y(), qIndices, qMatrix));
  result.c1(inverse_quantise_transform(qCoeffs.c1(), qIndices, qMatrix));
  result.c2(inverse_quantise_transform(qCoeffs.c2(), qIndices, qMatrix));
  return result;
}

// Quantise in-place transformed coefficients of a whole picture as slices
// Without LL (DC) subband prediction
// Uses a quantisation matrix
const Picture quantise_transform_np(const Picture& transform,
                                    const Array2D& qIndices,
                                    const Array1D& qMatrix) {
  Picture result(transform.format());
  result.y(quantise_transform_np(transform.y(), qIndices, qMatrix));
  result.c1(quantise_transform_np(transform.c1(), qIndices, qMatrix));
  result.c2(quantise_transform_np(transform.c2(), qIndices, qMatrix));
  return result;
}

// Quantise in-place transformed coefficients (without LL subband prediction)
const Picture quantise_transform_np(const Picture& transform,
                                    const int qIndex,
                                    const Array1D& qMatrix) {
  Picture result(transform.format());
  result.y(quantise_transform_np(transform.y(), qIndex, qMatrix));
  result.c1(quantise_transform_np(transform.c1(), qIndex, qMatrix));
  result.c2(quantise_transform_np(transform.c2(), qIndex, qMatrix));
  return result;
}

const Picture inverse_quantise_transform_np(const Picture& qCoeffs,
                                            const Array2D& qIndices,
                                            const Array1D& qMatrix) {
  Picture result(qCoeffs.format());
  result.y(inverse_quantise_transform_np(qCoeffs.y(), qIndices, qMatrix));
  result.c1(inverse_quantise_transform_np(qCoeffs.c1(), qIndices, qMatrix));
  result.c2(inverse_quantise_transform_np(qCoeffs.c2(), qIndices, qMatrix));
  return result;
}
