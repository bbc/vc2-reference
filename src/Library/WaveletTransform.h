/*********************************************************************/
/* WaveletTransform.h                                                */
/* Author: Tim Borer                                                 */
/* This version 7th July 2011                                        */
/*                                                                   */
/* Declares wavelet transform stuff                                  */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

#ifndef WAVELETTRANSFORM_1MARCH10
#define WAVELETTRANSFORM_1MARCH10

#include <iosfwd>
#include "Arrays.h"
#include "Picture.h"

// Define enumeration for different ypes of wavelet kernel
// Kernels are: Deslauriers-Dubuc (9,7)
//              LeGall (5,3)
//              Deslauriers-Dubuc (13,7)
//              Haar with no shift
//              Haar with single shift per level
//              Fidelity filter
//              Daubechies (9,7) integer approximation
//              NullKernel (does nothing, for test purposes)
enum WaveletKernel {DD97, LeGall, DD137, Haar0, Haar1, Fidelity, Daub97, NullKernel};

std::ostream& operator<<(std::ostream& os, WaveletKernel kernel);

std::istream& operator>>(std::istream& strm, WaveletKernel& kernel);

// Return the size of a padded array given image size and wavelet depth.
const int paddedSize(int size, int depth);

// Check if a wavelet transform is possible with the given wavelet depth
const bool waveletTransformIsPossible(const int waveletDepth, const int lengthLuma, const int lengthChroma);

// Checks if a wavelet transform is possible with the given slice size for the given wavelet depth 
const int sliceSizeIsValid(const int waveletDepth, const int lengthLuma, const int lengthChroma, const int nSize);

// Suggest a wavelet depth that can encode the image
const int suggestWaveletDepth(const int lumaWidth, const int lumaHeight, const int chromaWidth, const int chromaHeight);
const int suggestWaveletDepth(const int lumaWidth, const int lumaHeight, const int chromaWidth, const int chromaHeight, int startingDepth);

// Calculate suggested slice size parameters for the user
const int suggestSliceSize(const int waveletDepth, const int lengthLuma, const int lengthChroma);
const int suggestSliceSize(const int waveletDepth, const int lengthLuma, const int lengthChroma, int startingSliceSize);

//Forward wavelet transform, including padding if necessary
const Array2D waveletTransform(const Array2D& picture, WaveletKernel kernel, int depth);

//Inverse wavelet transform, removes padding if necessary.
// "shape" give size of unpadded image
const Array2D inverseWaveletTransform(const Array2D& transform,
                                      WaveletKernel kernel,
                                      int depth,
                                      Shape2D shape);

// Return the default quantisation matrix for a given wavelet kernel and depth
const Array1D quantMatrix(WaveletKernel kernel, int depth);

// Convert a Array2D containing an in place wavelet transform into a 1D array of subbands
const BlockVector split_into_subbands(const Array2D& picture, const char waveletDepth);

// Converts a 1D array of subbands back to a single single Array2D corresponding
// to an in-place wavelet transform
const Array2D merge_subbands(const BlockVector& subbands);

//Forward wavelet transform, including padding if necessary
const Picture waveletTransform(const Picture& picture, enum WaveletKernel kernel, int depth);

//Inverse wavelet transform, removes padding if necessary.
// "format" specifies format of unpadded image
const Picture inverseWaveletTransform(const Picture& transform,
                                      enum WaveletKernel kernel,
                                      int depth,
                                      PictureFormat format);

#endif //WAVELETTRANSFORM_1MARCH10
