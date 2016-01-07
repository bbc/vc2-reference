/*********************************************************************/
/* Quantisation.h                                                    */
/* Author: Tim Borer,  BBC Research                                  */
/* This version 7th July 2011                                        */
/*                                                                   */
/* Declares stuff related to quantisation                            */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

#ifndef QUANTISATION_14MAY10
#define QUANTISATION_14MAY10

#include "Arrays.h"
#include "Picture.h"

const int adjust_quant_index(const int qIndex, const int qMatrix);

const Array1D adjust_quant_indices(const Array1D& qIndices, const int qMatrix);

const Array2D adjust_quant_indices(const Array2D& qIndices, const int qMatrix);

// Quantise according to parameter q
const int quant(int value, int q);

// Inverse quantise a value according to quantisation index q
const int scale(int value, int q);

// Quantise all the coefficients in a block using
// the same quantiser index
const Array2D quantise_block(const ConstView2D& block, int q);

// Quantise a block of coefficients using an array of quantisers
// The block to be quantised may either be the transform of the whole picture
// or a subband. In the former case this function will quantise slices, in the 
// latter case this function will quantise codeblocks
const Array2D quantise_block(const ConstView2D& block, const Array2D& qIndices);

// Inverse quantise all the coefficients in a block using
// the same quantiser index
const Array2D inverse_quantise_block(const ConstView2D& block, int q);

// Inverse quantise a block of quantised coefficients using an array of quantisers.
// The block to be inverse quantised may correspond either to the transform of the whole picture
// or to a subband. In the former case this function will inverse quantise slices, in the 
// latter case this function will inverse quantise codeblocks
const Array2D inverse_quantise_block(const ConstView2D& block, const Array2D& qIndices);

/***** Predictive Quantisation for Simple, Main and Low Delay Profiles *****/

// Predict LL subband coefficient, at position [y][x],
// from its neighbours above and to the left.
const int predictDC(const Array2D& llSubband, int y, int x);

// Quantise in-place transformed coefficients (using LL subband prediction)
const Array2D quantise_transform(const Array2D& coefficients,
                                 const Array2D& qIndices,
                                 const Array1D& qMatrix);

// Inverse quantise in-place transformed coefficients (using LL subband prediction)
const Array2D inverse_quantise_transform(const Array2D& qCoeffs,
                                         const Array2D& qIndices,
                                         const Array1D& qMatrix);

/***** Non-predictive Quantisation for High Quality Profile *****/

// Quantise in-place transformed coefficients (without LL subband prediction)
const Array2D quantise_transform_np(const Array2D& coefficients,
                                    const int qIndex,
                                    const Array1D& qMatrix);

const Array2D quantise_transform_np(const Array2D& coefficients,
                                    const Array2D& qIndices,
                                    const Array1D& qMatrix);

// Inverse quantise in-place transformed coefficients (without LL subband prediction)
const Array2D inverse_quantise_transform_np(const Array2D& qCoeffs,
                                             const Array2D& qIndices,
                                             const Array1D& qMatrix);

// Quantise in-place transformed coefficients (using LL subband prediction)
const Picture quantise_transform(const Picture& coefficients,
                                 const int qIndex,
                                 const Array1D& qMatrix);

const Picture quantise_transform(const Picture& coefficients,
                                 const Array2D& qIndices,
                                 const Array1D& qMatrix);

// Inverse quantise in-place transformed coefficients (using LL subband prediction)
const Picture inverse_quantise_transform(const Picture& qCoeffs,
                                         const int qIndex,
                                         const Array1D& qMatrix);

const Picture inverse_quantise_transform(const Picture& qCoeffs,
                                         const Array2D& qIndices,
                                         const Array1D& qMatrix);

// Quantise in-place transformed coefficients (without LL subband prediction)
const Picture quantise_transform_np(const Picture& coefficients,
                                    const int qIndex,
                                    const Array1D& qMatrix);

const Picture quantise_transform_np(const Picture& coefficients,
                                    const Array2D& qIndices,
                                    const Array1D& qMatrix);

// Inverse quantise in-place transformed coefficients (without LL subband prediction)
const Picture inverse_quantise_transform_np(const Picture& qCoeffs,
                                            const int qIndex,
                                            const Array1D& qMatrix);

const Picture inverse_quantise_transform_np(const Picture& qCoeffs,
                                            const Array2D& qIndices,
                                            const Array1D& qMatrix);

#endif //QUANTISATION_14MAY10
