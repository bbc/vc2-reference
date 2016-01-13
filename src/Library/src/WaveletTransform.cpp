/************************************************************************/
/* WaveletTransform.cpp                                                 */
/* Author: Tim Borer                                                    */
/* This version 2nd November 2011                                       */
/*                                                                      */
/* Defines wavelet transform stuff                                      */
/* 8th May 2011: Added transform of Picture                             */
/* 15th September 2011: quantMatrix now calculates values for any depth */
/* 2nd November 2011: added <cfloat> header                             */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file      */
/************************************************************************/

#include "WaveletTransform.h"

#include <iostream>
#include <string>
#include <stdexcept> // For invalid_argument
#include <cfloat> // For FLT_MAX in quantMatrix

std::ostream& operator<<(std::ostream& os, WaveletKernel kernel) {
  const char* s;
  switch (kernel) {
    case DD97:
      s = "Deslauriers-Dubuc (9,7) (\"DD97\")";
      break;
    case LeGall:
      s = "LeGall (5,3) (\"LeGall\")";
      break;
    case DD137:
      s = "Deslauriers-Dubuc (13,7) (\"DD137\")";
      break;
    case Haar0:
      s = "Haar (no shift) (\"Haar0\")";
      break;
    case Haar1:
      s = "Haar (one bit shift) (\"Haar1\")";
      break;
    case Fidelity:
      s = "Fidelity (\"Fidelity\")";
      break;
    case Daub97:
      s = "Daubechies (9,7) (\"Daub97\")";
      break;
    case NullKernel:
      s = "NullKernel (\"NullKernel\")";
      break;
    default:
      s = "Unknown wavelet kernel!";
      break;
  }
  return os<<s;
}

std::istream& operator>>(std::istream& strm, WaveletKernel& kernel) {
        std::string text;
        strm >> text;
        if (text == "DD97") kernel = DD97;
        else if (text == "LeGall") kernel = LeGall;
        else if (text == "DD137") kernel = DD137;
        else if (text == "Haar0") kernel = Haar0;
        else if (text == "Haar1") kernel = Haar1;
        else if (text == "Fidelity") kernel = Fidelity;
        else if (text == "Daub97") kernel = Daub97;
        else if (text == "NullKernel") kernel = NullKernel;
        else throw std::invalid_argument("invalid wavelet kernel");
        return strm;
}

#include "Utils.h"

const int paddedSize(int size, int depth) {
  const int cell = utils::pow(2, depth);
  return cell*((size+cell-1)/cell);
}

const Array2D waveletPad(const Array2D& picture, int depth) {
  const Index pictureHeight = picture.shape()[0];
  const Index pictureWidth = picture.shape()[1];
  const Index paddedHeight = paddedSize(pictureHeight, depth);
  const Index paddedWidth = paddedSize(pictureWidth, depth);
  const Shape2D paddedShape = {{paddedHeight, paddedWidth}};
  Array2D padded(paddedShape);
  for (int line=0; line<paddedHeight; ++line) {
    for (int pixel=0; pixel<paddedWidth; ++pixel) {
      const int picLine = (line<pictureHeight)?line:(pictureHeight-1);
      const int picPixel = (pixel<pictureWidth)?pixel:(pictureWidth-1);
      padded[line][pixel] = picture[picLine][picPixel];
    }
  }
  return padded;
}

// Forward declarations of functions to implement a single wavelet level
void waveletLevelDD97(View2D&, unsigned int shift);
void inverseWaveletLevelDD97(View2D&, unsigned int shift);
void waveletLevelLeGall(View2D&, unsigned int shift);
void inverseWaveletLevelLeGall(View2D&, unsigned int shift);
void waveletLevelDD137(View2D&, unsigned int shift);
void inverseWaveletLevelDD137(View2D&, unsigned int shift);
void waveletLevelHaar(View2D&, unsigned int shift);
void inverseWaveletLevelHaar(View2D&, unsigned int shift);
void waveletLevelFidelity(View2D&, unsigned int shift);
void inverseWaveletLevelFidelity(View2D&, unsigned int shift);
void waveletLevelDaub97(View2D&, unsigned int shift);
void inverseWaveletLevelDaub97(View2D&, unsigned int shift);

void waveletLevel(View2D& p, WaveletKernel kernel) {
  switch(kernel) {
    case DD97:
      // DD97 uses 1 accuracy bit (shift=1)
      waveletLevelDD97(p, 1);
      break;
    case LeGall:
      // LeGall uses 1 accuracy bit (shift=1)
      waveletLevelLeGall(p, 1);
      break;
    case DD137:
      // DD137 uses 1 accuracy bit (shift=1)
      waveletLevelDD137(p, 1);
      break;
    case Haar0:
      // Haar0 uses no accuracy bit (shift=0)
      waveletLevelHaar(p, 0);
      break;
    case Haar1:
      // Haar1 uses 1 accuracy bit (shift=1)
      waveletLevelHaar(p, 1);
      break;
    case Fidelity:
      // Fidelity uses 1 accuracy bit (shift=0)
      waveletLevelFidelity(p, 0);
      break;
    case Daub97:
      // Daub97 uses 1 accuracy bit (shift=1)
      waveletLevelDaub97(p, 1);
      break;
    case NullKernel:
      // Null Kernel does nothing (for testing)
      break;
    default:
      throw std::invalid_argument("invalid wavelet kernel");
  }
}

const Array2D waveletTransform(const Array2D& picture, WaveletKernel kernel, int depth) {

  Array2D transform = waveletPad(picture, depth);

  // Iterate over levels
  // Note: Level numbers go from zero for high frequencies to depth for
  // the lowest ("DC") frequencies. This is the opposite way round to
  // the level definitions in the VC-2 specification.
  for (int level=0; level<depth; ++level) {
    // Create a subsampled view of (padded)picture (include only low frequency samples)
    const Index height = transform.shape()[0];
    const Index width = transform.shape()[1];
    const Index stride = utils::pow(2, level);
    View2D view =
      transform[indices[Range(0,height,stride)][Range(0,width,stride)]];
    // Do one level of in place wavelet transform
    waveletLevel(view, kernel);
  }
  return transform;
}

void inverseWaveletLevel(View2D& p, WaveletKernel kernel) {
  switch(kernel) {
    case DD97:
      // DD97 uses 1 accuracy bit (shift=1)
      inverseWaveletLevelDD97(p, 1);
      break;
    case LeGall:
      // LeGall uses 1 accuracy bit (shift=1)
      inverseWaveletLevelLeGall(p, 1);
      break;
    case DD137:
      // DD137 uses 1 accuracy bit (shift=1)
      inverseWaveletLevelDD137(p, 1);
      break;
    case Haar0:
      // Haar0 uses no accuracy bit (shift=0)
      inverseWaveletLevelHaar(p, 0);
      break;
    case Haar1:
      // Haar1 uses 1 accuracy bit (shift=1)
      inverseWaveletLevelHaar(p, 1);
      break;
    case Fidelity:
      // Fidelity uses 1 accuracy bit (shift=0)
      inverseWaveletLevelFidelity(p, 0);
      break;
    case Daub97:
      // Daub97 uses 1 accuracy bit (shift=1)
      inverseWaveletLevelDaub97(p, 1);
      break;
    case NullKernel:
      // Null Kernel does nothing (for testing)
      break;
    default:
      throw std::invalid_argument("invalid wavelet kernel");
  }
}

const Array2D inverseWaveletTransform(const Array2D& transform,
                                      WaveletKernel kernel,
                                      int depth,
                                      Shape2D shape) {
  Array2D picture = transform;
  // Iterate over levels
  // Note: Level numbers go from zero for high frequencies to depth-1 for
  // the lowest frequencies. This is the opposite way round to the level 
  // definitions in the VC-2 specification.
  for (int level=depth-1; level>=0; --level) {
    // Create a subsampled view of (padded)picture (include only low frequency samples)
    const Index height = picture.shape()[0];
    const Index width = picture.shape()[1];
    const Index stride = utils::pow(2, level);
    View2D view =
      picture[indices[Range(0,height,stride)][Range(0,width,stride)]];
    // Do one level of in place wavelet transform
    inverseWaveletLevel(view, kernel);
  }
  picture.resize(shape); // remove wavelet padding
  return picture;
}

// Return the quantisation matrix for a given wavelet kernel and depth
const Array1D quantMatrix(WaveletKernel kernel, int depth) {
  using std::vector;
  using std::min;
  if (depth < 0) throw std::domain_error("wavelet depth may not be < 0");
  Array1D qMatrix(extents[3*depth+1]);
  if (depth == 0) return (qMatrix[0]=0, qMatrix);
  float alpha, beta;
  int shift;
  switch(kernel) {
    case DD97:
      alpha = 1.280868846f;
      beta = 0.820572875f;
      shift = 1;
      break;
    case LeGall:
      alpha = 1.224744871f;
      beta = 0.847791248f;
      shift = 1;
      break;
    case DD137:
      alpha = 1.280868846f;
      beta = 0.809253958f;
      shift = 1;
      break;
    case Haar0:
      alpha = 1.414213562f;
      beta = 0.707106871f;
      shift = 0;
      break;
    case Haar1:
      alpha = 1.414213562f;
      beta = 0.707106871f;
      shift = 1;
      break;
    case Fidelity:
      alpha = 0.682408629f;
      beta = 1.367856979f;
      shift = 0;
      break;
    case Daub97:
      alpha = 1.139917028f;
      beta = 0.887168005f;
      shift = 1;
      break;
    case NullKernel: // Null Kernel does nothing (for testing)
      alpha = 1.0f;
      beta = 1.0f;
      shift = 0;
      break;
    default:
      throw std::invalid_argument("invalid wavelet kernel");
  }
  const float a2 = alpha*alpha;
  const float ab = alpha*beta;
  const float b2 = beta*beta;
  vector<float> LLGain(depth+1), LHGain(depth+1), HHGain(depth+1); //Allow space for (unused) zero level
  float minGain = FLT_MAX;
  for (int level=depth; level>0; --level) {
    const float scale = pow(a2, depth-level)/pow(2.0f, shift*(depth-level+1));
    LLGain[level] = scale*a2;
    LHGain[level] = scale*ab;
    HHGain[level] = scale*b2;
    minGain = min(min(min(LLGain[level], LHGain[level]), HHGain[level]), minGain);
  }
  vector<int> LLQuant(depth+1), LHQuant(depth+1), HHQuant(depth+1); //Allow space for (unused) zero level
  for (int level=depth; level>0; --level) {
    LLQuant[level] = static_cast<int>(floor(4.0f*log(LLGain[level]/minGain)/log(2.0f)+0.5f));
    LHQuant[level] = static_cast<int>(floor(4.0f*log(LHGain[level]/minGain)/log(2.0f)+0.5f));
    HHQuant[level] = static_cast<int>(floor(4.0f*log(HHGain[level]/minGain)/log(2.0f)+0.5f));
  }
  int index = 0;
  qMatrix[index++] = LLQuant[1];
  for (int level=1; level<=depth; ++level) {
    qMatrix[index++] = LHQuant[level];
    qMatrix[index++] = LHQuant[level];
    qMatrix[index++] = HHQuant[level];
  }
  return qMatrix;
}

using utils::pow;

// Convert a Array2D containing an in place wavelet transform into a 1D array of subbands
const BlockVector split_into_subbands(const Array2D& picture, const char waveletDepth) {
  const Index pictureHeight = picture.shape()[0];
  const Index pictureWidth = picture.shape()[1];
  const int numberOfSubbands = 3*waveletDepth+1;
  // Define a 1D array of subbands, each subband is a Array2D
  BlockVector subbands(extents[numberOfSubbands]);
  Index stride, offset; // stride is subsampling factor, offset is subsampling phase
  stride = pow(2, waveletDepth);
  subbands[0] = // LL (Low horizontal, Low vertical) "DC" subband
    picture[indices[Range(0,pictureHeight,stride)][Range(0,pictureWidth,stride)]];
  for (char level=1, band=1; level<=waveletDepth; ++level) {
    stride = pow(2, waveletDepth+1-level);
    offset = stride/2;
    // Each subband is copied from a view ((i.e. subsampled version) of the whole picture.
    subbands[band++] = //HL subband (Hight horizontal, Low vertical)
      picture[indices[Range(0,pictureHeight,stride)][Range(offset,pictureWidth,stride)]];
    subbands[band++] = //LH subband (Low horizontal, High vertical)
      picture[indices[Range(offset,pictureHeight,stride)][Range(0,pictureWidth,stride)]];
    subbands[band++] = //HH subband (Hight horizontal, High vertical)
      picture[indices[Range(offset,pictureHeight,stride)][Range(offset,pictureWidth,stride)]];
  }
  return subbands;
}

// Converts a 1D array of subbands back to a single single Array2D corresponding
// to an in-place wavelet transform
const Array2D merge_subbands(const BlockVector& subbands) {
  // TO DO: Check numberOfSubbands==3*n+1
  const int numberOfSubbands = subbands.size();
  const char waveletDepth = (numberOfSubbands-1)/3;
  const Index pictureHeight = subbands[0].shape()[0]*pow(2, waveletDepth);
  const Index pictureWidth = subbands[0].shape()[1]*pow(2, waveletDepth);
  Array2D picture(extents[pictureHeight][pictureWidth]);
  Index stride, offset;
  stride = pow(2, waveletDepth);
  picture[indices[Range(0,pictureHeight,stride)][Range(0,pictureWidth,stride)]] =
    subbands[0]; // LL (Low horizontal, Low vertical) "DC" subband
  for (char level=1, band=1; level<=waveletDepth; ++level) {
    stride = pow(2, waveletDepth+1-level);
    offset = stride/2;
    picture[indices[Range(0,pictureHeight,stride)][Range(offset,pictureWidth,stride)]] =
      subbands[band++]; //HL subband (High horizontal, Low vertical);
    picture[indices[Range(offset,pictureHeight,stride)][Range(0,pictureWidth,stride)]] =
      subbands[band++]; //LH subband (Low horizontal, High vertical);
    picture[indices[Range(offset,pictureHeight,stride)][Range(offset,pictureWidth,stride)]] =
        subbands[band++]; //HH subband (High horizontal, High vertical);
  }
  return picture;
}

void waveletLevelDD97(View2D& p, unsigned int shift) {

  const Index height = p.shape()[0];
  const Index width = p.shape()[1];

  // Do shift to introduce accuracy bits
  if (shift) {
    for (int line=0; line<height; ++line) {
      for (int pixel=0; pixel<width; ++pixel) {
	    p[line][pixel] <<= shift;
      }
    }
  }

  // horizontal predict
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = ((pixel-2)>=0) ? (pixel-2) : 0;
      const int tap1 = pixel;
      const int tap2 = ((pixel+2)<width) ? (pixel+2) : (width-2);
      const int tap3 = ((pixel+4)<width) ? (pixel+4) : (width-2);
      p[line][pixel+1] -=
        (-p[line][tap0]+9*p[line][tap1]+9*p[line][tap2]-p[line][tap3]+8)>>4;
    }
  }

  // horizontal update
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = ((pixel-1)>=0) ? (pixel-1) : 1 ;
      const int tap1 = pixel+1;
      p[line][pixel] += (p[line][tap0]+p[line][tap1] + 2)>>2;
    }
  }

  // vertical predict
  for (int line=0; line<height; line+=2) {
    const int tap0 = ((line-2)>=0) ? (line-2) : 0;
    const int tap1 = line;
    const int tap2 = ((line+2)<height) ? (line+2) : (height-2);
    const int tap3 = ((line+4)<height) ? (line+4) : (height-2);
    for (int pixel=0; pixel<width; ++pixel) {
      p[line+1][pixel] -=
        (-p[tap0][pixel]+9*p[tap1][pixel]+9*p[tap2][pixel]-p[tap3][pixel]+8)>>4;
    }
  }

  // vertical update
  for (int line=0; line<height; line+=2) {
    const int tap0 = ((line-1)>=0) ? (line-1) : 1 ;
    const int tap1 = line+1;
    for (int pixel=0; pixel<width; ++pixel) {
      p[line][pixel] += (p[tap0][pixel]+p[tap1][pixel] + 2)>>2;
    }
  }
}


void inverseWaveletLevelDD97(View2D& p, unsigned int shift) {

  const Index height = p.shape()[0];
  const Index width = p.shape()[1];

  // vertical inverse update
  for (int line=0; line<height; line+=2) {
    const int tap0 = ((line-1)>=0) ? (line-1) : 1 ;
    const int tap1 = line+1;
    for (int pixel=0; pixel<width; ++pixel) {
      p[line][pixel] -= (p[tap0][pixel]+p[tap1][pixel] + 2)>>2;
    }
  }

  // vertical inverse predict
  for (int line=0; line<height; line+=2) {
    const int tap0 = ((line-2)>=0) ? (line-2) : 0;
    const int tap1 = line;
    const int tap2 = ((line+2)<height) ? (line+2) : (height-2);
    const int tap3 = ((line+4)<height) ? (line+4) : (height-2);
    for (int pixel=0; pixel<width; ++pixel) {
      p[line+1][pixel] +=
        (-p[tap0][pixel]+9*p[tap1][pixel]+9*p[tap2][pixel]-p[tap3][pixel]+8)>>4;
    }
  }

  // horizontal inverse update
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = ((pixel-1)>=0) ? (pixel-1) : 1 ;
      const int tap1 = pixel+1;
      p[line][pixel] -= (p[line][tap0]+p[line][tap1] + 2)>>2;
    }
  }

  // horizontal inverse predict
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = ((pixel-2)>=0) ? (pixel-2) : 0;
      const int tap1 = pixel;
      const int tap2 = ((pixel+2)<width) ? (pixel+2) : (width-2);
      const int tap3 = ((pixel+4)<width) ? (pixel+4) : (width-2);
      p[line][pixel+1] +=
        (-p[line][tap0]+9*p[line][tap1]+9*p[line][tap2]-p[line][tap3]+8)>>4;
    }
  }

  // Round & shift right "shift" bits (with rounding)
  if (shift) {
    View2D::element offset = utils::pow(2, shift-1);
    for (int pixel=0; pixel<width; ++pixel) {
      for (int line=0; line<height; ++line) {
		    p[line][pixel] += offset;
		    p[line][pixel] >>= shift;
      }
    }
  }
}

void waveletLevelLeGall(View2D& p, unsigned int shift) {

  const Index height = p.shape()[0];
  const Index width = p.shape()[1];

  // Do shift to introduce accuracy bits
  if (shift) {
    for (int line=0; line<height; ++line) {
      for (int pixel=0; pixel<width; ++pixel) {
	    p[line][pixel] <<= shift;
      }
    }
  }

  // horizontal LeGall (5,3): Predict
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = pixel;
      const int tap1 = ((pixel+2)<width) ? (pixel+2) : (width-2);
      p[line][pixel+1] -= (p[line][tap0]+p[line][tap1]+1)>>1;
    }
  }

  // horizontal LeGall (5,3): Update
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = ((pixel-1)>=0) ? (pixel-1) : 1 ;
      const int tap1 = pixel+1;
      p[line][pixel] += (p[line][tap0]+p[line][tap1] + 2)>>2;
    }
  }

  // vertical LeGall (5,3): Predict
  for (int line=0; line<height; line+=2) {
    const int tap0 = line;
    const int tap1 = ((line+2)<height) ? (line+2) : (height-2);
    for (int pixel=0; pixel<width; ++pixel) {
      p[line+1][pixel] -= (p[tap0][pixel]+p[tap1][pixel]+1)>>1;
    }
  }

  // vertical LeGall (5,3): Update
  for (int line=0; line<height; line+=2) {
    const int tap0 = ((line-1)>=0) ? (line-1) : 1 ;
    const int tap1 = line+1;
    for (int pixel=0; pixel<width; ++pixel) {
      p[line][pixel] += (p[tap0][pixel]+p[tap1][pixel] + 2)>>2;
    }
  }
}


void inverseWaveletLevelLeGall(View2D& p, unsigned int shift) {

  const Index height = p.shape()[0];
  const Index width = p.shape()[1];

  // vertical LeGall (5,3): Inverse Update
  for (int line=0; line<height; line+=2) {
    const int tap0 = ((line-1)>=0) ? (line-1) : 1 ;
    const int tap1 = line+1;
    for (int pixel=0; pixel<width; ++pixel) {
      p[line][pixel] -= (p[tap0][pixel]+p[tap1][pixel] + 2)>>2;
    }
  }

  // vertical LeGall (5,3): Inverse Predict
  for (int line=0; line<height; line+=2) {
    const int tap0 = line;
    const int tap1 = ((line+2)<height) ? (line+2) : (height-2);
    for (int pixel=0; pixel<width; ++pixel) {
      p[line+1][pixel] += (p[tap0][pixel]+p[tap1][pixel]+1)>>1;
    }
  }

  // horizontal LeGall (5,3): Inverse Update
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = ((pixel-1)>=0) ? (pixel-1) : 1 ;
      const int tap1 = pixel+1;
      p[line][pixel] -= (p[line][tap0]+p[line][tap1] + 2)>>2;
    }
  }

  // horizontal LeGall (5,3): Inverse Predict
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = pixel;
      const int tap1 = ((pixel+2)<width) ? (pixel+2) : (width-2);
      p[line][pixel+1] += (p[line][tap0]+p[line][tap1]+1)>>1;
    }
  }

  // Round & shift right "shift" bits (with rounding)
  if (shift) {
    View2D::element offset = utils::pow(2, shift-1);
    for (int pixel=0; pixel<width; ++pixel) {
      for (int line=0; line<height; ++line) {
		    p[line][pixel] += offset;
		    p[line][pixel] >>= shift;
      }
    }
  }
}

void waveletLevelDD137(View2D& p, unsigned int shift) {

  const Index height = p.shape()[0];
  const Index width = p.shape()[1];

  // Do shift to introduce accuracy bits
  if (shift) {
    for (int line=0; line<height; ++line) {
      for (int pixel=0; pixel<width; ++pixel) {
	    p[line][pixel] <<= shift;
      }
    }
  }

  // horizontal predict
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = ((pixel-2)>=0) ? (pixel-2) : 0;
      const int tap1 = pixel;
      const int tap2 = ((pixel+2)<width) ? (pixel+2) : (width-2);
      const int tap3 = ((pixel+4)<width) ? (pixel+4) : (width-2);
      p[line][pixel+1] -=
        (-p[line][tap0]+9*p[line][tap1]+9*p[line][tap2]-p[line][tap3]+8)>>4;
    }
  }

  // horizontal update
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = ((pixel-3)>=0) ? (pixel-3) : 1 ;
      const int tap1 = ((pixel-1)>=0) ? (pixel-1) : 1 ;
      const int tap2 = pixel+1;
      const int tap3 = ((pixel+3)<width) ? (pixel+3) : (width-1) ;
      p[line][pixel] +=
        (-p[line][tap0]+9*p[line][tap1]+9*p[line][tap2]-p[line][tap3]+16)>>5;
    }
  }

  // vertical predict
  for (int line=0; line<height; line+=2) {
    const int tap0 = ((line-2)>=0) ? (line-2) : 0;
    const int tap1 = line;
    const int tap2 = ((line+2)<height) ? (line+2) : (height-2);
    const int tap3 = ((line+4)<height) ? (line+4) : (height-2);
    for (int pixel=0; pixel<width; ++pixel) {
      p[line+1][pixel] -=
        (-p[tap0][pixel]+9*p[tap1][pixel]+9*p[tap2][pixel]-p[tap3][pixel]+8)>>4;
    }
  }

  // vertical update
  for (int line=0; line<height; line+=2) {
    const int tap0 = ((line-3)>=0) ? (line-3) : 1 ;
    const int tap1 = ((line-1)>=0) ? (line-1) : 1 ;
    const int tap2 = line+1;
    const int tap3 = ((line+3)<height) ? (line+3) : (height-1) ;
    for (int pixel=0; pixel<width; ++pixel) {
      p[line][pixel] +=
        (-p[tap0][pixel]+9*p[tap1][pixel]+9*p[tap2][pixel]-p[tap3][pixel]+16)>>5;
    }
  }
}


void inverseWaveletLevelDD137(View2D& p, unsigned int shift) {

  const Index height = p.shape()[0];
  const Index width = p.shape()[1];

  // vertical inverse update
  for (int line=0; line<height; line+=2) {
    const int tap0 = ((line-3)>=0) ? (line-3) : 1 ;
    const int tap1 = ((line-1)>=0) ? (line-1) : 1 ;
    const int tap2 = line+1;
    const int tap3 = ((line+3)<height) ? (line+3) : (height-1) ;
    for (int pixel=0; pixel<width; ++pixel) {
      p[line][pixel] -=
        (-p[tap0][pixel]+9*p[tap1][pixel]+9*p[tap2][pixel]-p[tap3][pixel]+16)>>5;
    }
  }

  // vertical inverse predict
  for (int line=0; line<height; line+=2) {
    const int tap0 = ((line-2)>=0) ? (line-2) : 0;
    const int tap1 = line;
    const int tap2 = ((line+2)<height) ? (line+2) : (height-2);
    const int tap3 = ((line+4)<height) ? (line+4) : (height-2);
    for (int pixel=0; pixel<width; ++pixel) {
      p[line+1][pixel] +=
        (-p[tap0][pixel]+9*p[tap1][pixel]+9*p[tap2][pixel]-p[tap3][pixel]+8)>>4;
    }
  }

  // horizontal inverse update
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = ((pixel-3)>=0) ? (pixel-3) : 1 ;
      const int tap1 = ((pixel-1)>=0) ? (pixel-1) : 1 ;
      const int tap2 = pixel+1;
      const int tap3 = ((pixel+3)<width) ? (pixel+3) : (width-1) ;
      p[line][pixel] -=
        (-p[line][tap0]+9*p[line][tap1]+9*p[line][tap2]-p[line][tap3]+16)>>5;
    }
  }

  // horizontal inverse predict
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = ((pixel-2)>=0) ? (pixel-2) : 0;
      const int tap1 = pixel;
      const int tap2 = ((pixel+2)<width) ? (pixel+2) : (width-2);
      const int tap3 = ((pixel+4)<width) ? (pixel+4) : (width-2);
      p[line][pixel+1] +=
        (-p[line][tap0]+9*p[line][tap1]+9*p[line][tap2]-p[line][tap3]+8)>>4;
    }
  }

  // Round & shift right "shift" bits (with rounding)
  if (shift) {
    View2D::element offset = utils::pow(2, shift-1);
    for (int pixel=0; pixel<width; ++pixel) {
      for (int line=0; line<height; ++line) {
		    p[line][pixel] += offset;
		    p[line][pixel] >>= shift;
      }
    }
  }
}

void waveletLevelHaar(View2D& p, unsigned int shift) {

  const Index height = p.shape()[0];
  const Index width = p.shape()[1];

  // Do shift to introduce accuracy bits
  if (shift) {
    for (int line=0; line<height; ++line) {
      for (int pixel=0; pixel<width; ++pixel) {
	    p[line][pixel] <<= shift;
      }
    }
  }

  // horizontal predict
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      p[line][pixel+1] -= p[line][pixel];
    }
  }

  // horizontal update
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      p[line][pixel] += ((p[line][pixel+1] + 1)>>1);
    }
  }

  // vertical predict
  for (int pixel=0; pixel<width; ++pixel) {
    for (int line=0; line<height; line+=2) {
      p[line+1][pixel] -= p[line][pixel];
    }
  }

  // vertical update
  for (int pixel=0; pixel<width; ++pixel) {
    for (int line=0; line<height; line+=2) {
      p[line][pixel] += ((p[line+1][pixel] + 1)>>1);
    }
  }

}


void inverseWaveletLevelHaar(View2D& p, unsigned int shift) {

  const Index height = p.shape()[0];
  const Index width = p.shape()[1];

  // vertical Haar: Inverse Update
  for (int line=0; line<height; line+=2) {
    for (int pixel=0; pixel<width; ++pixel) {
      p[line][pixel] -= ((p[line+1][pixel] + 1)>>1);
    }
  }

  // vertical Haar: Inverse Predict
  for (int line=0; line<height; line+=2) {
    for (int pixel=0; pixel<width; ++pixel) {
      p[line+1][pixel] += p[line][pixel];
    }
  }

  // horizontal Haar: Inverse Update
  for (int pixel=0; pixel<width; pixel+=2) {
    for (int line=0; line<height; ++line) {
      p[line][pixel] -= ((p[line][pixel+1] + 1)>>1);
    }
  }

  // horizontal Haar: Inverse Predict
  for (int pixel=0; pixel<width; pixel+=2) {
    for (int line=0; line<height; ++line) {
      p[line][pixel+1] += p[line][pixel];
	  }
  }

  // Round & shift right "shift" bits (with rounding)
  if (shift) {
    View2D::element offset = utils::pow(2, shift-1);
    for (int pixel=0; pixel<width; ++pixel) {
      for (int line=0; line<height; ++line) {
		    p[line][pixel] += offset;
		    p[line][pixel] >>= shift;
      }
    }
  }
}

void waveletLevelFidelity(View2D& p, unsigned int shift) {

  const Index height = p.shape()[0];
  const Index width = p.shape()[1];

  // Do shift to introduce accuracy bits
  if (shift) {
    for (int line=0; line<height; ++line) {
      for (int pixel=0; pixel<width; ++pixel) {
	    p[line][pixel] <<= shift;
      }
    }
  }

  // horizontal type 1
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = ((pixel-7)>=0) ? (pixel-7) : 1;
      const int tap1 = ((pixel-5)>=0) ? (pixel-5) : 1;
      const int tap2 = ((pixel-3)>=0) ? (pixel-3) : 1;
      const int tap3 = ((pixel-1)>=0) ? (pixel-1) : 1;
      const int tap4 = pixel+1;
      const int tap5 = ((pixel+3)<width) ? (pixel+3) : (width-1);
      const int tap6 = ((pixel+5)<width) ? (pixel+5) : (width-1);
      const int tap7 = ((pixel+7)<width) ? (pixel+7) : (width-1);
      p[line][pixel] +=
        (-8*p[line][tap0]+21*p[line][tap1]-46*p[line][tap2]+161*p[line][tap3]
         +161*p[line][tap4]-46*p[line][tap5]+21*p[line][tap6]-8*p[line][tap7]+128)>>8;
    }
  }

  // horizontal type 4
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = ((pixel-6)>=0) ? (pixel-6) : 0;
      const int tap1 = ((pixel-4)>=0) ? (pixel-4) : 0;
      const int tap2 = ((pixel-2)>=0) ? (pixel-2) : 0;
      const int tap3 = pixel;
      const int tap4 = ((pixel+2)<width) ? (pixel+2) : (width-2);
      const int tap5 = ((pixel+4)<width) ? (pixel+4) : (width-2);
      const int tap6 = ((pixel+6)<width) ? (pixel+6) : (width-2);
      const int tap7 = ((pixel+8)<width) ? (pixel+8) : (width-2);
      p[line][pixel+1] -=
        (-2*p[line][tap0]+10*p[line][tap1]-25*p[line][tap2]+81*p[line][tap3]
         +81*p[line][tap4]-25*p[line][tap5]+10*p[line][tap6]-2*p[line][tap7]+128)>>8;
    }
  }

  // vertical type 1
  for (int line=0; line<height; line+=2) {    
    const int tap0 = ((line-7)>=0) ? (line-7) : 1;
    const int tap1 = ((line-5)>=0) ? (line-5) : 1;
    const int tap2 = ((line-3)>=0) ? (line-3) : 1;
    const int tap3 = ((line-1)>=0) ? (line-1) : 1;
    const int tap4 = line+1;
    const int tap5 = ((line+3)<height) ? (line+3) : (height-1);
    const int tap6 = ((line+5)<height) ? (line+5) : (height-1);
    const int tap7 = ((line+7)<height) ? (line+7) : (height-1);
    for (int pixel=0; pixel<width; ++pixel) {
      p[line][pixel] +=
        (-8*p[tap0][pixel]+21*p[tap1][pixel]-46*p[tap2][pixel]+161*p[tap3][pixel]
         +161*p[tap4][pixel]-46*p[tap5][pixel]+21*p[tap6][pixel]-8*p[tap7][pixel]+128)>>8;
    }
  }

  // vertical type 4
  for (int line=0; line<height; line+=2) {
    const int tap0 = ((line-6)>=0) ? (line-6) : 0;
    const int tap1 = ((line-4)>=0) ? (line-4) : 0;
    const int tap2 = ((line-2)>=0) ? (line-2) : 0;
    const int tap3 = line;
    const int tap4 = ((line+2)<height) ? (line+2) : (height-2);
    const int tap5 = ((line+4)<height) ? (line+4) : (height-2);
    const int tap6 = ((line+6)<height) ? (line+6) : (height-2);
    const int tap7 = ((line+8)<height) ? (line+8) : (height-2);
    for (int pixel=0; pixel<width; ++pixel) {
      p[line+1][pixel] -=
        (-2*p[tap0][pixel]+10*p[tap1][pixel]-25*p[tap2][pixel]+81*p[tap3][pixel]
         +81*p[tap4][pixel]-25*p[tap5][pixel]+10*p[tap6][pixel]-2*p[tap7][pixel]+128)>>8;
    }
  }

}


void inverseWaveletLevelFidelity(View2D& p, unsigned int shift) {

  const Index height = p.shape()[0];
  const Index width = p.shape()[1];

  // vertical type 3
  for (int line=0; line<height; line+=2) {
    const int tap0 = ((line-6)>=0) ? (line-6) : 0;
    const int tap1 = ((line-4)>=0) ? (line-4) : 0;
    const int tap2 = ((line-2)>=0) ? (line-2) : 0;
    const int tap3 = line;
    const int tap4 = ((line+2)<height) ? (line+2) : (height-2);
    const int tap5 = ((line+4)<height) ? (line+4) : (height-2);
    const int tap6 = ((line+6)<height) ? (line+6) : (height-2);
    const int tap7 = ((line+8)<height) ? (line+8) : (height-2);
    for (int pixel=0; pixel<width; ++pixel) {
      p[line+1][pixel] +=
        (-2*p[tap0][pixel]+10*p[tap1][pixel]-25*p[tap2][pixel]+81*p[tap3][pixel]
         +81*p[tap4][pixel]-25*p[tap5][pixel]+10*p[tap6][pixel]-2*p[tap7][pixel]+128)>>8;
    }
  }

  // vertical type 2
  for (int line=0; line<height; line+=2) {
    const int tap0 = ((line-7)>=0) ? (line-7) : 1;
    const int tap1 = ((line-5)>=0) ? (line-5) : 1;
    const int tap2 = ((line-3)>=0) ? (line-3) : 1;
    const int tap3 = ((line-1)>=0) ? (line-1) : 1;
    const int tap4 = line+1;
    const int tap5 = ((line+3)<height) ? (line+3) : (height-1);
    const int tap6 = ((line+5)<height) ? (line+5) : (height-1);
    const int tap7 = ((line+7)<height) ? (line+7) : (height-1);
    for (int pixel=0; pixel<width; ++pixel) {
      p[line][pixel] -=
        (-8*p[tap0][pixel]+21*p[tap1][pixel]-46*p[tap2][pixel]+161*p[tap3][pixel]
         +161*p[tap4][pixel]-46*p[tap5][pixel]+21*p[tap6][pixel]-8*p[tap7][pixel]+128)>>8;
    }
  }

  // horizontal type 3
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = ((pixel-6)>=0) ? (pixel-6) : 0;
      const int tap1 = ((pixel-4)>=0) ? (pixel-4) : 0;
      const int tap2 = ((pixel-2)>=0) ? (pixel-2) : 0;
      const int tap3 = pixel;
      const int tap4 = ((pixel+2)<width) ? (pixel+2) : (width-2);
      const int tap5 = ((pixel+4)<width) ? (pixel+4) : (width-2);
      const int tap6 = ((pixel+6)<width) ? (pixel+6) : (width-2);
      const int tap7 = ((pixel+8)<width) ? (pixel+8) : (width-2);
      p[line][pixel+1] +=
        (-2*p[line][tap0]+10*p[line][tap1]-25*p[line][tap2]+81*p[line][tap3]
         +81*p[line][tap4]-25*p[line][tap5]+10*p[line][tap6]-2*p[line][tap7]+128)>>8;

    }
  }

  // horizontal type 2
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = ((pixel-7)>=0) ? (pixel-7) : 1;
      const int tap1 = ((pixel-5)>=0) ? (pixel-5) : 1;
      const int tap2 = ((pixel-3)>=0) ? (pixel-3) : 1;
      const int tap3 = ((pixel-1)>=0) ? (pixel-1) : 1;
      const int tap4 = pixel+1;
      const int tap5 = ((pixel+3)<width) ? (pixel+3) : (width-1);
      const int tap6 = ((pixel+5)<width) ? (pixel+5) : (width-1);
      const int tap7 = ((pixel+7)<width) ? (pixel+7) : (width-1);
      p[line][pixel] -=
        (-8*p[line][tap0]+21*p[line][tap1]-46*p[line][tap2]+161*p[line][tap3]
         +161*p[line][tap4]-46*p[line][tap5]+21*p[line][tap6]-8*p[line][tap7]+128)>>8;
    }
  }

  // Round & shift right "shift" bits (with rounding)
  if (shift) {
    View2D::element offset = utils::pow(2, shift-1);
    for (int pixel=0; pixel<width; ++pixel) {
      for (int line=0; line<height; ++line) {
		    p[line][pixel] += offset;
		    p[line][pixel] >>= shift;
      }
    }
  }
}

void waveletLevelDaub97(View2D& p, unsigned int shift) {

  const Index height = p.shape()[0];
  const Index width = p.shape()[1];

  // Do shift to introduce accuracy bits
  if (shift) {
    for (int line=0; line<height; ++line) {
      for (int pixel=0; pixel<width; ++pixel) {
	    p[line][pixel] <<= shift;
      }
    }
  }

  // horizontal type 4
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = pixel;
      const int tap1 = ((pixel+2)<width) ? (pixel+2) : (width-2);
      p[line][pixel+1] -= (6497*p[line][tap0]+6497*p[line][tap1]+2048)>>12;
    }
  }

  // horizontal type 2
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = ((pixel-1)>=0) ? (pixel-1) : 1;
      const int tap1 = pixel+1;
      p[line][pixel] -= (217*p[line][tap0]+217*p[line][tap1]+2048)>>12;
    }
  }

  // horizontal type 3
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = pixel;
      const int tap1 = ((pixel+2)<width) ? (pixel+2) : (width-2);
      p[line][pixel+1] += (3616*p[line][tap0]+3616*p[line][tap1]+2048)>>12;
    }
  }

  // horizontal type 1
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = ((pixel-1)>=0) ? (pixel-1) : 1;
      const int tap1 = pixel+1;
      p[line][pixel] += (1817*p[line][tap0]+1817*p[line][tap1]+2048)>>12;
    }
  }

  // vertical type 4
  for (int line=0; line<height; line+=2) {
    const int tap0 = line;
    const int tap1 = ((line+2)<height) ? (line+2) : (height-2);
    for (int pixel=0; pixel<width; ++pixel) {
      p[line+1][pixel] -= (6497*p[tap0][pixel]+6497*p[tap1][pixel]+2048)>>12;
    }
  }

  // vertical type 2
  for (int line=0; line<height; line+=2) { 
    const int tap0 = ((line-1)>=0) ? (line-1) : 1;
    const int tap1 = line+1;
    for (int pixel=0; pixel<width; ++pixel) {
      p[line][pixel] -= (217*p[tap0][pixel]+217*p[tap1][pixel]+2048)>>12;
    }
  }

  // vertical type 3
  for (int line=0; line<height; line+=2) {
    const int tap0 = line;
    const int tap1 = ((line+2)<height) ? (line+2) : (height-2);
    for (int pixel=0; pixel<width; ++pixel) {
      p[line+1][pixel] += (3616*p[tap0][pixel]+3616*p[tap1][pixel]+2048)>>12;
    }
  }

  // vertical type 1
  for (int line=0; line<height; line+=2) {
    const int tap0 = ((line-1)>=0) ? (line-1) : 1;
    const int tap1 = line+1;
    for (int pixel=0; pixel<width; ++pixel) {
      p[line][pixel] += (1817*p[tap0][pixel]+1817*p[tap1][pixel]+2048)>>12;
    }
  }
}


void inverseWaveletLevelDaub97(View2D& p, unsigned int shift) {

  const Index height = p.shape()[0];
  const Index width = p.shape()[1];

  // vertical type 2
  for (int line=0; line<height; line+=2) {
    const int tap0 = ((line-1)>=0) ? (line-1) : 1;
    const int tap1 = line+1;
    for (int pixel=0; pixel<width; ++pixel) {
      p[line][pixel] -= (1817*p[tap0][pixel]+1817*p[tap1][pixel]+2048)>>12;
    }
  }

  // vertical type 4
  for (int line=0; line<height; line+=2) {
    const int tap0 = line;
    const int tap1 = ((line+2)<height) ? (line+2) : (height-2);
    for (int pixel=0; pixel<width; ++pixel) {
      p[line+1][pixel] -= (3616*p[tap0][pixel]+3616*p[tap1][pixel]+2048)>>12;
    }
  }

  // vertical type 1
  for (int line=0; line<height; line+=2) { 
    const int tap0 = ((line-1)>=0) ? (line-1) : 1;
    const int tap1 = line+1;
    for (int pixel=0; pixel<width; ++pixel) {
      p[line][pixel] += (217*p[tap0][pixel]+217*p[tap1][pixel]+2048)>>12;
    }
  }

  // vertical type 3
  for (int line=0; line<height; line+=2) {
    const int tap0 = line;
    const int tap1 = ((line+2)<height) ? (line+2) : (height-2);
    for (int pixel=0; pixel<width; ++pixel) {
      p[line+1][pixel] += (6497*p[tap0][pixel]+6497*p[tap1][pixel]+2048)>>12;
    }
  }

  // horizontal type 2
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = ((pixel-1)>=0) ? (pixel-1) : 1;
      const int tap1 = pixel+1;
      p[line][pixel] -= (1817*p[line][tap0]+1817*p[line][tap1]+2048)>>12;
    }
  }

  // horizontal type 4
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = pixel;
      const int tap1 = ((pixel+2)<width) ? (pixel+2) : (width-2);
      p[line][pixel+1] -= (3616*p[line][tap0]+3616*p[line][tap1]+2048)>>12;
    }
  }

  // horizontal type 1
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = ((pixel-1)>=0) ? (pixel-1) : 1;
      const int tap1 = pixel+1;
      p[line][pixel] += (217*p[line][tap0]+217*p[line][tap1]+2048)>>12;
    }
  }

  // horizontal type 3
  for (int line=0; line<height; ++line) {
    for (int pixel=0; pixel<width; pixel+=2) {
      const int tap0 = pixel;
      const int tap1 = ((pixel+2)<width) ? (pixel+2) : (width-2);
      p[line][pixel+1] += (6497*p[line][tap0]+6497*p[line][tap1]+2048)>>12;
    }
  }

  // Round & shift right "shift" bits (with rounding)
  if (shift) {
    View2D::element offset = utils::pow(2, shift-1);
    for (int pixel=0; pixel<width; ++pixel) {
      for (int line=0; line<height; ++line) {
		    p[line][pixel] += offset;
		    p[line][pixel] >>= shift;
      }
    }
  }
}

const Picture waveletTransform(const Picture& input, WaveletKernel kernel, int waveletDepth) {
  const int lumaHeight = paddedSize(input.format().lumaHeight(), waveletDepth);
  const int lumaWidth = paddedSize(input.format().lumaWidth(), waveletDepth);
  const int chromaHeight = paddedSize(input.format().chromaHeight(), waveletDepth);
  const int chromaWidth = paddedSize(input.format().chromaWidth(), waveletDepth);
  const ColourFormat uvFormat = input.format().chromaFormat();
  PictureFormat const transformFormat(lumaHeight, lumaWidth, chromaHeight, chromaWidth, uvFormat);
  Picture transform(transformFormat);
  transform.y(waveletTransform(input.y(), kernel, waveletDepth));
  transform.c1(waveletTransform(input.c1(), kernel, waveletDepth));
  transform.c2(waveletTransform(input.c2(), kernel, waveletDepth));
  return transform;
}

const Picture inverseWaveletTransform(const Picture& transform,
                                      WaveletKernel kernel,
                                      int depth,
                                      PictureFormat format) {
  Picture picture(format);
  const Shape2D lumaShape(format.lumaShape());
  const Shape2D chromaShape(format.chromaShape());
  picture.y(inverseWaveletTransform(transform.y(), kernel, depth, lumaShape));
  picture.c1(inverseWaveletTransform(transform.c1(), kernel, depth, chromaShape));
  picture.c2(inverseWaveletTransform(transform.c2(), kernel, depth, chromaShape));
  return picture;
}
