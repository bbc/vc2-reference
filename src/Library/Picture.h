/*********************************************************************/
/* Picture.h                                                         */
/* Author: Tim Borer                                                 */
/* This version 4th September 2012                                   */
/*                                                                   */
/* Declares stuff related to colour pictures.                        */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

#ifndef PICTURE_3MAY11
#define PICTURE_3MAY11

#include <iosfwd>

#include "Arrays.h"

enum ColourFormat {UNKNOWN, CF444, CF422, CF420, RGB}; //UNKOWN needed for PictureFormat default constructor

std::ostream& operator<<(std::ostream& os, ColourFormat format);

std::istream& operator>>(std::istream& strm, ColourFormat& format);

class PictureFormat {
  public:
    PictureFormat();  //Needed for Picture default constructor
    PictureFormat(int height, int width, ColourFormat);
    // Use this constructor for non-standard picture (e.g. wavelet transforms)
    PictureFormat(int lumaHeight, int lumaWidth, int chromaHeight, int chromaWidth, ColourFormat);
    // Constructor checks format has right number of samples (or
    //   returns default if height or width = 0 (i.e. unknown) 
    PictureFormat(int height, int width, ColourFormat, int imageSamples);
    operator const bool() const; //Test if valid PictureFormat
    const ColourFormat chromaFormat() const;
    const Shape2D lumaShape() const; 
    const Shape2D chromaShape()const;
    const int lumaHeight() const;
    const int lumaWidth() const;
    const int chromaHeight() const;
    const int chromaWidth() const;
    const int samples() const;
  private:
    void construct(int height, int width, ColourFormat);
    // Guess frame format from number of samples in the frame and the colour format.
    // The list of supported frame formats is in "FrameResolutions.h"
    void guessFormat(int imageSamples, ColourFormat);
    int yHeight;
    int yWidth;
    int uvHeight;
    int uvWidth;
    ColourFormat uvFormat;
};

class Picture {
public:
  Picture();  //Needed for arrays of pictures. Picture can be set by assignment (only)
  Picture(const PictureFormat&);
  Picture(const int height, const int width, const ColourFormat);
  Picture(const PictureFormat&, const Array2D& y, const Array2D& c1, const Array2D& c2);
  PictureFormat format() const;
  const Array2D& y() const;
  const Array2D& c1() const;
  const Array2D& c2() const;
  void y(const Array2D&);
  void c1(const Array2D&);
  void c2(const Array2D&);
protected:
  PictureFormat picFormat;
  Array2D luma, chroma1, chroma2;
  friend std::istream& operator >> (std::istream&, Picture&);
};

typedef boost::multi_array<Picture, 2> PictureArray;

// Get the shape of a PictureArray
const Shape2D shape(const PictureArray&);

const PictureArray split_into_blocks(const Picture& picture, int ySlices, int xSlices);

const Picture merge_blocks(const PictureArray& blocks);

// Clip a Picture to specified limits
// First function clips all components to the same values (good for RGB)
const Picture clip(const Picture& picture, const int min_value, const int max_value);
// Second function clips luma and chroma values separately (good for YUV)
const Picture clip(const Picture& picture,
                   const int luma_min, const int luma_max,
                   const int chroma_min, const int chroma_max);

//**** Picture IO declarations ****//

namespace pictureio {

  // Use same io manipulators as in ArrayIO.h
  using arrayio::left_justified;
  using arrayio::right_justified;
  using arrayio::offset_binary;
  using arrayio::signed_binary;
  using arrayio::unsigned_binary;
  using arrayio::text;

  //Use same wordWidth class as in ArrayIO.h
  using arrayio::wordWidth;

  using arrayio::ioFormat;

  // Use to set data format as an alternative to manipulators above or
  // use for different data formats for luma and chroma
  class format {
    public:
      format(ioFormat f): lumaFormat(f), chromaFormat(f) {};
      format(ioFormat lf, ioFormat cf): lumaFormat(lf), chromaFormat(cf) {}; 
      void operator () (std::ios_base& stream) const;
    private:
      const ioFormat lumaFormat;
      const ioFormat chromaFormat;
  };

  class bitDepth {
    public:
      bitDepth(int d): lumaBitDepth(d), chromaBitDepth(d) {};
      bitDepth(int ld, int cd): lumaBitDepth(ld), chromaBitDepth(cd) {}; 
      void operator () (std::ios_base& stream) const;
    private:
      const int lumaBitDepth;
      const int chromaBitDepth;
  };

  class offset {
    public:
      offset(int o): lumaOffset(o), chromaOffset(o) {};
      offset(int lo, int co): lumaOffset(lo), chromaOffset(co) {};
      void operator () (std::ios_base& stream) const;
    private:
      const int lumaOffset;
      const int chromaOffset;
  };

} // end namespace pictureio

// ostream io format, format manipulator
std::ostream& operator << (std::ostream& stream, pictureio::format);

// istream io format, format manipulator
std::istream& operator >> (std::istream& stream, pictureio::format);

// ostream bit Depth format manipulator
std::ostream& operator << (std::ostream& stream, pictureio::bitDepth);

// istream bit Depth format manipulator
std::istream& operator >> (std::istream& stream, pictureio::bitDepth);

// ostream offset format manipulator
std::ostream& operator << (std::ostream& stream, pictureio::offset);

// istream offset format manipulator
std::istream& operator >> (std::istream& stream, pictureio::offset);

std::istream& operator >> (std::istream& stream, Picture& array);

std::ostream& operator << (std::ostream& stream, const Picture& array);
  
#endif //PICTURE_3MAY11
