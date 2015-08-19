/*********************************************************************/
/* Frame.cpp                                                         */
/* Author: Tim Borer                                                 */
/* This version 26th August 2011                                     */
/*                                                                   */
/* Defines stuff related to colour frames.                           */
/* 6th May 2011: Initial version                                     */
/* 2nd June 2011: Added "frame" methods and some other implementaions*/
/* 5th August 2011: Corrected header declaration                     */
/* 26th Sept 2011: Added alternate constructor (see header)          */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

#include "Frame.h"

Frame::Frame(int height, int width, ColourFormat format, bool  interlaced, bool  topFieldFirst):
  Picture(PictureFormat(height, width, format)), intl(interlaced), tff(topFieldFirst) {
}

Frame::Frame(const PictureFormat& format, bool  interlaced, bool  topFieldFirst):
  Picture(format), intl(interlaced), tff(topFieldFirst) {
}

bool Frame::interlaced() const {
  return intl;
}

void Frame::interlaced(bool i) {
  intl = i;
}

bool Frame::topFieldFirst() const {
  return tff;
}

void Frame::topFieldFirst(bool t) {
  tff = t;
}

const Picture Frame::topField() const {
  // Get parameters of field
  const int height = format().lumaHeight()/2;
  const int width = format().lumaWidth();
  const ColourFormat chromaFormat = format().chromaFormat();
  // Construct interlaced field
  Picture picture(PictureFormat(height, width, chromaFormat));
  // Set components of interlaced field
  const int top = 0;
  const int yBottom = format().lumaHeight();
  const int uvBottom = format().chromaHeight();
  using boost::indices;
  picture.y(y()[indices[Range(top,yBottom,2)][Range()]]);
  picture.c1(c1()[indices[Range(top,uvBottom,2)][Range()]]);
  picture.c2(c2()[indices[Range(top,uvBottom,2)][Range()]]);
  return picture;
}

void Frame::topField(const Picture& f) {
  const int top = 0;
  const int yBottom = format().lumaHeight();
  const int uvBottom = format().chromaHeight();
  using boost::indices;
  luma[indices[Range(top,yBottom,2)][Range()]] = f.y();
  chroma1[indices[Range(top,uvBottom,2)][Range()]] = f.c1();
  chroma2[indices[Range(top,uvBottom,2)][Range()]] = f.c2();
}

const Picture Frame::bottomField() const {
  // Get parameters of field
  const int height = format().lumaHeight()/2;
  const int width = format().lumaWidth();
  const ColourFormat chromaFormat = format().chromaFormat();
  // Construct interlaced field
  Picture picture(PictureFormat(height, width, chromaFormat));
  // Set components of interlaced field
  const int top = 1;
  const int yBottom = format().lumaHeight();
  const int uvBottom = format().chromaHeight();
  using boost::indices;
  picture.y(y()[indices[Range(top,yBottom,2)][Range()]]);
  picture.c1(c1()[indices[Range(top,uvBottom,2)][Range()]]);
  picture.c2(c2()[indices[Range(top,uvBottom,2)][Range()]]);
  return picture;
}

void Frame::bottomField(const Picture& f) {
  const int top = 1;
  const int yBottom = format().lumaHeight();
  const int uvBottom = format().chromaHeight();
  using boost::indices;
  luma[indices[Range(top,yBottom,2)][Range()]] = f.y();
  chroma1[indices[Range(top,uvBottom,2)][Range()]] = f.c1();
  chroma2[indices[Range(top,uvBottom,2)][Range()]] = f.c2();
}

const Picture Frame::firstField() const {
  return (tff? topField() : bottomField());
}

void Frame::firstField(const Picture& f) {
  tff ? topField(f) : bottomField(f);
}

const Picture Frame::secondField() const {
  return (tff? bottomField() : topField());
}

void Frame::secondField(const Picture& f) {
  tff ? bottomField(f) : topField(f);
}

const Picture& Frame::frame() const {
  return *this;
}

void Frame::frame(const Picture& f) {
  Picture::operator=(f);
}
