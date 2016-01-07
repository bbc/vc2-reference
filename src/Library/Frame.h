/************************************************************************/
/* Frame.h                                                              */
/* Author: Tim Borer                                                    */
/* This version 26th September 2011                                     */
/*                                                                      */
/* Declares stuff related to colour frames.                             */
/* 6th May 2011: Initial version                                        */
/* 2nd June 2011: Added "frame" method                                  */
/* 26th Sept 2011: Added constructor from height, width & colour format */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file      */
/************************************************************************/

#ifndef FRAME_6MAY11
#define FRAME_6MAY11

#include "Picture.h"

class Frame: public Picture {
  public:
    Frame(int height, int width, ColourFormat, bool interlaced=false, bool topFieldFirst=true);
    Frame(const PictureFormat&, bool interlaced=false, bool topFieldFirst=true);
    bool interlaced() const;
    void interlaced(bool);
    bool topFieldFirst() const;
    void topFieldFirst(bool);
    const Picture topField() const;
    void topField(const Picture&);
    const Picture bottomField() const;
    void bottomField(const Picture&);
    const Picture firstField() const;
    void firstField(const Picture&);
    const Picture secondField() const;
    void secondField(const Picture&);
    const Picture& frame() const;
    void frame(const Picture&);
  private:
    bool intl;
    bool tff;
};
  
#endif //FRAME_6MAY11
