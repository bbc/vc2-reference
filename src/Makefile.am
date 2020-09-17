AUTOMAKE_OPTIONS = foreign

if ENABLE_FRAME_DECODER
OPT_SUBDIRS = DecodeFrame
else
OPT_SUBDIRS = 
endif

SUBDIRS = boost tclap Library DecodeStream EncodeStream $(OPT_SUBDIRS)

DISTCLEANFILES = vc2reference-stdint.h
