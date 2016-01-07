AUTOMAKE_OPTIONS = foreign

if ENABLE_OTHER_DECODERS
OPT_SUBDIRS = DecodeHQ DecodeLD
else
OPT_SUBDIRS = 
endif

SUBDIRS = boost tclap Library DecodeStream EncodeHQ-CBR EncodeHQ-ConstQ EncodeLD $(OPT_SUBDIRS)

DISTCLEANFILES = vc2reference-stdint.h
