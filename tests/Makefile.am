check_LTLIBRARIES = libgtest.la
libgtest_la_SOURCES = ../googletest/googletest/src/gtest-all.cc
libgtest_la_CPPFLAGS =  -I$(top_srcdir)/googletest/googletest/include \
                        -I$(top_srcdir)/googletest/googletest -pthread 
libgtest_la_LDFLAGS = -pthread

check_PROGRAMS = gtest
gtest_SOURCES = gtest.cpp 
gtest_LDADD =   libgtest.la \
                $(top_srcdir)/src/Library/libVC2.la \
                $(BOOST_LDFLAGS) \
                $(BOOST_SYSTEM_LIB)\
                $(BOOST_THREAD_LIB)
gtest_LDFLAGS = -pthread
gtest_CPPFLAGS =    -I$(top_srcdir)/src/ \
                    -I$(top_srcdir)/src/Library \
                    -I$(top_srcdir)/googletest/googletest/include \
                    -I$(top_srcdir)/googletest/googletest -pthread  \
                    $(VC2REFERENCE_CFLAGS) \
	                $(BOOST_CFLAGS)

TESTS = gtest