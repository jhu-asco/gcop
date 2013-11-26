CC=g++
CXX=g++
SUBDIRS = .

AM_CPPFLAGS = -Wall -Wno-deprecated -fno-exceptions -O3 -g -I. -fPIC
AM_CFLAGS = $(AM_CPPFLAGS) -fno-rtti

lib_LIBRARIES = libest.a

libest_a_LIBADD = *.o

libest_a_DEPENDENCIES = $(libest_a_LIBADD)

libest_a_SOURCES = \
	model.cc model.h \
	filter.cc filter.h \
	ukf.cc ukf.h \
	quat.cc quat.h \
	insmodel.cc insmodel.h \
	insukf.cc insukf.h




libestdir=$(includedir)/est
libest_HEADERS = model.h filter.h ukf.h quat.h insmodel.h insukf.h