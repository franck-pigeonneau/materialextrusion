include $(shell rheolef-config --libdir)/rheolef/rheolef.mk
LD = $(CXX)
CC = $(CXX)
CXXFLAGS  = $(INCLUDES_RHEOLEF)
LDLIBS    = $(LIBS_RHEOLEF)
default: calculdts
calculdts.o : calculdts.cc donneesdts.cc advscalar.cc
