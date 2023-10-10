include $(shell rheolef-config --libdir)/rheolef/rheolef.mk
LD = $(CXX)
CC = $(CXX)
CXXFLAGS  = $(INCLUDES_RHEOLEF)
LDLIBS    = $(LIBS_RHEOLEF)
default: gennewtonianfluid
gennewtonianfluid.o : gennewtonianfluid.cc donneesGNewFluid.cc viscopowerlaw.cc rateofstraintensor.cc calculviscosite.cc stokessolver.cc thermalfield.cc
