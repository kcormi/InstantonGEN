CXXFLAGS=$(shell root-config --cflags)
GLIBS=$(shell root-config --glibs)

SRCS:=$(wildcard src/*.cc)

all: InstantonGEN

InstantonGEN: $(SRCS)
	g++ $(CXXFLAGS) $^ -o $@ -I$(HOME)/local/include/ $(GLIBS) $(shell /afs/cern.ch/user/j/jinw/public/LHAPDF/bin/lhapdf-config --cflags --ldflags)
