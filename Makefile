CXXFLAGS=$(shell root-config --cflags)
GLIBS=$(shell root-config --glibs)

SRCS:=$(wildcard src/*.cc)

all: InstantonGEN

InstantonGEN: $(SRCS)
	g++ $(CXXFLAGS) $^ -o $@ -I$(HOME)/local/include/ $(GLIBS) $(shell /work/kcormier/instantons/hepmc/instanton_production/LHAPDF/bin/lhapdf-config --cflags --ldflags)

clean:
	rm ./InstantonGEN
