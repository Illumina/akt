CC=gcc
CXX=g++

OMP=-fopenmp

CXXFLAGS = -O2  $(OMP) -mpopcnt
CFLAGS = -O2  $(OMP)

all: akt

HTSDIR=htslib-1.5
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
IFLAGS = -I$(HTSDIR)  -I./
LFLAGS = -lz -lm

no_omp: CXXFLAGS = -O2 
no_omp: CFLAGS = -O2 
no_omp: LFLAGS = -lz -lm -lpthread
no_omp: all

default: CXXFLAGS = -O2  $(OMP)
default: CFLAGS = -O2  $(OMP)
default: all

release: CXXFLAGS = -O2  $(OMP)
release: CFLAGS = -O2  $(OMP)
release: LFLAGS = -lz -lm -static
release: all

debug: CXXFLAGS = -g -O1 -lz -lm -lpthread
debug: CFLAGS =  -g -O1  -lz -lm -lpthread
debug: all

profile: CXXFLAGS = -pg -O2 $(OMP)
profile: CFLAGS =  -pg -O2 $(OMP)
profile: all

##generates a version
GIT_HASH := $(shell git describe --abbrev=4 --always )
BCFTOOLS_VERSION=1.5
VERSION = 0.3.0
ifneq "$(wildcard .git)" ""
VERSION = $(shell git describe --always)
endif
version.h:
	echo '#define AKT_VERSION "$(VERSION)"' > $@
	echo '#define BCFTOOLS_VERSION "$(BCFTOOLS_VERSION)"' >> $@

OBJS= utils.o pedphase.o family.o reader.o vcfpca.o relatives.o kin.o pedigree.o  unrelated.o cluster.o
##bcftools code
filter.o: filter.c filter.h
	$(CC)  $(CFLAGS) $(IFLAGS) -c $<  
version.o: version.c filter.h version.h
	$(CC)  $(CFLAGS) $(IFLAGS) -c $<  
##akt code
cluster.o: cluster.cpp cluster.hpp
	$(CXX) $(CXXFLAGS)  -c cluster.cpp $(IFLAGS)
family.o: family.cpp family.hpp
	$(CXX) $(CXXFLAGS)  -c family.cpp $(IFLAGS)
relatives.o: relatives.cpp 
	$(CXX) $(CXXFLAGS)  -c relatives.cpp $(IFLAGS)
unrelated.o: unrelated.cpp 
	$(CXX) $(CXXFLAGS)  -c unrelated.cpp $(IFLAGS)
vcfpca.o: vcfpca.cpp RandomSVD.hpp
	$(CXX) $(CXXFLAGS) -c vcfpca.cpp $(IFLAGS)
kin.o: kin.cpp 
	$(CXX) $(CXXFLAGS)  -c kin.cpp $(IFLAGS)
pedigree.o: pedigree.cpp pedigree.h
	$(CXX) $(CXXFLAGS)  -c $< $(IFLAGS)
reader.o: reader.cpp 
	$(CXX) $(CXXFLAGS)  -c reader.cpp $(IFLAGS)
pedphase.o: pedphase.cpp pedphase.h utils.h
	$(CXX) $(CXXFLAGS) -c $< $(IFLAGS)
utils.o: utils.cpp utils.h
	$(CXX) $(CXXFLAGS) -c $< $(IFLAGS)
akt: akt.cpp version.h $(OBJS) $(HTSLIB)
	$(CXX) $(CXXFLAGS)   -o akt akt.cpp $(OBJS) $(IFLAGS) $(HTSLIB) $(LFLAGS) $(CXXFLAGS)
clean:
	rm *.o akt version.h
test: akt
	cd test/;bash -e test.sh
