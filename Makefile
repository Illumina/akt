CC=gcc
CXX=g++

OMP=-fopenmp

CXXFLAGS = -O3 -DNDEBUG $(OMP)
CFLAGS = -O3 -DNDEBUG $(OMP)

all: akt

HTSDIR=htslib-1.3
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
IFLAGS = -I$(HTSDIR)  -I./
LFLAGS = -lz -lm

no_omp: CXXFLAGS = -O3 -DNDEBUG 
no_omp: CFLAGS = -O3 -DNDEBUG
no_omp: LFLAGS = -lz -lm -lpthread
no_omp: all

default: CXXFLAGS = -O3 -DNDEBUG $(OMP)
default: CFLAGS = -O3 -DNDEBUG $(OMP)
default: all

debug: CXXFLAGS = -g -O1  $(OMP)
debug: CFLAGS =  -g -O1   $(OMP)
debug: all

profile: CXXFLAGS = -pg -O3 $(OMP)
profile: CFLAGS =  -pg -O3 $(OMP)
profile: all

##generates a version
GIT_HASH := $(shell git describe --abbrev=4 --always )
BCFTOOLS_VERSION=1.3
VERSION = 0.1.0
ifneq "$(wildcard .git)" ""
VERSION = $(shell git describe --always)
endif
version.h:
	echo '#define VERSION "$(VERSION)"' > $@
	echo '#define BCFTOOLS_VERSION "$(BCFTOOLS_VERSION)"' >> $@

##bcftools code
filter.o: filter.c filter.h
	$(CC)  $(CFLAGS) $(IFLAGS) -c $<  
version.o: version.c filter.h version.h
	$(CC)  $(CFLAGS) $(IFLAGS) -c $<  
##akt code
relatives.o: relatives.cpp 
	$(CXX) $(CXXFLAGS)  -c relatives.cpp $(IFLAGS)
vcfpca.o: vcfpca.cpp 
	$(CXX) $(CXXFLAGS) -c vcfpca.cpp $(IFLAGS)
cluster.o: cluster.cpp 
	$(CXX) $(CXXFLAGS)  -c cluster.cpp $(IFLAGS)
kin.o: kin.cpp  kin.hpp
	$(CXX) $(CXXFLAGS)  -c kin.cpp $(IFLAGS)
ibd.o: ibd.cpp 
	$(CXX) $(CXXFLAGS)  -c ibd.cpp $(IFLAGS)
grm.o: grm.cpp 
	$(CXX) $(CXXFLAGS)  -c $< $(IFLAGS)
pedigree.o: pedigree.cpp pedigree.h
	$(CXX) $(CXXFLAGS)  -c $< $(IFLAGS)
mendel.o: pedigree.h mendel.cpp 
	$(CXX) $(CXXFLAGS)  -c mendel.cpp $(IFLAGS)
stats.o: stats.cpp 
	$(CXX) $(CXXFLAGS)  -c stats.cpp $(IFLAGS)
reader.o: reader.cpp 
	$(CXX) $(CXXFLAGS)  -c reader.cpp $(IFLAGS)
ldplot.o: ldplot.cpp 
	$(CXX) $(CXXFLAGS)  -c ldplot.cpp $(IFLAGS)
admix.o: admix.cpp 
	$(CXX) $(CXXFLAGS) -c admix.cpp $(IFLAGS)
metafreq.o: metafreq.cpp 
	$(CXX) $(CXXFLAGS) -c metafreq.cpp $(IFLAGS)
akt: version.h akt.cpp metafreq.o admix.o ldplot.o reader.o vcfpca.o relatives.o kin.o ibd.o cluster.o stats.o pedigree.o mendel.o filter.o version.o grm.o $(HTSLIB)
	$(CXX) $(CXXFLAGS)   -o akt akt.cpp metafreq.o admix.o ldplot.o reader.o vcfpca.o relatives.o kin.o ibd.o cluster.o stats.o pedigree.o mendel.o filter.o version.o grm.o $(IFLAGS) $(HTSLIB) $(LFLAGS) $(CXXFLAGS)
clean:
	rm *.o akt version.h

test: akt
	cd test/;bash -e test.sh
