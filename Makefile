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

debug: CXXFLAGS = -g -O1
debug: CFLAGS =  -g -O1
debug: all

profile: CXXFLAGS = -pg -O2 $(OMP)
profile: CFLAGS =  -pg -O2 $(OMP)
profile: all

##generates a version
GIT_HASH := $(shell git describe --abbrev=4 --always )
BCFTOOLS_VERSION=1.5
VERSION = 0.2.0
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
family.o: family.cpp family.hpp
	$(CXX) $(CXXFLAGS)  -c family.cpp $(IFLAGS)
relatives.o: relatives.cpp 
	$(CXX) $(CXXFLAGS)  -c relatives.cpp $(IFLAGS)
unrelated.o: unrelated.cpp 
	$(CXX) $(CXXFLAGS)  -c unrelated.cpp $(IFLAGS)
vcfpca.o: vcfpca.cpp RandomSVD.hpp
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
circularBuffer.o: circularBuffer.cpp circularBuffer.hpp
	$(CXX) $(CXXFLAGS) -c  circularBuffer.cpp $(IFLAGS)	
tag.o: tag.cpp  tag.hpp 
	$(CXX) $(CXXFLAGS) -c tag.cpp $(IFLAGS)
prune.o: prune.cpp  prune.hpp 
	$(CXX) $(CXXFLAGS) -c prune.cpp $(IFLAGS)
tdt.o: tdt.cpp 
	$(CXX) $(CXXFLAGS) -c tdt.cpp $(IFLAGS)
pedphase.o: pedphase.cpp pedphase.h utils.h
	$(CXX) $(CXXFLAGS) -c $< $(IFLAGS)
utils.o: utils.cpp utils.h
	$(CXX) $(CXXFLAGS) -c $< $(IFLAGS)
akt: version.h akt.cpp utils.o tdt.o pedphase.o family.o admix.o ldplot.o reader.o vcfpca.o relatives.o kin.o ibd.o cluster.o stats.o pedigree.o mendel.o filter.o version.o grm.o metafreq.o tag.o prune.o circularBuffer.o unrelated.o $(HTSLIB)
	$(CXX) $(CXXFLAGS)   -o akt akt.cpp  utils.o pedphase.o family.o tdt.o metafreq.o admix.o ldplot.o reader.o vcfpca.o relatives.o kin.o ibd.o cluster.o stats.o pedigree.o mendel.o filter.o version.o prune.o grm.o tag.o circularBuffer.o unrelated.o $(IFLAGS) $(HTSLIB) $(LFLAGS) $(CXXFLAGS)
clean:
	rm *.o akt version.h

test: akt
	cd test/;bash -e test.sh
