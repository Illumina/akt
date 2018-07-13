CC=gcc
CXX=g++

OMP=-fopenmp

CXXFLAGS = -std=c++11 -O2  $(OMP) -mpopcnt
CFLAGS = -O2  $(OMP) -mpopcnt

.PHONY: all
all: akt

HTSDIR=htslib-1.6
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
IFLAGS = -I$(HTSDIR)  -I./
LFLAGS = -lz -lm  -lpthread

.PHONY: no_omp
no_omp: CXXFLAGS = -std=c++11 -O2 
no_omp: CFLAGS = -O2 
no_omp: all

.PHONY: release
release: CXXFLAGS = -std=c++11 -O2  $(OMP) -mpopcnt
release: CFLAGS = -O2  $(OMP) -mpopcnt
release: LFLAGS +=  -static
release: all

.PHONY: debug
debug: CXXFLAGS = -std=c++11 -g -O1 -Wall
debug: CFLAGS = -g -O1  -lz -lm -lpthread
debug: all

.PHONY: profile
profile: CXXFLAGS = -pg -O2 $(OMP)
profile: CFLAGS =  -pg -O2 $(OMP)
profile: all

##generates a version
GIT_HASH := $(shell git describe --abbrev=4 --always )
BCFTOOLS_VERSION=1.6
VERSION = 0.3.2
ifneq "$(wildcard .git)" ""
VERSION = $(shell git describe --always)
endif
version.hh:
	echo '#define AKT_VERSION "$(VERSION)"' > $@
	echo '#define BCFTOOLS_VERSION "$(BCFTOOLS_VERSION)"' >> $@

OBJS= utils.o pedphase.o family.o reader.o vcfpca.o relatives.o kin.o pedigree.o unrelated.o cluster.o HaplotypeBuffer.o Genotype.o
.cpp.o:
	$(CXX) $(CXXFLAGS) $(IFLAGS) -c -o $@ $<
.c.o:
	$(CC) $(CXXFLAGS) -c -o $@ $<

##akt code
cluster.o: cluster.cpp cluster.hh
family.o: family.cpp family.hh
relatives.o: relatives.cpp 
unrelated.o: unrelated.cpp 
vcfpca.o: vcfpca.cpp RandomSVD.hh
kin.o: kin.cpp 
pedigree.o: pedigree.cpp pedigree.hh
reader.o: reader.cpp 
pedphase.o: pedphase.cpp pedphase.hh utils.hh HaplotypeBuffer.o
utils.o: utils.cpp utils.hh
HaplotypeBuffer.o: HaplotypeBuffer.cpp HaplotypeBuffer.hh
akt: akt.cpp version.hh $(OBJS) $(HTSLIB)
	$(CXX) $(CXXFLAGS)   -o akt akt.cpp $(OBJS) $(IFLAGS) $(HTSLIB) $(LFLAGS) $(CXXFLAGS)

.PHONY: clean
clean:
	rm -f *.o akt version.hh

.PHONY: test
test: akt
	cd test/;bash -e test.sh
