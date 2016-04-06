#ifndef AKT_H
#define AKT_H

#define __STDC_LIMIT_MACROS

#include "version.h"
#include <stdint.h>

#include "htslib/hts.h"
#include <vector>
#include <string>
#include <getopt.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <stdlib.h> 
#include <iterator>
#include <omp.h>
#include <map>
#include <list>
#include <algorithm>
#include <set>
#include <stdexcept>     
#include "samples.hpp"

extern "C" {
#include "htslib/synced_bcf_reader.h"
#include "filter.h"
}

template <typename T> string to_string( T x ){ 
	return static_cast< std::ostringstream & >( ( std::ostringstream() << std::dec << x ) ).str(); 
}

        
using namespace std;

void umessage(const char type);

int pca_main(int argc,char **argv);

int kin_main(int argc,char **argv);

int relatives_main(int argc,char **argv);

int ibd_main(int argc,char **argv);

int cluster_main(int argc,char **argv);

int stats_main(int argc,char **argv);

int mendel_main(int argc,char **argv);

int ldplot_main(int argc,char **argv);

int admix_main(int argc,char **argv);


#endif
