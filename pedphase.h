//
// Created by O'Connell, Jared on 8/3/17.
//

#ifndef AKT_PEDPHASE_H
#define AKT_PEDPHASE_H
#include "akt.hpp"
#include "pedigree.h"
#include "utils.h"
#include "version.h"
#include <deque>
#include <iomanip>
#include <numeric>
#include <stdlib.h>
#include <bitset>


typedef struct _args
{
  int nthreads;
  bool regions_is_file;
  bool targets_is_file;
  char output_type;
  string exclude_chromosomes;
  const char *pedigree, *inputfile, *include, *regions, *targets, *outfile;
} args;


class PedPhaser
{

 public:
  PedPhaser(args &a);
  ~PedPhaser();

 private:
  int mendelPhase(int idx,int *gt_array,int *ps_array=NULL);  
  void setup_io(args &a);
  void setup_output(args &a);  
  int flushBuffer();
  bool chromosome_is_in_ignore_list(bcf1_t *record);
  bcf_srs_t *_bcf_reader;
  htsFile *_out_file;
  bcf_hdr_t *_out_header,*_in_header;
  sampleInfo *_pedigree;
  int _num_sample;
  deque<bcf1_t *> _line_buffer;
  int *_gt_array,*_gt_array_dup;
  int32_t *_ps_array,*_rps_array;
  vector<int> _chromosomes_to_ignore;//dont phase these chromosomes
  vector<bool> _sample_has_been_phased;
  void main();
};


#endif //AKT_PEDPHASE_H
