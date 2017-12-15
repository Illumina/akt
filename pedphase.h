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



//This algorithm to harmonise read-back phased variants with pedigree transmission.
//1. Perform phase-by-transmission, temporarily ignoring all read-back-phasing information.
//2. Count the number of times a read-back phased variant agrees with the pedigree phasing (per sample).
//3. If the count in step 2 equals <50% of pedigree resolved hets, flip all read-back phased alleles in the sample.
//4. Calculate the concordance of the flipped alleles with pedigree inheritance.
//5. If the value from 4 is 100%, move FORMAT/PS to FORMAT/RPS to indicate this phase set is fully in agreement with the pedigree.
class HaplotypeBuffer
{
public:
    HaplotypeBuffer(size_t num_sample,sampleInfo *pedigree);
    void push_back(int32_t *gt_array, int32_t *ps_array=nullptr);
    void clear();
    void phase();
    void align(HaplotypeBuffer & haps_to_align);
    Genotype get_genotype(size_t variant_index,size_t sample_index);
    int get_num_variant() {return _num_variant;};
    int get_num_sample() {return _num_sample;};    
    void update_bcf1_genotypes(int linenum,int32_t *gt_array, int32_t *ps_array,int32_t *rps_array);
    bool is_mendel_consistent(int linenum);
private:
    size_t _num_sample,_num_variant;
    vector< vector< Genotype > > _kid,_dad,_mum;
    sampleInfo *_pedigree;
    vector<int> _index_of_first_child;
    map< pair<int,int>,pair<int,int> > _phase_set_vote;
    vector< vector<bool> >_sample_was_mendel_phased;
    vector<bool> _line_is_mendel_consistent;
};

bool is_mendel_inconsistent(Genotype kid,Genotype dad,Genotype mum);
class PedPhaser
{

 public:
  PedPhaser(args &a);
  ~PedPhaser();

 private:
  int mendel_phase(int idx,int *gt_array,int *ps_array=NULL);  
  void setup_io(args &a);
  void setup_output(args &a);

  int flush_buffer();

  //Simply checks if this chromosome should be ignored and piped to output as is (eg. chrMT or chrY).
  bool chromosome_is_in_ignore_list(bcf1_t *record);
  bcf_srs_t *_bcf_reader;
  htsFile *_out_file;
  bcf_hdr_t *_out_header,*_in_header;
  sampleInfo *_pedigree;
  int _num_sample;
  deque<bcf1_t *> _line_buffer;
  int _num_gt,_num_ps;//stores length of ps/gt
  int *_gt_array,*_gt_array_dup;
  int32_t *_ps_array,*_rps_array;
  vector<int> _chromosomes_to_ignore;//dont phase these chromosomes
  vector<bool> _sample_has_been_phased;
  vector< pair<int,int> >  _parental_genotypes;
  void main();
};

#endif //AKT_PEDPHASE_H
