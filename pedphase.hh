//
// Created by O'Connell, Jared on 8/3/17.
//

#ifndef AKT_PEDPHASE_H
#define AKT_PEDPHASE_H
#include <deque>
#include <iomanip>
#include <numeric>
#include <stdlib.h>
#include <bitset>

#include "akt.hh"
#include "pedigree.hh"
#include "utils.hh"
#include "version.hh"
#include "HaplotypeBuffer.hh"
#include "Genotype.hh"

typedef struct _args
{
    int nthreads;
    bool regions_is_file;
    bool targets_is_file;
    char output_type;
    string exclude_chromosomes;
    const char *pedigree, *inputfile, *include, *regions, *targets, *outfile;
} args;


int phase_by_transmission(Genotype & kid_gt,Genotype & dad_gt,Genotype & mum_gt);
bool is_mendel_inconsistent(Genotype  kid,Genotype  dad,Genotype  mum);

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
    int _num_gt,_num_ps,_num_rps;//stores length of ps/gt
    int *_gt_array;
    int32_t *_ps_array,*_rps_array;
    vector<int> _chromosomes_to_ignore;//dont phase these chromosomes
    vector<bool> _sample_has_been_phased;
    vector< pair<int,int> >  _parental_genotypes;
    void main();
};

#endif //AKT_PEDPHASE_H
