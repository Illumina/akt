//
// Created by O'Connell, Jared on 8/3/17.
//

#ifndef AKT_PEDPHASE_H
#define AKT_PEDPHASE_H
#include "akt.hpp"
#include "pedigree.h"
#include "mendel.h"
#include "utils.h"
#include <deque>
#include <iomanip>
#include <numeric>
#include <stdlib.h>

typedef struct _args
{
    bool regions_is_file;
    bool targets_is_file;
    char output_type;
    char *pedigree, *inputfile, *include, *regions, *targets, *outfile;
} args;


class PedPhaser
{

public:
    PedPhaser(args &a);
    int phaseTrio(int idx,int *gt_arr);

private:
    int flushBuffer();
    bcf_srs_t *_sr;
    bcf_hdr_t *_hdr;
    htsFile *_out_fh;
    bcf_hdr_t *_out_hdr;
    sampleInfo *_ped;
    int _nsample;
    deque<bcf1_t *> _line_buffer;

};


#endif //AKT_PEDPHASE_H
