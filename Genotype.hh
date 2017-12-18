#ifndef  GENOTYPE_H
#define GENOTYPE_H

#include <string>
#include <iostream>

extern "C" {
#include "htslib/vcf.h"
}

#include "utils.hh"

class Genotype
{
public:
    Genotype(int g0,int g1);
    Genotype(int idx,int *gt_array,int *ps_array=NULL); 
    void setGenotype(int g0,int g1);   
    bool isHet();
    bool isMissing();
    bool isHaploid();
    bool is_phased();
    void setPhase(bool phased);
    int first();
    int second();
    int getGenotype(int idx);
    int swap(); 
    int update_bcf_gt_array(int *gt_array,int index,int *ps_array=nullptr);
    std::string print();
    int ps();
    
private:
    int _g0,_g1,_ps;
    bool _is_haploid,_is_phased;
};

#endif
