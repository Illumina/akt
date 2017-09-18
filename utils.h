//
// Created by O'Connell, Jared on 8/3/17.
//

#ifndef AKT_UTILS_H
#define AKT_UTILS_H

#include "akt.hpp"
#include <string>
#include <vector>
#include <sstream>

bool is_genotyped(int *gt,int idx);

int stringSplit(string & s,vector<string> & split);

class Genotype
{
public:
    Genotype(int g0,int g1);
    Genotype(int idx,int *gt_arr); 
    void setGenotype(int g0,int g1);   
    bool isHet();
    bool isMissing();
    bool isHaploid();
    int first();
    int second();
    int getGenotype(int idx);
    int swap(); 
    int update_bcf_gt_array(int *gt_array,int index);
    
private:
    int _g0,_g1;
    bool _is_haploid;    
};


#endif //AKT_UTILS_H
