#ifndef IBD_H
#define IBD_H

#include "akt.hpp"
#include "logs.hpp"
#include "reader.hpp"
#include <bitset>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <ctype.h>
#include <limits> 
#include <iterator>
#include <map>
#include <string>

///Size of bitset. Arbitrary but 128 works well in practice
#define L 128

using namespace std;

extern void read_pairs(ifstream &in, vector< pair<string, string> > &relpairs, map<string,int> &name_to_id);

extern void make_pair_list(vector< pair<string, string> > &relpairs, vector<string> names);

//main class for kinship calculations.
//stores the genotype bitset and relevant counters
class Kinship 
{
public:
    Kinship(int nsample);
    void estimateKinship(int j1,int j2,float & ibd0, float & ibd1, float & ibd2,float & ibd3,float &ks,int method) ;
    void estimateIBD(float & ibd0, float & ibd1, float & ibd2,float & ibd3,bool normalise=true) ;
    void addGenotypes(int *gt_arr,float p);
    void addGenotypes(int *gt_arr);
    float _n00,_n10,_n11,_n20,_n21,_n22;
    int _nsample,_markers,_bc;  
    vector<float> _af;//allele freqs
    vector< vector< vector< bitset<L> > > > _bits; ///[sample][site][type][val]
};

#endif
