#ifndef IBD_H
#define IBD_H

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

using namespace std;

extern void read_pairs(ifstream &in, vector< pair<string, string> > &relpairs, map<string,int> &name_to_id);

extern void make_pair_list(vector< pair<string, string> > &relpairs, vector<string> names);

//trying to keep some of the kinship stuff in this class
//i will expand this to hold genotype bitset as well - jared.
class Kinship 
{
public:
    Kinship(int nsample);
    void estimate_ibd(float & ibd0, float & ibd1, float & ibd2,float & ibd3,bool normalise=true) ;
    void update_n(float p);
    float _n00,_n10,_n11,_n20,_n21,_n22;
    int _nsample;  
    vector<float> _af;//allele freqs
};

#endif
