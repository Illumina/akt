#pragma once

#include "akt.hpp"
#include "utils.h"

//data structure for pedigree analysis.
class sampleInfo
{
public:
    //reads a pedigree froma  plink.fam file
    sampleInfo(string fname);
    sampleInfo(string fname, bcf_hdr_t *hdr);
    sampleInfo(bcf_hdr_t *hdr);

    vector <string> fid, id, dad, mum;
    vector<int> dadidx, mumidx;
    vector<int> sex, status;
    map <pair<int, int>, vector<int> > parent_map;
    int N, ntrio, nduo;

    int getStatus(int i)
    {
        assert(i < N);
        return status[i];
    }

    string *getID(int i)
    {
        assert(i < N);
        return &id[i];
    }

    int getDadIndex(int idx);
    int getMumIndex(int idx);
private:
    int readFromPlinkFam(string fname);
    void alignWithVcf(bcf_hdr_t *hdr);
    int buildIndex();

};

