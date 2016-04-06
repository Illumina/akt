#pragma once

#include "akt.hpp"

//data structure for pedigree analysis.
class sampleInfo{ 
 public:
  //reads a pedigree froma  plink.fam file
  sampleInfo(string fname);
  //note this function changes the samples the hdr will read
  sampleInfo(string fname,bcf_hdr_t *hdr);

  vector<string> fid,id,dad,mum; 
  vector<int> dadidx,mumidx;
  int N,ntrio,nduo;

 private:
  int readFromPlinkFam(string fname);
  int buildIndex();
};

