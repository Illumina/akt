//
// Created by O'Connell, Jared on 8/3/17.
//

#include "utils.h"

using namespace std;

bool is_genotyped(int *gt,int idx)
{
    if(idx>=0)
    {
        return( !bcf_gt_is_missing(gt[2*idx]) && !bcf_gt_is_missing(gt[2*idx+1]));
    }
    else
    {
        return(0);
    }
  //  return( !bcf_gt_is_missing(gt[2*idx]) && gt[2*idx]!=bcf_int32_vector_end && !bcf_gt_is_missing(gt[2*idx+1]) && gt[2*idx+1]!=bcf_int32_vector_end );

}


int stringSplit(string & s,vector<string> & split)
{
    split.clear();
    stringstream ss;
    ss << s;
    string tmp;
    while(ss>>tmp)
    {
        split.push_back(tmp);
    }
    return(split.size());
}

void getGenotype(int idx,int *gt_arr,int *ret)
{
    if(idx<0)
    {
        ret[0]=ret[1]=bcf_gt_missing;
    }
    else
    {
        ret[0] = bcf_gt_allele(gt_arr[idx * 2]);
        ret[1] = bcf_gt_allele(gt_arr[idx * 2 + 1]);
    }
    //cerr << bcf_gt_allele(gt_arr[idx * 2])<< " " << bcf_gt_allele(gt_arr[idx * 2+1]) << " -> " <<ret[0]<<" "<<ret[1]<<endl;
}
