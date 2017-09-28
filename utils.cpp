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

pair<int,int> getGenotype(int idx,int *gt_array)
{
    pair<int,int> ret;
    if(idx<0)
    {
        ret.first = ret.second = bcf_gt_missing;
    }
    else
    {
        ret.first = bcf_gt_allele(gt_array[idx * 2]);
        ret.second = bcf_gt_allele(gt_array[idx * 2 + 1]);
    }
    return(ret);
}

Genotype::Genotype(int g0,int g1)
{
    setGenotype(g0,g1);
}

void Genotype::setGenotype(int g0,int g1)
{
    _is_haploid=false;
    if(bcf_gt_is_missing(g0))
    {
	_g0 = -1;
    }
    else
    {
	_g0=bcf_gt_allele(g0);	    
    }
    
    if(bcf_gt_is_missing(g1))
    {
	_g1 = -1;
    }
    else if(g1==bcf_int32_vector_end)
    {
	_g1=bcf_int32_vector_end;
	_is_haploid=true;
    }
    else
    {
	_g1=bcf_gt_allele(g1);	    
    }

}

Genotype::Genotype(int idx,int *gt_array,int *ps_array)
{
    _is_phased = ps_array && ps_array[idx]!=bcf_int32_missing;
    if(idx<0)
    {
	_g0=-1;
	_g1=-1;
	_is_haploid=false;
    }
    else
    {
	setGenotype(gt_array[2*idx],gt_array[2*idx+1]);
    }
}

bool Genotype::isPhased()
{
    return( _is_phased || !(isHet() || isHaploid()) );
}

bool Genotype::isHet()
{
    return(_g0!=-1 && _g1!=-1 && !_is_haploid && _g0!=_g1);
}

bool Genotype::isMissing()
{
    return(_g0==-1||_g1==-1);
}

bool Genotype::isHaploid()
{
    return(_is_haploid);
}

int Genotype::first()
{
    return(_g0);
}

int Genotype::second()
{
    if(_is_haploid)
    {
	return(bcf_int32_vector_end);
    }
    else
    {
	return(_g1);
    }
}

int Genotype::getGenotype(int idx)
{
    if(idx==0)
    {
	return(_g0);
    }
    
    if(idx==1)
    {
	if(_is_haploid)
	{
	    return(bcf_int32_vector_end);
	}
	else
	{
	    return(_g1);
	}
    }
    cerr << "idx="<<idx<<endl;
    die("invalid genotype.");
    return(-1);
}

int Genotype::swap()
{
    if(!_is_haploid)
    {
	int tmp = _g0;
	_g0 = _g1;
	_g1 = tmp;
    }
    return(0);
}

int Genotype::update_bcf_gt_array(int *gt_array,int index)
{
    if(!isMissing())
    {
	gt_array[index * 2] = bcf_gt_phased(first());
	if(isHaploid())
	{
	    gt_array[index * 2 + 1] = bcf_int32_vector_end;
	}
	else
	{
	    gt_array[index * 2 + 1] = bcf_gt_phased(second());	    
	}
    }
    
    return(0);
}
