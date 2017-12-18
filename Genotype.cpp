#include "Genotype.hh"

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
    if(idx<0)
    {
	_g0=-1;
	_g1=-1;
	_is_haploid=false;
	_is_phased=false;
	_ps=bcf_int32_missing;
    }
    else
    {
	_is_phased = ps_array && ps_array[idx]!=bcf_int32_missing;
	_ps = ps_array ? ps_array[idx] : bcf_int32_missing;
	setGenotype(gt_array[2*idx],gt_array[2*idx+1]);
    }
}

int Genotype::ps()
{
    return _ps;
}

bool Genotype::is_phased()
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
    std::cerr << "idx="<<idx<<std::endl;
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

int Genotype::update_bcf_gt_array(int *gt_array,int index,int32_t *ps_array)
{
    if(!isMissing())
	gt_array[index * 2] = _is_phased ? bcf_gt_phased(first()) : bcf_gt_unphased(first());
    else
	gt_array[index * 2] = bcf_gt_missing;
    
    if(isHaploid())
	gt_array[index * 2 + 1] = bcf_int32_vector_end;
    else
	gt_array[index * 2 + 1] = _is_phased ? bcf_gt_phased(second()) : bcf_gt_unphased(second()); 
    
    if(ps_array)
	ps_array[index]=_ps;
    return(0);
}

void Genotype::setPhase(bool phased)
{
    _is_phased=phased;
}

std::string Genotype::print()
{
    if(isMissing()) return "./.";
    std::stringstream ss;
    ss<<first();
    if(isHaploid())
	return ss.str();
    
    if(_is_phased)
	ss<<"|";
    else
	ss<<"/";
    ss<<second();
    if(_ps!=bcf_int32_missing)
	ss<<":"<<_ps;
    return ss.str();
}
