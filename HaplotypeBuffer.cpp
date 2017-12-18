#include "HaplotypeBuffer.hh"

HaplotypeBuffer::HaplotypeBuffer(size_t num_sample,sampleInfo *pedigree) :
    _num_sample(num_sample),_pedigree(pedigree)
{
    _num_variant=0;
    _index_of_first_child.assign(_num_sample,-1);
    _kid_vote.resize(_num_sample);
    _dad_vote.resize(_num_sample);
    _mum_vote.resize(_num_sample);
}

Genotype HaplotypeBuffer::get_genotype(size_t variant_index,size_t sample_index)
{
    assert(variant_index<_num_variant);
    assert(sample_index<_num_sample);
    return(_kid[variant_index][sample_index]);
}

bool HaplotypeBuffer::is_mendel_consistent(size_t linenum)
{
    assert(linenum>=0 && linenum<_num_variant);
    return _line_is_mendel_consistent[linenum];
}

void HaplotypeBuffer::push_back(int32_t *gt_array, int32_t *ps_array)
{
    _kid.push_back(std::vector<Genotype>());
    _dad.push_back(std::vector<Genotype>());
    _mum.push_back(std::vector<Genotype>());
    _is_aligned_with_pedigree.push_back(std::vector<bool>(_num_sample,true));
    for(size_t kid_index=0;kid_index<_num_sample;kid_index++)
    {
	int dad_index = _pedigree->getDadIndex(kid_index);
	int mum_index = _pedigree->getMumIndex(kid_index);
	_kid.back().emplace_back(kid_index, gt_array, ps_array);
	_dad.back().emplace_back(dad_index, gt_array, ps_array);
	_mum.back().emplace_back(mum_index, gt_array, ps_array);

	if(dad_index!=-1)
	    if(_index_of_first_child[dad_index]==-1)
		_index_of_first_child[dad_index] = kid_index;
	if(mum_index!=-1)
	    if(_index_of_first_child[mum_index]==-1)
		_index_of_first_child[mum_index] = kid_index;	
    }
    _num_variant++;
    assert(_kid.size()==_num_variant);
    assert(_mum.size()==_num_variant);
    assert(_dad.size()==_num_variant);
}


void HaplotypeBuffer::copy_from_parents()
{
    for(size_t variant_index=0;variant_index<_num_variant;variant_index++)
    {    
	for(size_t dst_index=0;dst_index<_num_sample;dst_index++)
	{
	    int src_index = _index_of_first_child[dst_index];	    
	    int dad_index = _pedigree->getDadIndex(dst_index);
	    int mum_index = _pedigree->getMumIndex(dst_index);
	    if(dad_index==-1&&mum_index==-1&&src_index!=-1)
	    {		
		if(dst_index== _pedigree->getMumIndex(src_index))
		    _kid[variant_index][dst_index]=_mum[variant_index][src_index];
		else if(dst_index== _pedigree->getDadIndex(src_index))
		    _kid[variant_index][dst_index]=_dad[variant_index][src_index];
		else
		    die("invalid pedigree");
		
		_is_aligned_with_pedigree[variant_index][dst_index]=_is_aligned_with_pedigree[variant_index][src_index];
	    }   
	}	
    }    
}

void HaplotypeBuffer::phase()
{
    _line_is_mendel_consistent.assign(_num_variant,true);
    for(size_t variant_index=0;variant_index<_num_variant;variant_index++)
    {
	for(size_t sample_index=0;sample_index<_num_sample;sample_index++)
	{
	    int status = phase_by_transmission(_kid[variant_index][sample_index],_dad[variant_index][sample_index],_mum[variant_index][sample_index]);
	    if(status==-1) _line_is_mendel_consistent[variant_index]=false;
	}
    }
    copy_from_parents();
}

void HaplotypeBuffer::swap(int variant,int sample) { _kid[variant][sample].swap(); } 

void HaplotypeBuffer::setPhase(int variant,int sample,bool phase) { _kid[variant][sample].setPhase(phase);}

void HaplotypeBuffer::align(HaplotypeBuffer & haps_to_align)
{
    assert(haps_to_align.get_num_sample() == get_num_sample());
    assert(haps_to_align.get_num_variant() == get_num_variant());
    align_sample(kid(),haps_to_align.kid(),_kid_vote);
    align_sample(mum(),haps_to_align.mum(),_mum_vote); 
    align_sample(dad(),haps_to_align.dad(),_dad_vote);
    check_pedigree_aligned();
    copy_from_parents();    
}

void HaplotypeBuffer::align_sample(std::vector< std::vector< Genotype > > & dst,
				   std::vector< std::vector< Genotype > > & src,
				   std::vector< std::unordered_map<int,pair<int,int> > >  & phase_set_vote)
{
    for(size_t sample_index=0;sample_index<_num_sample;sample_index++)
    {
	for(size_t variant_index=0;variant_index<_num_variant;variant_index++)
	{
	    Genotype g = src[variant_index][sample_index];
	    if(g.ps()!=bcf_int32_missing && g.is_phased() && g.isHet() && dst[variant_index][sample_index].is_phased())
	    {
		assert(g.ps()!=bcf_int32_missing);
		if(!phase_set_vote[sample_index].count(g.ps()))
		    phase_set_vote[sample_index][g.ps()] = pair<int,int>(0,0);
		phase_set_vote[sample_index][g.ps()].second++;
		if(g.first() != dst[variant_index][sample_index].first())	
		    phase_set_vote[sample_index][g.ps()].first++;
	    }
	}      
	for(size_t variant_index=0;variant_index<_num_variant;variant_index++)
	{
	    Genotype g = src[variant_index][sample_index];	    
	    if(g.ps()!=bcf_int32_missing)
	    {
		bool flip=phase_set_vote[sample_index][g.ps()].first > phase_set_vote[sample_index][g.ps()].second/2;
		if(phase_set_vote[sample_index][g.ps()].second>0)
		{
		    if(flip)
		    {
			src[variant_index][sample_index].swap();
			src[variant_index][sample_index].setPhase(true);
		    }
		}
		dst[variant_index][sample_index] = src[variant_index][sample_index];
	    }
	}
    }    
}

bool HaplotypeBuffer::is_sample_phased(int variant,int sample,
				       std::vector< std::vector< Genotype > > & genotypes,
				       std::vector< std::unordered_map<int,pair<int,int> > >  & phase_set_vote)
{
    bool phased = genotypes[variant][sample].is_phased();
    int ps = genotypes[variant][sample].ps();
    if(ps!=bcf_int32_missing)
    {
	int a=phase_set_vote[sample][ps].first;
	int b=phase_set_vote[sample][ps].second;
//	std::cerr<<sample<<" "<<ps<<" "<<a<<","<<b<<std::endl;
	phased &= b>0;
	phased &= a==0 || a==b;
    }
    return(phased);
}
    
void HaplotypeBuffer::check_pedigree_aligned()
{
    for(size_t variant_index=0;variant_index<_num_variant;variant_index++)
    {
	for(size_t sample_index=0;sample_index<_num_sample;sample_index++)
	{
	    bool consistent = _kid[variant_index][sample_index].first() == _mum[variant_index][sample_index].first();
	    bool phased = is_sample_phased(variant_index,sample_index,_kid,_kid_vote);
	    phased &= is_sample_phased(variant_index,sample_index,_mum,_mum_vote);
	    _is_aligned_with_pedigree[variant_index][sample_index]= !phased || consistent;
	    consistent = _kid[variant_index][sample_index].second() == _dad[variant_index][sample_index].first();
	    phased = is_sample_phased(variant_index,sample_index,_kid,_kid_vote);
	    phased &= is_sample_phased(variant_index,sample_index,_dad,_dad_vote);
	    _is_aligned_with_pedigree[variant_index][sample_index] = _is_aligned_with_pedigree[variant_index][sample_index] && (!phased || consistent);
	}
    }
}

bool HaplotypeBuffer::is_phase_set_aligned_with_pedigree(int sample,int phaseset)
{
    assert(phaseset!=bcf_int32_missing);
    bool ret = _kid_vote[sample][phaseset].second>0;
    ret &= _kid_vote[sample][phaseset].first==0 || _kid_vote[sample][phaseset].first==_kid_vote[sample][phaseset].second;
    for(size_t variant_index=0;variant_index<_num_variant;variant_index++)
	if(_kid[variant_index][sample].ps()==phaseset)
	    ret &= _is_aligned_with_pedigree[variant_index][sample];
    return(ret);
}

void HaplotypeBuffer::update_bcf1_genotypes(size_t linenum,int32_t *gt_array, int32_t *ps_array,int32_t *rps_array)
{
    assert(linenum>=0 && linenum<_num_variant);
    for(size_t i=0;i<_num_sample;i++)
    {
	ps_array[i]=rps_array[i]=bcf_int32_missing;
	if(_kid[linenum][i].ps()!=bcf_int32_missing)
	{
	    if(is_phase_set_aligned_with_pedigree(i,_kid[linenum][i].ps()))
		_kid[linenum][i].update_bcf_gt_array(gt_array,i,rps_array);
	    else
		_kid[linenum][i].update_bcf_gt_array(gt_array,i,ps_array);
	}
	else
	{
	    _kid[linenum][i].update_bcf_gt_array(gt_array,i,nullptr);	    
	}
    }
}
