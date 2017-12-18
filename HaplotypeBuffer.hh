#ifndef AKT_HAPLOTYPEBUFFER_H
#define AKT_HAPLOTYPEBUFFER_H

#include "akt.hh"
#include "pedigree.hh"
#include "pedphase.hh"
#include "Genotype.hh"

//This algorithm to harmonise read-back phased variants with pedigree transmission.
//1. Perform phase-by-transmission, temporarily ignoring all read-back-phasing information.
//2. Count the number of times a read-back phased variant agrees with the pedigree phasing (per sample).
//3. If the count in step 2 equals <50% of pedigree resolved hets, flip all read-back phased alleles in the sample.
//4. Calculate the concordance of the flipped alleles with pedigree inheritance.
//5. If the value from 4 is 100%, move FORMAT/PS to FORMAT/RPS to indicate this phase set is fully in agreement with the pedigree.
class HaplotypeBuffer
{
public:
    HaplotypeBuffer(size_t num_sample,sampleInfo *pedigree);
    void push_back(int32_t *gt_array, int32_t *ps_array=nullptr);
    void clear();
    void phase();
    void align(HaplotypeBuffer & haps_to_align);
    Genotype get_genotype(size_t variant_index,size_t sample_index);
    int get_num_variant() {return _num_variant;};
    int get_num_sample() {return _num_sample;};    
    void update_bcf1_genotypes(size_t linenum,int32_t *gt_array, int32_t *ps_array,int32_t *rps_array);
    bool is_mendel_consistent(size_t linenum);
    void swap(int variant,int sample);
    void setPhase(int variant,int sample,bool phase);
    void copy_from_parents();
    void align_sample(vector< vector< Genotype > > & dst,vector< vector< Genotype > > & src,vector< map<int,pair<int,int> > >  & phase_set_vote    );
				   
    bool is_phase_set_aligned_with_pedigree(int sample,int phaseset);
    vector< vector< Genotype > > & kid() {return _kid;};
    vector< vector< Genotype > > & dad() {return _dad;};
    vector< vector< Genotype > > & mum() {return _mum;};
    bool is_sample_phased(int variant,int sample,vector< vector< Genotype > > & genotypes,vector< map<int,pair<int,int> > >  & phase_set_vote);    
private:
    
    void check_pedigree_aligned();    
    size_t _num_sample,_num_variant;
    vector< vector< Genotype > > _kid,_dad,_mum;
    sampleInfo *_pedigree;
    vector<int> _index_of_first_child;
    vector< vector<bool> >_sample_was_mendel_phased;
    vector<bool> _line_is_mendel_consistent;
    vector< vector<bool> > _is_aligned_with_pedigree;
    vector< map<int,pair<int,int> > > _kid_vote,_dad_vote,_mum_vote;
};

#endif //
