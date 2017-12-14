#include "gtest/gtest.h"
#include "akt.hpp"
#include "pedigree.h"
#include "pedphase.h"

TEST(pedphase,HaplotypeBufferPhaseByTransmission)
{
    bcf_srs_t *bcf_reader = bcf_sr_init(); 
    bcf_sr_add_reader(bcf_reader,"test/pedphase/test16.vcf.gz");
    bcf_hdr_t *in_header = bcf_reader->readers[0].header;
    sampleInfo *pedigree = new sampleInfo(in_header);
    int  num_sample = bcf_hdr_nsamples(in_header);
    HaplotypeBuffer hap(num_sample,pedigree);
    bcf1_t *line;
    int *gt=nullptr,ngt=0;
    while (bcf_sr_next_line(bcf_reader))
    {
	line = bcf_sr_get_line(bcf_reader, 0);
	ASSERT_EQ(bcf_get_genotypes(in_header, line, &gt, &ngt),2*num_sample);
    }
    hap.phase();
    for(int i=0;i<hap.get_num_variants();i++)
	ASSERT_FALSE(is_mendel_inconsistent(hap.get_genotype(i,2),hap.get_genotype(i,0),hap.get_genotype(i,1)));

    free(gt);    
}
