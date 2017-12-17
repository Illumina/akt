#include "gtest/gtest.h"
#include "akt.hh"
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
    HaplotypeBuffer hap_ps(num_sample,pedigree);    
    bcf1_t *line;
    int *gt=nullptr,ngt=0;
    int *ps=nullptr,nps=0;    
    while (bcf_sr_next_line(bcf_reader))
    {
	line = bcf_sr_get_line(bcf_reader, 0);
	bcf_unpack(line, BCF_UN_ALL);
	ASSERT_EQ(bcf_get_genotypes(in_header, line, &gt, &ngt),2*num_sample);
	int status = bcf_get_format_int32(in_header, line,"PS", &ps, &nps);
	hap.push_back(gt);
	if(status<=0)
	    hap_ps.push_back(gt,nullptr);
	else if(status==num_sample)
	    hap_ps.push_back(gt,ps);
	else
	    die("invalid PS = "+std::to_string(status));
    }
    hap.phase();
    for(int i=0;i<hap.get_num_variant();i++)
    {
	for(int j=0;j<4;j++)  std::cerr << hap.get_genotype(i,j).print() <<"\t";
	std::cerr<<std::endl;
	ASSERT_FALSE(is_mendel_inconsistent(hap.get_genotype(i,2),hap.get_genotype(i,0),hap.get_genotype(i,1)));
    }
    std::cerr<<std::endl;
    hap.align(hap_ps);
    for(int i=0;i<hap.get_num_variant();i++)
    {
	for(int j=0;j<4;j++)  std::cerr << hap.get_genotype(i,j).print() <<"\t";
	std::cerr<<std::endl;
	ASSERT_FALSE(is_mendel_inconsistent(hap.get_genotype(i,2),hap.get_genotype(i,0),hap.get_genotype(i,1)));
    }
    
    free(gt);    
}
