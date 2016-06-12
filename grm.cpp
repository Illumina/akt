#include "akt.hpp"

//#define DEBUG 

static void usage()
{
    cerr << "\nAbout: Calculate genetic relation matrix from a BCF/VCF" << endl;	
    cerr << "Usage: akt grm [options] <in.bcf>" << endl;
    cerr << "Expects input.bcf to contain genotypes." << endl;
}

int grm_main(int argc, char** argv)
{
    int c;
    if(argc<3) usage();
    static struct option loptions[] =    {
	{"regions",1,0,'r'},
	{"regions-file",1,0,'R'},
	{"samples",1,0,'s'},
	{"samples-file",1,0,'S'},
	{0,0,0,0}
    };
    string regions="",input="",samples="";
    int regions_is_file=1;
    int samples_is_file=1;
    while ((c = getopt_long(argc, argv, "r:R:s:S",loptions,NULL)) >= 0 )
    {
	switch (c)
	{	    
	case 'r': regions = (optarg); regions_is_file=0; break;
	case 'R': regions = (optarg); break;
	case 's': samples = (optarg); samples_is_file=0; break;
	case 'S': samples = (optarg); break;	    
	}
    }
    optind++;
    if(optind>=argc) die("no input provided");
    bcf_srs_t *sr =  bcf_sr_init() ; ///htslib synced reader.
    if(!regions.empty())
    {
	if( bcf_sr_set_regions(sr, regions.c_str(), regions_is_file)<0 )
	{
	    die("problem setting regions");
	}
    }
    
    if(!(bcf_sr_add_reader (sr, argv[optind])))
    {
	die("problem opening"+(string)argv[optind]);
    }
    int nsample = bcf_hdr_nsamples(sr->readers[0].header);
    cerr <<nsample<<" samples in "<<argv[optind]<<endl;
    float lookup[3][3];
    vector<float> grm(nsample*nsample,0.0);
//    grm.assign(nsample*(nsample - 1)/2 + nsample,0.0);

    bcf1_t *line;
    int *gt_arr=NULL;
    int ngt=0;
    int nsnp=0;
    while(bcf_sr_next_line (sr))  ///read file
    {
	line =  bcf_sr_get_line(sr, 0);
	if(line->n_allele==2) {
	    ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt);  
	    assert(ngt==2*nsample);
	    int ac=0,an=0;
	    for(int i=0;i<nsample;i++)
	    {
		if(bcf_gt_is_missing(gt_arr[2*i])||bcf_gt_is_missing(gt_arr[2*i+1]))
		{
		    gt_arr[i] = 3;
		}
		else
		{
		    gt_arr[i]=bcf_gt_allele(gt_arr[2*i])+bcf_gt_allele(gt_arr[2*i+1]);
		    an+=2;
		    ac+=gt_arr[i];
		}
	    }
	    float af = (float)ac/(float)an;//allele frequency.
//	    if(af>0.5) af = 1-af;
#ifdef DEBUG
	    cerr << line->rid<<":"<<line->pos+1<<" ac="<<ac<<" an="<<an<<" af="<<af<<endl;
#endif
	    if(ac>0 && ac<an)
	    {
		for(int g0=0;g0<4;g0++) 
		{
		    for(int g1=0;g1<4;g1++)
		    {
			if(g0<3&&g1<3)
			{
			    float m = 2*af;
			    float v = 2*af*(1-af);
			    lookup[g0][g1]=((float)g0-m)*((float)g1-m)/v;
			}
			else
			{
			    lookup[g0][g1]=0;
			}
		    }
		}
//#pragma omp parallel for	
		    for(int j1=0;j1<nsample;j1++) 
		    {
			int offset=j1*nsample - j1*(j1+1)/2;
			for(int j2=j1;j2<nsample;j2++) 
			{
			    int idx = j1*nsample+j2;//(offset+j2);
#ifdef DEBUG
			    assert( idx < grm.size());
			    assert( gt_arr[j1]>=0 && gt_arr[j1]<4);
			    assert( gt_arr[j2]>=0 && gt_arr[j2]<4);
#endif
			    grm[idx]+=lookup[gt_arr[j1]][gt_arr[j2]];
			}
		    }
		    nsnp++;
	    }
	}
    }
    for(int j1=0;j1<nsample;j1++) 
    {
	int offset=j1*nsample - j1*(j1+1)/2;
	for(int j2=j1;j2<nsample;j2++) 
	{
	    int idx = j1*nsample+j2;// offset+j2;
	    assert( idx < grm.size());
	    grm[idx]/=nsnp;
	    cout<<idx<<" "<<sr->readers[0].header->samples[j1]<<"\t"<<sr->readers[0].header->samples[j2]<<"\t"<<grm[idx]<<endl;
	}
    }
    cerr << nsnp << " markers used"<<endl;
    cerr << "done."<<endl;
    return(0);
}
