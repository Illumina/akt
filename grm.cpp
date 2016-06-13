#include "akt.hpp"

//#define DEBUG 

static void usage()
{
    cerr << "\nAbout: Calculate genetic relationship matrix from a BCF/VCF" << endl;	
    cerr << "Usage: akt grm [options] <in.bcf>" << endl;
    cerr<< "    -r, --regions <region>        restrict to comma-separated list of regions"<<endl;
    cerr<< "    -R, --regions-file <file>     restrict to regions listed in a file"<<endl;
    cerr<< "    -s, --samples [^]<list>       comma separated list of samples to include (or exclude with \"^\" prefix)"<<endl;
    cerr<< "    -S, --samples-file [^]<file>  file of samples to include (or exclude with \"^\" prefix)"<<endl;
    cerr<< "    -n                            number of threads to use for calculations."<<endl;
    exit(1);
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
	{"nthreads",1,0,'n'},
	{0,0,0,0}
    };
    string regions="",input="",samples="";
    int regions_is_file=1;
    int samples_is_file=1;
    int nthreads=1;
    while ((c = getopt_long(argc, argv, "r:R:s:S:n:",loptions,NULL)) >= 0 )
    {
	switch (c)
	{	    
	case 'r': regions = (optarg); regions_is_file=0; break;
	case 'R': regions = (optarg); break;
	case 's': samples = (optarg); samples_is_file=0; break;
	case 'S': samples = (optarg); break;	    
	case 'n': nthreads = atoi(optarg); break;
	}
    }
    optind++;
    if(optind>=argc) 
    {
	die("no input provided");
    }

    if(nthreads < 1)
    { 
	nthreads = 1; 
    }
    omp_set_num_threads(nthreads);
#pragma omp parallel
    {
	if(omp_get_thread_num() == 0)
	{
	    if( omp_get_num_threads() != 1)
	    {
		cerr << "Using " << omp_get_num_threads() << " threads" << endl;
		nthreads = omp_get_num_threads();
	    }
	}
    }



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
    bcf_hdr_t *hdr=sr->readers[0].header;
    if( !samples.empty() )
    {
	bcf_hdr_set_samples(hdr, samples.c_str(), samples_is_file);
    }
    int nsample = bcf_hdr_nsamples(hdr);
    cerr <<nsample<<" samples in "<<argv[optind]<<endl;
    vector<float> grm(nsample*(nsample - 1)/2 + nsample,0.0);
    vector<int> nmiss(nsample*(nsample - 1)/2 + nsample,0);
    bcf1_t *line;
    int ngt=2*nsample;
    int *gt_arr=new int[ngt];
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

#ifdef DEBUG
	    cerr << line->rid<<":"<<line->pos+1<<" ac="<<ac<<" an="<<an<<" af="<<af<<endl;
#endif
	    float lookup[4][4];		
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
#pragma omp parallel for	
		for(int j1=0;j1<nsample;j1++) 
		{
		    int offset=j1*nsample - j1*(j1+1)/2;
		    for(int j2=j1;j2<nsample;j2++) 
		    {
			int idx = offset+j2;
#ifdef DEBUG
			assert( idx>=0 && idx<nmiss.size());
			assert( idx>=0 && idx<grm.size());
			assert( gt_arr[j1]>=0 && gt_arr[j1]<4);
			assert( gt_arr[j2]>=0 && gt_arr[j2]<4);
#endif
			int g0=gt_arr[j1];
			int g1=gt_arr[j2];
			if(g0<3&&g1<3)
			{
			    grm[idx]+=lookup[g0][g1];
			}
			else
			{
			    nmiss[idx]++;
			}
		    }
		}
		nsnp++;
	    }
	}
    }

    /////code for lower triangular matrix
    // for(int j1=0;j1<nsample;j1++) 
    // {
    // 	for(int j2=0;j2<=j1;j2++) 
    // 	{
    // 	    int idx = j2*nsample - j2*(j2+1)/2 + j1;
    for(int j1=0;j1<nsample;j1++) 
    {
	int offset=j1*nsample - j1*(j1+1)/2;
	for(int j2=j1;j2<nsample;j2++) 
	{
	    int idx = offset+j2;
#ifdef DEBUG
	    assert(idx>=0 && idx < nmiss.size());
	    assert(idx>=0 && idx < grm.size());
#endif
	    grm[idx]/=(nsnp-nmiss[idx]);
	    cout<<sr->readers[0].header->samples[j1]<<"\t"<<sr->readers[0].header->samples[j2]<<"\t"<<grm[idx]<<endl;
	}
    }
    cerr << nsnp << " markers used"<<endl;
    cerr << "done."<<endl;
    return(0);
}
