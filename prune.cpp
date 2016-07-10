#include "prune.hpp"
using namespace Eigen;

//#define DEBUG

/**
 * @name    usage
 * @brief   print out options
 *
 * List of input options
 *
 */
static void usage()
{ 
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   akt prune - LD prunes variants\n");
    fprintf(stderr, "Usage:   ./akt prune input.bcf\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -o, --out                   site-only vcf containg pruned variants\n");
    fprintf(stderr, "    -w, --window                window size within to consider possible tag SNPs (default 250)\n");
    fprintf(stderr, "    -b, --buffer                the buffer/step size to move forward\n");
    fprintf(stderr, "    -c, --correlation           minimum r^2 value to consider a SNP tagged (default 0.5)\n");
    fprintf(stderr, "    -n, --nthread               number of threads to use\n");
    exit(1);
}



int prune_main(int argc,char **argv) 
{
    int c;
    float r2_thresh=0.5;
    int w = 1000;
    int b = 1000;  
    string outfile="-";
    int nthread=1;
    string mode = "v";
    if(argc<3) usage();
    static struct option loptions[] =    {
	{"output",1,0,'o'},	
	{"output-type",1,0,'O'},	
	{"correlation",1,0,'c'},	
	{"window",1,0,'w'},	
	{"buffer",1,0,'b'},	
	{"nthread",1,0,'n'},	
	{0,0,0,0}
    };

    while ((c = getopt_long(argc, argv, "b:n:o:w:c:O:",loptions,NULL)) >= 0) 
    {  
	switch (c)      {
	case 'o': outfile = optarg; break;
	case 'O': mode = optarg; break;
	case 'c': r2_thresh = atof(optarg); break;
	case 'w': w = atoi(optarg); break;
	case 'b': b = atoi(optarg); break;
	case 'n': nthread = atoi(optarg); break;
	default: cerr << "Unknown argument:"+(string)optarg+"\n" << endl; exit(1);
	}
    }

    if(w<=0)
    {
	die("-w must be greater than 0");
    }
    if(r2_thresh<=0)
    {
	die("-c must be greater than 0");
    }

    optind++;
    if(optind==argc)
    {
	die("no input file provided");
    }
    circularBuffer cb(argv[optind],"",w,b);
    cb.setOutput(outfile,("w"+mode));
    int chunk=0,nread=0;
    int nsnp=w + 2*b;
    MatrixXf R2 = MatrixXf::Zero(nsnp,nsnp);    
    vector<vector<int > > G;
    //create a 0-sample header for writing.
    bcf_hdr_t *hdr_in = cb.getHeader();
    omp_set_num_threads(nthread);

    vector<float> frq,s;
    int nsample = cb.getNumberOfSamples();
    int npruned=0;
    while(cb.next())
    {
	cerr <<"Read "<<cb.getNread()<<" markers\r";
	int nsnp=cb.getNumberOfSNPs();
	int start=b;
	if(chunk==0) start=0;
	int stop =w+b;
	if(nsnp<(w+2*b)) stop=nsnp;
	meanvar(cb,frq,s);
	for(int i=start;i<stop;i++)
	{
	    if(!cb.isFiltered(i))
	    {
		int *x=cb.getGT(i);
		float x_maf=cb.getMAF(i);
#pragma omp parallel for
		for(int j=max(0,i-b);j<min(nsnp,i+b);j++)
		{
		    int *y=cb.getGT(j);
		    if(!cb.isFiltered(j)&&i!=j)
		    {
			float r2 = correlation(x,y,frq[i],frq[j],s[i],s[j],nsample);
			if(r2>=r2_thresh)
			{
//			    cerr << i << " " << j << " "<< r2<<endl;
			    if( x_maf>cb.getMAF(j))
			    {
				cb.setFilter(i);
			    }
			    else
			    {
				cb.setFilter(j);
			    }
			    npruned++;
			}
		    }		    
		}
	    }
	}
	chunk++;
	nread+=w;
    }
    cerr << npruned<<"/"<<cb.getNread()<<" variants pruned"<<endl;
    return(0);
}
