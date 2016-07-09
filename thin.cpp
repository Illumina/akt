#include "thin.hpp"
#include <time.h>

using namespace Eigen;

//#define DEBUG

/**
 * @name    usage
 * @brief   print out options
 *
 * List of input options
 *
 */
static void usage(){ 
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   akt thin - selects a set of tag SNPs\n");
    fprintf(stderr, "Usage:   ./akt thin input.bcf\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -o, --out                   site only vcf with tag SNPs\n");
    fprintf(stderr, "    -w, --window                window size within to consider possible tag SNPs (default 250)\n");
    fprintf(stderr, "    -c, --correlation           minimum r^2 value to consider a SNP tagged (default 0.8)\n");
    fprintf(stderr, "    -k, --nkeep                 number of tag SNPs to selection (default 20000)\n");
    fprintf(stderr, "    -n, --nthread               number of threads to use\n");
    exit(1);
}



int selectTagSNPs(vector<vector<int > > & R2,vector<float> maf,int window,vector<int> & keep,int nkeep,float r2_thresh)
{
    keep.clear();
    assert(window>0);
    assert(nkeep>0);
    int ntag=0,nselect=0;
    int b = window/2;
    int nsnp=R2.size();
    vector<bool> tagged(nsnp,false);    
    keep.clear();
    while( nselect<nkeep && ntag<nsnp)
    {
	if(nselect%100==0)
	{
	    cerr <<"Selected "<<nselect<<" tag markers\r";
	}
//1. find the best tag SNP. if there is a tie, take the more common.
	int maxidx=0;
	int maxscore = -1;
	ntag=0;
	for(size_t i=0;i<R2.size();i++)
	{
	    if(!tagged[i])
	    {
		int score = R2[i].size();
		if(score>=maxscore)
		{
		    for(vector<int>::iterator it=R2[i].begin();it!=R2[i].end();it++)
		    {
			assert(*it>=0);
			assert(*it<nsnp);
			if(tagged[*it])
			{
			    score--;
			}
		    }
		    assert(score>=0);
		    if(score==maxscore && maf[i]>maf[maxidx])
		    {
			maxscore = score;
			maxidx=i;
		    }
		    else if(score > maxscore)
		    {
			maxscore = score;
			maxidx=i;
		    }
//		cerr << i <<" "<< score<<" " << R2[i].size()<<endl;
		}
	    }
	    else
	    {
		ntag++;
	    }
	}

//2. update list of tagged SNPs (will no longer be counted in scores).
#ifdef DEBUG
	cerr << "maxidx="<<maxidx<<" maxscore="<<maxscore<<" ntag="<<ntag<<endl;
#endif
	assert(ntag==nsnp || maxscore>=1);
	if(maxscore>=1)
	{
	    if(tagged[maxidx])
	    {
		die("duplicated tag SNP selected");
	    }

	    for(vector<int>::iterator it=R2[maxidx].begin();it!=R2[maxidx].end();it++)
	    {
		tagged[*it] = true;
	    }
	    assert(tagged[maxidx]);
	    keep.push_back(maxidx);
	    nselect++;
	}
    }
    sort(keep.begin(),keep.end());
    cerr <<endl<<"Finshed tag selection"<<endl;
    return(0);
}

int thin_main(int argc,char **argv) 
{
    int nkeep=20000;
    int c;
    float r2_thresh=0.5;
    int w = 250;
    char *outfile=NULL;
    int nthread=1;
    if(argc<3) usage();
    static struct option loptions[] =    {
	{"output",1,0,'o'},	
	{"correlation",1,0,'c'},	
	{"window",1,0,'w'},	
	{"nkeep",1,0,'k'},	
	{"nthread",1,0,'n'},	
	{0,0,0,0}
    };

    while ((c = getopt_long(argc, argv, "k:n:o:w:c:",loptions,NULL)) >= 0) 
    {  
	switch (c)      {
	case 'o': outfile = optarg; break;
	case 'c': r2_thresh = atof(optarg); break;
	case 'w': w = atoi(optarg); break;
	case 'k': nkeep = atoi(optarg); break;
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
    if(nkeep<=0)
    {
	die("-c must be greater than 0");
    }

    optind++;
    circularBuffer cb(argv[optind],"",w,w);
    int chunk=0;
    int nsnp=3*w;
    MatrixXf R2 = MatrixXf::Zero(nsnp,nsnp);    
    vector<vector<int > > G;
    G.reserve(1e6);//might be pretty big. ~7million SNPs with MAF>=5% in 1000G.
    int count=0;
    int nread;
    vector<float> maf;
    maf.reserve(1e6);
    vector<bcf1_t *> lines;
    lines.reserve(1e6);
    //create a 0-sample header for writing.
    bcf_hdr_t *hdr_in = cb.getHeader();
    bcf_hdr_t *hdr_out = bcf_hdr_subset(cb.getHeader(),0,NULL,NULL);
    bcf_hdr_add_sample(hdr_out, NULL);
    omp_set_num_threads(nthread);

    clock_t wall0 = clock();
    while(cb.next()) 
    {
	cerr <<"Read "<<count<<" markers\r";
	int nsnp=cb.getNumberOfSNPs();
	correlationMatrix(cb,R2);
	int start=w;
	int stop =2*w;
	if(chunk==0)
	{
	    start=0;
	}
	if(nsnp<(3*w)) 
	{
	    stop=nsnp;
	}
	for(int i=start;i<stop;i++)
	{
	    G.push_back(vector<int>());
	    maf.push_back(cb.getMAF(i));
	    lines.push_back(bcf_dup(cb.getLine(i)));
	    bcf_subset(hdr_in,lines.back(),0,NULL);
	    for(int j=max(0,i-w);j<min(nsnp,i+w);j++)
	    {
		if(R2(i,j)>=r2_thresh)
		{
		    G.back().push_back(count+j);
		}
	    }
	    assert(G.back().size()>0);
	}
	count+=w;
	chunk++;
    }
    cerr << "Read "<<G.size() << " SNPs for pruning"<<endl;
    cerr << "Reading took "<<(float)(clock()-wall0)/CLOCKS_PER_SEC<<" seconds"<<endl;

#ifdef DEBUG
    for(size_t i=0;i<G.size();i++)
    {
	cerr << i << ": " << G[i].size() << ": ";
	for(vector<int>::iterator it=G[i].begin();it!=G[i].end();it++)
	{
	    cerr << *it << " ";
	}
	cerr << endl;
    }
#endif

    for(size_t i=0;i<G.size();i++)
    {
	if(G[i].size()==0)
	{
	    cerr << "SNP "<<i<<" did not tag anything (including itself)"<<endl;
	    exit(1);
	}
    }

    wall0 = clock();
    vector<int> keep;
    selectTagSNPs(G,maf,w,keep,nkeep,r2_thresh);
    cerr << "tagging took "<<(float)(clock()-wall0)/CLOCKS_PER_SEC<<" seconds"<<endl;
    
    htsFile *fout;
    if(outfile==NULL)
    {
	fout=hts_open("-", "wv");    
    }
    else
    {
	fout=hts_open(outfile, "wz");    
    }
    
    bcf_hdr_write(fout,hdr_out);
    cerr << "Kept "<<keep.size()<<"/"<<lines.size()<<" markers after tag SNP selection"<<endl;
    int last_pos=-1;
    int prev_rid=-1;
    for(vector<int>::iterator it=keep.begin();it!=keep.end();it++)
    {
#ifdef DEBUG
	cerr << *it << "/"<<lines.size()<<" "<<lines[*it]->pos<<endl;
#endif
	if(lines[*it]->rid==prev_rid && lines[*it]->pos <= last_pos)
	{
	    die("duplicated positions");
	}
	prev_rid = lines[*it]->rid;
	last_pos = lines[*it]->pos;
	bcf_write1(fout,hdr_out,lines[*it]);
    }

    //cleanup
    for(vector<bcf1_t *>::iterator it=lines.begin();it!=lines.end();it++)
    {
	bcf_destroy(*it);
    }
    hts_close(fout);
    bcf_hdr_destroy(hdr_out);
    return(0);
}
