#include "prune.hpp"

using namespace Eigen;

//#define DEBUG

circularBuffer::circularBuffer(const string & fname,const string & filter,int window,int overlap)
{
    _write=false;
    _sr =  bcf_sr_init() ; 
    cerr<<"Input: "<<fname<<endl;
    if(bcf_sr_add_reader(_sr,fname.c_str())!=1)
    {
	die("problem opening "+fname);
    }
    _filter=NULL;
    if(!filter.empty())   
    {
	_filter = filter_init(_hdr_in, filter.c_str());
    }
    _hdr_in=_sr->readers[0].header;
    _nsample =  bcf_hdr_nsamples(_hdr_in);
    _window=window;
    _overlap=overlap;
    _bufsize=window+2*_overlap;
    _offset=0;
    _nsnp=0;
    _line = new bcf1_t*[_bufsize];
    _gt = (int **)malloc(_bufsize*sizeof(int *));
    _ngt=2*_nsample;  
    for(int i=0;i<_bufsize;i++)
    {
	_line[i] = bcf_init1();
	_gt[i] = (int *)malloc(_ngt*sizeof(int));
    }

//    _hdr_out = bcf_hdr_subset(_hdr_in,0,NULL,NULL);
//    bcf_hdr_add_sample(_hdr_out, NULL);     

    _maf.assign(_bufsize,-1);
    if(_write)
    {
	_hdr_out=_hdr_in;
	_fout = hts_open("-", "wb");
	_keep.assign(_bufsize,true);
	bcf_hdr_write(_fout,_hdr_out);
    }
    _nread=_nwrote=0;
}

circularBuffer::~circularBuffer()
{
    if(_write)
    {
	flush(_nsnp);
	hts_close(_fout);
	bcf_hdr_destroy(_hdr_out);
    }
    bcf_sr_destroy(_sr);

    for(int i=0;i<_bufsize;i++)
    {
	free(_gt[i]);
	bcf_destroy(_line[i]);
    }    
    free(_gt);
    delete[] _line;
}

int circularBuffer::flush(int nflush)
{
    int count=0;
    int prev_pos=0;
    while(_nsnp>0 && count<nflush)
    {
	int idx = _offset;
	assert(idx<_bufsize);

	if(_write&&_keep[idx])
	{
	    bcf_write1(_fout, _hdr_out, _line[idx]);
	    _nwrote++;
	}
	assert(_line[idx]->pos>=prev_pos);
	prev_pos=_line[idx]->pos;
	_offset++;
	_offset%=_bufsize;
	_nsnp--;
	count++;
    }
    return(count);
}

int circularBuffer::next()
{
    flush(_overlap);
    assert(_nsnp<_bufsize);
    int count=0;
    while(_nsnp<_bufsize && bcf_sr_next_line(_sr))
    {
	int idx = (_offset+_nsnp)%_bufsize;
	bcf_copy(_line[idx],bcf_sr_get_line(_sr,0));
	bcf_unpack(_line[idx],BCF_UN_ALL);
	if(_line[idx]->n_allele==2)
	{
	    if(_write)
	    {	    
		_keep[idx]=true;
	    }
	    int *test=NULL,ntest=0;
	    if(bcf_get_genotypes(_hdr_in, _line[idx], &_gt[idx], &_ngt)!=2*_nsample)
	    {
		die("problem reading genotypes ngt != 2*n");
	    }
	    _maf[idx]=0.0;
	    for(int i=0;i<_nsample;i++) 
	    {
		if(bcf_gt_is_missing(_gt[idx][2*i]) || bcf_gt_is_missing(_gt[idx][2*i+1]))
		{
		    _gt[idx][i] = 3;
		}
		else
		{
		    _gt[idx][i] = bcf_gt_allele(_gt[idx][2*i]) + bcf_gt_allele(_gt[idx][2*i+1]);
		    _maf[idx]+=_gt[idx][i];
		}
	    }
	    _maf[idx]/=2.;
	    _nsnp++;
	    count++;
	}
	_nread++;	    
    }	
    return(count);
}

int meanvar(circularBuffer & cb, vector<float> & mean, vector<float> & sd)
{
    int L = cb.getNumberOfSNPs();
    int N = cb.getNumberOfSamples();
    mean.assign(L,0.0);
    sd.assign(L,0.0);

    for(int i=0;i<L;i++)
    {
	int ac=0,n=0;
	int *x=cb.getGT(i);
	for(int k=0;k<N;k++)
	{
	    assert(x[k]>=0 && x[k]<=3);
	    if(x[k]<3)
	    {
		ac+=x[k];
		n++;
	    }
	}
	if(n==0) die("no data");
	mean[i] = (float)ac/(float)n;	    
	if(ac>0 && ac<(2*n) )
	{
	    for(int k=0;k<N;k++)
	    {
		if(x[k]<3)
		{
		    sd[i] += pow(x[k] - mean[i],2.0);
		}
	    }	
	}
	else
	{
	    sd[i]=0.0;
	}
//	sd[i] /= (n-1); //we want the sum-of-squares not the variance.
	sd[i] = sqrt(sd[i]);
//	cerr <<i << " "<< mean[i] << " "<< sd[i]<< " n="<<n<<endl;
	assert(sd[i]>0);
	assert(mean[i]>=0&&mean[i]<=2.);
    }
    return(0);
}

float correlation(int *x,int *y, float mx, float my, float sx,float sy,int N)
{
    float lookup[4][4];
    for(int g1=0;g1<4;g1++)
    {
	for(int g2=0;g2<4;g2++)
	{
	    if(g1<3&&g2<3)
	    {
		lookup[g1][g2]=((float)g1 - mx)*((float)g2 - my) ;
	    }
	    else
	    {
		lookup[g1][g2]=0.0;
	    }
	}
    }
    float r2=0.0;
    for(int k=0;k<N;k++)
    {
	if(x[k]<3&&y[k]<3)
	{
	    r2+=lookup[x[k]][y[k]];
	}
    }
    r2 /=  (sx*sy);
    r2  *=   r2;
    return(r2);
}

int correlationMatrix(circularBuffer & cb,MatrixXf & R2)
{
    int L = cb.getNumberOfSNPs();
    assert(R2.rows() == R2.cols());    
    assert(L<=R2.rows());
    int N = cb.getNumberOfSamples();
    assert(L>0);
    assert(N>0);
    vector<float> frq,s;
    meanvar(cb,frq,s);
#pragma omp parallel
    for(int i=0;i<L;i++)
    {
	int *x=cb.getGT(i);
	for(int j=i;j<L;j++)
	{
	    int *y=cb.getGT(j);
	    R2(i,j)=correlation(x,y,frq[i],frq[j],s[i],s[j],N);
	    R2(j,i) =  R2(i,j);
	}
    }    
    return(0);
}

int prune(vector<list<int > > & R2,vector<float> maf,int window,vector<int> & keep,int nkeep)
{
    float r2_thresh = 0.8;
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

//1. find the best tag SNP. if there is a tie, take the more common.
	int maxidx=0;
	int maxscore = -1;
	ntag=0;
	for(size_t i=0;i<R2.size();i++)
	{
	    if(!tagged[i])
	    {
		int score = R2[i].size();
		for(list<int>::iterator it=R2[i].begin();it!=R2[i].end();it++)
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

	    for(list<int>::iterator it=R2[maxidx].begin();it!=R2[maxidx].end();it++)
	    {
		tagged[*it] = true;
	    }
	    assert(tagged[maxidx]);
	    keep.push_back(maxidx);
	    nselect++;
	}
    }
    sort(keep.begin(),keep.end());
    
    return(0);
}

int prune_main(int argc,char **argv) {
    assert(argc==3);
    int w = 250;
    int b = 250;
    float r2_thresh=0.8;
    circularBuffer cb(argv[2],"",w,b);
    int chunk=0;
    int nsnp=w+2*b;
    MatrixXf R2 = MatrixXf::Zero(nsnp,nsnp);    
    vector<list<int > > G;
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
    omp_set_num_threads(24);
    while(cb.next()) 
    {
	cerr << "count="<<count<<endl;
	int nsnp=cb.getNumberOfSNPs();
	correlationMatrix(cb,R2);
	int start=b;
	int stop =b+w;
	if(chunk==0)
	{
	    start=0;
	}
	if(nsnp<(w+2*b)) 
	{
	    stop=nsnp;
	}
	for(int i=start;i<stop;i++)
	{
	    G.push_back(list<int>());
	    maf.push_back(cb.getMAF(i));
	    lines.push_back(bcf_dup(cb.getLine(i)));
	    bcf_subset(hdr_in,lines.back(),0,NULL);
	    for(int j=max(0,i-b);j<min(nsnp,i+b);j++)
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
    cerr << "buffered "<<G.size() << " SNPs for pruning"<<endl;


#ifdef DEBUG
    for(size_t i=0;i<G.size();i++)
    {
	cerr << i << ": " << G[i].size() << ": ";
	for(list<int>::iterator it=G[i].begin();it!=G[i].end();it++)
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

    vector<int> keep;
    prune(G,maf,2*b,keep,1000);
    
    htsFile *fout = hts_open("-", "wv");    
    bcf_hdr_write(fout,hdr_out);
    cerr << "keeping "<<keep.size()<<"/"<<lines.size()<<endl;
    int last_pos=-1;

    for(vector<int>::iterator it=keep.begin();it!=keep.end();it++)
    {
#ifdef DEBUG
	cerr << *it << "/"<<lines.size()<<" "<<lines[*it]->pos<<endl;
#endif
	if(lines[*it]->pos <= last_pos)
	{
	    die("duplicated positions");
	}
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
