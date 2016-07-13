#include "circularBuffer.hpp"
#include "Eigen/Dense"

using namespace Eigen;

void circularBuffer::setOutput(const string & fname,const string & mode)
{
    _write=true;
    _hdr_out=_hdr_in;
    _fout = hts_open(fname.c_str(),mode.c_str());
    _keep.assign(_bufsize,true);
    bcf_hdr_write(_fout,_hdr_out);
}

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


    _maf.assign(_bufsize,-1);

    _nread=_nwrote=0;
}

circularBuffer::~circularBuffer()
{
    if(_write)
    {
	flush(_nsnp);
	hts_close(_fout);
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
    int prev_rid=_line[_offset]->rid;
    while(_nsnp>0 && count<nflush)
    {
	int idx = _offset;
	assert(idx<_bufsize);

	if(_write&&_keep[idx])
	{
	    bcf_write1(_fout, _hdr_out, _line[idx]);
	    _nwrote++;
	}
	if(_line[idx]->rid==prev_rid)
	{
	    assert(_line[idx]->pos>=prev_pos);
	}
	prev_rid=_line[idx]->rid;
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
#pragma omp parallel for
    for(int i=0;i<L;i++)
    {
	int *x=cb.getGT(i);
	for(int j=i;j<L;j++)
	{
	    int *y=cb.getGT(j);
	    if(s[i]>0 && s[j]>0)
	    {
		R2(i,j)=correlation(x,y,frq[i],frq[j],s[i],s[j],N);
		R2(j,i) =  R2(i,j);
	    }
	}
    }    
    return(0);
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
//	assert(sd[i]>0);
	assert(mean[i]>=0&&mean[i]<=2.);
    }
    return(0);
}


