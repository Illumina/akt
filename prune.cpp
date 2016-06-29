#include "prune.hpp"

using namespace Eigen;


circularBuffer::circularBuffer(const string & fname,const string & filter,int window,int overlap)
{
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
    _hdr_out=_hdr_in;
//    _hdr_out = bcf_hdr_subset(_hdr_in,0,NULL,NULL);
//    bcf_hdr_add_sample(_hdr_out, NULL);     
    _keep.assign(_bufsize,true);
    _maf.assign(_bufsize,-1);
    _fout = hts_open("-", "wb");
    bcf_hdr_write(_fout,_hdr_out);
    _nread=_nwrote=0;
}

circularBuffer::~circularBuffer()
{
    flush(_nsnp);
    hts_close(_fout);
    bcf_sr_destroy(_sr);
    for(int i=0;i<_bufsize;i++)
    {
	free(_gt[i]);
//	bcf_destroy(_line[i]);
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
#ifdef DEBUG
	cerr <<"Write: "<<count<< "/"<<nflush<<" "<<idx <<" "<< _line[idx]->pos+1<<endl;
#endif
	if(_keep[idx])
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
	_line[idx] = bcf_dup(bcf_sr_get_line(_sr,0));
	if(	_line[idx]->n_allele==2)
	{
#ifdef DEBUG
	    cerr<<"Read: idx="<<idx<<" nsnp="<<_nsnp<<" bufsize="<<_bufsize<<endl;
	    cerr<<_line[idx]->pos+1<<endl;
#endif
	    _keep[idx]=true;
	    assert(bcf_get_genotypes(_hdr_in, _line[idx], &_gt[idx], &_ngt)==2*_nsample);
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
	float n = 0;
	int *x=cb.getGT(i);
	for(int k=0;k<N;k++)
	{
	    assert(x[k]>=0 && x[k]<=3);
	    if(x[k]<3)
	    {
		mean[i]+=x[k];
		n++;
	    }
	}	
	mean[i] /=n;	    
	for(int k=0;k<N;k++)
	{
	    if(x[k]<3)
	    {
		sd[i] += pow(x[k] - mean[i],2.0);
	    }
	}	
//	sd[i] /= (n-1); //we want the sum-of-squares not the variance.
	sd[i] = sqrt(sd[i]);
//	cerr <<i << " "<< mean[i] << " "<< sd[i]<<endl;
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
    cerr << R2 << endl <<endl;;
    exit(1);
    return(0);
}

int prune(circularBuffer & cb,MatrixXf & R2,int buf,int nkeep)
{
    float r2_thresh = 0.8;
    assert(buf < 2*R2.rows());
    int npruned = 0;
    vector<int> score;
    int start = buf;
    int stop = R2.rows() - buf;
    int width = stop - start;
    assert(nkeep<=width);
    list<int> tags;
    int ntag=0;

    while( ntag<nkeep )
    {
//1. calculate score (number of snps tagged) per SNP
	score.assign(width,0);
	for(int i=start;i<stop;i++)
	{
	    if(!cb.isFiltered(i))
	    {
		for(int j=i-buf;j<(i+buf);j++)
		{
		    if(!cb.isFiltered(j))
		    {
//		score[i-start] += R2(i,j);
			if(R2(i,j)>=r2_thresh)
			{
			    score[i-start]++;
			}
		    }
		}
	    }
	    cerr << i << " "<< score[i-start]<<endl;
	}

//2. find the best tag SNP. if there is a tie, take the more common.
	int maxidx=0;
	for(int i=start;i<stop;i++)
	{
	    if(score[i-start]>score[maxidx])
	    {
		maxidx=i-start;
	    }
	    else if(score[i-start]==score[maxidx] && cb.getMAF(i)>cb.getMAF(maxidx))
	    {
		maxidx=i-start;		    
	    }
	}

//3. set all tagged SNPs to filtered.
	cerr << "maxidx="<<maxidx;
	for(int i=0;i<cb.getNumberOfSNPs();i++)
	{
	    if(R2(i,maxidx+start)>=r2_thresh)
	    {
		cb.setFilter(i);
		cerr << " "<<i;
	    }
	}
	cerr << endl;
	tags.push_back(maxidx);
	ntag++;
    }
    exit(1);
    for(list<int>::iterator it=tags.begin();it!=tags.end();it++)
    {
	cb.setKeep(*it);
    }
    return(0);
}

int prune_main(int argc,char **argv) {
    assert(argc==3);
    int w = 10;
    int b = 10;

    circularBuffer cb(argv[2],"",w,b);
    int chunk=0;
    int nsnp=w+2*b;
    MatrixXf R2 = MatrixXf::Zero(nsnp,nsnp);    
    while(cb.next()) 
    {
	correlationMatrix(cb,R2);
	prune(cb,R2,w,10);
	chunk++;
    }
    cerr << "Retained "<<cb.getNwrote()<<"/"<<cb.getNread()<<" after pruning"<<endl;
    return(0);
}
