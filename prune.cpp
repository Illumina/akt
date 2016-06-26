#include "prune.hpp"

circularBuffer::circularBuffer(const string & fname,const string & filter,int window,int overlap)
{
    _sr =  bcf_sr_init() ; 
    cerr<<"Input: "<<fname<<endl;
    if(bcf_sr_add_reader(_sr,fname.c_str())!=1)
    {
	die("problem opening "+fname);
    }
    filter_t *_filter=NULL;
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
    _fout = hts_open("-", "wb");
    bcf_hdr_write(_fout,_hdr_out);
}

circularBuffer::~circularBuffer()
{
    flush(_nsnp);
    hts_close(_fout);
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
    while(_nsnp>0 && count<nflush)
    {
	int idx = (count+_offset)%_bufsize;
	assert(idx<_bufsize);
//	cerr << _line[idx]->pos+1<<endl;
	if(_keep[idx])
	{
	    bcf_write1(_fout, _hdr_out, _line[idx]);
	}
	_offset++;
	_offset%=_bufsize;
	_nsnp--;
    }
    return(count);
}

int circularBuffer::next()
{
    flush(_nsnp - _overlap);
    assert(_nsnp<_bufsize);
    while(_nsnp<_bufsize && bcf_sr_next_line(_sr))
    {
	int idx = (_offset+_nsnp)%_bufsize;
	_line[idx] = bcf_sr_get_line(_sr,0);
//	cerr<<"idx="<<idx<<" nsnp="<<_nsnp<<" bufsize="<<_bufsize<<endl;
//	cerr<<_line[idx]->pos+1<<endl;
	_keep[idx]=true;
	assert(bcf_get_genotypes(_hdr_in, _line[idx], &_gt[idx], &_ngt)==2*_nsample);
	for(int i=0;i<2*_nsample;i++) 
	{
	    if(bcf_gt_is_missing(_gt[idx][2*i]) || bcf_gt_is_missing(_gt[idx][2*i+1]))
	    {
		_gt[idx][i] = 3;
	    }
	    else
	    {
		_gt[idx][i] = bcf_gt_allele(_gt[idx][2*i]) + bcf_gt_allele(_gt[idx][2*i+1]);
	    }
	}
	_nsnp++;
    }
    return(_nsnp);
}


int prune_main(int argc,char **argv) {
    assert(argc==3);
    circularBuffer cb(argv[2],"",500,500);
    while(cb.next()) 
    {
	
    }
    return(0);
}
