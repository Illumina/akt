#include prune.hpp

circularBuffer::circularBuffer(const string & fname,const string & filter,int window,int buffer)
{
    if(bcf_sr_add_reader(_sr,fname.c_str())!=1)
    {
	die("problem opening "+fname);
    }
    filter_t *_filter=NULL;
    if(!filter.empty())   
    {
	_filter = filter_init(_hdr_in, filter);
    }
    _hdr_in=sr->readers[0].header;
    _nsample =  bcf_hdr_nsamples(_hdr_in);
    assert(window>(2*buffer));
    _window=window;
    _bufsize=window+2*buffer;
    _offset=0;
    _nsnp=0;
    _line = (bcf1_t *)malloc(_bufsize*sizeof(bcf1_t));
    _gt = (int **)malloc(_bufsize*sizeof(int *));
    _ngt=2*_nsample;  
    bcf1_t *line_ptr=_line;
    for(int i=0;i<_bufsize;i++)
    {
	line_ptr = bcf_init1();
	_gt[i] = (int *)malloc(_ngt*sizeof(int));
	_line_ptr++;
    }

    _hdr_out = bcf_hdr_subset(_hdr_in,0,NULL,NULL);
    bcf_hdr_add_sample(_hdr_out, NULL);     
    _fout = hts_open("-", "wv");
}

circularBuffer::~circularBuffer
{
    flush(_nsnp);
    hts_close(_fout);
    bcf_sr_destroy(sr);
    bcf1_t *line_ptr=_line;
    for(int i=0;i<_bufsize;i++)
    {
	free(_gt[i]);
	bcf_destroy(line_ptr);
	line_ptr++;
    }    
    free(_gt);
    free(_line);
}

int circularBuffer::next()
{
    flush(_nsnp - _overlap);

    bcf1_t *ptr = line+_offset+_nsnp;
    while(_nsnp<_bufsize && bcf_sr_next_line(_sr))
    {
	int idx = _offset+_nsnp;
	ptr = bcf_sr_get_line(_sr,0);
	assert(bcf_get_genotypes(_hdr_in, ptr, &_gt[idx], &_ngt)==2*_nsample);
	for(int i=0;i<2*_nsample;i++) 
	{
	    if(bcf_gt_is_missing(_gt[2*i]) || bcf_gt_is_missing(_gt[2*i+1]))
	    {
		_gt[idx][i] = 3;
	    }
	    else
	    {
		_gt[idx][i] = bcf_gt_allele(_gt[2*i]) + bcf_gt_allele(_gt[2*i+1]);
	    }
	}
	ptr++;
	_nsnp++;
    }
    return(_nsnp);
}


int prune_main(int argc,char **argv) {
    assert(argc==3);
    circularBuffer cb(argv[2],"",1000,100);
    while(cb.next()) 
    {
	
    }
    return(0);
}
