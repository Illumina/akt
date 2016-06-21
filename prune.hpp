#include "akt.hpp"

class circularBuffer
{
    circularBuffer(const string & fname,const string & filter,int window,int buffer);
    ~circularBuffer();
    int next();
    int flush(int nflush);
    int **_gt,_ngt;
    bcf_srs_t *_sr;
    bcf1_t *_line;
    bcf_hdr_t *_hdr_in,*_hdr_out; 
    filter_t *_filter=NULL;
    int _window,_bufsize,_nsnp,_offset,_nsample;
    htsFile *_fout;
};

