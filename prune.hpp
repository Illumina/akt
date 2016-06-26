#include "akt.hpp"

class circularBuffer
{
public:
    circularBuffer(const string & fname,const string & filter,int window,int buffer);
    ~circularBuffer();
    int next();
    int flush(int nflush);
private:
    int **_gt,_ngt;
    bcf_srs_t *_sr;
    bcf1_t **_line;
    bcf_hdr_t *_hdr_in,*_hdr_out; 
    filter_t *_filter=NULL;
    int _nsnp,_offset,_nsample,_overlap,_window,_bufsize;
    htsFile *_fout;
    vector <bool> _keep;
};

