#include "akt.hpp"
#include "Eigen/Dense"

class circularBuffer
{
public:
    circularBuffer(const string & fname,const string & filter,int window,int buffer);
    ~circularBuffer();
    int next();
    int flush(int nflush);
    int getNumberOfSNPs() {return (_nsnp) ; }
    int *getGT(int i) { return(_gt[(_offset+i)%_bufsize]); }
    float getMAF(int i) { return(_maf[(_offset+i)%_bufsize]); }
    int  getNumberOfSamples() { return( _nsample ); }
    bool isFull() { return( _nsnp==_bufsize); }
    int setFilter(int i) { _keep[(_offset+i)%_bufsize]=false; }
    int setKeep(int i) { _keep[(_offset+i)%_bufsize]=true; }
    bool isFiltered(int i) { return(!_keep[(_offset+i)%_bufsize]); }
    int getNread() {return(_nread);}
    int getNwrote() {return(_nwrote);};

private:
    int **_gt,_ngt;
    bcf_srs_t *_sr;
    bcf1_t **_line;
    bcf_hdr_t *_hdr_in,*_hdr_out; 
    filter_t *_filter;
    int _nsnp,_offset,_nsample,_overlap,_window,_bufsize;
    htsFile *_fout;
    vector <bool> _keep;
    int _nread,_nwrote;
    vector<float> _maf;
};

