#include "akt.hpp"
#include "Eigen/Dense"

class circularBuffer
{
public:
    circularBuffer(const string & fname,const string & filter,int window,int buffer);
    ~circularBuffer();
    int next();
    int flush(int nflush);
    
    //get stuff
    bcf_hdr_t *getHeader() {return ( _hdr_in); }
    int getNumberOfSNPs() {return (_nsnp) ; }
    int *getGT(int i) { return(_gt[(_offset+i)%_bufsize]); }
    float getMAF(int i) { return(_maf[(_offset+i)%_bufsize]); }
    int  getNumberOfSamples() { return( _nsample ); }
    int getNread() {return(_nread);}
    int getNwrote() {return(_nwrote);};
    bcf1_t *getLine(int i) {return( _line[(_offset+i)%_bufsize]) ; }

//set stuff
    int setFilter(int i) { _keep[(_offset+i)%_bufsize]=false; }
    int setKeep(int i) { _keep[(_offset+i)%_bufsize]=true; }
    void setOutput(const string & fname, const string & mode);
//check stuff
    bool isFull() { return( _nsnp==_bufsize); }
    bool isFiltered(int i) { return(!_keep[(_offset+i)%_bufsize]); }


private:
    int **_gt,_ngt;
    bcf_srs_t *_sr;
    bcf1_t **_line;
    bcf_hdr_t *_hdr_in,*_hdr_out; 
    filter_t *_filter;
    int _nsnp,_offset,_nsample,_overlap,_window,_bufsize;
    htsFile *_fout;
    vector <bool> _keep;
    bool _write;//write output?
    int _nread,_nwrote;
    vector<float> _maf;
};


int correlationMatrix(circularBuffer & cb,Eigen::MatrixXf & R2);
float correlation(int *x,int *y, float mx, float my, float sx,float sy,int N);
int meanvar(circularBuffer & cb, vector<float> & mean, vector<float> & sd);
