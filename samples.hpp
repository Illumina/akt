#ifndef __SAMPLES_H__
#define __SAMPLES_H__

using namespace std;

class sample_args{
public:

char *sample_names;
int sample_is_file;
bool subsample;

	sample_args(){
		sample_names = NULL;
		sample_is_file = 0;
		subsample = false;
	}
	

};


#endif
