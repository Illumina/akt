#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <iterator>
#include "reader.hh"

using namespace std;

/**
 * @name    readMatrix
 * @brief   read text file into std vector
 *
 * @param [in] in  		file containing tab separated data
 * @param [in] data   		holds stuff to cluster
 * @param [in] labels   	holds other stuff
 * @param [in] dims			which columns should it read
 */
void readMatrix(ifstream &in, vector<vector <float> > &data, vector< vector<string> > &labels, string dims){
	
	if(!in.is_open()){
		cout << "Failed to open file." << endl;
		exit(1);
	}
	
	stringstream iss(dims);
	vector<int> dim(2,1);
	int d = -1;

	//split string on "-"
	if(dims != ""){
		string item;
		int ct = 0;
		while( getline(iss, item, '-') && ct < 2){ dim[ct++] = atoi(item.c_str()); }
		d = (dim[1] - dim[0]) + 1;
	}		

	string line = ""; 	
	while(getline(in,line)) // loop through the file
	{
		stringstream is(line);
		istream_iterator<string> begin(is);
		istream_iterator<string> end;
        vector<string> tokens(begin, end);

		//require all rows to be same width
		if( d<0 ){ 
			d = tokens.size(); 
		} else {
			if(dims == "" && tokens.size() != (size_t)d){
				cerr << "ragged array elements at " << data.size() << endl; exit(1);
			}
		}
		if( (int)tokens.size()-dim[0]+1 < d ){ cerr << "too few array elements at " << data.size() << endl; exit(1); }
		
		//stuff before the data cols
		vector<string> tmps;
		for(int i=0; i<dim[0]-1; ++i){ 
				tmps.push_back( tokens[i] );
		} 
		//data cols
		vector<float> tmp;
		for(int i=dim[0]-1; i<dim[0]-1+d; ++i){ tmp.push_back( atof(tokens[i].c_str()) ); }
		//stuff after data cols
		for(int i=dim[0]-1+d; i<(int)tokens.size(); ++i){ 
				tmps.push_back( tokens[i] );
		} 
		
		data.push_back(tmp);
		labels.push_back(tmps);
	}
	
}



