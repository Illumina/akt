#ifndef READER_H
#define READER_H

#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include <iterator>
#include <stdexcept>     
#include <stdlib.h> 

using namespace std;

//utility for reading colum based text files with mixed string and float data
void readMatrix(ifstream &in, vector<vector <float> > &data, vector< vector<string> > &labels, string dims);

class range{
public:
	string chrom;
	int begin;
	int end;
	
	range() : begin(0), end(0), chrom(""){}
	range(string _chrom, int _begin, int _end) : chrom(_chrom), begin(_begin), end(_end) {}
	range(string regions){
		
		int minpos, maxpos;

		stringstream ss(regions);
		string item;
		getline(ss, chrom, ':'); getline(ss, item, ':');

		if(item!=""){
			stringstream ss2(item);
			string lt, rt;
			getline(ss2, lt, '-'); getline(ss2, rt, '-');	
			
			
			try {
				begin = atoi(lt.c_str()); //c++11 stoi !
			} catch (const invalid_argument& ia) {
				cerr << "Parallel input errorn: converting " << lt << " to int" << endl; exit(1);
			}		
			try {
				end = atoi(rt.c_str());
			} catch (const invalid_argument& ia) {
				cerr << "Parallel input errorn: converting " << rt << " to int" << endl; exit(1);
			}
		} else {
			cerr << "Parallel input error: no start-end found in " << regions << endl; exit(1);
		}

	}
	range(const range& other) : begin(other.begin), end(other.end), chrom(other.chrom){}
	
	vector<range> split(int n){
		
		vector<range> ranges(n);
		int rinc = (end - begin)/n;
		int st = begin;
		
		for(int i=0; i<n;++i){
			range r(chrom, st, st+rinc);
			st += rinc+1;
			ranges[i] = r;
		}
		
		return ranges;
	}	
		
	friend void swap(range& first, range& second)
    {
        swap(first.begin, second.begin);
        swap(first.end, second.end);
        swap(first.chrom, second.chrom);
    }	
    range& operator=(range other)
	{
		swap(*this, other); 
		return *this;
	}

};

#endif
