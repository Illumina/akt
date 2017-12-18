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

//class for reading colum based text files with mixed string and float data
void readMatrix(ifstream &in, vector<vector <float> > &data, vector< vector<string> > &labels, string dims);

#endif
