#ifndef IBD_H
#define IBD_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <ctype.h>
#include <limits> 
#include <iterator>
#include <map>
#include <string>

using namespace std;

extern void read_pairs(ifstream &in, vector< pair<string, string> > &relpairs, map<string,int> &name_to_id);

extern void make_pair_list(vector< pair<string, string> > &relpairs, vector<string> names);

#endif
