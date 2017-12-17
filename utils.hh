//
// Created by O'Connell, Jared on 8/3/17.
//

#ifndef AKT_UTILS_H
#define AKT_UTILS_H

#include "akt.hh"
#include <string>
#include <vector>
#include <sstream>

bool is_genotyped(int *gt,int idx);

//splits s by any whitespace
int stringSplit(string & s,vector<string> & split);

//splits input by char split
int stringSplit(const string &input, const char split, vector<string> &out);

void die(const string &s);
void umessage(const char type);
#endif //AKT_UTILS_H
