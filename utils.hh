//
// Created by O'Connell, Jared on 8/3/17.
//

#ifndef AKT_UTILS_H
#define AKT_UTILS_H

#include <iostream>
#include <string>
#include <vector>
#include <sstream>


extern "C" {
#include "htslib/vcf.h"
}

bool is_genotyped(int *gt,int idx);

//splits s by any whitespace
int stringSplit(std::string & s,std::vector<std::string> & split);

//splits input by char split
int stringSplit(const std::string &input, const char split, std::vector<std::string> &out);

void die(const std::string &s);
void umessage(const char type);
#endif //AKT_UTILS_H
