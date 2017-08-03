//
// Created by O'Connell, Jared on 8/3/17.
//

#ifndef AKT_UTILS_H
#define AKT_UTILS_H

#include "akt.hpp"
#include <string>
#include <vector>
#include <sstream>

bool is_genotyped(int *gt,int idx);

int stringSplit(string & s,vector<string> & split);

void getGenotype(int idx,int *gt_arr,int *ret);

#endif //AKT_UTILS_H
