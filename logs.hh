#ifndef LOGS_H
#define LOGS_H

#include <stdint.h>
#include <stdlib.h> 
#include <ctype.h>
#include <limits> 

using namespace std;

/************************
*	Constants			*
************************/
const double logz = -numeric_limits<double>::infinity();
const double logten = log(10);

//log(a + b) = log(a ( 1 + b/a) ) = log(a) + log(1 + b/a)
template<typename T> void log_sum(T &in, T add){

	if(add == in && add == logz){ in = logz;}
	else{
		T a = max(in, add);
		T sum = exp(in-a) + exp(add-a);
		in = log(sum) + a;
	}
	
}

#endif
