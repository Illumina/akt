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
const double logtwo = log(2);
const double logh = log(0.5);
const double logq = log(0.25);
const double logten = log(10);


const double Trans[3][3][3] = 
{
	{
		{1,0,0},
		{0.5,0.5,0},
		{0,1,0}
	},
	{
		{0.5,0.5,0},
		{0.25,0.5,0.25},
		{0,0.5,0.5}
	},
	{
		{0,1,0},
		{0,0.5,0.5},
		{0,0,1}
	}
};

const double LogTrans[3][3][3] = 
{
	{
		{0,logz,logz},
		{logh,logh,logz},
		{logz,0,logz}
	},
	{
		{logh,logh,logz},
		{logq,logh,logq},
		{logz,logh,logh}
	},
	{
		{logz,0,logz},
		{logz,logh,logh},
		{logz,logz,0}
	}
};

template<typename T> void log_sum(T &in, T add){

	if(add == in && add == logz){ in = logz;}
	else{
		T a = max(in, add);
		T sum = exp(in-a) + exp(add-a);
		in = log(sum) + a;
	}
	
}

#endif
