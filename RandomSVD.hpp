// Randomised SVD implementation
// Copyright (c) 2016 Illumina, Inc.
// See LICENSE
// Auther: Jared O'Connell <joconnell@illumina.com>
//
//This is an implementation of the algorithm described in Halko 2011:
//http://arxiv.org/pdf/0909.4061.pdf

#pragma once

#include "math.h"
#include <stdlib.h>    
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

class RandomSVD {

public:
 
    RandomSVD(Eigen::MatrixXf & mat,int e,int q=3)	{
	int r = e;
	if(r>mat.rows())
	{
	    r=mat.rows();
	}
	if(q<1) 
	{
	    q=1;
	}

	Eigen::MatrixXf R;
	rnorm(R,mat.cols(),r);
	Eigen::MatrixXf Y  = mat * R;
	orthonormalize(Y);
	Eigen::MatrixXf Ystar;
	for(int i=0;i<q;i++)
	{
	    Ystar=mat.transpose() * Y;
	    orthonormalize(Ystar);
	    Y=mat * Ystar;
	    orthonormalize(Y);
	}	    
	Eigen::MatrixXf B = Y.transpose() * mat;
	Eigen::JacobiSVD<Eigen::MatrixXf> svd(B, Eigen::ComputeThinU | Eigen::ComputeThinV);
	_U = Y * svd.matrixU();
	_S = svd.singularValues();
	_V = svd.matrixV();
    }
    
    Eigen::MatrixXf matrixU() const{
	return _U;
    }
    
    Eigen::VectorXf singularValues() const {
	return _S;
    }
    
    Eigen::MatrixXf matrixV() const {
	return _V;
    }
    
private:
    Eigen::MatrixXf _U;
    Eigen::VectorXf _S;
    Eigen::MatrixXf _V;
    inline void rnorm(Eigen::MatrixXf & X,int nrow, int ncol) {
	X.resize(nrow,ncol);
	float pi = 3.141592653589793238462643383279502884;
	float rmax = (float)RAND_MAX;
	for(int i=0;i<nrow;i++)
	{
	    for(int j=0;j<ncol;j+=2)
	    {
		float v = ((float)std::rand() + 1.) / (rmax+2.);
		float u = ((float)std::rand() + 2.) / (rmax+2.);
		float c = sqrt(-2. * log(v));
		X(i,j) = c * cos(2. * pi * u);
		if(j<ncol-1)
		{
		    X(i,j+1) = c * sin(2. * pi * u);
		}
	    }	
	}
    }
    
    inline void orthonormalize(Eigen::MatrixXf & mat) {
	Eigen::MatrixXf  thinQ;
	thinQ.setIdentity(mat.rows(), mat.cols());
	mat = mat.householderQr().householderQ()*thinQ;
    }
};


    
