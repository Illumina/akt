#pragma once

#include "math.h"
#include <stdlib.h>    
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

using namespace std;

class RandomSVD {

public:
 
    RandomSVD(Eigen::MatrixXf & mat,int e,int q)
	{
	    int r = e;
	    if(r>mat.rows())
	    {
		r=mat.rows();
	    }

	    Eigen::MatrixXf R;
	    rnorm(R,mat.cols(),mat.rows());
	    Eigen::MatrixXf Y  = mat * R;
	    orthonormalize(Y);
	    for(int i=0;i<q;i++)
	    {
		Y=mat.transpose() * Y;                                                                                                                                                                    
		orthonormalize(Y);                                                                                                                                                                      
		Y=mat * Y;                                                                                                                                                                                
		orthonormalize(Y);                                                                                                                                                                      
	    }                                                                                                                                                                                           
	    
	    Eigen::MatrixXf B = Y.transpose() * mat;                                                                                                                                                          
	    Eigen::JacobiSVD<Eigen::MatrixXf> svdOfC(B, Eigen::ComputeThinU | Eigen::ComputeThinV);                                                                                                         
	    _U = Y * svdOfC.matrixU();                                                                                                                                                           
	    _S = svdOfC.singularValues();                                                                                                                                                        
	    _V =  svdOfC.matrixV();         
	}
    
    Eigen::MatrixXf matrixU() const
	{
	    return _U;
	}
    
    Eigen::VectorXf singularValues() const
	{
	    return _S;
	}
    
    Eigen::MatrixXf matrixV() const
	{
	    return _V;
	}
    
private:
    Eigen::MatrixXf _U;
    Eigen::VectorXf _S;
    Eigen::MatrixXf _V;
    inline void rnorm(Eigen::MatrixXf & X,int nrow, int ncol)
	{
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

    inline void orthonormalize(Eigen::MatrixXf & mat)
	{
	    Eigen::MatrixXf  thinQ;
	    thinQ.setIdentity(mat.rows(), mat.cols());
	    mat = mat.householderQr().householderQ()*thinQ;
	}
       
};

    
