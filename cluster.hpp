#ifndef CLUSTER_H
#define CLUSTER_H

#include <Eigen/Dense>
#include <limits> 
#include <stdlib.h> 
#include <vector>
#include <limits> 
#include <algorithm>

using namespace std;
using namespace Eigen;

//compare vector by first component
extern bool cf_idx(vector<int> i, vector<int> j);
//compare pair by first value
extern bool cf_pair(pair<int, float> i, pair<int, float> j);

//let's try to use Eigen here so I don't have to write many loops
class Cluster
{
	public:
	
		//input
		MatrixXf P; //the data  N x d (this way for easier interface with vector< vector >
		int d; //dimension of the data
		int N; //number of data points
		int K; //number of clusters
		float eps;	//target accuracy for EM
		float dc; //cutoff for density cluster
		
		//output
		MatrixXf centres;
		VectorXi assignment;
		VectorXi sizes;
		
		//em stuff
		VectorXf alpha; 
		vector< MatrixXf > covariance; 
		MatrixXf prob; 
		MatrixXf weight;

		//density stuff
		VectorXi rho;
		VectorXf delta;
		
		//silhouette
		VectorXf sil;
		bool silset;
		
		//empty constructor
		Cluster(){
			d = 0;
			N = 0;
			K = 0;
			eps = 1e-4;
		}
		//no data constructor
		Cluster(int K_, float eps_=1e-4){
			K = K_;
			eps = eps_;
		}
		//init
		void initP(const Ref<MatrixXf> P_){
			P = P_;
			d = P.cols();
			N = P.rows();
		}
		//init
		void initK(int K_){
			K = K_;

			centres.resize(K, d);
			assignment = VectorXi::Zero(N);
			sizes = VectorXi::Zero(K); 
			
			alpha.resize(K);
			for(int k=0; k<K; ++k){ covariance.push_back( MatrixXf::Zero(d,d) ); }
			weight = MatrixXf::Zero(K,N);
			prob.resize(K,N);
			
			
			sil = VectorXf::Zero(N); 
			silset = false;
		}
		//full constructor
		Cluster(const Ref<MatrixXf> P_, int K_, float eps_=1e-4){
			eps = eps_;
			initP(P_);
			initK(K_);
		}
		//copy-constructor
		Cluster(const Cluster& other) : d(other.d), N(other.N), K(other.K), eps(other.eps)
		{
			P = other.P;  
		}
		//copy-swap assignment
		friend void swap(Cluster& first, Cluster& second) 
		{
			swap(first.P, second.P);
			swap(first.d, second.d);
			swap(first.N, second.N);
			swap(first.K, second.K);
			swap(first.eps, second.eps);
		}
		Cluster& operator=(Cluster other) 
		{
			swap(*this, other); 
			return *this;
		}
		//move constructors are stupid
				
		//output functions
		void clustered_data_dump(vector<vector<string> > &labels);

		//basic init functions
		void initialiseCentres();
		void assignCentres(const Ref<MatrixXf> centres_){ centres = centres_; }
		void clusterAssign();
		
		//em init functions
		float EMweights();
		void EMassign();
	
		//density init functions
		void localDensity(float dc_);
		void minDistance();
		void densityPlot();
		void densityCluster(int min_rho, float min_delta);

		//clustering
		void kMeans(int max_it);
		void kppMeans(int max_it);
		void EMcluster(int max_its);

		//cluster QA
		void silhouette();

};

#endif
