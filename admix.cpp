#include <Eigen/Dense>
#include <Eigen/QR>
#include "akt.hpp"
#include "reader.hpp"
#include "logs.hpp"

using namespace std;
using namespace Eigen;

/**
 * @name    usage
 * @brief   print out options
 *
 * List of input options
 *
 */
static void usage(){ 
    cerr << "Approximate admixture fractions" << endl;
    cerr << "Usage:   ./akt admix input.txt -C centres -c 2-3 " << endl;
    umessage('c');
    cerr << "\t -C --cfile:			initial guess for cluster centres in text file" << endl;
    exit(1);
}

int admix_main(int argc,char **argv) {

    int c;

    if(argc<3) usage();
    static struct option loptions[] =    {
        {"cfile",1,0,'C'},
        {"cols",1,0,'c'},
        {0,0,0,0}
    };
    
    string cfile = "";
    string dims="";
   
    
    while ((c = getopt_long(argc, argv, "c:C:",loptions,NULL)) >= 0) {  
        switch (c)
        {
        case 'C': cfile = (optarg); break;
        case 'c': dims = (optarg); break;
        case '?': usage();
        default: cerr << "Unknown argument:"+(string)optarg << endl; exit(1);
        }
    }

	optind++;
    string input = argv[optind];
    cerr <<"Input: " << input << endl; 
	
	vector< vector<float> > data;
	vector< vector<string> > labels;
	ifstream in_file(input.c_str());
	readMatrix(in_file, data, labels, dims);	//read PCA projections
	in_file.close();

	int N = data.size();
	size_t d = data[0].size();
	
	//vector to Eigen
	MatrixXf P(d,N);
	for (int i=0; i<N; i++){ P.col(i) = VectorXf::Map(&data[i][0],d); }
	
	MatrixXf mu;

	ifstream in_cfile(cfile.c_str());
	vector< vector<float> > mu_data;
	vector< vector<string> > mu_lab;
	readMatrix(in_cfile, mu_data, mu_lab, "1-" + to_string(d) ); //read transformation matrix
	in_cfile.close();

	ifstream in_cfile2(cfile.c_str());
	vector< vector<float> > prop_data;
	vector< vector<string> > prop_lab; 
	readMatrix(in_cfile2, prop_data, prop_lab, to_string(d+2) + "-" + to_string(2*d+1) ); //read transformation matrix
	in_cfile2.close();
	
	if(mu_data.size() < d || mu_data.size() > d+1){ 
		cerr << "For " << d << " dimensional data need " << d << " or " << d+1 << " input vectors" << endl; 
		exit(1); 
	}	

	int K = d;
	int ct;
	MatrixXf L(d,K);
	VectorXf sub = VectorXf::Zero(d);
	VectorXf perp = VectorXf::Random(d);
	int subi=-1;
	
	if(mu_data.size() == d){	//specified admixture fractions but not zero vector
		mu.resize(d,K);
		float av_len = 0;
		//find centre of data
		for (size_t i=0; i<d; i++){ 
			mu.col(i) = VectorXf::Map(&mu_data[i][0],d); 
			sub += ( mu.col(i) );
			av_len += mu.col(i).norm();
		}
		sub /= d; //centre
		av_len /= d;
		//sides of tetrahedron
		MatrixXf sides(d,d);
		size_t ns = 0;
		for(size_t i=0; i<d; ++i){
			for(size_t j=i+1; j<d; ++j){
				sides.col(ns++) = mu.col(i) - mu.col(j);
				if(ns == d){break;}
			}	
			if(ns == d){break;}
		}
		sides.col(d-1) = perp;
		
		//find vector orthogonal to all sides
		ColPivHouseholderQR<MatrixXf> sQR(sides);
		MatrixXf Qm = sQR.matrixQ();
		perp = Qm.col(d-1);

		//zero vector is in perp direction distance of av_len
		sub += ( perp * av_len );
		for (int i=0; i<K; i++){ L.col(i) = VectorXf::Map(&prop_data[i][0],d); }
		for (int i=0; i<K; i++){ mu.col(i) -= sub; }
		
	} else {	//know the subtraction vector (one with smallest norm)
		subi = -1;
		float min_norm = -logz; //inf
		for(size_t k=0; k<prop_data.size(); ++k){
			float norm = 0; 
			for(size_t i=0; i<prop_data[k].size(); ++i){ norm += prop_data[k][i]*prop_data[k][i]; }
			if( norm < min_norm ){ min_norm = norm; subi = k; }
		}
		cerr << "Subtract minimum norm vector: line " << subi << " norm = " << min_norm << endl;

		mu.resize(d,K);
		ct = 0; 
		for (int i=0; i<K+1; i++){
			if(i!=subi){ mu.col(ct++) = VectorXf::Map(&mu_data[i][0],d); } 
		}
		
		if(subi!=-1){
			sub = VectorXf::Map(&mu_data[subi][0],d);
			for (int i=0; i<K; i++){ mu.col(i) -= sub; }
		}
		ct = 0; for (int i=0; i<K+1; i++){ if(i!=subi){ L.col(ct++) = VectorXf::Map(&prop_data[i][0],d); } }
	}
	cerr << "Transforming " << N << " " << d << "-vectors into " << d << " population admixtures." << endl; 

	//exact SVD of data
	JacobiSVD<MatrixXf> svd(mu, ComputeFullU | ComputeFullV);
    float  tol=1.e-5; //lower than this, set to zero
    VectorXf sing_inv( svd.singularValues().rows() );

	//calculate inverse of pca to admix transformation
    for ( int i=0; i<sing_inv.rows(); ++i) {
       if ( svd.singularValues()(i) > tol ){
          sing_inv(i)=1.0/svd.singularValues()(i);
       } else { 
		   sing_inv(i)=0; 
	   }
    }
    MatrixXf Tinv = L*( svd.matrixV()*sing_inv.asDiagonal()*svd.matrixU().transpose());
	
	//output
	cout << "SampleID\t"; 
	if( mu_lab.size() == d ){
		for (int i=0; i<K; i++){cout << mu_lab[i][0] << "\t"; } cout << endl;
	} else {
		for (int i=0; i<K+1; i++){ if(i!=subi){ cout << mu_lab[i][0] << "\t"; } } cout << endl;
	}
	for(int i=0; i<N; i++){
		
		VectorXf ad = Tinv * (P.col(i) - sub);	//apply pca to admix transformation
		cout << labels[i][0] << "\t" << ad.transpose() << endl;
	
	}
	
    return(0);
}

