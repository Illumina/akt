#include "akt.hpp"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

static void usage(){
	cerr << "Calculate correlation across a population" << endl;	
	cerr << "Usage:" << endl;
	cerr << "./akt ldplot input_filename.vcf" << endl;
	cerr << "Expects input_filename.vcf to contain hard genotypes" << endl;
	umessage('r');
	umessage('R');
	umessage('s');
	umessage('S');
	exit(1);
}

	

int ldplot_main(int argc, char* argv[])
{
	int c;
    
    if(argc<3) usage();
    static struct option loptions[] =    {
        {"regions",1,0,'r'},	
		{"regions-file",1,0,'R'},
		{"samples",1,0,'s'},
		{"samples-file",1,0,'S'},
        {0,0,0,0}
    };

  string regions = "";
  bool regions_is_file = false;
  bool used_r = false;
  bool used_R = false;

  sample_args sargs;

    while ((c = getopt_long(argc, argv, "r:R:s:S:",loptions,NULL)) >= 0) {  
        switch (c)
        {
		case 'r': regions = (optarg); used_r = true; break;
		case 'R': regions = (optarg); used_R = true; regions_is_file = true; break;
        case 's': sargs.sample_names = (optarg); sargs.subsample = true; break;
        case 'S': sargs.sample_names = (optarg); sargs.subsample = true; sargs.sample_is_file = 1; break;
        case '?': usage();
        default: cerr << "Unknown argument:"+(string)optarg+"\n" << endl; exit(1);
        }
    }
    if( used_r && used_R ){ cerr << "-r and -R cannot be used simultaneously" << endl; exit(1); }

	optind++;
    string filename = argv[optind];
	
	//Setup htslib reader
	bcf_srs_t *sr =  bcf_sr_init() ; ///htslib synced reader.
	if(regions != ""){
		sr->require_index = 1;
		if ( bcf_sr_set_regions(sr, regions.c_str(), regions_is_file)<0 ){
			cerr << "Failed to read the regions: " <<  regions << endl; exit(1);
		}
	}
	if(!(bcf_sr_add_reader (sr, filename.c_str() ))){ cerr << "Problem opening " << filename << endl; exit(1); }

    if(sargs.subsample){ bcf_hdr_set_samples(sr->readers[0].header, sargs.sample_names, sargs.sample_is_file); }
		
	int N = bcf_hdr_nsamples(sr->readers[0].header);	///number of samples	
	cerr << N << " samples in " << filename << endl;
	
	int ngt_arr=N*2;
	bcf1_t *line;///bcf/vcf line structure.
	vector< int* > gt;
	int M=0;

	//read the input
	while(bcf_sr_next_line (sr)) { 
		
		line =  bcf_sr_get_line(sr, 0);
		
		//exclude haploid sites & multi allelics
		if( line->n_allele == 2){	
			int *gt_arr=(int *)malloc(N*2*sizeof(int));
			int ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr);   
			if(ngt > 0){
				gt.push_back( gt_arr );
				++M;
			}
		}
			
	}
	bcf_sr_destroy(sr);	
	cerr << M << " variants to correlate\nComputing " << M*M << " elements" << endl;
	
	float nm = 1.0/N;
	MatrixXf D(M,N);
	VectorXf mu = VectorXf::Zero(M); //these are the allele freqs
	VectorXf sig = VectorXf::Zero(M); //these are the allele vars
	for(int i=0; i<M; ++i){
		for(int j=0; j<N; ++j){
			if(gt[i][2*j] == bcf_gt_missing || gt[i][2*j+1] == bcf_gt_missing ){
				//cerr << "Found a missing site." << endl; exit(1); 
				gt[i][2*j] = 0;
				gt[i][2*j+1] = 0;
			}
			
			D(i,j) = (float)(bcf_gt_allele(gt[i][2*j]) + bcf_gt_allele(gt[i][2*j+1]));
			mu(i) += D(i,j);
				
		} 
		mu(i) *= nm;
		for(int j=0; j<N; ++j){ sig(i) += ( D(i,j) - mu(i) )*( D(i,j) - mu(i) ); } 
		sig(i) *= nm; 
	}	
	cerr << "Calculating Sigma." << endl;
	MatrixXf Sigma = MatrixXf::Zero(M,M);
	for(int i=0; i<M; ++i){
		for(int j=0; j<M; ++j){
			if( sig(i) != 0 && sig(j) != 0 ){ 
				for(int n=0; n<N; ++n){
					Sigma(i,j) += ( D(i,n) - mu(i) )*( D(j,n) - mu(j) );
				}
				Sigma(i,j) *= nm;
				Sigma(i,j) /= sqrt( sig(i) * sig(j) );
			}
		}
	}
	while(!gt.empty()) free(gt.back()), gt.pop_back(); 

	cout << Sigma << endl;
	return (0);
}

