/**
 * @file   vcfpca.cpp
 * @Author Rudy Arthur (rudy.d.arthur@gmail.com) and Jared O'Connell (joconnell@illumina.com)
 * @date   May, 2015
 * @brief  PCA tool.
 *
 * A simple tool to read vcf/bcf files and perform principle component analysis.
 * Calculate principle components from scratch or use predetermined weights given in
 * a site only vcf.
 */
 
//because the one big multiplication it isn't worth the set up
#define EIGEN_DONT_PARALLELIZE
 
#include "akt.hpp"
#include "Eigen/Dense"
#include "RedSVD.hpp"
#include "reader.hpp"

using namespace Eigen;

/**
 * @name    usage
 * @brief   print out options
 *
 * List of input options
 *
 */
static void usage(){ 

    cerr << "Performs principal component analysis on a vcf/bcf" << endl;
    cerr << "Usage:" << endl;
    cerr << "./akt pca input.bcf" << endl;
	umessage('T');
	umessage('o');
	umessage('O');
	umessage('r');
	umessage('R');
	umessage('s');
	umessage('S');
	umessage('h');
	umessage('m');
    cerr << "\t -w --weight:			VCF with weights for PCA" << endl;
    cerr << "\t -N --npca:			first N principle components" << endl;
    cerr << "\t -a --alg:			exact SVD (slow)" << endl;
    cerr << "\t -C --covdef:			definition of SVD matrix: 0=(G-mu) 1=(G-mu)/sqrt(p(1-p)) 2=diag-G(2-G) default(1)" << endl;
    cerr << "\t -e --extra:			extra vectors for Red SVD" << endl;
    cerr << "\t -F --svfile:			File containing singular values" << endl;
    exit(1);
}
        
        
/**
 * @name    pca
 * @brief   Use weights in vcf2 to project vcf1
 *
 * @param [in] vcf1  vcf file to project
 * @param [in] vcf2  site only vcf containing PCA weights
 *
 */
void pca(string vcf1,string vcf2, bool don, int maxn, sample_args sargs) {
	
	int Nsamples;
	int Npca=0;  

	vector<string> names;
	
	int n0=0,n1=0;
    
	bcf_srs_t *sr =  bcf_sr_init() ; ///htslib synced reader.
	sr->require_index = 1;
	sr->collapse=COLLAPSE_NONE;

	///set the regions
	if ( bcf_sr_set_regions(sr, vcf2.c_str(), 1)<0 ){
		cerr << "Failed to read the regions: " <<  vcf2 << endl; exit(1);
	}
	///input file
	if(!(bcf_sr_add_reader (sr, vcf1.c_str() ))){
		cerr << "Problem opening " + vcf1 << endl; exit(1);
	}
	///sites file
	if(!(bcf_sr_add_reader (sr, vcf2.c_str() ))){
		cerr << "Problem opening " + vcf2 << endl; exit(1);
	}
	///subsample
    if(sargs.subsample){ bcf_hdr_set_samples(sr->readers[0].header, sargs.sample_names, sargs.sample_is_file); }
		
	int ret;
	bcf1_t *line0, *line1;
	int *gt_arr=NULL;int ngt=0,ngt_arr=0;
	float *wts=NULL;int nwts=0;

	int N = bcf_hdr_nsamples(sr->readers[0].header); ///number of samples;
	if(N<=0) {
	  cerr<<"ERROR: no samples found in "+vcf1<<endl;
	  exit(1);
	}
	cerr << N << " samples" << endl;
	Nsamples = N;
	for(int i=0; i<N; ++i){ 
		string tmp = sr->readers[0].header->samples[i]; 
		names.push_back(tmp);
	}	
    
	int nPC=0; ///number of pcs;
	vector< vector<float> > PC(N);///principal components.

	vector<float> gs(N);///genotypes, defaults to 0/0 - these should be easily called sites.	
	float *af_ptr=NULL,af;///read AF
	int nval = 0;
		
	while(bcf_sr_next_line (sr)) { //read

	if(bcf_sr_has_line(sr,1)){ ++n1; }	//in sites file	
	
	if(bcf_sr_has_line(sr,0) &&  (bcf_sr_has_line(sr,1))){	//in both files
		++n0;
		line1 =  bcf_sr_get_line(sr, 1);
		ret =  bcf_get_info_float(sr->readers[1].header, line1, "AF", &af_ptr, &nval); ///get allele frequency
		af = af_ptr[0];

		line0 =  bcf_sr_get_line(sr, 0);
		ngt = bcf_get_genotypes(sr->readers[0].header, line0, &gt_arr, &ngt_arr);	//get genotypes

		if(ngt==2*N){	///only diploid sites
			for(int i=0; i<2*N; i+=2){	///all samples in vcf
				if(gt_arr[i]!=bcf_gt_missing && gt_arr[i+1]!=bcf_gt_missing){
				  gs[i/2] = (float)(bcf_gt_allele(gt_arr[i]) + bcf_gt_allele(gt_arr[i+1]));
				} else { 
				  gs[i/2] = 2*af;	///assign missing sites to expected af
				}
			}
		}
		
		if(ret!=1) {
			cerr << "WARNING: no AF field at "<<line1->pos+1<<endl;
		}  else {	
			for(int n=0; n<N; ++n){
				
				if(!(gs[n]>=0 && gs[n]<=2)) {
				  cerr << "ERROR at " << bcf_hdr_id2name(sr->readers[0].header,line0->rid) <<":"<<line0->pos+1 
				  <<" g = "<<gs[n]<<" af="<<af<<endl;
				  exit(1);
				}
				//read PCA loadings
				ret =  bcf_get_info_float(sr->readers[1].header , line1, "WEIGHT", &wts, &nwts);
				if(ret<=0) {
				  cerr << bcf_hdr_id2name(sr->readers[1].header,line1->rid)<<":"<<line1->pos+1 << endl;
				  cerr << "no weights" << endl; 
				  exit(1);
				}
			
				  
				if(!isnan(wts[0])) {///sometimes you get nan weights due to monormoprhic sites in 1000g
				  gs[n] -= 2*af; //zero mean
				  if(nPC==0) {	//set correct number of Principle components
					if( don ){
						if(maxn > 0){ nPC = min( maxn, nwts ); }
						else { nPC = nwts;  }
					} else {
						nPC = nwts;
					}
					if(nPC==0) { cerr << "No principle components found in file" << endl; exit(1); }
					cerr << "Using " << nPC << " PCs from input file." << endl;
					Npca = nPC;
				  }
				  if(PC[n].size() != (unsigned)nPC){	//initialise vector
					PC[n].assign(nPC,0.);
				  }
				  for(int i=0;i<nPC;i++){	//projection
					PC[n][i] += wts[i]*gs[n];
				  }
				  if( isnan( PC[n][0] ) ){
					cerr << "nan value found. something went wrong." << endl; 
					exit(1);
				  }
				}
			}
		}
	}
		
	}
	bcf_sr_destroy(sr);
	free(gt_arr);
	free(wts);

	cerr << n0 << "/" << n1 << " of sites were in "<< vcf1 << endl;

	if(n0==0){
		cerr << "No intersecting SNPs found.  Check chromosome prefix matches on sites and input file." << endl; exit(1);
	}
	///print projections to stdout
	for(int n=0; n<Nsamples; ++n){
		cout << names[n] << "\t";
		for(int i=0;i<Npca;i++){ 
			cout << PC[n][i] << "\t";
		}
		cout << "\n";
	}

}




/**
 * @name    DatatoMatrix
 * @brief   Turn a std vector into an Eigen Matrix
 *
 * Eigen matrices have no push operator, since number of markers is variable
 * we have to fill the matrix ourselves.
 * @param [in] A   Matrix to get written
 * @param [in] G   vector containing vcf info
 * @param [in] AF  vector of allele frequencies * 2
 * @param [in] N   number of samples
 * @param [in] M   number of markers
 *
 */
void DatatoMatrix(Ref<MatrixXf> A, vector<float> &G, vector<float> AF, int N, int M, int covn){
	int ct = 0;
	for(int i=0; i<M; ++i){

		for(int j=0; j<N; ++j){
				A(j,i) = G[ct++];
		} 
		for(int j=0; j<N; ++j){
				A(j,i) -= AF[i];		///mean subtracted, required for PCA
				if(covn ==1){
					A(j,i) /= sqrt(0.5 * AF[i] * ( 1.0 - AF[i] * 0.5 ));	//divide by std dev
				}
		}
	} 
}
/**
 * @name    DatatoSymmMatrix
 * @brief   Turn a std vector into an Eigen Matrix subtract bias correction
 * http://www.ncbi.nlm.nih.gov/pubmed/26482676
 *
 * Eigen matrices have no push operator, since number of markers is variable
 * we have to fill the matrix ourselves.
 * @param [in] A   Matrix to get written
 * @param [in] G   vector containing vcf info
 * @param [in] AF  vector of allele frequencies * 2
 * @param [in] N   number of samples
 * @param [in] M   number of markers
 *
 */
void DatatoSymmMatrix(Ref<MatrixXf> A, vector<float> &G, vector<float> AF, int N, int M){

float norm = 0;
for(int i=0; i<M; ++i){ norm += 0.5 * AF[i] * ( 1.0 - AF[i] * 0.5 ); } norm /= 4;

for(int j1=0; j1<N; ++j1){
for(int j2=j1; j2<N; ++j2){
	A(j1,j2) = 0;
	if( j1 == j2){
		for(int i=0; i<M; ++i){
			A(j1,j2) += (G[i*N + j1] - AF[i]) * (G[i*N + j1] - AF[i]) - G[i*N + j1]*(2 - G[i*N + j1]);
		}
	} else {
		for(int i=0; i<M; ++i){
			A(j1,j2) += (G[i*N + j1] - AF[i]) * (G[i*N + j2] - AF[i]);
		}
	}

	A(j1,j2) /= norm;
	A(j2,j1) = A(j1,j2);
	
}}

}

/**
 * @name    calcpca
 * @brief   Calculate principle components and project onto leading ones
 *
 * @param [in] input_name   vcf containing data
 * @param [in] o   			switch to output weights
 * @param [in] output_name  if(o) file to output weights
 * @param [in] m   			minimum allele frequency required
 * @param [in] k   			thinning factor (skip k markers)
 * @param [in] a   			use Jacobi SVD
 * @param [in] npca   		number of principle components required
 * @param [in] extra   		number of extra vectors for RedSVD
 * @param [in] regions   	which variants to use
 * @param [in] pfile   		intersecting variant list
 *
 */
void calcpca(string input_name, bool o, string outf, string output_name, float m, int k, bool a, int npca, int extra,
string regions, bool regions_is_file, string pfilename, sample_args sargs, int covn, string svfilename) {
	
	cerr << "Reading data..." << endl;
	  
  	int nkept=0,nline=0,npanel=0;

	vector<string> names;
  
	bcf_srs_t *sr =  bcf_sr_init() ; ///htslib synced reader.
	sr->require_index = 1;

	//set the regions
	if(regions != ""){
		if ( bcf_sr_set_regions(sr, regions.c_str(), regions_is_file)<0 ){
			cerr << "Failed to read the regions: " <<  regions << endl; exit(1);
		}
		if(regions_is_file){ pfilename = regions; }
	}
	//input file
	if(!(bcf_sr_add_reader (sr, input_name.c_str() ))){
		string tmp = input_name;
		cerr << "Problem opening " + tmp << endl; exit(1);
	}
	//sites file
	if(pfilename != ""){
		if(!(bcf_sr_add_reader (sr, pfilename.c_str() ))){
			string tmp = pfilename;
			cerr << "Problem opening " + pfilename << endl; exit(1);
		}
	}
	//subset samples
    if(sargs.subsample){ bcf_hdr_set_samples(sr->readers[0].header, sargs.sample_names, sargs.sample_is_file); }
		
	int N = bcf_hdr_nsamples(sr->readers[0].header);	///number of samples
	if(N<=0) {
	  cerr<<"ERROR: no samples found in "+input_name<<endl;
	  exit(1);
	}
	int *gt_arr=(int *)malloc(N*2*sizeof(int)),ngt=N*2,ngt_arr=N*2;
	
	  cerr << N << " samples" << endl;
	  for(int i=0; i<N; ++i){ 
		string tmp = sr->readers[0].header->samples[i]; 
		names.push_back(tmp);
	  }
    

	vector<float> G; G.reserve(50000*N); ///genotypes stored here temporarily. 
	vector<float> AF;
	vector<int> sites;
	
	int count=0;
	bcf1_t *line;///bcf/vcf line structure.
	
	while(bcf_sr_next_line (sr)) { //read
		
		bool read = ( pfilename == "" ) ? true : (bcf_sr_has_line(sr,0) && bcf_sr_has_line(sr,1));
		if( read ){	//present in sites file and sample file.
			
			line =  bcf_sr_get_line(sr, 0);
			ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr);    
			if(ngt < 0){ cerr << "Bad genotypes at " << 
				 bcf_hdr_id2name(sr->readers[0].header,line->rid) << ":" << line->pos+1 << endl; exit(1); }
			
			int mac = 0,nmiss=0;
			for(int i=0;i<2*N;i++){	//calc allele count
				if(gt_arr[i]!=bcf_gt_missing){
					mac += bcf_gt_allele(gt_arr[i]);
				} else {
					++nmiss;
				}
			}
			float frq = (float)mac / (float)(2*N-nmiss);	///allele frequency

			//minor allele freq
			if(mac > (2*N-nmiss)/2) mac = (2*N-nmiss)-mac;
			if(mac > (2*N-nmiss)*m) ++count;

			//keep every k of these sites
			if(count%k==0 && mac > (2*N-nmiss)*m ){ //remember, 0%k == 0
				
				float mu = 0;	///actual mean =/= frq because default = 2*frq
				sites.push_back(nline);

				for(int i=0;i<N;i++){
					if( gt_arr[2*i] < 0 || gt_arr[2*i+1] < 0 ){
						cerr << "Fix Ploidy on " << line->rid << ":" << line->pos+1 << " sample " 
						<< sr->readers[0].header->samples[i] << endl; exit(1);
					}
					if(gt_arr[2*i]!=bcf_gt_missing && gt_arr[2*i+1]!=bcf_gt_missing){
						G.push_back((float)(bcf_gt_allele(gt_arr[2*i])+bcf_gt_allele(gt_arr[2*i+1])));
					} else {
						G.push_back(2*frq);      ///if missing push the "expected genotype" based on allele frequency.
					}
					mu += G.back();
				}
				AF.push_back(mu/(float)N );
				++nkept;
			}
		} //end sites file check
		if( bcf_sr_has_line(sr,0) ){ ++nline; } //lines in sample file
		if( pfilename != "" && bcf_sr_has_line(sr,1) ) { ++npanel; };
	}  
	
	bcf_sr_destroy(sr);	
	free(gt_arr);	
	
	if( pfilename != "" ){ cerr << nkept << "/"<<npanel<<" of study markers were in the sites file"<<endl; }
	else{
			cerr << "Kept " << nkept << " markers out of " << nline << endl;
	}

	if( nkept == 0 ){ 
		cerr << "ERROR: no intersecting SNPs found.  Check chromosome prefix matches on sites and input file." << endl; exit(1);
	}
	int M = nkept;
	
	int vsize = M;
	//vector to Eigen
	MatrixXf A(N, vsize); ///rows = samples, cols = markers
	if(covn >= 2){
		vsize = N;
		DatatoSymmMatrix(A, G, AF, N, M);
	} else {
		DatatoMatrix(A, G, AF, N, M, covn);
	}
	
    MatrixXf P = MatrixXf::Zero(N,npca);
    MatrixXf V(vsize, npca);
    
    ///Do the SVD of A
    bool out_sv = false; 
    ofstream out_file;
    if(svfilename != ""){ out_sv = true; out_file.open ( svfilename.c_str() ); }
		
    if(a){ //Jacobi algorithm
		JacobiSVD<MatrixXf> svd(A, ComputeThinU | ComputeThinV);
		for(int j=0; j<npca; ++j){ 
			P.col(j).noalias() = svd.matrixU().col(j) * svd.singularValues()(j) ;
			if(out_sv){ out_file << svd.singularValues()(j) << "\n";  }
		}
		V.noalias() = svd.matrixV().block(0,0,vsize,npca);
    } else {	//RedSVD algorithm
		int e = min(  min(N,vsize)-npca  , extra);
		RedSVD::RedSVD<MatrixXf> svd(A, npca + e);
		for(int j=0; j<npca; ++j){ 
			P.col(j).noalias() = svd.matrixU().col(j) * svd.singularValues()(j) ;
			if(out_sv){ out_file << svd.singularValues()(j) << "\n"; }
		}
		V.noalias() = svd.matrixV().block(0,0,vsize,npca);
	} 

	if(out_sv){ out_file.close(); }
	
	if(o){ 	//output sites file
		
		cerr <<"Printing coefficients to " << output_name << endl; 
		  
		bcf_srs_t *reader =  bcf_sr_init() ; ///htslib synced reader.
		reader->require_index = 1;

		if(regions != ""){
			if ( bcf_sr_set_regions(reader, regions.c_str(), regions_is_file)<0 ){
				cerr << "Failed to read the regions: " <<  regions << endl; exit(1);
			}
		}
		if(!(bcf_sr_add_reader (reader, input_name.c_str() ))){
			string tmp = input_name;
			cerr << "Problem opening " + tmp << endl; exit(1);
		}
		if(pfilename != ""){
			if(!(bcf_sr_add_reader (reader, pfilename.c_str() ))){
				string tmp = pfilename;
				cerr << "Problem opening " + pfilename << endl; exit(1);
			}
		}
		
		bcf_hdr_t *hdr = bcf_hdr_dup(reader->readers[0].header);
		//new header
		bcf_hdr_append(hdr, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate allele frequency\">");
		bcf_hdr_append(hdr, "##INFO=<ID=WEIGHT,Number=20,Type=Float,Description=\"PCA loading\">");

		bcf_hdr_t *new_hdr = bcf_hdr_subset(hdr,0,NULL,NULL); ///creates a new subsetted header (with 0 samples) from src_header
		bcf_hdr_add_sample(new_hdr, NULL);      /// update internal structures
		bcf1_t *rec = bcf_init1() ;

		htsFile *out_fh  = hts_open(output_name.c_str(), outf.c_str());
		bcf_hdr_write(out_fh, new_hdr);

		bcf1_t *line;///bcf/vcf line structure.
	
		int idx = 0;
		nline = 0;

		float *weights = (float *)malloc(npca*sizeof(float));
		while(bcf_sr_next_line (reader)) { 
			bool read = ( pfilename == "" ) ? true : (bcf_sr_has_line(reader,0) && bcf_sr_has_line(reader,1));
			if( read ){	//present in sites file and sample file.
			
				line =  bcf_sr_get_line(reader, 0); 
				if( (unsigned)idx < sites.size() && sites[idx] == nline ){	//tracks sites which were included in PCA
					
					//copy line essentials
					rec->rid = line->rid;
					rec->pos = line->pos;
					rec->qual = line->qual;

					bcf_update_id(new_hdr, rec, line->d.id);
					string ref = line->d.allele[0];
					string alt = line->d.allele[1];
					string alleles = ref + "," + alt;
					bcf_update_alleles_str(new_hdr, rec, alleles.c_str());
				
					float taf = AF[idx]*0.5;
					bcf_update_info_float(new_hdr, rec, "AF", &taf, 1);	
					
					for(int i=0; i<npca; ++i){ weights[i] = V(idx, i); }
					bcf_update_info_float(new_hdr,rec,"WEIGHT",weights,npca);
					
					bcf_unpack(rec, BCF_UN_ALL);								
					bcf_write1(out_fh, new_hdr, rec) ;
					bcf_clear1(line) ;
					++idx;
				}
			}
			if( bcf_sr_has_line(reader,0) ){ ++nline; } //lines in sample file
		}
		free(weights);

		bcf_destroy(rec);
		bcf_sr_destroy(reader);
		hts_close(out_fh);
		bcf_hdr_destroy(hdr);
		bcf_hdr_destroy(new_hdr);
		
	} 
	///print projections to stdout
	for(int j=0; j<N; ++j){ cout << names[j] << "\t" << P.row(j) << endl; }
	

}



        
int pca_main(int argc,char **argv) {
    
    int c;
    if(argc<3) usage();
    static struct option loptions[] =    {
        {"out",1,0,'o'},	
        {"outf",1,0,'O'},	
        {"weight",1,0,'w'},
        {"region",1,0,'r'},
		{"regions-file",1,0,'R'},
		{"targets-file",1,0,'T'},
        {"maf",1,0,'m'},
        {"thin",1,0,'h'},
        {"npca",1,0,'N'},
        {"alg",0,0,'a'},
        {"covdef",0,0,'C'},
        {"extra",1,0,'e'},
        {"samples",1,0,'s'},
        {"samples-file",1,0,'S'},
        {0,0,0,0}
    };
    float  m=.05;
    int thin=400;
    int n=20; bool don = false;
    bool a=false;
    int e=200;
    bool o=false; string out_filename="";
    bool w=false; string weight_filename;
    string outf = "w";
    
    sample_args sargs;

    string pfilename = "";

    string regions = "";
    bool regions_is_file = false;
    bool used_r = false;
    bool used_R = false;

	int covn = 1;
	
	string svfilename = "";
	
    while ((c = getopt_long(argc, argv, "o:O:w:m:h:N:ae:r:R:s:S:T:C:F:",loptions,NULL)) >= 0) {  
        switch (c)
        {
		case 'o': o = true; out_filename = optarg; break;
        case 'O': outf += (string)(optarg); break;
        case 'w': w = true; weight_filename = (string)optarg; break;
        case 'm': m = atof(optarg); break;

        case 'h': thin = atoi(optarg); break;
        case 'N': don = true; n = atoi(optarg); break;
        case 'C': covn = atoi(optarg); break;
        case 'a': a = true; break;
        case 'e': e = atoi(optarg); break;

        case 'r': regions = (optarg); used_r = true; break;
		case 'R': regions = (optarg); used_R = true; regions_is_file = true; break;
		case 'T': pfilename = (optarg);  break;    
		case 'F': svfilename = (optarg);  break;    

        case 's': sargs.sample_names = (optarg); sargs.subsample = true; break;
        case 'S': sargs.sample_names = (optarg); sargs.subsample = true; sargs.sample_is_file = 1; break;
        case '?': usage();
        default: 
	  if(optarg!=NULL) {cerr << "Unknown argument:"+(string)optarg << endl; exit(1);}
	  else {cerr << "Unknown argument:"; exit(1);}
        }
    }

    if(optind>=argc-1) {cerr<<"No input .bcf/.vcf provided!"<<endl;exit(1);}
    if(!pfilename.empty()||used_R) thin=1;
    if( used_r && used_R ){ cerr << "-r and -R cannot be used simultaneously" << endl; exit(1); }
	

	optind++;
    string input = argv[optind];
    cerr <<"Input: " << input << endl; 
	if(w){ 
		cerr << "Using file " << weight_filename << " for PCA weights" << endl; 
		pca(input,weight_filename, don, n, sargs);
	}
	else{
		cerr << "MAF lower bound: " << m << "\nThin: "<< thin <<" \nNumber principle components: "<<n<<endl;
		calcpca(input,o,outf,out_filename,m,thin,a,n,e,regions,regions_is_file,pfilename,sargs, covn, svfilename);
	}


    return(0);
}

