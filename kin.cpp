/**
 * @file   kin.cpp
 * @Author Rudy Arthur (rudy.d.arthur@gmail.com)
 * @brief  Kinship calculator.
 *
 * A simple tool to read vcf/bcf files and calculate IBD and kinship.
 * Calculate allele frequencies from data or from input file.
 */
 
#include "akt.hpp"
#include "logs.hpp"
#include "kin.hpp"
#include "reader.hpp"
#include <bitset>
#include <iomanip>

///Size of bitset. Arbitrary but 128 works well in practice
#define L 128

using namespace std;

/**
 * @name    read_pair
 * @brief   read input file containing sample pair names
 *
 * @param [in] in file 		handle of input list
 * @param [in] relpairs  	vector of pairs found
 * @param [in] name_to_id  	samples found in VCF
 *
 */
void read_pairs(ifstream &in, vector< pair<string, string> > &relpairs, map<string,int> &name_to_id )
{
	
    if(!in.is_open())
    {
	cout << "Failed to open file." << endl;
	exit(1);
    }
    
    string line = ""; 	
    while(getline(in,line)) // loop through the file
    {
	stringstream is(line);
	istream_iterator<string> begin(is);
	istream_iterator<string> end;
	vector<string> tokens(begin, end);
      
	if ( name_to_id.find(tokens[0]) == name_to_id.end() ) 
	{
	    cerr << tokens[0] << " not found in input" << endl; exit(1);
	} 
	else if ( name_to_id.find(tokens[1]) == name_to_id.end() ) 
	{
	    cerr << tokens[1] << " not found in input" << endl; exit(1);
	} 
		
	pair<string, string> tmps = make_pair(tokens[0], tokens[1]);
	relpairs.push_back(tmps);
    }	
}

/**
 * @name    make_pair_list
 * @brief   make a list of all possible sample pairs
 *
 * @param [in] relpairs  	vector of pairs
 * @param [in] names  		samples found in VCF
 *
 */
void make_pair_list(vector< pair<string, string> > &relpairs, vector<string> names)
{
    for(size_t j1=0; j1<names.size(); ++j1)
    {
	for(size_t j2=j1+1; j2<names.size(); ++j2)
	{
	    relpairs.push_back( make_pair(names[j1], names[j2]) );
	}
    }
}

/**
 * @name    usage
 * @brief   print out options
 *
 * List of input options
 *
 */
static void usage()
{
    cerr << "\nAbout: Calculate IBD stats from a VCF" << endl;	
    cerr << "Usage: akt kin [options] <in.bcf>" << endl;
    cerr << "Expects input.bcf to contain genotypes." << endl;

    cerr << "\nKinship calculation options:"<<endl;
    cerr << "\t -k --minkin:			threshold for relatedness output (none)" << endl;
//    cerr << "\t -u --unnorm:			If present don't normalize" << endl;
//  cerr << "\t -c --calc:			calculate AF from data" << endl;
    cerr << "\t -F --freq-file:                a file containing population allele frequencies to use in kinship calculation"<<endl;
    cerr << "\t -M --method:			type of estimator. 0:plink (default) 1:king-robust 2:genetic relationship matrix" << endl;
    umessage('a');
    umessage('n');
    cerr << "\nSite filtering options:"<<endl;  
    umessage('R');
//  umessage('T');
    umessage('r');
    umessage('m');
    umessage('h');
    cerr << "\nSample filtering options:"<<endl;
    umessage('s');
    umessage('S');
//  umessage('f');

    exit(1);
}

Kinship::Kinship(int nsample)
{
    _nsample=nsample;
    _n00=0;
    _n10=0;
    _n11=0;
    _n20=0;
    _n21=0;
    _n22=0;
    _af.resize(50000);//shouldnt need more markers than this.
}

void Kinship::update_n(float p) 
{
    float q = 1-p;
    _n00 += 2*p*p*q*q;
						
    _n11 += 2*p*q; ///2ppq + 2qqp = 2pq(p+q) = 2pq == Hij
    _n10 += 4*p*q*(p*p + q*q); ///4pppq + 4qqqp = 4pq(pp + qq)
						
    _n20 += p*p*p*p + q*q*q*q + 4*q*q*p*p; ///pppp + qqqq + 4ppqq 
    _n21 += p*p + q*q; ///ppp + qqq + ppq + pqq = pp(p+q) + qq(q+p) = pp + qq
    _n22 += 1;
    _af.push_back(p);    
}

void Kinship::estimate_ibd(float & ibd0, float & ibd1, float & ibd2,float & ibd3,bool normalise) 
{
    ///method of moments
    ibd0 /= _n00;	
    ibd1 = (ibd1 - ibd0*_n10)/_n11;
    ibd2 = (ibd2 - ibd0*_n20 - ibd1*_n21)/_n22;
    ibd3 = _n22 - ibd3;

    ///_normalize i_n [0,1]
    if(normalise)
    {
	if(ibd0 > 1)//very unrelated, project to 100
	{
	    ibd0 = 1;  ibd1 = 0; ibd2 = 0; 
	}
	if(ibd1 < 0)
	{ 
	    ibd1 = 0; 
	}
	if(ibd2 < 0)
	{ 
	    ibd2 = 0; 
	}
		
	float sum = ibd0 + ibd1 + ibd2;
	ibd0/=sum;
	ibd1/=sum;
	ibd2/=sum;
    }
}

int kin_main(int argc, char* argv[])
{
	
    int c;
  
    if(argc<3) usage();
    static struct option loptions[] =    {
	{"regions",1,0,'r'},	
	{"nthreads",1,0,'n'},
	{"freq-file",1,0,'F'},	
	{"regions-file",1,0,'R'},
	//    {"targets-file",1,0,'T'},
	{"aftag",1,0,'a'},
	{"minkin",1,0,'k'},
	{"thin",1,0,'h'},
//	{"unnorm",1,0,'u'},
	{"maf",1,0,'m'},
	{"method",1,0,'M'},
//    {"pairfile",1,0,'f'},
	{"samples",1,0,'s'},
	{"samples-file",1,0,'S'},
	{0,0,0,0}
    };
    int method=0;
    string regions = "";
    bool norm = true;
    float min_kin = 0; bool tk = false;
    int thin = 1;
    int nthreads = -1;
    float min_freq = 0;
    string pairfile="";
    string af_tag = "AF";
    sample_args sargs;
    bool regions_is_file = false;
    bool used_r = false;
    bool used_R = false;
    bool target = false;
    string frq_file="";  
    while ((c = getopt_long(argc, argv, "r:R:M:F:k:h:n:m:f:a:s:S:c",loptions,NULL)) >= 0) 
    {  
	switch (c)
	{
	case 'r': regions = (optarg); used_r = true; break;
	case 'F': frq_file=optarg;break;
	case 'M': method=atoi(optarg);break;
	case 'R': regions = (optarg); used_R = true; regions_is_file = true; break;
	    //      case 'T': target = true; pfilename = (optarg); break;
	case 'k': tk = true; min_kin = atof(optarg); break;
	case 'h': thin = atoi(optarg); break;
//	case 'u': norm = false; break;
	case 'm': min_freq = atof(optarg); break;
	case 'n': nthreads = atoi(optarg); break;
//      case 'f': pairfile = (optarg); break;
	case 'a': af_tag = string(optarg); break;
	case 's': sargs.sample_names = (optarg); sargs.subsample = true; break;
	case 'S': sargs.sample_names = (optarg); sargs.subsample = true; sargs.sample_is_file = 1; break;
	case '?': usage();
	default: cerr << "Unknown argument:"+(string)optarg+"\n" << endl; exit(1);
	}
    }

    if(method<0 || method>1) {
	cerr << "ERROR: method must be one of 0/1"<<endl;
	exit(1);
    }
    if( used_r && used_R )
    { 
	cerr << "-r and -R cannot be used simultaneously" << endl; exit(1); 
    }
    //  if( !used_R && !target ){ cerr << "Must use one of -R or -T" << endl; exit(1); }

    assert(min_freq>=0 && min_freq<=1);
    if((used_R&&target) )  
    {
	cerr<<"ERROR: -h/-m and -R/-T are incompatible"<<endl;
	exit(1);
    }
    if(!frq_file.empty() && !regions.empty()) 
    {
	cerr<<"ERROR: -F and -R/-r are incompatible!"<<endl;
	exit(1);
    } 
    else if(regions.empty())
    {
	regions=frq_file;
	regions_is_file=true;
    }

    if(frq_file=="")  
    {
	cerr<<"No frequency VCF provided (-F). Allele frequencies will be estimated from the data."<<endl;
    }
    else 
    {
	cerr <<"Taking allele frequencies from "<<frq_file<<" using INFO/"<<af_tag<<endl;
    }
    if(nthreads < 1)
    { 
	nthreads = 1; 
    }
    omp_set_num_threads(nthreads);
#pragma omp parallel
    {
	if(omp_get_thread_num() == 0)
	{
	    if( omp_get_num_threads() != 1)
	    {
		cerr << "Using " << omp_get_num_threads() << " threads" << endl;
		nthreads = omp_get_num_threads();
	    }
	}
    }

    optind++;
    string filename = argv[optind];	///input VCF
			
    int Nsamples;
    vector<string> names;			///sample names
    map<string,int> name_to_id;   ///sample ids
	  
    int sites=0,markers=0,num_sites=0,num_study=0;
  
    bcf_srs_t *sr =  bcf_sr_init() ; ///htslib synced reader.
    sr->collapse = COLLAPSE_NONE;		///require matching ALTs
    sr->require_index = 1;			///require indexed VCF

    ///subset regions
    if(!regions.empty())
    {
	if ( bcf_sr_set_regions(sr, regions.c_str(), regions_is_file)<0 ){
	    cerr << "Failed to read the regions: " <<  regions << endl; return 0;
	}
    }

    ///open input VCF
    if(!(bcf_sr_add_reader (sr, filename.c_str() )))
    { 
	cerr << "Problem opening " << filename << endl; 
	cerr << "Input file not found." << endl;
	bcf_sr_destroy(sr);	
	return 0;
    }
    if(bcf_hdr_nsamples(sr->readers[0].header)==0) 
    {
	cerr <<"ERROR: no samples in "<<filename<<endl;
	exit(1);
    }
    ///Open file of allele freqs
    if(frq_file!="" && !(bcf_sr_add_reader (sr, frq_file.c_str() )))
    { 
	cerr << "Problem opening " << frq_file << endl; 
	cerr << "Sites file not found." << endl;
	bcf_sr_destroy(sr);	
	return 0;
    }
    ///subsample input vcf
    if(sargs.subsample)
    { 
	bcf_hdr_set_samples(sr->readers[0].header, sargs.sample_names, sargs.sample_is_file); 
    }
  
    int N = bcf_hdr_nsamples(sr->readers[0].header);	///number of samples
  
    cerr << N << " samples" << endl;
    Nsamples = N;
    for(int i=0; i<N; ++i)
    { 
	string tmp = sr->readers[0].header->samples[i]; 
	names.push_back(tmp);
	name_to_id[tmp] = i;
    }
    Kinship K(Nsamples);
    vector< vector< vector< bitset<L> > > > bits(N); ///[sample][site][type][val]

    int count=0;

    bcf1_t *line, *line2;///bcf/vcf line structure.
	
    int *gt_arr=(int *)malloc(N*2*sizeof(int)),ngt=N*2,ngt_arr=N*2;
    float *af_ptr=(float *)malloc(1*sizeof(float)); int nval = 1;

    int bc = 0;
    bool use_frq = !frq_file.empty();
    cerr << "Reading genotypes...";
    while(bcf_sr_next_line (sr))  ///read file
    {
	if(bcf_sr_has_line(sr,0) && (!use_frq||bcf_sr_has_line(sr,1)) )  ///present in the study file (and frequency file)
	{
	    int nmiss=0;
	    int npres=0;
	    int sum = 0;	///AC

	    line =  bcf_sr_get_line(sr, 0);
//	    cerr<<line->pos+1<<endl;
	    if(line->n_allele == 2 && (count++)%thin==0 )		///bi-allelic
	    {
		ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr);  
		if(ngt < 0)
		{ 
		    cerr << "Bad genotypes at " << line->pos+1 << endl; exit(1); 
		}
			
		///htslib -> int
		for(int i=0;i<2*N;i++)
		{ 
		    if(gt_arr[i]==bcf_gt_missing || (bcf_gt_allele(gt_arr[i])<0) || (bcf_gt_allele(gt_arr[i])>2)  )
		    {
			gt_arr[i] = -1;
			++nmiss;
		    } 
		    else 
		    {
			sum += bcf_gt_allele(gt_arr[i]);
			++npres;
		    }
		}
		float p;
		if(frq_file=="")	///calculate AF from data
		{
		    p = (float)sum / (float)(npres);	///allele frequency
		} 
		else 
		{
		    assert(bcf_sr_has_line(sr,1));	///present in sites file.      
		    line2 =  bcf_sr_get_line(sr, 1);
		    num_sites++;
		    ++sites;
		    int ret = bcf_get_info_float(sr->readers[1].header, line2, af_tag.c_str(), &af_ptr, &nval);
		    if( ret<0 || nval != 1 )
		    { 
			cerr << af_tag << " read error at " << line2->rid << ":" << line->pos+1 << endl; exit(1); 
		    }
		    p = af_ptr[0];
		}
		if( (p < 0.5) ? ( p > min_freq ) : (1-p > min_freq) )///min af	  
		{
		    K.update_n(p);
		    ///for each site record truth table
		    ///g=0  g=1  g=2  g=missing
		    for(int i=0; i<N; ++i)
		    { 
			if(bc == 0)
			{ 
			    bits[i].push_back( vector< bitset<L> >(4, bitset<L>() ) ); 
			}
			if(gt_arr[2*i] != -1 && gt_arr[2*i+1] != -1)
			{
			    bits[i].back()[ bcf_gt_allele(gt_arr[2*i]) + bcf_gt_allele(gt_arr[2*i+1]) ][bc] = 1;
			} 
			else 
			{
			    bits[i].back()[ 3 ][bc] = 1;
			}
		    }
		    ///chunks of size L
		    bc = (bc+1)%(L);
			  
		    ++markers;
		}	
	    } //thin			
	    ++num_study;
	} //in study
    }//reader
    cerr << "done."<<endl;
    free(gt_arr);
    free(af_ptr);
    bcf_sr_destroy(sr);	

    if(frq_file.empty()) 
    {
	cerr << "Using "<<markers<<" markers for calculations"<<endl;
    }
    else
    {
	cerr << "Kept " << markers << " markers out of " << sites << " in panel." << endl;
	cerr << num_study << "/"<<num_sites<<" of study markers were in the sites file"<<endl;
    }
    cerr << "Calculating kinship values...";

#pragma omp parallel for	
	for(int j1=0;j1<Nsamples;j1++) 
	{
	    for(int j2=j1;j2<Nsamples;j2++) 
	    {
		float ibd0 = 0;
		float ibd1 = 0;
		float ibd2 = 0;
		float ibd3 = 0;
		float ks=-1;
		for(size_t i=0; i<bits[j1].size(); ++i)
		{
		    //opposite homozygotes. NAA,aa
		    ibd0 += (bits[j1][i][0] & bits[j2][i][2]).count() + (bits[j1][i][2] & bits[j2][i][0]).count();	
		    //same genotype.  NAA,AA + Naa,aa
		    ibd2 += (bits[j1][i][0] & bits[j2][i][0]).count() + (bits[j1][i][1] & bits[j2][i][1]).count() + (bits[j1][i][2] & bits[j2][i][2]).count();	
		    //missing in both.
		    ibd3 += (bits[j1][i][3] | bits[j2][i][3]).count();	
		}
		///consistent with IBD1 N - (NAA,AA + Naa,aa)
		ibd1 = markers-ibd3-ibd0-ibd2;

		if(method==0) 
		{
		    K.estimate_ibd(ibd0,ibd1,ibd2,ibd3);
		    ks = 0.5 * ibd2 + 0.25 * ibd1;
		}
		if(method==1)//king
		{
		    int Nhet_1=0,Nhet_2=0,Nhet_12=0;
		    for(size_t i=0; i<bits[j1].size(); ++i)
		    {
			Nhet_1 += bits[j1][i][1].count(); //NAa^i
			Nhet_2 += bits[j2][i][1].count(); //NAa^j
			Nhet_12 += (bits[j1][i][1] & bits[j2][i][1]).count(); //NAa,Aa
		    }
		    int minhet=min(Nhet_1,Nhet_2);
		    ks = (Nhet_12 - 2*ibd0)/(2*minhet) + 0.5 - 0.25*(Nhet_1+Nhet_2)/minhet;
		    K.estimate_ibd(ibd0,ibd1,ibd2,ibd3);
		}
	    
		if( !tk || ks > min_kin )
		{
#pragma omp critical
		    {		    
			cout  << names[j1] << "\t" << names[j2] << "\t" << left << " " << setprecision(5) << fixed << ibd0  << left << " " << setprecision(5) << fixed << ibd1  << left << " " << setprecision(5) << fixed << ibd2  << left << " " << setprecision(5) << fixed << ks << " " << setprecision(0) <<ibd3 << "\n";
		    }
		}
	    }
	}
    cerr << "done."<<endl;
    return 0;
}
