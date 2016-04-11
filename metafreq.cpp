#include "akt.hpp"
#include "logs.hpp"
#include <math.h>       /* erfc */

using namespace std;

static void usage(){
  cerr << "Compare AFs between two cohorts" << endl;	
  cerr << "Usage:" << endl;
  cerr << "akt metafreq a.sites.vcf.gz b.sites.vcf.gz -Oz -o meta.sites.vcf.gz" << endl;
  cerr << "Expects input.bcf to contain genotypes." << endl;
  exit(1);
}

double normalCDF(double x)
{
   return 0.5 * erfc(-x * M_SQRT1_2);
}

void cumulative_sum( vector<double> &log_sum){
	
	log_sum[0] = 0;
	log_sum[1] = 0;
	for(size_t i=2; i<=log_sum.size(); ++i){
		log_sum[i] = log_sum[i-1] + log(i);
	}
	
}

inline double log_binomialP(double k, unsigned long n, double p, vector<double> &log_cum){
    return log_cum[n] - log_cum[k] - log_cum[n-k] + k*log(p) + (n-k)*log(1.0-p);
} 

inline double binomialP(double k, unsigned long n, double p, vector<double> &log_cum){
	if(k > n || k < 0){return 0;}
    return exp( log_binomialP(k, n, p, log_cum) );
} 

double binomialCDF_leq(double k, unsigned long n, double p, vector<double> &log_cum){
	if(k > n || k < 0){return 0;}
	double sum = logz;
	for(int i=0; i<=k; ++i){ log_sum(sum, log_binomialP(i,n,p,log_cum) ); }
    return exp(sum);
} 
double binomialCDF_geq(double k, unsigned long n, double p, vector<double> &log_cum){
	if(k > n || k < 0){return 0;}
	double sum = logz;
	for(int i=k; i<=n; ++i){ log_sum(sum, log_binomialP(i,n,p,log_cum) ); }
    return exp(sum);
} 

int metafreq_main(int argc, char* argv[])
{

  /*vector<double> log_c(300,0);
  cumulative_sum(log_c);

  cout << binomialCDF_geq( 2*(100*0.2) - 7,100,0.2, log_c) << endl; 
  cout << binomialCDF_leq(7,100,0.2, log_c) << endl; 
  
  cout << binomialP( 147, 192, 0.994792, log_c) << endl;
  cout << binomialCDF_geq( 147, 192, 0.994792, log_c) << endl;
  
  exit(1);*/
		
  /*cout << normalCDF(0) << " " << log10( 2*normalCDF(0) )<< endl; 
  cout << normalCDF(-1) << " " << log10( 2*normalCDF(-1) )<< endl; 
  cout << normalCDF(-2) << " " << log10( 2*normalCDF(-2) )<< endl; 
  cout << normalCDF(-3) << " " << log10( 2*normalCDF(-3) )<< endl; 
  cout << normalCDF(-4) << " " << log10( 2*normalCDF(-4) )<< endl; 
  cout << normalCDF(-5) << " " << log10( 2*normalCDF(-5) )<< endl; 
  exit(1);*/
  			
  int c;
  
  if(argc<3) usage();
  static struct option loptions[] =    {
    {"regions",1,0,'r'},	
	{"regions-file",1,0,'R'},
	{"targets-file",1,0,'T'},
    {"output",1,0,'o'},
    {"outputfmt",1,0,'O'},
    {"aftag",1,0,'a'},
    {0,0,0,0}
  };

  string regions = "";
  string pfilename = "";

  bool regions_is_file = false;
  bool used_r = false;
  bool used_R = false;
  bool target = false;
  string output = "out.vcf";
  string outputfmt = "w";
  string af_tag = "";

  while ((c = getopt_long(argc, argv, "r:R:T:o:O:a:",loptions,NULL)) >= 0) {  
    switch (c)
      {
      case 'r': regions = (optarg); used_r = true; break;
      case 'R': regions = (optarg); pfilename = (optarg); used_R = true; regions_is_file = true; break;
      case 'T': target = true; pfilename = (optarg); break;
      case 'o': output = (string)(optarg); break;
      case 'O': outputfmt += (string)(optarg); break;
      case 'a': af_tag = string(optarg)+"_"; break;
      case '?': usage();
      default: cerr << "Unknown argument:"+(string)optarg+"\n" << endl; exit(1);
      }
  }

  if( used_r && used_R ){ cerr << "-r and -R cannot be used simultaneously" << endl; exit(1); }

  optind++;
  string filename1 = argv[optind];
  optind++;
  string filename2 = argv[optind];
			
  int Nsamples;
  
  bcf_srs_t *sr =  bcf_sr_init() ; ///htslib synced reader.
  sr->collapse = COLLAPSE_NONE;
  sr->require_index = 1;
  
  if(regions != ""){
		if ( bcf_sr_set_regions(sr, regions.c_str(), regions_is_file)<0 ){
			cerr << "Failed to read the regions: " <<  regions << endl; return 0;
		}
  }
	
  if(!(bcf_sr_add_reader (sr, filename1.c_str() ))){ 
    cerr << "Problem opening " << filename1 << endl; 
    cerr << "Input file not found." << endl;
    bcf_sr_destroy(sr);	
    return 0;
  }
  if(!(bcf_sr_add_reader (sr, filename2.c_str() ))){ 
	cerr << "Problem opening " << filename2 << endl; 
	cerr << "Sites file not found." << endl;
	bcf_sr_destroy(sr);	
	return 0;
  }
  
  //Setup output
	bcf_hdr_t *hdr = bcf_hdr_dup(sr->readers[0].header);
	bcf_hdr_t *new_hdr = bcf_hdr_subset(hdr,0,NULL,NULL); ///creates a new subsetted header (with 0 samples) from src_header
	bcf_hdr_add_sample(new_hdr, NULL);      /// update internal structures

	bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "AF1,Number=A,Type=Float,Description=\"Allele frequency file1\">").c_str());
	bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "AF2,Number=A,Type=Float,Description=\"Allele frequency file2\">").c_str());
	bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "AF,Number=A,Type=Float,Description=\"Allele frequency combined\">").c_str());
	bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "QB,Number=A,Type=Float,Description=\"Q-score for different distributions - binomial variance\">").c_str());
	bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "QS,Number=A,Type=Float,Description=\"Q-score for different distributions - sample variance\">").c_str());
	bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "QE,Number=A,Type=Float,Description=\"Q-score exact binomial test\">").c_str());

	htsFile *out_fh  = hts_open( output.c_str(), outputfmt.c_str() );
	bcf_hdr_write(out_fh, new_hdr);
	
  
  int Nsamples1 = bcf_hdr_nsamples(sr->readers[0].header);	///number of samples in pedigree
  int Nsamples2 = bcf_hdr_nsamples(sr->readers[1].header);	///number of samples in pedigree
  cerr << Nsamples1 << " samples in " << filename1 << endl;
  cerr << Nsamples2 << " samples in " << filename2 << endl;

  vector<double> log_cum(2*max(Nsamples1, Nsamples2),0);
  cumulative_sum( log_cum );


  bcf1_t *line1, *line2;///bcf/vcf line structure.
  bcf1_t *rec = bcf_init1() ;

  int *gt_arr1=(int *)malloc(Nsamples1*2*sizeof(int)),ngt1=Nsamples1*2,ngt_arr1=Nsamples1*2;
  int *gt_arr2=(int *)malloc(Nsamples2*2*sizeof(int)),ngt2=Nsamples2*2,ngt_arr2=Nsamples2*2;
  
  int total = 0;
  int kept = 0;
  
  				int zt = 0;
				int bt = 0;
				
  while(bcf_sr_next_line (sr)) { 
	
	if( bcf_sr_has_line(sr,0) && bcf_sr_has_line(sr,1) ){ //in both
			line1 =  bcf_sr_get_line(sr, 0);
			line2 =  bcf_sr_get_line(sr, 1);
			if( line1->n_allele == 2 && line2->n_allele == 2 ){

				ngt1 = bcf_get_genotypes(sr->readers[0].header, line1, &gt_arr1, &ngt_arr1);  
				if(ngt1 < 0){ cerr << "Bad genotypes at " << bcf_hdr_id2name(sr->readers[0].header,line1->rid)<<":"<<line1->pos+1 << " in " << filename1 << endl; exit(1); }
				ngt2 = bcf_get_genotypes(sr->readers[1].header, line2, &gt_arr2, &ngt_arr2);  
				if(ngt2 < 0){ cerr << "Bad genotypes at " << bcf_hdr_id2name(sr->readers[1].header,line2->rid)<<":"<<line2->pos+1 << " in " << filename2 << endl; exit(1); }
					
				int nmiss1 = 0, npres1 = 0;
				double sum1 = 0, sumsq1 = 0;
				for(int i=0;i<Nsamples1;++i){ 
					if(
					(gt_arr1[2*i]==bcf_gt_missing || (bcf_gt_allele(gt_arr1[2*i])<0) || (bcf_gt_allele(gt_arr1[2*i])>2) ) ||
					(gt_arr1[2*i+1]==bcf_gt_missing || (bcf_gt_allele(gt_arr1[2*i+1])<0) || (bcf_gt_allele(gt_arr1[2*i+1])>2) ) 
					){
						++nmiss1;
					} else {
						double gt = (bcf_gt_allele(gt_arr1[2*i]) + bcf_gt_allele(gt_arr1[2*i+1]))/2.0;
						sum1 += gt;
						sumsq1 += gt*gt;
						++npres1;
					}
				}
				
				int nmiss2 = 0, npres2 = 0;
				double sum2 = 0, sumsq2 = 0;
				for(int i=0;i<Nsamples2;++i){ 
					if(
					(gt_arr2[2*i]==bcf_gt_missing || (bcf_gt_allele(gt_arr2[2*i])<0) || (bcf_gt_allele(gt_arr2[2*i])>2) ) ||
					(gt_arr2[2*i+1]==bcf_gt_missing || (bcf_gt_allele(gt_arr2[2*i+1])<0) || (bcf_gt_allele(gt_arr2[2*i+1])>2) ) 
					){
						++nmiss2;
					} else {
						double gt = (bcf_gt_allele(gt_arr2[2*i]) + bcf_gt_allele(gt_arr2[2*i+1]))/2.0;
						sum2 += gt;
						sumsq2 += gt*gt;
						++npres2;
					}
				}
				
				double p1 = (double)(sum1)/(double)(npres1);
				double p2 = (double)(sum2)/(double)(npres2);
				
				float QB=-1;
				float QS=-1;
				float QE=-1;
				float pav = (npres1*p1 + npres2*p2)/(double)(npres1 + npres2);
				float p1f = p1;
				float p2f = p2;
					

				//z-test should work in this case			
				if(npres1 * p1 > 5 && npres2 * p2 > 5 && npres1 * (1 - p1) > 5 && npres2 * (1 - p2) > 5){

					double bvar1 = p1*(1-p1);
					double bvar2 = p2*(1-p2);
					double zscore = (p1 - p2)/sqrt( (bvar1/npres1) + (bvar2/npres2) );
					
					double svar1 = ((double)(sumsq1)/(double)(npres1)) - p1*p1;
					double svar2 = ((double)(sumsq2)/(double)(npres2)) - p2*p2;
					double szscore = (p1 - p2)/sqrt( (svar1/npres1) + (svar2/npres2) );
					
					QB = -log10( 2*normalCDF( -fabs(zscore) ) );
					QS = -log10( 2*normalCDF( -fabs(szscore) ) );
					++zt;
				} //else {	//calculate binomial test
				double p = 0; 
				
				if( p1 == 0 ){ p = p2; }
				else if (p2 == 0){ p = p1; }
				else if (p1 == 1){ p = p2; }
				else if (p2 == 1){ p = p1; }
				else if(npres1 == npres2){ p = (p1 > p2) ? p1 : p2; } 
				else{ p = (npres1 > npres2) ? p1 : p2; }
				if(p != 0 && p != 1){ 
				
					int expectation, trials, successes;
					if( (p==p1) ){
						expectation = 2*round(sum1);
						trials = 2*npres2;
						successes = 2*sum2;
					} else {
						expectation = 2*round(sum2);
						trials = 2*npres1;
						successes = 2*sum1;
					}
		
					//cout << expectation << " " << trials << " " << successes << " " << p << " :: ";
					
					double pval;
					if(successes > expectation){ //got more
						pval = binomialCDF_geq( successes,trials,p, log_cum); //one tail
						pval += binomialCDF_leq( expectation - (successes-expectation),trials,p, log_cum); //2nd tail 
					} else if(successes < expectation) { //got less 
						pval = binomialCDF_geq(expectation + (expectation - successes),trials,p, log_cum); 
						pval += binomialCDF_leq(successes,trials,p, log_cum);
					} else {	//got exactly what we expected
						pval = binomialP(successes,trials,p, log_cum);
					}
					QE = -log10(pval);
					cout << bcf_hdr_id2name(sr->readers[0].header,line1->rid)<<":"<<line1->pos+1 << " " << QE << endl;
				}
				++bt;
			
				//copy line record 'essentials'
				rec->rid = line1->rid;
				rec->pos = line1->pos;
				rec->qual = line1->qual;
				rec->rlen = line1->rlen;

				bcf_update_id(new_hdr, rec, line1->d.id);
				string ref = line1->d.allele[0];
				string alt = line1->d.allele[1];
				string alleles = ref + "," + alt;
				bcf_update_alleles_str(new_hdr, rec, alleles.c_str());
		
				//add AF annotations
				bcf_update_info_float(new_hdr, rec, (af_tag + "AF1").c_str(), &p1f, 1);	
				bcf_update_info_float(new_hdr, rec, (af_tag + "AF2").c_str(), &p2f, 1);	
				bcf_update_info_float(new_hdr, rec, (af_tag + "AF").c_str(), &pav, 1);	
				if(QB > 0)bcf_update_info_float(new_hdr, rec, (af_tag + "QB").c_str(), &QB, 1);	
				if(QS > 0)bcf_update_info_float(new_hdr, rec, (af_tag + "QS").c_str(), &QS, 1);	
				if(QE > 0)bcf_update_info_float(new_hdr, rec, (af_tag + "QE").c_str(), &QE, 1);	

				bcf_unpack(rec, BCF_UN_ALL);				
				bcf_write1(out_fh, new_hdr, rec);
				bcf_clear1(rec);	
					
				++kept;
			}//biallelic
		} //in both
		++total;
  }//reader
  
  cerr << zt << " " << bt << endl;
  free(gt_arr1);
  free(gt_arr2);
  bcf_sr_destroy(sr);	
  hts_close(out_fh);
  bcf_hdr_destroy(new_hdr);
  bcf_destroy1(rec);
	
  cerr << "Kept " << kept << " markers out of " << total << endl;
 
  return 0;
}
