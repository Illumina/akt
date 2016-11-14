#include "akt.hpp"
#include "logs.hpp"
#include <math.h>       /* erfc */
#include <iomanip>

using namespace std;

/**
 * @name    usage
 * @brief   print out options
 *
 * List of input options
 *
 */
static void usage(){
  cerr << "Compare AFs between two cohorts" << endl;	
  cerr << "Usage:" << endl;
  cerr << "akt metafreq a.vcf.gz b.vcf.gz -Oz -o meta.sites.vcf.gz" << endl;
  umessage('r');
  umessage('R');
  umessage('T');
  umessage('o');
  umessage('O');
  umessage('a');
  exit(1);
}

//log(1!), log(2!), log(3!), ...
void update_cumulative_sum( vector<double> &log_sum, int top){
	
	if(top <= 1 || log_sum.size() < 2 ){ log_sum.resize(2); log_sum[0] = 0; log_sum[1] = 0; }
	
	for(int i=log_sum.size()-1; i<top; ++i){
		log_sum.push_back( log_sum.back() + log(i+1) );
	}
	
}



//contingency table
//a b
//c d
//precomputed log sums for hypergeometric distribution
//eps for finite precision comparison
double fisher_prob(int a, int b, int c, int d, vector<double> &log_cum, double eps=1e-8){

    //row sums
	int r1 = a+b;
	int r2 = c+d;
	//col sums
	int c1 = a+c;
	int c2 = b+d;
	//lowest marginal sum
	int least = min( min(r1,r2), min(c1,c2));
	int t = (least == r1) ? 0 : (least == r2) ? 1 : (least == c1 ) ? 2 : (least == c2) ? 3 : 4;
	if(t == 4){cerr << "failed to find min marginal in fisher table" << endl; exit(1); }

	//hypergeometric numerator
	double num = log_cum[r1] + log_cum[r2] + log_cum[c1] + log_cum[c2] - log_cum[r1+r2];
	//least 'extreme' probability
	double cutoff = num - (log_cum[a]+log_cum[b]+log_cum[c]+log_cum[d]) + eps;

	//Go through all the allowed tables computing their probability
	int i=0;	
	double prob = logz;
	if(t == 0){ //stretch row1
		while(i <= r1){
			double term = num - (log_cum[i] + log_cum[ r1-i ] + log_cum[c1-i] + log_cum[ r2-c1+i]);
			if(term <= cutoff){ log_sum( prob, term ); }
			++i;
		}      
	} else if(t == 1){ //stretch row2
		while(i <= r2){
			double term = num - (log_cum[i] + log_cum[ r2-i ] + log_cum[c2-i] + log_cum[ r1-c2+i]);
			if(term <= cutoff){ log_sum( prob, term ); }
			++i;	          
		}
	} else if(t == 2){ //stretch col1
		while(i <= c1){
			double term = num - (log_cum[i] + log_cum[ r1-i ] + log_cum[c1-i] + log_cum[ r2-c1+i]);
			if(term <= cutoff){ log_sum( prob, term ); }
			++i;
		
		}
	} else if(t == 3){ //stretch col2
		while(i <= c2){
			double term = num - (log_cum[i] + log_cum[ r2-i ] + log_cum[c2-i] + log_cum[ r1-c2+i]);
			if(term <= cutoff){ log_sum( prob, term );  }
			++i;	
		}
	} 
	return exp(prob);
	
}

int metafreq_main(int argc, char* argv[])
{

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
  string filename1 = argv[optind];	//file1
  optind++;
  string filename2 = argv[optind];	//file2
			
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
	bcf_hdr_merge(hdr, sr->readers[1].header);
	bcf_hdr_t *new_hdr = bcf_hdr_subset(hdr,0,NULL,NULL); ///creates a new subsetted header (with 0 samples) from src_header
	bcf_hdr_add_sample(new_hdr, NULL);      /// update internal structures

	bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "AF1,Number=A,Type=Float,Description=\"Allele frequency file1\">").c_str());
	bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "AC1,Number=A,Type=Integer,Description=\"Allele count file1\">").c_str());
	bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "AN1,Number=A,Type=Integer,Description=\"Allele number file1\">").c_str());
	bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "AF2,Number=A,Type=Float,Description=\"Allele frequency file2\">").c_str());
	bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "AC2,Number=A,Type=Integer,Description=\"Allele count file2\">").c_str());
	bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "AN2,Number=A,Type=Integer,Description=\"Allele number file2\">").c_str());
	bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "AF,Number=A,Type=Float,Description=\"Allele frequency combined\">").c_str());
	bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "QF,Number=A,Type=Float,Description=\"Q-score from fisher's exact test\">").c_str());
	bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "QX,Number=A,Type=Float,Description=\"Q-score from chi squared test\">").c_str());

	htsFile *out_fh  = hts_open( output.c_str(), outputfmt.c_str() );
	bcf_hdr_write(out_fh, new_hdr);
	bcf_hdr_destroy(hdr);

  
  int Nsamples1 = bcf_hdr_nsamples(sr->readers[0].header);	///number of samples in pedigree
  int Nsamples2 = bcf_hdr_nsamples(sr->readers[1].header);	///number of samples in pedigree
  cerr << Nsamples1 << " samples in " << filename1 << endl;
  cerr << Nsamples2 << " samples in " << filename2 << endl;

  //largest factorial required is sum of all elements in  contingency table
  vector<double> log_cum;

  bcf1_t *line1, *line2, *line_copy;///bcf/vcf line structure.

  int *gt_arr1=(int *)malloc(Nsamples1*2*sizeof(int)),ngt1=Nsamples1*2,ngt_arr1=Nsamples1*2;
  int *gt_arr2=(int *)malloc(Nsamples2*2*sizeof(int)),ngt2=Nsamples2*2,ngt_arr2=Nsamples2*2;
  int *ac_ptr=(int *)malloc(1*sizeof(int)); int nval = 1;
  int *an_ptr=(int *)malloc(1*sizeof(int)); 

  bcf1_t *rec = bcf_init1() ;

  int total = 0;
  int found1 = 0;
  int found2 = 0;
  int kept = 0;
  
  int ma1 = 0;
  int ma2 = 0;
  int both = 0;
  
  int np1 = 0;
  int np2 = 0;
  
  while(bcf_sr_next_line (sr)) { 
	  
	++total;
	
	int nmiss1 = 0, npres1 = 0, sum1 = 0;
	int nmiss2 = 0, npres2 = 0, sum2 = 0;

	bool masite = false;
	
	if( bcf_sr_has_line(sr,0) ){ //in first
		line1 =  bcf_sr_get_line(sr, 0);
		++found1;
		if( line1->n_allele == 2 ){ //biallelic

			ngt1 = bcf_get_genotypes(sr->readers[0].header, line1, &gt_arr1, &ngt_arr1);  
			if(ngt1 < 0){ //didn't read genotypes, look for AC, AN
				int ret = bcf_get_info_int32(sr->readers[0].header, line1, "AC", &ac_ptr, &nval);
				if(ret < 0){
					cerr << "No AC at " << bcf_hdr_id2name(sr->readers[0].header,line1->rid)<<":"<<line1->pos+1 << " in " << filename1 << endl; exit(1); 
				}
				ret = bcf_get_info_int32(sr->readers[0].header, line1, "AN", &an_ptr, &nval);
				if(ret < 0){
					cerr << "No AN at " << bcf_hdr_id2name(sr->readers[0].header,line1->rid)<<":"<<line1->pos+1 << " in " << filename1 << endl; exit(1); 
				}
				npres1 = an_ptr[0];
				sum1 = ac_ptr[0];
				nmiss1 = 0;
			} else { //read genotypes, calculate AC, AN
				for(int i=0;i<2*Nsamples1;++i){ 
					if(gt_arr1[i]==bcf_gt_missing || (bcf_gt_allele(gt_arr1[i])<0) || (bcf_gt_allele(gt_arr1[2])>2)){
						++nmiss1;
					} else {
						sum1 += bcf_gt_allele(gt_arr1[i]);
						++npres1;
					}
				}
			}
			//cout << filename1 << " " << bcf_hdr_id2name(sr->readers[0].header,line1->rid)<<":"<<line1->pos+1 << " "
			//<< npres1 << " " << sum1 << endl;
		} else { //not biallelic
			++ma1;
			masite = true;
		}
	}
	
	if( bcf_sr_has_line(sr,1) ){ //in second
		line2 =  bcf_sr_get_line(sr, 1);
		++found2;
		if( line2->n_allele == 2 ){//biallelic

			ngt2 = bcf_get_genotypes(sr->readers[1].header, line2, &gt_arr2, &ngt_arr2);  
			if(ngt2 < 0){ //didn't read genotypes, look for AC, AN

				int ret = bcf_get_info_int32(sr->readers[1].header, line2, "AC", &ac_ptr, &nval);
				if(ret < 0){
					cerr << "No AC at " << bcf_hdr_id2name(sr->readers[1].header,line2->rid)<<":"<<line2->pos+1 << " in " << filename2 << endl; exit(1); 
				}
				ret = bcf_get_info_int32(sr->readers[1].header, line2, "AN", &an_ptr, &nval);
				if(ret < 0){
					cerr << "No AN at " << bcf_hdr_id2name(sr->readers[1].header,line2->rid)<<":"<<line2->pos+1 << " in " << filename2 << endl; exit(1); 
				}
				npres2 = an_ptr[0];
				sum2 = ac_ptr[0];
				nmiss2 = 0;
			} else { //read genotypes, calculate AC, AN
				for(int i=0;i<2*Nsamples2;++i){ 
					if(gt_arr2[i]==bcf_gt_missing || (bcf_gt_allele(gt_arr2[i])<0) || (bcf_gt_allele(gt_arr2[i])>2)){
						++nmiss2;
					} else {
						sum2 += bcf_gt_allele(gt_arr2[i]);
						++npres2;
					}
				}
			}
		} else { //multiallelic
			++ma2;
			masite = true;
		}
	}
	
	int uh = 0; //which header to use
	if( bcf_sr_has_line(sr,0) && bcf_sr_has_line(sr,1) && !masite ){ //in both and not multiallelic
		++both; line_copy =  bcf_sr_get_line(sr, 0); uh = 0; 
	} 
	else if( bcf_sr_has_line(sr,0) && !bcf_sr_has_line(sr,1) && !masite ){ //in first and not multiallelic
		npres2 = np2; sum2 = 0; line_copy =  bcf_sr_get_line(sr, 0); uh = 0;
	}
	else if( !bcf_sr_has_line(sr,0) && bcf_sr_has_line(sr,1) && !masite ){ //in second and not multiallelic
		npres1 = np1; sum1 = 0; line_copy =  bcf_sr_get_line(sr, 1); uh = 1;
	} 
	else {
		masite = true;
	}
	if(!masite && npres1>0 && npres2>0){	//can test site
		
		//AF
		double p1 = (double)(sum1)/(double)(npres1);
		double p2 = (double)(sum2)/(double)(npres2);
		
		float QF=-1;
		float QX=-1;
		float pav = (npres1*p1 + npres2*p2)/(double)(npres1 + npres2);	//average AF
		float p1f = p1;
		float p2f = p2;

		if(sum1 == 0 && sum2 == 0){ QF = 0; QX = 0; }
		else{
			
			//Fisher
			//      	file1 		file2
			//novar 	n1-sum1 	n2-sum2		| n1+n2 - sum1 - sum2
			//var		sum1  		sum2		| sum1+sum2
			//			----------------------------
			//			n1			n2			| n1+n2	
			update_cumulative_sum( log_cum , npres1+npres2); 							
			double pval = fisher_prob( (npres1-sum1),(npres2-sum2), sum1, sum2, log_cum);
			QF = min(100.0, -log10(pval) );
			
			
			//Chisquared test, with Yates's correction
			double N = npres1 + npres2;
			double a = (npres1-sum1);
			double b = (npres2-sum2);
			double c = (sum1);
			double d = (sum2);
			double c1 = npres1;
			double c2 = npres2;
			double r2 = sum1+sum2;
			double r1 = N - r2;
			double chi = max(0.0,  fabs( a*d - b*c ) - N*0.5  ) * sqrt( N/(r1*r2*c1*c2) ) ;

			//chisq CDF is
			//\gamma(k/2, x/2) / \Gamma(k/2) for k dof
			//2x2 table has 1 dof
			//\Gamma(1/2) = sqrt(pi)
			//\gamma(1/2,x/2) = sqrt(pi) erf( sqrt(x/2) )
			double cpval = 1 - erf( chi*M_SQRT1_2 );
			QX = min(100.0, -log10(cpval) );
		}
		
		//copy line record 'essentials'
		rec->rid = bcf_hdr_name2id( new_hdr, bcf_hdr_id2name(sr->readers[uh].header, line_copy->rid) ); 	///two headers may have different representation of chr in hash table
		rec->pos = line_copy->pos;
		rec->qual = line_copy->qual;
		rec->rlen = line_copy->rlen;

		bcf_update_id(new_hdr, rec, line_copy->d.id);
		string ref = line_copy->d.allele[0];
		string alt = line_copy->d.allele[1];
		string alleles = ref + "," + alt;
		bcf_update_alleles_str(new_hdr, rec, alleles.c_str());

		//add AF annotations
		bcf_update_info_float(new_hdr, rec, (af_tag + "AF1").c_str(), &p1f, 1);	
		bcf_update_info_int32(new_hdr, rec, (af_tag + "AC1").c_str(), &sum1, 1);	
		bcf_update_info_int32(new_hdr, rec, (af_tag + "AN1").c_str(), &npres1, 1);	
		bcf_update_info_float(new_hdr, rec, (af_tag + "AF2").c_str(), &p2f, 1);	
		bcf_update_info_int32(new_hdr, rec, (af_tag + "AC2").c_str(), &sum2, 1);	
		bcf_update_info_int32(new_hdr, rec, (af_tag + "AN2").c_str(), &npres2, 1);	
		bcf_update_info_float(new_hdr, rec, (af_tag + "AF").c_str(), &pav, 1);	

		//cout << "Writing " << line_copy->pos+1 << " " << npres1 << " " << npres2 << " " << sum1 << " " << sum2 << " " << p1f << " " << p2f << " " << QF << " " << QX << endl;

		if(QF <= 0){QF = 0;} //-0 annoying
		bcf_update_info_float(new_hdr, rec, (af_tag + "QF").c_str(), &QF, 1);	
		if(QX <= 0){QX = 0;} //-0 annoying
		bcf_update_info_float(new_hdr, rec, (af_tag + "QX").c_str(), &QX, 1);	
				
		bcf_unpack(rec, BCF_UN_ALL);		
		bcf_write1(out_fh, new_hdr, rec);
		bcf_clear1(rec);	

		++kept;
		
		
		np1 = npres1;
		np2 = npres2;
		
		
	}
  }//reader
  
  bcf_destroy1(rec);
  free(gt_arr1);
  free(gt_arr2);
  free(ac_ptr);
  free(an_ptr);
  bcf_sr_destroy(sr);	
  hts_close(out_fh);
  bcf_hdr_destroy(new_hdr);
  
  cerr << "Tested " << kept << " sites out of " << total << endl;
  cerr << "Found " << found1 << " sites in " << filename1 << " and skipped " << ma1 << " multiallelic ones" << endl;
  cerr << "Found " << found2 << " sites in " << filename2 << " and skipped " << ma2 << " multiallelic ones" << endl;
  cerr << both << " biallelic sites in both files." << endl;

 
  return 0;
}
