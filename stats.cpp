#include "akt.hpp"
#include <Eigen/Dense>

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
	cerr << "Calculate correlation across a population" << endl;	
	cerr << "Usage:" << endl;
	cerr << "./akt stats input_filename.vcf" << endl;
	cerr << "Expects input_filename.vcf to contain hard genotypes" << endl;
	cerr << "\t -F --flank:			size of left and right flanking region (1000bp)" << endl;
	cerr << "\t -b --block:			-F argument is number of markers instead of number of base pairs" << endl;
	cerr << "\t -x --afonly:			calculate allele freq only" << endl;
	cerr << "\t -c --output_cor:		output sitewise correlation (false)" << endl;
	cerr << "\t -C --output_cormin:		output sitewise correlation greater than (0)" << endl;
	umessage('a');
	umessage('o');
	umessage('O');
	umessage('r');
	umessage('R');
	umessage('s');
	umessage('S');
	exit(1);
}

	
//read VCF record and calculate AF and VAR
static void process_line(bcf_srs_t *sr, bcf1_t *line, 
					int N, int *gt_arr, int ngt, int ngt_arr,
					list< vector<int> > &G, list< float > &af, list< float > &sig, list< int > &position)
{
	 
	ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr);  
	for(int i=0;i<ngt_arr;i++){ 
		if(gt_arr[i]==bcf_gt_missing || (bcf_gt_allele(gt_arr[i])<0) || (bcf_gt_allele(gt_arr[i])>2)  ){
			gt_arr[i]=bcf_gt_phased(0);//inserting homref instead of missing.
		}
	}

	float p = 0;	//mean
	vector<int> tmp(N); 
	for(int i=0; i<N; ++i){ tmp[i] = bcf_gt_allele(gt_arr[2*i]) + bcf_gt_allele(gt_arr[2*i+1]); p += (float)tmp[i]; } 
	p /= N; 
	float s = 0;	//variance
	for(int i=0; i<N; ++i){ s += ( (float)tmp[i] - p )*( (float)tmp[i] - p ); } 
	s /= N;
        
	G.push_back(tmp);
	af.push_back( p );
	sig.push_back( s );
	position.push_back( line->pos+1 );
			
}

//calculate correlation from neighbouring site data and save in bcf record 
static void process_output(
bcf_hdr_t *new_hdr,
bcf1_t *line, bcf1_t *rec, 
int N, 
list< vector<int> >::iterator G_it, list<float>::iterator af_it, list<float>::iterator sig_it, list<int>::iterator calc_it,
list< vector<int> > &G, list< float > &af, list< float > &sig, list< int > &position,
string af_tag, bool xsig, bool output_cor, 
htsFile *out_fh, float ocmin
){
	//copy line record 'essentials'
	rec->rid = line->rid;
	rec->pos = line->pos;
	rec->qual = line->qual;
	rec->rlen = line->rlen;

	bcf_update_id(new_hdr, rec, line->d.id);
	string ref = line->d.allele[0];
	string alt = line->d.allele[1];
	string alleles = ref + "," + alt;
	bcf_update_alleles_str(new_hdr, rec, alleles.c_str());
	
	//update allele freq
	float taf = *af_it/2.0; //af[sites]/2.0;
	bcf_update_info_float(new_hdr, rec, (af_tag + "AF").c_str(), &taf, 1);	
	float tvar = *sig_it/4.0; //sig[sites]/4.0; //var( c x ) = c^2 var(x)
	bcf_update_info_float(new_hdr, rec, (af_tag + "SIG").c_str(), &tvar, 1);	
	
	if( !xsig && !position.empty() && *sig_it > 0){	//calculate cor = true && we have data && we have some variants

		vector<float> cor;	//correlation values
		vector<int> corp;	//correlation positions

		float LDscore = 0;	
		
		list<int>::iterator fin = --position.end();	//stop here
		
		list<int>::iterator it=position.begin(); 
		list< vector<int> >::iterator git=G.begin(); 
		list<float>::iterator ait = af.begin(); 
		list<float>::iterator sit = sig.begin();

		//calculate correlation
		for (; it != fin; ++it, ++git, ++ait, ++sit){
			if( *sit > 0 ){ 
				cor.push_back(0);
				corp.push_back( *it );

				for(int j=0; j<N; ++j){
					cor.back() += ( (*git)[j] -(*ait) )*( (*G_it)[j] - (*af_it) );
				}
				cor.back() /= (float)N;
				cor.back() /= sqrt( (*sit)*(*sig_it) );
				
				float r2 = cor.back()*cor.back();
				if(N > 2){ r2 -= (1.0-r2)/(N-2.0); } //unbiased estimator
				LDscore += r2;
			}
		}
		bcf_update_info_float(new_hdr, rec, (af_tag + "LD").c_str(), &LDscore, 1);	

		//save raw correlation
		if(output_cor){
			vector<float> cor_red;
			vector<int> corp_red;
			for(size_t i=0; i<cor.size(); ++i){
				if( fabs(cor[i]) >= ocmin ){
					cor_red.push_back(cor[i]);
					corp_red.push_back(corp[i]);
				}
			}
			bcf_update_info_float(new_hdr,rec,(af_tag + "COR").c_str(),&cor_red[0],cor_red.size() );
			bcf_update_info_int32(new_hdr,rec,(af_tag + "CORP").c_str(),&corp_red[0],corp_red.size() );
		}

		
	}	
	bcf_unpack(rec, BCF_UN_ALL);				
	bcf_write1(out_fh, new_hdr, rec);
	bcf_clear1(rec);
}
int stats_main(int argc, char* argv[])
{
	int c;
    
    if(argc<3) usage();
    static struct option loptions[] =    {
        {"afonly",1,0,'x'},
        {"output_cor",1,0,'c'},
        {"output_cormin",1,0,'C'},
        {"flank",1,0,'F'},
        {"block",0,0,'b'},
        {"aftag",1,0,'a'},
        {"regions",1,0,'r'},	
		{"regions-file",1,0,'R'},
		{"samples",1,0,'s'},
		{"samples-file",1,0,'S'},
		{"output",1,0,'o'},
        {"outputfmt",1,0,'O'},
        {0,0,0,0}
    };

    bool xsig = false;
    bool output_cor = false;
    int flank = 1000;
	bool do_block = false; 
  
  string output = "out.vcf";
  string outputfmt = "w";
  string regions = "";
  bool regions_is_file = false;

  sample_args sargs;
  string af_tag = "";
  float ocmin = 0;
  
    while ((c = getopt_long(argc, argv, "cC:xF:a:o:O:s:S:r:R:b",loptions,NULL)) >= 0) {  
        switch (c)
        {
        case 'F': flank = atoi(optarg); break;
        case 'b': do_block = true; break;
        case 'x': xsig = true; break;
        case 'c': output_cor = true; break;
        case 'C': ocmin = atof(optarg); break;
        case 'o': output = (string)(optarg); break;
        case 'O': outputfmt += (string)(optarg); break;
        case 'a': af_tag = string(optarg)+"_"; break;
		case 'r': regions = (optarg); break;
		case 'R': regions = (optarg); regions_is_file = true; break;
        case 's': sargs.sample_names = (optarg); sargs.subsample = true; break;
        case 'S': sargs.sample_names = (optarg); sargs.subsample = true; sargs.sample_is_file = 1; break;
        case '?': usage();
        default: cerr << "Unknown argument:"+(string)optarg+"\n" << endl; exit(1);
        }
    }

    if( !xsig ){ 
		if(do_block){ cerr << "Calculating LD within " << flank << " variants" << endl; }
		else{ cerr << "Calculating LD within " << flank << "bp" << endl; }
	}
	optind++;
    string filename = argv[optind];
	
	//reads same file twice, F sites out of sync
	//Setup htslib reader OUTER
	bcf_srs_t *sr =  bcf_sr_init() ; ///htslib synced reader.
	sr->collapse = COLLAPSE_ANY;
	if(regions != ""){
		sr->require_index = 1;
		if ( bcf_sr_set_regions(sr, regions.c_str(), regions_is_file)<0 ){
			cerr << "Failed to read the regions: " <<  regions << endl; exit(1);
		}
	}
	if(!(bcf_sr_add_reader (sr, filename.c_str() ))){ cerr << "Problem opening " << filename << endl; exit(1); }
	//Setup htslib reader INNER
	bcf_srs_t *reader =  bcf_sr_init() ; ///htslib synced reader.
	reader->collapse = COLLAPSE_ANY;
	if(regions != ""){
		reader->require_index = 1;
		if ( bcf_sr_set_regions(reader, regions.c_str(), regions_is_file)<0 ){
			cerr << "Failed to read the regions: " <<  regions << endl; exit(1);
		}
	}
	if(!(bcf_sr_add_reader (reader, filename.c_str() ))){ cerr << "Problem opening " << filename << endl; exit(1); }
	//Setup output
	bcf_hdr_t *hdr = bcf_hdr_dup(reader->readers[0].header);
	bcf_hdr_t *new_hdr = bcf_hdr_subset(hdr,0,NULL,NULL); ///creates a new subsetted header (with 0 samples) from src_header
	bcf_hdr_add_sample(new_hdr, NULL);      /// update internal structures

	bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "AF,Number=A,Type=Float,Description=\"Allele frequency\">").c_str());
	bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "SIG,Number=A,Type=Float,Description=\"Variance\">").c_str());
	if(!xsig){ 
		bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "LD,Number=A,Type=Float,Description=\"LD Score\">").c_str());
		if(output_cor){
			bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "COR,Number=.,Type=Float,Description=\"Correlation\">").c_str()); 
			bcf_hdr_append(new_hdr, ("##INFO=<ID=" + af_tag + "CORP,Number=.,Type=Float,Description=\"Correlation Position\">").c_str()); 
		}
	}
	htsFile *out_fh  = hts_open( output.c_str(), outputfmt.c_str() );
	bcf_hdr_write(out_fh, new_hdr);
	
    if(sargs.subsample){ 
		bcf_hdr_set_samples(sr->readers[0].header, sargs.sample_names, sargs.sample_is_file); 
	}
		
	int N = bcf_hdr_nsamples(sr->readers[0].header);	///number of samples	
	cerr << N << " samples in " << filename << endl;
	
	bcf1_t *line;///bcf/vcf line structure.
	bcf1_t *rec = bcf_init1() ;

	int *gt_arr=(int *)malloc(N*2*sizeof(int)),ngt=N*2,ngt_arr=N*2;	

	//should contain just enough data in front and behind
	list< vector<int> > G; 
	list< float > af;
	list< float > sig;
	list<int> position;

	list< vector<int> >::iterator G_it;
	list<float>::iterator af_it;
	list<float>::iterator sig_it;
	list<int>::iterator calc_it;

	int chr_id = -1;
	
	while(bcf_sr_next_line (sr) ) { 
			
		line =  bcf_sr_get_line(sr, 0);
		
		if( line->n_allele == 2 && bcf_is_snp(line) ){ //biallelic

		if( xsig ){	//AF only case
			
			float p = 0;
			ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr); 
			for(int i=0; i<2*N; ++i){ p += (float)(bcf_gt_allele(gt_arr[i])); } 
			p /= (2*N); 
			//copy line record 'essentials'
			rec->rid = line->rid;
			rec->pos = line->pos;
			rec->qual = line->qual;
			rec->rlen = line->rlen;

			bcf_update_id(new_hdr, rec, line->d.id);
			string ref = line->d.allele[0];
			string alt = line->d.allele[1];
			string alleles = ref + "," + alt;
			bcf_update_alleles_str(new_hdr, rec, alleles.c_str());
	
			bcf_update_info_float(new_hdr, rec, (af_tag + "AF").c_str(), &p, 1);	

			bcf_unpack(rec, BCF_UN_ALL);				
			bcf_write1(out_fh, new_hdr, rec);
			bcf_clear1(rec);	
		} else {	

			if( chr_id != line->rid ){ //new chromosome
				
				if( chr_id != -1 ){ //not the first time
					
					int last_pos = position.back();

					//Finish off the right edge of previous chromosome
					while( bcf_sr_next_line (reader) ){

						line =  bcf_sr_get_line(reader, 0);
						
						if( line->n_allele == 2 && bcf_is_snp(line) ){ //biallelic
			
						process_output( new_hdr, line, rec, N, G_it, af_it, sig_it, calc_it, 
						G, af, sig, position, af_tag, xsig, output_cor, out_fh, ocmin); 
						
						++calc_it; ++af_it; ++sig_it; ++G_it;
						while( calc_it != position.end() && 
						((do_block && distance( position.begin(), calc_it ) > flank ) ||
						(!do_block && *calc_it - position.front() > flank ))
						){ //within flank
							G.pop_front();
							af.pop_front();
							sig.pop_front();
							position.pop_front();
						}
						if(line->pos+1 >= last_pos){ break; } //don't read past the end of the chromosome
						}
					}
					position.clear(); //throw away everything
					af.clear(); sig.clear(); G.clear();
				}
				line =  bcf_sr_get_line(sr, 0);

				chr_id = line->rid; 
				process_line( sr, line, N, gt_arr, ngt, ngt_arr, G, af, sig, position );
				
				calc_it = position.begin(); af_it = af.begin(); sig_it = sig.begin(); G_it = G.begin();
	
			} else {
				process_line( sr, line, N, gt_arr, ngt, ngt_arr, G, af, sig, position );
			}
			//read enough to calculate next site
				while( 
					(
					(do_block && distance( calc_it, position.end() ) > flank ) ||
					(!do_block && position.back() - *calc_it > flank) 
					) &&
					bcf_sr_next_line (reader)
				){

					line =  bcf_sr_get_line(reader, 0);
					if( line->n_allele == 2 && bcf_is_snp(line) ){ //biallelic
					

					process_output( new_hdr, line, rec, N, G_it, af_it, sig_it, calc_it, 
					G, af, sig, position, af_tag, xsig, output_cor, out_fh, ocmin); 
								
					++calc_it; ++af_it; ++sig_it; ++G_it;
					while( calc_it != position.end() && 
					((do_block && distance( position.begin(), calc_it ) > flank ) ||
					(!do_block && *calc_it - position.front() > flank ))
					){ 
						G.pop_front();
						af.pop_front();
						sig.pop_front();
						position.pop_front();		
					}
					}
				}
			
		} //xsig false
		} //biallelic
	}


	//Done most of the reading now finish off the right edge of last chromosome
	while( !xsig && bcf_sr_next_line (reader) && !position.empty() ){

		line =  bcf_sr_get_line(reader, 0);
				
		if( line->n_allele == 2 && bcf_is_snp(line) ){ //biallelic
					

		process_output( new_hdr, line, rec, N, G_it, af_it, sig_it, calc_it, 
		G, af, sig, position, af_tag, xsig, output_cor, out_fh, ocmin); 

		++calc_it; ++af_it; ++sig_it; ++G_it;
		while( calc_it != position.end() && 
		((do_block && distance( position.begin(), calc_it ) > flank ) ||
		(!do_block && *calc_it - position.front() > flank ))
		){ 
			G.pop_front();
			af.pop_front();
			sig.pop_front();
			position.pop_front();			
		}
		}
	}
    
    
	free(gt_arr);
	bcf_sr_destroy(sr);	
	bcf_sr_destroy(reader);	
	hts_close(out_fh);
	bcf_hdr_destroy(hdr);
	bcf_hdr_destroy(new_hdr);
	bcf_destroy1(rec);
	
	return (0);

}

