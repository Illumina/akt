#define __STDC_LIMIT_MACROS

#include "akt.hpp"
#include "logs.hpp"
#include "kin.hpp"

using namespace std;

float calc_snp_logodds(vector< vector<int> > &G, vector< float > &af, vector< float > &ld, int ip, int j1, int j2, vector< vector<float> > &GQ){
	
		//log allele freq
		float l1af = log(1 - af[ip]);
		float laf = log(af[ip]);
		
		//IBD
		float gl1[3]; for(int g=0; g<3; ++g){ gl1[g] = ( g == G[ip][j1] ) ? 1.0-GQ[ip][j1] : GQ[ip][j1]/2.0; } 
		float gl2[3]; for(int g=0; g<3; ++g){ gl2[g] = ( g == G[ip][j2] ) ? 1.0-GQ[ip][j2] : GQ[ip][j2]/2.0; } 

		//no IBD -> HWE
		float den = log( 
				gl1[0] * (1 - af[ip]) * (1 - af[ip])
			+	gl1[1] * (1 - af[ip]) * af[ip] * 2
			+	gl1[2] * af[ip] * af[ip]
			);

		float num = log(
			gl2[0] * ( gl1[0]*(1 - af[ip]) + gl1[1]*af[ip] )
		+	gl2[1] * ( gl1[0]*(0.5*(1 - af[ip])) + gl1[1]*0.5 + gl1[2]*0.5*af[ip] )
		+	gl2[2] * ( gl1[1]*(1 - af[ip]) + gl1[2]*af[ip] )
		);
		
		return (num - den)/ld[ip];
}

static void usage(){
  cerr << "Find IBD regions from a VCF" << endl;	
  cerr << "Usage:" << endl;
  cerr << "./akt ibd in.bcf -p sites.bcf" << endl;
  cerr << "Expects input_filename.vcf to contain hard genotypes" << endl;
  umessage('T');
  umessage('r');
  umessage('R');
  umessage('s');
  umessage('S');
  umessage('h');
  umessage('a');
  umessage('f');
  umessage('n');
  umessage('m');
  cerr << "\t -e --error: 			If no GQ then this is error probability(1e-3)" << endl;	 
  cerr << "\t -M --thresh: 			likelihood output threshold (default 50)" << endl;
  cerr << "\t -L --length: 			length output threshold (default 100000)" << endl;
  cerr << "\t -w --wsize: 			window size (default 20)" << endl;
  cerr << "\t -x --maxerr: 			stop when encountering a window with this many ( 0/0 , 1/1 ) pairs (default 2)" << endl;
  cerr << "\t -l --lsize: 			try to join long regions closer than this (default 100000)" << endl;
  exit(1);
}

class region{
public:
	int begin;
	int end;
	int matches;
	float score;	

	region() : begin(0), end(0), matches(0), score(logz){}
	region(const region& other) : begin(other.begin), end(other.end), matches(other.matches), score(other.score){}
		
	friend void swap(region& first, region& second)
    {
        swap(first.begin, second.begin);
        swap(first.end, second.end);
        swap(first.matches, second.matches);
        swap(first.score, second.score);
    }	
    region& operator=(region other)
	{
		swap(*this, other); 
		return *this;
	}
	
	region& operator+=(const region &rhs) {
		end = rhs.end;
		matches += rhs.matches;
		score +=rhs.score;
		return *this;
    }
    const region operator+(const region &other) const {
		return region(*this) += other;
    }
    
    friend bool operator<(const region& l, const region& r){ return l.score < r.score; }
    friend bool operator>(const region& l, const region& r){ return r<l; }
    friend bool operator<=(const region& l, const region& r){ return !(l > r); }
    friend bool operator>=(const region& l, const region& r){ return !(l < r); }
};

bool chrom_compare(const region& l, const region& r){ return l.begin < r.begin; }

int ibd_main(int argc, char* argv[])
{
		
	int c;
  
	if(argc<3) usage();
	static struct option loptions[] =    {
	{"regions",1,0,'r'},	
	{"nthreads",1,0,'n'},	
	{"regions-file",1,0,'R'},
	{"targets-file",1,0,'T'},
	{"thin",1,0,'h'},
	{"wsize",1,0,'w'},
	{"lsize",1,0,'l'},
	{"aftag",1,0,'a'},
	{"pairfile",1,0,'f'},
	{"error",1,0,'e'},
	{"min_freq",1,0,'m'},
	{"thresh",1,0,'M'},
	{"length",1,0,'L'},
	{"maxerr",1,0,'x'},
	{"samples",1,0,'s'},
    {"samples-file",1,0,'S'},
	{0,0,0,0}
	};

	string pfilename = "";
	bool norm = false;
	int thin = 1;
	int nthreads = -1;
	int wsize = 20;
	int lsize = 100000;
	string af_tag = "";
	float thresh = 50;
	float min_len = 100000;
	string pairfile="";
	float error = 1e-3;

	float min_freq = 0;
    sample_args sargs;

    bool get_regions = false; string regions = "";
    bool regions_is_file = false;
    bool used_r = false;
    bool used_R = false;

	int maxerr = 2;
	
	while ((c = getopt_long(argc, argv, "r:R:h:n:w:l:a:T:L:f:e:M:m:x:s:S:",loptions,NULL)) >= 0) {  
		switch (c) {
			case 'r': get_regions = true; regions = (optarg); used_r = true; break;
			case 'R': get_regions = true; regions = (optarg); used_R = true; regions_is_file = true; break;
			case 'T': pfilename = (optarg); break;
			case 'h': thin = atoi(optarg); break;
			case 'n': nthreads = atoi(optarg); break;
			case 'M': thresh = atof(optarg); break;
			case 'L': min_len = atoi(optarg); break;
			case 'm': min_freq = atof(optarg); break;
			case 'e': error = atof(optarg); break;
			case 'w': wsize = atoi(optarg); break;
			case 'l': lsize = atoi(optarg); break;
			case 'a': af_tag = string(optarg) + "_"; break;
			case 'f': pairfile = (optarg); break;
			case 'x': maxerr = atoi(optarg); break;
			case 's': sargs.sample_names = (optarg); sargs.subsample = true; break;
			case 'S': sargs.sample_names = (optarg); sargs.subsample = true; sargs.sample_is_file = 1; break;
			case '?': usage();
			default: cerr << "Unknown argument:"+(string)optarg+"\n" << endl; exit(1);
		}
	}
	if(!get_regions){ cerr << "-r argument is mandatory" << endl; exit(1); }
	if(optind>=argc-1) { cerr<< "No input .bcf/.vcf provided!"<<endl;exit(1);}
    if(pfilename.empty()){ cerr << "No sites file provided" << endl;exit(1); }
    
	omp_set_num_threads(nthreads); 
	#pragma omp parallel
	{
		if(omp_get_thread_num() == 0){
			if( omp_get_num_threads() != 1){
				cerr << "Using " << omp_get_num_threads() << " threads" << endl;
			}
		}
	}
	
	optind++;
	string filename = argv[optind];
					
	bcf_srs_t *sr =  bcf_sr_init() ; ///htslib synced reader.
	sr->collapse = COLLAPSE_NONE;
	sr->require_index = 1;
	if(get_regions){
	if ( bcf_sr_set_regions(sr, regions.c_str(), false)<0 ){ 
	cerr << "Failed to read the regions: " <<  regions << endl; 
	return 0;
	}}
	if(!(bcf_sr_add_reader (sr, filename.c_str() ))){ 
	cerr << "Problem opening " << filename << endl; 
	bcf_sr_destroy(sr);	
	return 0;
	}
	if(!(bcf_sr_add_reader (sr, pfilename.c_str() ))){ 
	cerr << "Problem opening " << pfilename << endl; 
	bcf_sr_destroy(sr);	
	return 0;
	}
	//sub_sample_header(&sr->readers[0].header, sargs);
    if(sargs.subsample){ bcf_hdr_set_samples(sr->readers[0].header, sargs.sample_names, sargs.sample_is_file); }
		
	int N = bcf_hdr_nsamples(sr->readers[0].header);	///number of samples 
	cerr << N << " samples" << endl;
	vector<string> names;
	map<string,int> name_to_id;
	set<string> unique_names;

	for(int i=0; i<N; ++i){ 
	string tmp = sr->readers[0].header->samples[i]; 
	names.push_back(tmp);
	name_to_id[tmp] = i;
	}

	vector< pair<string, string> > relpairs;
	if( pairfile != "" ){
		ifstream in(pairfile.c_str());
		read_pairs(in, relpairs, name_to_id);
		in.close();
	} else {
		make_pair_list(relpairs, names);
	}

	//Dont't waste time reading all samples
	for(int r=0; r<relpairs.size(); ++r){ //all sample pairs
		unique_names.insert( relpairs[r].first  );
		unique_names.insert( relpairs[r].second  );
    }	  
	set<string>::iterator ufin = --unique_names.end();
	set<string>::iterator uit  = unique_names.begin();
	string sample_names="";
	for (; uit!=ufin; ++uit){
		sample_names += (*uit) + ",";
	} sample_names += (*uit);

    int sret = bcf_hdr_set_samples(sr->readers[0].header, sample_names.c_str(), 0);
	N = bcf_hdr_nsamples(sr->readers[0].header);	///number of samples 
	
	name_to_id.clear();
	for(int i=0; i<N; ++i){ 
		string tmp = sr->readers[0].header->samples[i]; 
		name_to_id[tmp] = i;
	}
	
	vector< vector<int> > G; 

	int count=0,sites=0,markers=0;
	bcf1_t *line, *line2;///bcf/vcf line structure.

	vector< float > af;
	vector< float > ld;
	int *gt_arr=(int *)malloc(N*2*sizeof(int)),ngt=N*2,ngt_arr=N*2;
	int *gq_arr=(int *)malloc(N*sizeof(int)),ngq=N,ngq_arr=N;
	float *af_ptr=(float *)malloc(1*sizeof(float)); int nval = 1;
	float *ld_ptr=(float *)malloc(1*sizeof(float)); 

	vector<int> position;

	vector< vector<float> > GQ;
    
	int num_sites=0,num_study=0;
	
	
	while(bcf_sr_next_line (sr)) { 	
	if( bcf_sr_has_line(sr,1) ){	//present in sites file.
			
	line2 =  bcf_sr_get_line(sr, 1);
	if(line2->n_allele == 2){	//bi-allelic
	int nmiss=0;
	int sum = 0;
	num_sites++;
	if(bcf_sr_has_line(sr,0)) { //present in the study file
	line =  bcf_sr_get_line(sr, 0);
	ngt = bcf_get_genotypes(sr->readers[0].header, line, &gt_arr, &ngt_arr);    
	for(int i=0;i<2*N;i++){ 
		if(gt_arr[i]==bcf_gt_missing || (bcf_gt_allele(gt_arr[i])<0) || (bcf_gt_allele(gt_arr[i])>2)  ){
			gt_arr[i]=bcf_gt_phased(0);//inserting homref instead of missing.
		}
		sum += bcf_gt_allele(gt_arr[i]);
		++nmiss; 
	}
	++num_study;

	//only take one of every 'thin' of these sites
	if( (count++)%thin==0 ){ 

		float p;
		int ret = bcf_get_info_float(sr->readers[1].header, line2, (af_tag + "AF").c_str(), &af_ptr, &nval);
		if( nval != 1 ){ cerr << (af_tag + "AF") << " read error at " << line2->rid << ":" << line->pos+1 << endl; exit(1); }
		p = af_ptr[0];
		
		if( (p < 0.5) ? ( p > min_freq ) : (1-p > min_freq) ){ 	
						
			vector<int> tmp(N); 
			for(int i=0; i<N; ++i){ tmp[i] = bcf_gt_allele(gt_arr[2*i]) + bcf_gt_allele(gt_arr[2*i+1]); } 
			G.push_back(tmp);							
			af.push_back( p );
			position.push_back( line->pos+1 );
			
			ngq = bcf_get_format_float(sr->readers[0].header, line, "GQ", &gq_arr, &ngq_arr);
			vector<float> tmpq(N);
			if( ngq < 0 ){
				for(int i=0; i<N; ++i){ tmpq[i] = error; } 
			} else {
				for(int i=0; i<N; ++i){ tmpq[i] = exp(-0.1*gq_arr[i]*logten); } 
			}
			GQ.push_back( tmpq );
			
			float tl;
			int ret = bcf_get_info_float(sr->readers[1].header, line2, (af_tag + "LD").c_str(), &ld_ptr, &nval);
			if( ret<=-1 ){ 
				tl = 1; 
			} else {
				tl = ld_ptr[0];
			}
			ld.push_back(tl);
			
			++markers;
		}

	}//thin
	
	}//in both	
	}//biallelic
	}//in intersection
		if( bcf_sr_has_line(sr,1) ){ ++sites; }
	}//reader


	free(gt_arr);
	free(gq_arr);
	free(af_ptr);

	cerr << "Kept " << markers << " markers out of " << sites << " in panel." << endl;
	cerr << num_study << "/"<<num_sites<<" of study markers were in the sites file"<<endl;

	bcf_sr_destroy(sr);	
  
  vector< list< region > > all_extended_match_region(relpairs.size());
  //vector< vector<float> > all_mscore(relpairs.size());
  cerr << "computing IBD sharing for " << relpairs.size() << " sample pairs." << endl;
  
  #pragma omp parallel for	
  for(int r=0; r<relpairs.size(); ++r){ //all sample pairs
  
	  int j1 = name_to_id[ relpairs[r].first  ]; //sample ids
	  int j2 = name_to_id[ relpairs[r].second ];	  
	 
	  list<region> matches;
	  ///count number of IBS>0 matches between j1 and j2 in each window
	  for(int i=0; i<G.size(); i+=wsize){	
		  region rg;
		  rg.begin = i;
		  rg.end = i+wsize;
		  for(int ip=i; ip<std::min(i+wsize,(int)G.size()); ++ip){ 
			  rg.matches += (int)( abs(G[ip][j2] - G[ip][j1]) < 2 ); // if 0/0 && 1/1 == 2 else < 2
		  } 
		  if(i+wsize < (int)G.size() 
		  && position[i+wsize]-position[i] < lsize){ matches.push_back(rg); } //chop off last window and long windows
	  }
	  	 
	  ///contiguous regions with IBS>0 
	  for (list<region>::iterator it=matches.begin(); it!=matches.end(); ++it){

			if( it->matches == wsize ){	//found a perfect window
				region rg( *it );
				
				//window score
				float score = 0;
				for(int wp=0; wp<wsize; ++wp){ score += calc_snp_logodds(G, af, ld, (it->begin) + wp, j1, j2, GQ); }
				rg.score = score;

				//done this window, drop it from list
				//123(it = 4)567...
				//123(it = 5)67...
				it = matches.erase( it );

				while( it !=matches.end() && it->matches == wsize ){ //while we can add perfect windows
					
					//window score
					score = 0; 
					for(int wp=0; wp<wsize; ++wp){ score += calc_snp_logodds(G, af, ld, (it->begin) + wp, j1, j2, GQ); }
					it->score = score;
					
					rg += (*it);	//accumulate
					it = matches.erase( it ); //done this window, drop it from list
					
				}
				//save perfect region
				//123(it = 4)56
				//123(it = 4)(rg)56
				matches.insert(it, rg);
				--it; //undoes increment (necessary??)

			} //end perfect window
	  }
  
	  float initm = 0;
	  list<region> best_regions;
	  
	  //utility iterator, can't use prev and next
	  list<region>::iterator tmpit;	  	  
	  do{ 
		  //no matches!
		  if( matches.empty() ) break;

		  //most likely IBD region
		  list<region>::iterator mit = max_element(matches.begin(), matches.end());
		  initm = mit->score;
		  
		  float best_score = mit->score;
		  region best_region = (*mit);
		  region init_region = best_region;
		  
		  //to the right of best region
		  list<region>::iterator rit = matches.erase(mit);
		  list<region>::iterator init_rit = rit;
		
		  //to the left of best region
		  tmpit = rit; --tmpit; 
		  list<region>::iterator lit = (rit != matches.begin()) ? tmpit : rit;
		  list<region>::iterator init_lit = lit;


		  region cur_region = best_region;
		  
		  //extend left from max	  
		  for(; lit != matches.begin(); --lit){ 
			  if( lit->matches <= wsize - maxerr ){ break; } //too many mismatches
			  else {
				  if(lit->score == logz){ //calculate window score if necessary
					  float score = 0; 
					  for(int wp=0; wp<wsize; ++wp){ score += calc_snp_logodds(G, af, ld, (lit->begin) + wp, j1, j2, GQ); }
					  lit->score = score;
				  }
				  //extend current region left
				  region tmp_region = (*lit);
				  tmp_region += cur_region;
				  cur_region = tmp_region;

				  if( cur_region.score > best_region.score){ best_region = cur_region; } //save region with highest score encountered
			  }
		  } 
		  
		  //extend right from max
		  cur_region = best_region;	  
		  for(; rit != matches.end(); ++rit){ 
			  if( rit->matches <= wsize - maxerr ){ break; } //too many mismatches
			  else {
				  if(rit->score == logz){ //calculate window score if necessary
					  float score = 0; 
					  for(int wp=0; wp<wsize; ++wp){ score += calc_snp_logodds(G, af, ld, (rit->begin) + wp, j1, j2, GQ); }
					  rit->score = score;
				  }
				  //extend current region right
				  cur_region += (*rit);
				  if( cur_region.score > best_region.score){best_region = cur_region; } //save region with highest score encountered
			  }
		  
		  } 

		 if( best_region.end != init_region.end ){ //addeed to RHS
			list<region>::iterator it = init_rit; 
			//delete used windows
			//123(4)567
			//123(5)67...
			while(it != matches.end() && it->begin < best_region.end){it = matches.erase(it);} 
		 }
		 
		 if( best_region.begin != init_region.begin ){ //addeed to LHS
			list<region>::iterator it = init_lit; 
			//delete used windows
			//12(3)4567
			//12(4)567 
			//1(2)4567...
			while(it != matches.end() && it->end > best_region.begin){ 
				if(it == matches.begin() ){ break; }
				else{ it = matches.erase(it); --it; }
			}
		 }
		 //save best region
		 best_regions.push_back( best_region );

	} while (initm > 0); //starting region had positive likelihood

	//attempt to join neighbouring regions 
	if( !best_regions.empty() ){
		
	//sort regions by position
	best_regions.sort(chrom_compare);
	
	bool merged;
	do{
		merged = false;

		tmpit = best_regions.end(); --tmpit;
		for( list<region>::iterator mit = best_regions.begin(); mit != tmpit; ++mit){
		
			list<region>::iterator rit = mit; ++rit; //region to the right
			
			if( position[ rit->begin ] - position[ mit->end ] < lsize){ //if regions are close
				
				//calculate score from gap
				float gap_score = 0;
				for(int wp=mit->end; wp<rit->begin; ++wp){ gap_score += calc_snp_logodds(G, af, ld, wp, j1, j2, GQ); }
				if(rit->score + gap_score > 0 ){ //closing gap increases score of combined segment
					*mit += *rit; //add regions
					mit->score += gap_score; //add gap
					rit = best_regions.erase(rit); //remove merged region
					tmpit = best_regions.end(); --tmpit;
					merged = true;
				}
				
			}
		}
		
	} while(merged); //stuff is still getting merged
	
	all_extended_match_region[r] = best_regions; //save best regions
	
	} //match list empty
	
  } //end pair list
  
  for(int r=0; r<all_extended_match_region.size(); ++r){
		for(list<region>::iterator it=all_extended_match_region[r].begin(); 
		it != all_extended_match_region[r].end(); ++it){
			if(it->score >= thresh && position[ it->end ] - position[ it->begin ] >= min_len){
				
			cout << relpairs[r].first << "\t" << relpairs[r].second << "\t";
			cout << position[ it->begin ] << "\t" << position[ it->end ]
			<< "\t" << it->score << "\t" << it->matches << "\t" << position[ it->end ] - position[ it->begin ] << endl;
			
			}
		} 
  }
  
  
  return 0;
}
