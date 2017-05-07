#include "akt.hpp"
#include "pedigree.h"
#include <iomanip>  
using namespace std;

typedef struct _args
{
    bool regions_is_file;
    bool targets_is_file;
    char *pedigree,*inputfile,*include,*regions,*targets,*outfile;  
} args;

/**
 * @name    usage
 * @brief   print out options
 *
 * List of input options
 *
 */
static void usage()
{ 
  fprintf(stderr, "\n");
  fprintf(stderr, "About:   akt mendel - profiles duo/trios\n");
  fprintf(stderr, "Usage:   ./akt mendel input.bcf -p pedigree.fam\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -p, --pedigree              pedigree information in plink .fam format\n");
  fprintf(stderr, "    -o, --out                   output a site only vcf.gz annotated with site specific error rates\n");
  fprintf(stderr, "    -i, --include               variant filters to apply eg. -i 'TYPE==snp && QUAL>=10 && DP<100000 && HWE<10' \n");  
  fprintf(stderr, "    -t, --targets [^]<region>   Set regions. Exclude regions with \"^\" prefix\n");
  fprintf(stderr, "    -T, --targets-file <file>   restrict to targets listed in a file\n");
  fprintf(stderr, "    -r, --regions <region>      restrict to comma-separated list of regions\n");
  fprintf(stderr, "    -R, --regions-file <file>   restrict to regions listed in a file\n");
  exit(1);
}


/**
 * @name    duoMendel
 * @brief   returns true if duo genotypes are inconsistent with Mendel inheritance.
 *
 * @param [in] kid child  alternate allele count (3==missing)
 * @param [in] par parent alternate allele count (3==missing)
 *
 */
bool duoMendel(int kid,int par)
{
  if( (kid==0&&par==2) || (kid==2&&par==0) )
    return(true);
  else
    return(false);  
}

/**
 * @name    trioMendel
 * @brief   returns trio if duo genotypes are inconsistent with Mendel inheritance.
 *
 * @param [in] kid child  alternate allele count (3==missing)
 * @param [in] mum alternate allele count (3==missing)
 * @param [in] dad alternate allele count (3==missing)
 *
 */
bool trioMendel(int kid,int mum,int dad)
{
  if(kid>=3||(dad>=3&&mum>=3))//cant evaluate
    return(false);
  if(mum>=3)
    return(duoMendel(kid,dad));
  if(dad>=3)
    return(duoMendel(kid,mum));

  bool ret = true;
  if(mum==0&&dad==0)
    if(kid==0) ret=false;
  if((mum==0&&dad==1)|(mum==1&&dad==0))
    if(kid!=2) ret=false;
  if((mum==0&&dad==2)|(mum==2&&dad==0))
    if(kid==1) ret=false;
  if(mum==1&&dad==1)
    ret=false;
  if((mum==1&&dad==2)|(mum==2&&dad==1))
    if(kid!=0) ret=false;      
  if(mum==2&&dad==2)
    if(kid==2) ret=false;

  return(ret);
}



inline int genotype(int *gt,int i)
{
  if(!bcf_gt_is_missing(gt[i*2])&&!bcf_gt_is_missing(gt[2*i+1])&&gt[i*2]>=0&&gt[i*2+1]>=0)
  {
    return( bcf_gt_allele(gt[i*2]) + bcf_gt_allele(gt[i*2+1]) );
  }
  else
    return(3);
}

bcf1_t *bcf_copy_sans_format(bcf1_t *dst, bcf1_t *src)
{
  //    bcf1_sync(src);

    bcf_clear(dst);
    dst->rid  = src->rid;
    dst->pos  = src->pos;
    dst->rlen = src->rlen;
    dst->qual = src->qual;
    dst->n_info = src->n_info; dst->n_allele = src->n_allele;
    dst->n_fmt = 0; dst->n_sample = 0;

    dst->shared.m = dst->shared.l = src->shared.l;
    dst->shared.s = (char*) malloc(dst->shared.l);
    memcpy(dst->shared.s,src->shared.s,dst->shared.l);
    return dst;
}


int mendel(args & a)
{
  //open a file.
  bcf_srs_t *sr =  bcf_sr_init() ; 

  if(a.targets!=NULL){
    if ( bcf_sr_set_targets(sr, a.targets, a.targets_is_file ,0)<0 ) {
      cerr << "ERROR: Failed to set targets " <<a.targets<<endl;;
      exit(1);
    }
  }
  if(a.regions!=NULL){
    if ( bcf_sr_set_regions(sr, a.regions, a.regions_is_file)<0 ) {
      cerr << "ERROR: Failed to read the regions: " <<a.regions<<endl;;
      exit(1);
    }
  }

  if(bcf_sr_add_reader(sr,a.inputfile)!=1) exit(1);

  bcf_hdr_t *hdr=sr->readers[0].header;
  filter_t *filter=NULL;
  if(a.include!=NULL)   
    filter = filter_init(hdr, a.include);

  //stuff for writing out an annotate vcf
  int32_t *duo_counts=NULL;
  int32_t *trio_counts=NULL;
  htsFile *out_fh=NULL;
  bcf_hdr_t *out_hdr=NULL;
  if(a.outfile!=NULL) {
    out_hdr = bcf_hdr_subset(hdr,0,NULL,NULL); ///creates a new subsetted header (with 0 samples) from src_header
    bcf_hdr_add_sample(out_hdr, NULL);      /// update internal structures		
    out_fh  = hts_open(a.outfile, "wz");
    bcf_hdr_append(out_hdr, "##INFO=<ID=DUO,Number=9,Type=Integer,Description=\"parent-child duo genotype counts\">");
    bcf_hdr_append(out_hdr, "##INFO=<ID=TRIO,Number=27,Type=Integer,Description=\"mother-father-child trio genotype counts\">");
    bcf_hdr_write(out_fh, out_hdr);
    duo_counts=new int32_t[9];
    trio_counts=new int32_t[27];
  }

  //note: this calls set_samples on the hdr
  sampleInfo ped(a.pedigree,hdr);

  int nsnp=0;

  bcf1_t *line,*outline=bcf_init1();;
  int  nsample = ped.N;
  int ngt = 2*nsample;
  int *gt_arr=(int *)malloc(ngt * sizeof(int));
  vector<int> gt(nsample);

  vector <vector <vector <vector<int> > > > inheritance_counts(nsample, vector< vector< vector<int> > >(3,vector< vector<int> >(3,vector<int>(3,0))));
  bool diploid_warn=false;

  cerr << "Reading input from "<<a.inputfile<<endl;
  while(bcf_sr_next_line (sr)) { 
    line =  bcf_sr_get_line(sr, 0);
    bcf_unpack(line,BCF_UN_ALL);
    if(a.include==NULL||filter_test(filter,line,NULL)) {
      if( line[0].n_allele==2) {
	if(a.outfile!=NULL) {//WRITE SOME OUTPUT.
	  for(int i=0;i<9;i++)    duo_counts[i]=0;	  
	  for(int i=0;i<27;i++)   trio_counts[i]=0;
	}
	int ret = bcf_get_genotypes(hdr, line, &gt_arr, &ngt);
	bool diploid = ret==2*nsample;
	for(int i=0;i<2*nsample;i++)  {
	  diploid = diploid && gt_arr[i]>=0;
	  if(!diploid) 
	    break;
	}
	if(diploid) {
	  for(int i=0;i<nsample;i++) 
	    gt[i] =  genotype(gt_arr,i);

	  for(int i=0;i<nsample;i++)  {
	    int dad=3,mum=3;
	    if(ped.dadidx[i]!=-1)
	      dad = ped.dadidx[i];
	    if(ped.mumidx[i]!=-1) 
	      mum = ped.mumidx[i];
	  
	    if(ped.dadidx[i]!=-1&&ped.mumidx[i]!=-1) {
	      if(gt[i]<3&&gt[dad]<3&&gt[mum]<3) {
		inheritance_counts[i][gt[dad]][gt[mum]][gt[i]]++;
		if(a.outfile)		trio_counts[gt[dad]*9+gt[mum]*3+gt[i]]++;
	      }
	    }
	    else if(ped.dadidx[i]!=-1) {
	      if(gt[i]<3&&gt[dad]<3){
		inheritance_counts[i][gt[dad]][0][gt[i]]++;
		if(a.outfile)		duo_counts[gt[dad]*3+gt[i]]++;
	      }
	    }
	    else if(ped.mumidx[i]!=-1) {
	      if(gt[i]<3&&gt[mum]<3){
		inheritance_counts[i][0][gt[mum]][gt[i]]++;
		if(a.outfile)		duo_counts[gt[mum]*3+gt[i]]++;
	      }
	    }
	  }
	  nsnp++;
	}
	else {
	  if(!diploid_warn) {
	    cerr << "\tWARNING: found non-diploid site ("<<bcf_hdr_id2name(hdr,line->rid)<<":"<<line->pos+1<<"). Not counting these sites.\n\tThis warning will not be repeated for subsequent sites.\n\tSex chromosomes will not be handled correctly, it is safest to ignore non-autosomal sites with -t ^MT,X,Y  (for example)\n"<<endl;
	    diploid_warn=true;
	  }
	}
	if(a.outfile!=NULL) {	  //WRITE SOME OUTPUT.
	  bcf_copy_sans_format(outline,line);
	  bcf_update_info_int32(out_hdr, outline, "DUO", duo_counts, 9);
	  bcf_update_info_int32(out_hdr, outline, "TRIO", trio_counts,27);
	  bcf_write1(out_fh, out_hdr, outline) ;
	}
      }
    }
  }
 
  if(a.outfile)  {
    hts_close(out_fh);
    delete[] duo_counts;
    delete[] trio_counts;
  }

  free(gt_arr);
  if(a.include!=NULL)  filter_destroy(filter);
  bcf_sr_destroy(sr);	
  if(nsnp>0)
    cerr << nsnp << " variants checked"<<endl;
  else
    die("no variants were checked!");

  cout << "PED_ID\tCHILD_ID\tDAD_ID\tMUM_ID\tDAD_GT\tMUM_GT\tCHILD_RR\tCHILD_RA\tCHILD_AA\tNERROR\tERROR_RATE\tHET_RATE"<<endl;

  string look[] = {"RR","RA","AA"};
  for(int i=0;i<nsample;i++)  {
    //    cerr << i<<"\t"<<ped.id[i]<<"\t"<<"\t"<<ped.dadidx[i]<<"\t"<<ped.mumidx[i]<<endl;
    if(ped.dadidx[i]!=-1||ped.mumidx[i]!=-1 ) {      

      if(ped.dadidx[i]!=-1&&ped.mumidx[i]!=-1) {      
	for(int j=0;j<3;j++){
	  for(int k=0;k<3;k++) {
	    int num=0,den=0;
	    cout  <<  i<<"\t"<<ped.id[i] << "\t" << ped.id[ped.dadidx[i]]<< "\t" << ped.id[ped.mumidx[i]]<<"\t"<<look[j]<<"\t"<<look[k];
	    for(int l=0;l<3;l++)	{
	      cout  <<"\t"<<inheritance_counts[i][j][k][l];
	      if(trioMendel(l,j,k))
		num+=inheritance_counts[i][j][k][l];
	      den+=inheritance_counts[i][j][k][l];
	    }
	    cout <<"\t"<<num<<"\t"<<std::setw(10)<<(float)num/(float)den<<"\t"<<std::setw(10)<<(float)inheritance_counts[i][j][k][1]/(float)den<<endl;
	  }
	}
      }
      else if(ped.dadidx[i]!=-1) {
	for(int j=0;j<3;j++){
	  int num=0,den=0;
	  cout  <<  i<<"\t"<<ped.id[i] << "\t" << ped.id[ped.dadidx[i]]<< "\t.\t"<<look[j]<<"\t.";
	  for(int l=0;l<3;l++)	{
	    cout  <<"\t"<<inheritance_counts[i][j][0][l];
	    if(duoMendel(l,j)) 
	      num+=inheritance_counts[i][j][0][l];
	    den+=inheritance_counts[i][j][0][l];
	  }
	  cout <<"\t"<<num<<"\t"<<std::setw(10)<<(float)num/(float)den<<"\t"<<std::setw(10)<<(float)inheritance_counts[i][j][0][1]/(float)den<<endl;
	}
      }
      else if(ped.mumidx[i]!=-1) {
	for(int k=0;k<3;k++){
	  int num=0,den=0;
	  cout  <<  i<<"\t"<<ped.id[i] << "\t.\t" << ped.id[ped.mumidx[i]]<< "\t.\t"<<look[k];
	  for(int l=0;l<3;l++)	{
	    cout  <<"\t"<<inheritance_counts[i][0][k][l];
	    if(duoMendel(l,k)) 
	      num+=inheritance_counts[i][0][k][l];
	    den+=inheritance_counts[i][0][k][l];
	  }
	  cout <<"\t"<<num<<"\t"<<std::setw(10)<<(float)num/(float)den<<"\t"<<std::setw(10)<<(float)inheritance_counts[i][0][k][1]/(float)den<<endl;
	}
      }    

    }
  }  

  return(0);
}


int mendel_main(int argc,char **argv) {
  int c;
  args a;
  if(argc<3) usage();
  static struct option loptions[] =    {
    {"out",1,0,'o'},	
    {"pedigree",1,0,'p'},
    {"include",1,0,'i'},
    {"targets",required_argument,NULL,'t'},
    {"targets-file",required_argument,NULL,'T'},
    {"regions-file",required_argument,NULL,'R'},
    {"regions",required_argument,NULL,'r'},
    {0,0,0,0}
  };
  a.regions_is_file=false;
  a.targets_is_file=false;
  a.targets=a.pedigree=a.inputfile=a.include=a.regions=a.outfile=NULL;

  while ((c = getopt_long(argc, argv, "o:p:i:t:T:r:R:",loptions,NULL)) >= 0) {  
    switch (c)      {
    case 'o': a.outfile = optarg; break;
    case 'p': a.pedigree = optarg; break;
    case 'i': a.include = optarg; break;
    case 't': a.targets = optarg; break;
    case 'T': a.targets = optarg; break;
    case 'r': a.regions = optarg; break;
    case 'R': a.regions = optarg; a.regions_is_file=true; break;
    }
  }
  optind++;
  a.inputfile=argv[optind];
  if(a.inputfile==NULL) die("no input provided");
  if(a.pedigree==NULL) die("the -p argument is required");
  mendel(a);
  return(0);
}
