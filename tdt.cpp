#include "akt.hpp"
#include "pedigree.h"
#include <iomanip>  
using namespace std;

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
  fprintf(stderr, "About:   akt tdt - simple transmission counting for affected/unaffected\n");
  fprintf(stderr, "Usage:   ./akt tdt input.bcf -p pedigree.fam\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "    -p, --pedigree              pedigree information in plink .fam format\n");
  exit(1);
}

string GENOTYPE[] = {"RR","RA","AA","."};
    
inline int genotype(int *gt,int i)
{
  if(!bcf_gt_is_missing(gt[i*2])&&!bcf_gt_is_missing(gt[2*i+1])&&gt[i*2]>=0&&gt[i*2+1]>=0)
  {
    return( bcf_gt_allele(gt[i*2]) + bcf_gt_allele(gt[i*2+1]) );
  }
  else
    return(3);
}

int tdt(char *input_file,char *pedigree_file)
{
  //open a file.
  bcf_srs_t *sr =  bcf_sr_init() ; 

  if(bcf_sr_add_reader(sr,input_file)!=1) exit(1);

  bcf_hdr_t *hdr=sr->readers[0].header;

  //stuff for writing out an annotate vcf
  int32_t *trio_counts=NULL;

  //note: this calls set_samples on the hdr
  sampleInfo ped(pedigree_file,hdr);

  int nsnp=0;

  bcf1_t *line,*outline=bcf_init1();;
  int  nsample = ped.N;
  int ngt = 2*nsample;
  int *gt_arr=(int *)malloc(ngt * sizeof(int));
  vector<int> gt(nsample);

  vector <vector <vector <vector<int> > > > inheritance_counts(nsample, vector< vector< vector<int> > >(3,vector< vector<int> >(3,vector<int>(3,0))));
  bool diploid_warn=false;

  cerr << "Reading input from "<<input_file<<endl;
  cout << "CHROM\tPOS\tREF\tALT\tDAD\tMUM\tDAD_GT\tMUM_GT\tRR_0\tRA_0\tAA_0\tMISS_0\tRR_1\tRA_1\tAA_1\tMISS_1\tRR_2\tRA_2\tAA_2\tMISS_2"<<endl;
  while(bcf_sr_next_line (sr)) 
  { 
      line =  bcf_sr_get_line(sr, 0);
      bcf_unpack(line,BCF_UN_ALL);
      if( line[0].n_allele==2) 
      {
	  int ret = bcf_get_genotypes(hdr, line, &gt_arr, &ngt);
	  bool diploid = ret==2*nsample;
	  for(int i=0;i<2*nsample;i++)  {
	      diploid = diploid && gt_arr[i]>=0;
	      if(!diploid) 
		  break;
	  }
	  if(diploid) 
	  {
	      for(int i=0;i<nsample;i++) 
	      {
		  gt[i] =  genotype(gt_arr,i);
	      }	      
	      //loop over parents.
	      //check for hets on either parent
	      //tabulate genotypes in affected/unaffected children	      
	      //print:
	      //dad mum affected_RR affected_RA affected_AA  unaffected_RR unaffected_RA unaffected_AA
	      //RR  RA  0           0           0                        0             1             0
	      for( map<pair<int,int>, vector<int>  >::iterator it1=ped.parent_map.begin();it1!=ped.parent_map.end();it1++)
	      {
		  vector<int> freq(12,0);
		  int dad_gt=gt[it1->first.first];
		  int mum_gt=gt[it1->first.second];
		  if(dad_gt<3 && mum_gt<3 && (dad_gt>0||mum_gt>0))
		  {		      
		      cout << bcf_hdr_id2name(hdr,line->rid)<<"\t"<<line->pos+1<<"\t"<<line->d.allele[0]<<"\t"<<line->d.allele[1]<<"\t"<<*ped.getID(it1->first.first)<<"\t"<<*ped.getID(it1->first.second)<<"\t"<< GENOTYPE[dad_gt]<<"\t"<<GENOTYPE[mum_gt];
		      for(size_t i=0;i<it1->second.size();i++)
		      {
			  int child_gt=gt[it1->second[i]];
			  int child_status=ped.getStatus(it1->second[i]);
			  if(child_status>=0 && child_status<3)
			  {
			      freq[child_status*4 + child_gt]++;
			  } 
		      }
		      for(size_t i=0;i<freq.size();i++)
		      {
			  cout << "\t"<<freq[i];
		      }
		      cout << endl;
		  }
	      }	      
	  }
      }
  }
  return(0);
}


int tdt_main(int argc,char **argv) {
  int c;
  
  if(argc<3) usage();
  static struct option loptions[] =    {
    {"pedigree",1,0,'p'},
    {0,0,0,0}
  };
  char *inputfile=NULL,*pedigree=NULL;

  while ((c = getopt_long(argc, argv, "p:",loptions,NULL)) >= 0) 
  { 
      switch (c)      
      {
      case 'p': pedigree = optarg; break;
      }
  }
  optind++;
  inputfile=argv[optind];
  if(inputfile==NULL) die("no input provided");
  if(pedigree==NULL) die("the -p argument is required");
  tdt(inputfile,pedigree);
  return(0);
}
