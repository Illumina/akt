#include "akt.hpp"

using namespace std;

void umessage(const char type){
	switch(type){
		case 'r': cerr << "\t -r --regions:			chromosome region" << endl; break;
		case 'R': cerr << "\t -R --regions-file:		list of regions, file" << endl; break;
		case 'T': cerr << "\t -T --targets-file:		intersecting VCF" << endl; break;
		case 's': cerr << "\t -s --samples:			list of samples" << endl; break; 
		case 'S': cerr << "\t -S --samples-file:		list of samples, file" << endl; break;
		case 'n': cerr << "\t -n --nthreads: 		num threads" << endl; break;
		case 'a': cerr << "\t -a --aftag:			allele frequency tag" << endl; break;	
		case 'c': cerr << "\t -c --cols:			column range to read" << endl; break;	
		case 'o': cerr << "\t -o --output:			output vcf" << endl; break;	
		case 'O': cerr << "\t -O --outputfmt:		output vcf format" << endl; break;	
		case 'f': cerr << "\t -f --pairfile:			file containing sample pairs" << endl; break;
		case 'h': cerr << "\t -h --thin:			keep every t variants" << endl; break;
		case 'm': cerr << "\t -m --maf:			minimum MAF" << endl; break;
		default: cerr << "No default message exists for " << type << endl; break; 
	}
}

int main(int argc, char **argv) {


  if(argc < 2) {
    cerr << "\nProgram:\takt (Ancestry and Kinship Tools)" << endl;
    cerr << "Version:\t" << VERSION <<endl;
    cerr << "Contact:\trarthur@illumina.com\n" << endl;
    cerr << "Copyright (c) 2016, Illumina, Inc. All rights reserved. See LICENSE for further details.\n"<<endl;
    cerr << "Usage:\takt <command> [options]\n" << endl;
    cerr << "\tpca                      principal component analysis" << endl;
    cerr << "\tkin                      detect average IBD sharing" << endl;
    cerr << "\trelatives                discover pedigrees" << endl;
    cerr << "\tibd                      detect segments shared IBD" << endl;
    cerr << "\tmendel                   profile Mendelian inhertiance and inconsistencies in known pedigrees" << endl;
    cerr << "\tcluster                  perform cluster analyses" << endl;
    cerr << "\tLDplot                   output correlation matrix" << endl;
    cerr << "\tstats                    calculate AF and LD metrics" << endl<<endl;
    return(1);
  }
  else if(((string)argv[1]) == "pca") {
    pca_main(argc, argv);
  } else if(((string)argv[1]) == "ibd") {
    ibd_main(argc, argv);
  } else if(((string)argv[1]) == "kin") {
    kin_main(argc, argv);
  } else if(((string)argv[1]) == "relatives") {
    relatives_main(argc, argv);
  } else if(((string)argv[1]) == "cluster") {
    cluster_main(argc, argv);
  } else if(((string)argv[1]) == "stats") {
    stats_main(argc, argv);
  } else if(((string)argv[1]) == "mendel") {
    mendel_main(argc, argv); 
  } else if(((string)argv[1]) == "LDplot") {
    ldplot_main(argc, argv); 
  } else if(((string)argv[1]) == "admix") {
    admix_main(argc, argv); 
  } 
  else {
    cerr << "Invalid command: " << argv[1] << endl;
  }
}
