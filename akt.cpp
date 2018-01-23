/**
 * @file   akt.cpp
 * @brief  controller for ancestry and kinship toolkit.
 *
 * parse command options and call correct functions
 */

#include "akt.hh"

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
    cerr << "\nProgram:\takt (Ancestry and Kinship Tools)" << endl;
    cerr << "Version:\t" << AKT_VERSION << endl;
    cerr << "Copyright (c) 2018, Illumina, Inc. All rights reserved. See LICENSE for further details.\n" << endl;
    cerr << "Usage:\takt <command> [options]\n" << endl;
    cerr << "\tpca                      principal component analysis" << endl;
    cerr << "\tkin                      calculate kinship coefficients" << endl;
    cerr << "\trelatives                discover pedigrees" << endl;
    cerr << "\tunrelated                generate a list of unrelated individuals" << endl;
    cerr << "\tpedphase                 Mendelian transmission phasing for duos/trios" << endl;    
    cerr << endl;
    exit(1);
}

int main(int argc, char **argv)
{


    if (argc < 2)
    {
        usage();
        return (1);
    } else if (((string) argv[1]) == "pca")
    {
        pca_main(argc, argv);
    } else if (((string) argv[1]) == "ibd")
    {
        die("ibd is deprecated");
//    ibd_main(argc, argv);
    } else if (((string) argv[1]) == "kin")
    {
        kin_main(argc, argv);
    } else if (((string) argv[1]) == "relatives")
    {
        relatives_main(argc, argv);
    } else if (((string) argv[1]) == "cluster")
    {
	die("cluster is deprecated");
//        cluster_main(argc, argv);
    } else if (((string) argv[1]) == "stats")
    {
	die("stats is deprecated");
//        stats_main(argc, argv);
    }
    else if (((string) argv[1]) == "mendel")
    {
	die("mendel command is deprecated. See the bcftools +mendelian plugin for equivalent functionality.");
        //mendel_main(argc, argv);
    } else if (((string) argv[1]) == "LDplot")
    {
	die("ldplot is deprecated");
//        ldplot_main(argc, argv);
    } else if (((string) argv[1]) == "admix")
    {
	die("admix is deprecated");
	//      admix_main(argc, argv);
	
    } else if (((string) argv[1]) == "metafreq")
    {
	die("metafreq is deprecated");
//        metafreq_main(argc, argv);
    } else if (((string) argv[1]) == "prune")
    {
	die("prune is deprecated");
//        prune_main(argc, argv);
    } else if (((string) argv[1]) == "unrelated")
    {
        unrelated_main(argc, argv);
    } else if (((string) argv[1]) == "tag")
    {
	die("tag is deprecated");
//        tag_main(argc, argv);
    }
    else if (((string) argv[1]) == "pedphase")
    {
        pedphase_main(argc, argv);
    }

        // else if(((string)argv[1]) == "grm") {
        //   grm_main(argc, argv);
        // }
    else
    {
        cerr << "Invalid command: " << argv[1] << endl;
        usage();
    }
}
