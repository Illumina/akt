/**
 * @file   akt.cpp
 * @brief  controller for ancestry and kinship toolkit.
 *
 * parse command options and call correct functions
 */

#include "akt.hpp"

using namespace std;

/**
 * @name    umessage
 * @brief   common command options
 *
 * Enforcing consistency on input option description
 *
 */
void umessage(const char type)
{
    switch (type)
    {
        case 'r':
            cerr << "\t -r --regions:			chromosome region" << endl;
            break;
        case 'R':
            cerr << "\t -R --regions-file:		restrict to regions listed in a file" << endl;
            break;
        case 'T':
            cerr << "\t -T --targets-file:		similar to -R but streams rather than index-jumps" << endl;
            break;
        case 't':
            cerr << "\t -t --targets:		        similar to -r but streams rather than index-jumps" << endl;
            break;
        case 's':
            cerr << "\t -s --samples:			list of samples" << endl;
            break;
        case 'S':
            cerr << "\t -S --samples-file:		list of samples, file" << endl;
            break;
        case 'n':
            cerr << "\t -n --nthreads: 		num threads" << endl;
            break;
        case 'a':
            cerr << "\t -a --aftag:			allele frequency tag (default AF)" << endl;
            break;
        case 'c':
            cerr << "\t -c --cols:			column range to read" << endl;
            break;
        case 'o':
            cerr << "\t -o --output:			output vcf" << endl;
            break;
        case 'O':
            cerr << "\t -O --outputfmt:		output vcf format" << endl;
            break;
        case 'f':
            cerr << "\t -f --pairfile:			file containing sample pairs to perform calculations on" << endl;
            break;
        case 'h':
            cerr << "\t -h --thin:			keep every h variants" << endl;
            break;
        case 'm':
            cerr << "\t -m --maf:			minimum MAF" << endl;
            break;
        default:
            cerr << "No default message exists for " << type << endl;
            break;
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
    cerr << "\nProgram:\takt (Ancestry and Kinship Tools)" << endl;
    cerr << "Version:\t" << VERSION << endl;
    cerr << "Copyright (c) 2016, Illumina, Inc. All rights reserved. See LICENSE for further details.\n" << endl;
    cerr << "Usage:\takt <command> [options]\n" << endl;
    cerr << "\tpca                      principal component analysis" << endl;
    cerr << "\tkin                      detect average IBD sharing" << endl;
    cerr << "\trelatives                discover pedigrees" << endl;
    cerr << "\tunrelated                generate a list of unrelated individuals" << endl;
//    cerr << "\tibd                      detect segments shared IBD" << endl;
    cerr << "\tmendel                   profile Mendelian inhertiance and inconsistencies in known pedigrees" << endl;
    cerr << "\tcluster                  perform cluster analyses" << endl;
    cerr << "\tLDplot                   output correlation matrix" << endl;
    cerr << "\tstats                    calculate AF and LD metrics" << endl;
    cerr << "\tprune                    perorms LD pruning of variants" << endl;
    cerr << "\ttag                      selects a set of K tagging variants" << endl;
    cerr << "\tmetafreq                 examine two files for AF differences" << endl;
    cerr << "\ttdt                      basic parent-child transmission counting for use in TDT" << endl;
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
        cluster_main(argc, argv);
    } else if (((string) argv[1]) == "stats")
    {
        stats_main(argc, argv);
    } else if (((string) argv[1]) == "mendel")
    {
        mendel_main(argc, argv);
    } else if (((string) argv[1]) == "LDplot")
    {
        ldplot_main(argc, argv);
    } else if (((string) argv[1]) == "admix")
    {
        admix_main(argc, argv);
    } else if (((string) argv[1]) == "metafreq")
    {
        metafreq_main(argc, argv);
    } else if (((string) argv[1]) == "prune")
    {
        prune_main(argc, argv);
    } else if (((string) argv[1]) == "unrelated")
    {
        unrelated_main(argc, argv);
    } else if (((string) argv[1]) == "tag")
    {
        tag_main(argc, argv);
    } else if (((string) argv[1]) == "tdt")
    {
        tdt_main(argc, argv);
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

/**
 * @name    die
 * @brief   exit with error message
 *
 * Exit "gracefully"
 *
 * @param [in] s Error string.
 */
void die(const string &s)
{
    cerr << "ERROR: " << s << "\nExiting..." << endl;
    exit(1);
}
