//
// Created by O'Connell, Jared on 8/3/17.
//

#include "utils.hh"

using namespace std;

bool is_genotyped(int *gt,int idx)
{
    if(idx>=0)
    {
        return( !bcf_gt_is_missing(gt[2*idx]) && !bcf_gt_is_missing(gt[2*idx+1]));
    }
    else
    {
        return(0);
    }
}


int stringSplit(std::string & s,std::vector<std::string> & split)
{
    split.clear();
    std::stringstream ss;
    ss << s;
    std::string tmp;
    while(ss>>tmp)
    {
        split.push_back(tmp);
    }
    return(split.size());
}


int stringSplit(const std::string &input, const char split, std::vector<std::string> &out)
{
    istringstream ss(input);
    out.clear();
    std::string tmp;
    int count = 0;
    while (std::getline(ss, tmp, split))
    {
        out.push_back(tmp);
        count++;
    }
    return (count);
}


pair<int,int> getGenotype(int idx,int *gt_array)
{
    pair<int,int> ret;
    if(idx<0)
    {
        ret.first = ret.second = bcf_gt_missing;
    }
    else
    {
        ret.first = bcf_gt_allele(gt_array[idx * 2]);
        ret.second = bcf_gt_allele(gt_array[idx * 2 + 1]);
    }
    return(ret);
}


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
            std::cerr << "\t -r --regions:			chromosome region" << endl;
            break;
        case 'R':
            std::cerr << "\t -R --regions-file:		restrict to regions listed in a file" << endl;
            break;
        case 'T':
            std::cerr << "\t -T --targets-file:		similar to -R but streams rather than index-jumps" << endl;
            break;
        case 't':
            std::cerr << "\t -t --targets:		        similar to -r but streams rather than index-jumps" << endl;
            break;
        case 's':
            std::cerr << "\t -s --samples:			list of samples" << endl;
            break;
        case 'S':
            std::cerr << "\t -S --samples-file:		list of samples, file" << endl;
            break;
        case 'n':
            std::cerr << "\t -n --nthreads: 		num threads" << endl;
            break;
        case 'a':
            std::cerr << "\t -a --aftag:			allele frequency tag (default AF)" << endl;
            break;
        case 'c':
            std::cerr << "\t -c --cols:			column range to read" << endl;
            break;
        case 'o':
            std::cerr << "\t -o --output:			output vcf" << endl;
            break;
        case 'O':
            std::cerr << "\t -O --outputfmt:		output vcf format" << endl;
            break;
        case 'f':
            std::cerr << "\t -f --pairfile:			file containing sample pairs to perform calculations on" << endl;
            break;
        case 'h':
            std::cerr << "\t -h --thin:			keep every h variants" << endl;
            break;
        case 'm':
            std::cerr << "\t -m --maf:			minimum MAF" << endl;
            break;
        default:
            std::cerr << "No default message exists for " << type << endl;
            break;
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
void die(const std::string &s)
{
    std::cerr << "ERROR: " << s << "\nExiting..." << endl;
    exit(1);
}
