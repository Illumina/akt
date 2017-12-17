//
// Created by O'Connell, Jared on 8/3/17.
//

#include "utils.h"

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


int stringSplit(string & s,vector<string> & split)
{
    split.clear();
    stringstream ss;
    ss << s;
    string tmp;
    while(ss>>tmp)
    {
        split.push_back(tmp);
    }
    return(split.size());
}


int stringSplit(const string &input, const char split, vector<string> &out)
{
    istringstream ss(input);
    out.clear();
    string tmp;
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

Genotype::Genotype(int g0,int g1)
{
    setGenotype(g0,g1);
}

void Genotype::setGenotype(int g0,int g1)
{
    _is_haploid=false;
    if(bcf_gt_is_missing(g0))
    {
	_g0 = -1;
    }
    else
    {
	_g0=bcf_gt_allele(g0);	    
    }
    
    if(bcf_gt_is_missing(g1))
    {
	_g1 = -1;
    }
    else if(g1==bcf_int32_vector_end)
    {
	_g1=bcf_int32_vector_end;
	_is_haploid=true;
    }
    else
    {
	_g1=bcf_gt_allele(g1);	    
    }

}

Genotype::Genotype(int idx,int *gt_array,int *ps_array)
{
    if(idx<0)
    {
	_g0=-1;
	_g1=-1;
	_is_haploid=false;
	_is_phased=false;
    }
    else
    {
	_is_phased = ps_array && ps_array[idx]!=bcf_int32_missing;
	_ps = ps_array ? ps_array[idx] : bcf_int32_missing;
	setGenotype(gt_array[2*idx],gt_array[2*idx+1]);
    }
}

int Genotype::ps()
{
    return _ps;
}

bool Genotype::is_phased()
{
    return( _is_phased || !(isHet() || isHaploid()) );
}

bool Genotype::isHet()
{
    return(_g0!=-1 && _g1!=-1 && !_is_haploid && _g0!=_g1);
}

bool Genotype::isMissing()
{
    return(_g0==-1||_g1==-1);
}

bool Genotype::isHaploid()
{
    return(_is_haploid);
}

int Genotype::first()
{
    return(_g0);
}

int Genotype::second()
{
    if(_is_haploid)
    {
	return(bcf_int32_vector_end);
    }
    else
    {
	return(_g1);
    }
}

int Genotype::getGenotype(int idx)
{
    if(idx==0)
    {
	return(_g0);
    }
    
    if(idx==1)
    {
	if(_is_haploid)
	{
	    return(bcf_int32_vector_end);
	}
	else
	{
	    return(_g1);
	}
    }
    cerr << "idx="<<idx<<endl;
    die("invalid genotype.");
    return(-1);
}

int Genotype::swap()
{
    if(!_is_haploid)
    {
	int tmp = _g0;
	_g0 = _g1;
	_g1 = tmp;
    }
    return(0);
}

int Genotype::update_bcf_gt_array(int *gt_array,int index,int32_t *ps_array)
{
    if(!isMissing())
    {
	gt_array[index * 2] = _is_phased ? bcf_gt_phased(first()) : bcf_gt_unphased(first());
	if(isHaploid())
	    gt_array[index * 2 + 1] = bcf_int32_vector_end;
	else
	    gt_array[index * 2 + 1] = _is_phased ? bcf_gt_phased(second()) : bcf_gt_unphased(second()); 
    }
    if(ps_array)
	ps_array[index]=_ps;
    return(0);
}

void Genotype::setPhase(bool phased)
{
    _is_phased=phased;
}

string Genotype::print()
{
    if(isMissing()) return "./.";
    std::stringstream ss;
    ss<<first();
    if(isHaploid())
	return ss.str();
    
    if(_is_phased)
	ss<<"|";
    else
	ss<<"/";
    ss<<second();
    if(_ps!=bcf_int32_missing)
	ss<<":"<<_ps;
    return ss.str();
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
