#include "pedphase.h"

using namespace std;
//#define DEBUG
#define NUM_LEAVES 8

static void usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   akt pedphase - simple Mendel inheritance phasing of duos/trios\n");
    fprintf(stderr, "Usage:   ./akt pedphase input.vcf.gz -p pedigree.fam\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -p, --pedigree                 pedigree information in plink .fam format\n");
    fprintf(stderr, "    -o, --output-file <file>       output file name [stdout]\n");
    fprintf(stderr, "    -O, --output-type <b|u|z|v>    b: compressed BCF, u: uncompressed BCF, z: compressed VCF, v: uncompressed VCF [v]\n");
    fprintf(stderr, "    -@, --threads                  number of compression/decompression threads to use\n");
    exit(1);
}


static int bcf_int32_count_missing(int32_t *array,int n)
{
    int num_missing=0;
    for(int i=0;i<n;i++)
    {
	if(array[i]==bcf_int32_missing)
	{
	    num_missing++;
	}
    }
    return(num_missing);
}


static void bcf_int32_set_missing(int32_t *array,int n)
{
    int num_missing=0;
    for(int i=0;i<n;i++)
    {
	array[i]=bcf_int32_missing;
    }
}

//performs simple duo/trio phasing using mendelian inheritance.
//returns
//-1: mendelian inconsistent
//0:  unphaseable
//1:  phased
int PedPhaser::mendelPhase(int kid_index,int *gt_array,int *ps_array)
{
    int pedigree_size = 3;//we might make this dynamic later
    int dad_index = _pedigree->getDadIndex(kid_index);
    int mum_index = _pedigree->getMumIndex(kid_index);
    
    Genotype kid_gt(kid_index,gt_array,ps_array);
    Genotype dad_gt(dad_index,gt_array,ps_array);
    Genotype mum_gt(mum_index,gt_array,ps_array);
    
    if( (dad_gt.isMissing() && mum_gt.isMissing()) || kid_gt.isMissing() )
    {
        return(0);
    }

    //phasetree is perfect binary tree (stored as array) that enumerates every genotype configuration in the pedigree
    //the leaf of each tree is 0 if the genotype configuration is inconsistent with inheritance and 1 otherwise.
    //redundant leaves (eg. where a sample is homozygous) are also 0
    vector<bool> phasetree(pow(2,pedigree_size),0);
    assert(phasetree.size()<=NUM_LEAVES);
    //this loop enumerates the 2**n possible phase configurations and checks which are compatible with inheritance
    for(size_t i=0;i<phasetree.size();i++)
    {
	bitset<NUM_LEAVES> leaf((int)i);
	bool kid_branch = leaf[0];
	bool dad_branch = leaf[1];
	bool mum_branch = leaf[2];

//	if( (kid_gt.isHaploid()||kid_gt.isHet()||!kid_branch) && (dad_gt.isHaploid()||dad_gt.isHet()||!dad_branch) && (mum_gt.isHaploid()||mum_gt.isHet()||!mum_branch) )
	if( (!kid_gt.isPhased()||!kid_branch) && (!dad_gt.isPhased()||!dad_branch) && (!mum_gt.isPhased()||!mum_branch) )
	{
	    bool is_inheritance_consistent = dad_gt.isMissing()||dad_gt.getGenotype(dad_branch)==kid_gt.getGenotype((kid_branch+1)%2);
	    is_inheritance_consistent &= mum_gt.isMissing()||mum_gt.getGenotype(mum_branch)==kid_gt.getGenotype(kid_branch);
	    if(is_inheritance_consistent)
	    {
		phasetree[i] = 1;
	    }
	}
    }

    int sum = accumulate(phasetree.begin(), phasetree.end(), 0);
    if (sum > 1)//multiple solutions - cannot phase
    {
        return(0);
    }
    else if (sum == 1)//found a unique solution for phasing. update the genotyping array.
    {
        int leaf = find(phasetree.begin(), phasetree.end(), 1) - phasetree.begin();
        if ((leaf >> 0) % 2)
        {
	    kid_gt.swap();
        }
        if ((leaf >> 1) % 2)
        {
	    dad_gt.swap();
        }
        if (leaf >> 2)
        {
	    mum_gt.swap();
        }

	kid_gt.update_bcf_gt_array(gt_array,kid_index);
	dad_gt.update_bcf_gt_array(gt_array,dad_index);
	mum_gt.update_bcf_gt_array(gt_array,mum_index);

        return(1);
    }
    else//inconsistent with mendelian inheritance.
    {
        return(-1);
    }
}

int PedPhaser::flushBuffer()
{
    if(_line_buffer.empty())
    {
        return(0);
    }

    int ngt=0,nps=0;
    int count=0;
    vector< map< int,pair<int,int> > > flip(_num_sample); // the key stores PS and the pair stores <number of inconsistencies,number of phased variants>
    int kid_gt[2],dad_gt[2],mum_gt[2];

    //this loops over the buffered vcf rows and performs mendel phasing
    //the number of inconsistencies between the mendel phased genotypes and the PS phased genotypes is stored in flip
    for(deque<bcf1_t *>::iterator it1=_line_buffer.begin();it1!=_line_buffer.end();it1++)
    {
        bcf1_t *line = *it1;
        assert(bcf_get_genotypes(_out_header, line, &_gt_array, &ngt)==2*_num_sample);
        _gt_array_dup = (int *)realloc(_gt_array_dup,ngt*sizeof(int));
        memcpy(_gt_array_dup,_gt_array,ngt*sizeof(int));


        for(int i=0;i<(2*_num_sample);i++)//wipe any existing phase information.
        {
            _gt_array_dup[i] = bcf_gt_unphased(bcf_gt_allele(_gt_array_dup[i]));
        }

        for(int i=_num_sample-1;i>=0;i--)//mendel phase each child
        {
            int phase=mendelPhase(i,_gt_array_dup);
            count++;
        }

        if(bcf_get_format_int32(_out_header, line, "PS", &_ps_array, &nps)>0)
        {
            for(int i=0;i<_num_sample;i++)
            {
		//if:
		//1. sample has a phase set
		//2. sample is het
		//3. sample is also mendel phased
		//then count the consistency with the phase set
                if(_ps_array[i]!=bcf_int32_missing && bcf_gt_allele(_gt_array[2*i])!=bcf_gt_allele(_gt_array[2*i+1]) && bcf_gt_is_phased(_gt_array_dup[2*i+1]))
                {
                    if(!flip[i].count(_ps_array[i]))
                    {
                        flip[i][_ps_array[i]] = pair<int, int>(0, 0);
                    }
                    flip[i][_ps_array[i]].second++;
                    if(bcf_gt_allele(_gt_array_dup[2*i])!=bcf_gt_allele(_gt_array[2*i]))
                    {
                        flip[i][_ps_array[i]].first++;
                    }
                }
            }
        }
    }

    //now iterate through the buffer again and flip phase sets that are inconsistent with pedigrees
    for(deque<bcf1_t *>::iterator it1=_line_buffer.begin();it1!=_line_buffer.end();it1++)
    {
        bcf1_t *line = *it1;
        assert(bcf_get_genotypes(_out_header, line, &_gt_array, &ngt)==2*_num_sample);
        _gt_array_dup = (int *)realloc(_gt_array_dup,ngt*sizeof(int));
        memcpy(_gt_array_dup,_gt_array,ngt*sizeof(int));
        bcf_int32_set_missing(_rps_array,_num_sample);
        int num_ps_values = bcf_get_format_int32(_out_header, line, "PS", &_ps_array, &nps);
        vector<bool> sample_has_been_flipped(_num_sample,0);

        for (int i=_num_sample-1;i>=0;i--)
        {
            int phase = mendelPhase(i,_gt_array_dup);
            if(phase == -1)
            {
                bcf_update_info_flag(_out_header, line, "MENDELCONFLICT", NULL, 1);
            }
            else
            {
		int mum = _pedigree->getMumIndex(i);
                int dad = _pedigree->getDadIndex(i);
                vector<int> indices(1,i);
                indices.push_back(dad);
                indices.push_back(mum);
                for(vector<int>::iterator index=indices.begin();index!=indices.end();index++)
                {
		    bool sample_not_missing = (*index)!= -1;
 		    bool sample_has_ps = sample_not_missing && num_ps_values>0 && _ps_array[*index] != bcf_int32_missing;
		    if(sample_not_missing)
		    {
			if (sample_has_ps)
			{
			    int num_hets_in_ps=flip[*index][_ps_array[*index]].second;
			    int num_inconsistencies_with_ps = flip[*index][_ps_array[*index]].first;
			    bool phase_set_is_inconsistent_with_pedigree = num_inconsistencies_with_ps > num_hets_in_ps / 2;
			    if (!sample_has_been_flipped[*index] && phase_set_is_inconsistent_with_pedigree)
			    {
				int tmp = _gt_array[2 * (*index)];
				_gt_array[2 * (*index)    ] = bcf_gt_phased(bcf_gt_allele(_gt_array[2 * (*index) + 1]));
				_gt_array[2 * (*index) + 1] = bcf_gt_phased(bcf_gt_allele(tmp));
				sample_has_been_flipped[*index] = true;
			    }
			    if(num_hets_in_ps>0)//this phase set also had pedigree information
			    {
				if(num_inconsistencies_with_ps==0 || num_inconsistencies_with_ps==num_hets_in_ps)
				{//we have a 100% agreement hence move PS -> RPS
				    _rps_array[i] = _ps_array[i];
				}
			    }				
			}
			else if(phase==1) //was mendel phased and not-missing
			{
			    _gt_array[2*(*index)]=_gt_array_dup[2*(*index)];
			    _gt_array[2*(*index)+1]=_gt_array_dup[2*(*index)+1];
			}
		    }		    
                }
            }
        }
	//do a final pass of phase-by-transmission to propagate read-back phasing throughout the pedigree
        for (int i=_num_sample-1;i>=0;i--)
        {
            int phase = mendelPhase(i,_gt_array,_rps_array);
	}	

        bcf_update_genotypes(_out_header, line, _gt_array, ngt);
	if(bcf_int32_count_missing(_rps_array,_num_sample)<_num_sample)
	{
	    if(memcmp(_rps_array,_ps_array,_num_sample)==0)
	    {//REMOVE FORMAT/PS
		bcf_update_format_int32(_out_header,line,"PS",NULL,0);
	    }
	    assert(bcf_update_format_int32(_out_header,line,"RPS",_rps_array,_num_sample)==0);
	}
    } //end for(deque<bcf1_t *>::iterator it1=_line_buffer.begin();it1!=_line_buffer.end();it1++)

    while(!_line_buffer.empty())//flushes out the deque
    {
        bcf1_t *tmp_line = _line_buffer.front();
        _line_buffer.pop_front();
        bcf_write(_out_file,_out_header,tmp_line);
        bcf_destroy(tmp_line);
    }

    return(0);
}

void PedPhaser::setup_io(args &a)
{
        //open a file.
    _bcf_reader = bcf_sr_init();

    if (a.targets != NULL)
    {
        if (bcf_sr_set_targets(_bcf_reader, a.targets, a.targets_is_file, 0) < 0)
        {
            cerr << "ERROR: Failed to set targets " << a.targets << endl;
            exit(1);
        }
    }

    if (a.regions != NULL)
    {
        if (bcf_sr_set_regions(_bcf_reader, a.regions, a.regions_is_file) < 0)
        {
            cerr << "ERROR: Failed to read the regions: " << a.regions << endl;;
            exit(1);
        }
    }

    if (bcf_sr_add_reader(_bcf_reader, a.inputfile) != 1)
    {
        cerr << "ERROR: problem opening " << a.inputfile << endl;
        exit(1);
    }

    _in_header = _bcf_reader->readers[0].header;
    _out_header = bcf_hdr_dup(_in_header);
    char output_type[] = "wv";
    output_type[1] = a.output_type;
    _out_file = hts_open(a.outfile, output_type);
    if(a.nthreads>0)
    {
        bcf_sr_set_threads(_bcf_reader,a.nthreads);
        hts_set_threads(_out_file,a.nthreads);
    }

    bcf_hdr_append(_out_header,"##FORMAT=<ID=RPS,Number=1,Type=Integer,Description=\"Read back phase set.\">");
    bcf_hdr_append(_out_header,"##INFO=<ID=MENDELCONFLICT,Number=0,Type=Flag,Description=\"this variant has at least one Mendelian inconsistency\">");
    bcf_hdr_write(_out_file, _out_header);

    if (a.pedigree == NULL)
    {
        _pedigree = new sampleInfo(_out_header);
    }
    else
    {
        _pedigree = new sampleInfo(a.pedigree, _out_header);
    }
    if(_pedigree->N<=0)
    {
	die("no pedigree detected");
    }
}

PedPhaser::PedPhaser(args &a)
{
    setup_io(a);
    
    bcf1_t *line;
    _num_sample = bcf_hdr_nsamples(_out_header);
    int ngt=0,nps=0;
    _gt_array=NULL;
    _ps_array=NULL;
    _gt_array=NULL;
    _gt_array_dup=NULL;        
    _rps_array=(int32_t *)malloc(_num_sample * sizeof(int32_t));

    vector<int> gt(_num_sample);

    bool diploid_warn = false;
    int nsnp=0;
    vector<int32_t> phase_set(_num_sample,bcf_int32_missing);//stores the current phase set for each sample
    int prev_rid = -1;
    cerr << "Reading input from " << a.inputfile << endl;
    int kid_gt[2],dad_gt[2],mum_gt[2];
    while (bcf_sr_next_line(_bcf_reader))
    {
        line = bcf_sr_get_line(_bcf_reader, 0);
        bcf_unpack(line, BCF_UN_ALL);

        if(bcf_get_format_int32(_in_header, line, "PS",&_ps_array,&nps)>0)
        {// variant has samples with phase set (PS) set. we need to buffer these.
            bcf1_t *new_line = bcf_dup(line);
            _line_buffer.push_back(new_line);
        }
        else
        {//no phase set. phase+flush the deque and then perform standard phase-by-transmission on the current line
            flushBuffer();
            int ret = bcf_get_genotypes(_in_header, line, &_gt_array, &ngt);
            bool diploid = ret > _num_sample; //are there any non-haploid genotypes??
            if (diploid)
            {
                for (int i=_num_sample-1;i>=0;i--)
                {
                     int phase = mendelPhase(i,_gt_array);
                     if(phase == -1)
                     {
                         bcf_update_info_flag(_out_header, line, "MENDELCONFLICT", NULL, 1);
                     }
                }
            }
            bcf_update_genotypes(_out_header, line, _gt_array, ret);
            bcf_write1(_out_file, _out_header, line);
        }
    }
    flushBuffer();
}

PedPhaser::~PedPhaser()
{
    delete _pedigree;
    hts_close(_out_file);
    free(_ps_array);
    free(_gt_array);
    free(_gt_array_dup);    
    free(_rps_array);    
    bcf_sr_destroy(_bcf_reader);
    bcf_hdr_destroy(_out_header);
}

int pedphase_main(int argc, char **argv)
{
    int c;
    args arguments;
    arguments.output_type = 'v';
    if (argc < 3) usage();
    static struct option loptions[] = {
            {"out",      1, 0,                        'o'},
            {"output-type",      1, 0,                'O'},
            {"pedigree", 1, 0,                        'p'},
            {"threads",      required_argument, NULL, '@'},
            {"targets",      required_argument, NULL, 't'},
            {"targets-file", required_argument, NULL, 'T'},
            {"regions-file", required_argument, NULL, 'R'},
            {"regions",      required_argument, NULL, 'r'},
            {0,          0, 0,                        0}
    };
    arguments.regions_is_file = false;
    arguments.targets_is_file = false;
    arguments.targets = arguments.pedigree = arguments.inputfile = arguments.include = arguments.regions = NULL;
    arguments.outfile="-";
    arguments.nthreads = 0;

    while ((c = getopt_long(argc, argv, "o:p:t:T:r:R:O:@:", loptions, NULL)) >= 0)
    {
        switch (c)
        {
            case 'o':
                arguments.outfile = optarg;
                break;
            case 'O':
                arguments.output_type = optarg[0];
                break;
            case 'p':
                arguments.pedigree = optarg;
                break;
            case 'i':
                arguments.include = optarg;
                break;
            case 't':
                arguments.targets = optarg;
                break;
            case 'T':
                arguments.targets = optarg;
                break;
            case 'r':
                arguments.regions = optarg;
                break;
            case '@':
                arguments.nthreads = atoi(optarg);
                break;
            case 'R':
                arguments.regions = optarg;
                arguments.regions_is_file = true;
                break;
            default:
                die("unknown argument");
        }
    }
    optind++;
    arguments.inputfile = argv[optind];
    if (arguments.inputfile == NULL)
    {
        die("no input provided");
    }

    cerr << "Output file: " << arguments.outfile<<endl;
    PedPhaser p(arguments);
    return (0);
}

