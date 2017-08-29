#include "pedphase.h"

using namespace std;
//#define DEBUG

static void usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About:   akt pedphase - simple Mendel inheritance phasing of duos/trios\n");
    fprintf(stderr, "Usage:   ./akt pedphase input.vcf.gz -p pedigree.fam\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "    -p, --pedigree              pedigree information in plink .fam format\n");
    fprintf(stderr, "    -o, --out                   output a site only vcf.gz annotated with site specific error rates\n");
    fprintf(stderr, "    -@, --threads               number of compression/decompression threads to use\n");
    exit(1);
}

void swap(int & a, int & b)
{
    int tmp = a;
    a=b;
    b=tmp;
}

//performs simple duo/trio phasing using mendelian inheritance.
//returns
//-1: mendelian inconsistent
//0:  unphaseable
//1:  phased
int PedPhaser::mendelPhase(int kid_idx,int *gt_array)
{
    int pedigree_size = 3;//we might make this dynamic later
    int dad_idx = _ped->getDadIndex(kid_idx);
    int mum_idx = _ped->getMumIndex(kid_idx);
    
    Genotype kid_gt(kid_idx,gt_array);
    Genotype dad_gt(dad_idx,gt_array);
    Genotype mum_gt(mum_idx,gt_array);
    
    if( (dad_gt.isMissing() && mum_gt.isMissing()) || kid_gt.isMissing() )
    {
        return(0);
    }

    //phasetree is perfect binary tree (stored as array) that enumerates every genotype configuration in the pedigree
    //the leaf of each tree is 0 if the genotype configuration is inconsistent with inheritance and 1 otherwise.
    //redundant leaves (eg. where a sample is homozygous) are also 0
    vector<bool> phasetree(pow(2,pedigree_size),0);
    //this loop enumerates the 2**n possible phase configurations and checks which are compatible with inheritance
    for(size_t i=0;i<phasetree.size();i++)
    {
	bitset<8> leaf((int)i);
	bool kid_branch = leaf[0];
	bool dad_branch = leaf[1];
	bool mum_branch = leaf[2];

	if( (kid_gt.isHaploid()||kid_gt.isHet()||!kid_branch) && (dad_gt.isHaploid()||dad_gt.isHet()||!dad_branch) && (mum_gt.isHaploid()||mum_gt.isHet()||!mum_branch) )
	{
	    bool is_inheritance_consistent = dad_gt.isMissing()||dad_gt.getGenotype(dad_branch)==kid_gt.getGenotype(kid_branch);
	    is_inheritance_consistent &= mum_gt.isMissing()||mum_gt.getGenotype(mum_branch)==kid_gt.getGenotype((kid_branch+1)%2);
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

	if(!kid_gt.isMissing())
	{
	    gt_array[kid_idx * 2] = bcf_gt_phased(kid_gt.first());
	    if(kid_gt.isHaploid())
	    {
		gt_array[kid_idx * 2 + 1] = bcf_int32_vector_end;
	    }
	    else
	    {
		gt_array[kid_idx * 2 + 1] = bcf_gt_phased(kid_gt.second());	    
	    }
	}
	
        if(!dad_gt.isMissing())
        {
	    gt_array[dad_idx * 2 ] = bcf_gt_phased(dad_gt.first());	    
            if(dad_gt.isHaploid())
            {
                gt_array[dad_idx * 2 + 1] =  bcf_int32_vector_end;				
            }
            else
            {
                gt_array[dad_idx * 2 + 1] = bcf_gt_phased(dad_gt.second());		
            }
        }
	
        if(!mum_gt.isMissing())
        {
	    gt_array[mum_idx * 2 ] = bcf_gt_phased(mum_gt.first());	    
            if(mum_gt.isHaploid())
            {
                gt_array[mum_idx * 2 + 1] =  bcf_int32_vector_end;		
            }
            else
            {
                gt_array[mum_idx * 2 + 1] = bcf_gt_phased(mum_gt.second());		
            }
        }

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
    int *gt_array=NULL,*ps_array=NULL,*gt_array_dup=NULL;

    int count=0;
    vector< map< int,pair<int,int> > > flip(_nsample);
    int kid_gt[2],dad_gt[2],mum_gt[2];
    for(deque<bcf1_t *>::iterator it1=_line_buffer.begin();it1!=_line_buffer.end();it1++)
    {
        bcf1_t *line = *it1;
        assert(bcf_get_genotypes(_hdr, line, &gt_array, &ngt)==2*_nsample);
        gt_array_dup = (int *)realloc(gt_array_dup,ngt*sizeof(int));
        memcpy(gt_array_dup,gt_array,ngt*sizeof(int));

        //wipe any existing phase information.
        for(int i=0;i<(2*_nsample);i++)
        {
            gt_array_dup[i] = bcf_gt_unphased(bcf_gt_allele(gt_array_dup[i]));
        }

        for(int i=0;i<_nsample;i++)
        {
            int phase=mendelPhase(i,gt_array_dup);
            count++;
        }

        if(bcf_get_format_int32(_hdr, line, "PS", &ps_array, &nps)>0)
        {
            for(int i=0;i<_nsample;i++)
            {
                if(ps_array[i]!=bcf_int32_missing && bcf_gt_allele(gt_array[2*i])!=bcf_gt_allele(gt_array[2*i+1]) && bcf_gt_is_phased(gt_array_dup[2*i+1]))
                {
                    if(!flip[i].count(ps_array[i]))
                    {
                        flip[i][ps_array[i]] = pair<int, int>(0, 0);
                    }
                    flip[i][ps_array[i]].second++;
                    if(bcf_gt_allele(gt_array_dup[2*i])!=bcf_gt_allele(gt_array[2*i]))
                    {
                        flip[i][ps_array[i]].first++;
                    }
                }
            }
        }
    }

    for(deque<bcf1_t *>::iterator it1=_line_buffer.begin();it1!=_line_buffer.end();it1++)
    {
        bcf1_t *line = *it1;

        assert(bcf_get_genotypes(_hdr, line, &gt_array, &ngt)==2*_nsample);
        gt_array_dup = (int *)realloc(gt_array_dup,ngt*sizeof(int));
        memcpy(gt_array_dup,gt_array,ngt*sizeof(int));

        int has_ps = bcf_get_format_int32(_hdr, line, "PS", &ps_array, &nps);
        vector<bool> flipped(_nsample,0);
        for (int i = 0; i < _nsample; i++)
        {
            int phase = mendelPhase(i,gt_array_dup);

            if(phase == -1)
            {
                bcf_update_info_flag(_hdr, line, "MENDELCONFLICT", NULL, 1);
            }
            else
            {
                int mum = _ped->getMumIndex(i);
                int dad = _ped->getDadIndex(i);
                vector<int> indices(1,i);
                indices.push_back(dad);
                indices.push_back(mum);
                for(vector<int>::iterator idx=indices.begin();idx!=indices.end();idx++)
                {
                    if ((*idx)!= -1 && has_ps>0 && ps_array[*idx] != bcf_int32_missing )
                    {
                        if (!flipped[*idx] && flip[*idx][ps_array[*idx]].first > flip[*idx][ps_array[*idx]].second / 2)
                        {
                            int tmp = gt_array[2 * (*idx)];
                            gt_array[2 * (*idx)] = bcf_gt_phased(bcf_gt_allele(gt_array[2 * (*idx) + 1]));
                            gt_array[2 * (*idx) + 1] = bcf_gt_phased(bcf_gt_allele(tmp));
                            flipped[*idx] = true;
                        }
                    }
                    else if(phase==1 && (*idx)!=-1)
                    {
                        gt_array[2*(*idx)]=gt_array_dup[2*(*idx)];
                        gt_array[2*(*idx)+1]=gt_array_dup[2*(*idx)+1];
                    }
                }
            }
        }
        bcf_update_genotypes(_hdr, line, gt_array, ngt);
    }

    while(!_line_buffer.empty())//flushes out the deque
    {
        bcf1_t *tmp_line = _line_buffer.front();
        _line_buffer.pop_front();
        bcf_write(_out_fh,_out_hdr,tmp_line);
        bcf_destroy(tmp_line);
    }

    free(gt_array);
    free(ps_array);
    return(0);
}

PedPhaser::PedPhaser(args &a)
{
    //open a file.
    _sr = bcf_sr_init();

    if (a.targets != NULL)
    {
        if (bcf_sr_set_targets(_sr, a.targets, a.targets_is_file, 0) < 0)
        {
            cerr << "ERROR: Failed to set targets " << a.targets << endl;
            exit(1);
        }
    }

    if (a.regions != NULL)
    {
        if (bcf_sr_set_regions(_sr, a.regions, a.regions_is_file) < 0)
        {
            cerr << "ERROR: Failed to read the regions: " << a.regions << endl;;
            exit(1);
        }
    }

    if (bcf_sr_add_reader(_sr, a.inputfile) != 1)
    {
        cerr << "ERROR: problem opening " << a.inputfile << endl;
        exit(1);
    }

    _hdr = _sr->readers[0].header;
    _out_hdr = bcf_hdr_dup(_hdr);
    char output_type[] = "wv";
    output_type[1] = a.output_type;
    _out_fh = hts_open(a.outfile, output_type);
    if(a.nthreads>0)
    {
        bcf_sr_set_threads(_sr,a.nthreads);
        hts_set_threads(_out_fh,a.nthreads);
    }


    bcf_hdr_append(_out_hdr,"##INFO=<ID=MENDELCONFLICT,Number=0,Type=Flag,Description=\"this variant has at least one Mendelian inconsistency\">");
    bcf_hdr_write(_out_fh, _out_hdr);

    if (a.pedigree == NULL)
    {
        _ped = new sampleInfo(_out_hdr);
    }
    else
    {
        _ped = new sampleInfo(a.pedigree, _out_hdr);
    }


    bcf1_t *line = bcf_init1();
    _nsample = bcf_hdr_nsamples(_out_hdr);
    int ngt=0,nps=0;
    int *gt_array=NULL,*ps_array=NULL;

    vector<int> gt(_nsample);

    bool diploid_warn = false;
    int nsnp=0;
    vector<int32_t> phase_set(_nsample,bcf_int32_missing);//stores the current phase set for each sample
    int prev_rid = -1;
    cerr << "Reading input from " << a.inputfile << endl;
    int kid_gt[2],dad_gt[2],mum_gt[2];
    while (bcf_sr_next_line(_sr))
    {
        line = bcf_sr_get_line(_sr, 0);
        bcf_unpack(line, BCF_UN_ALL);

        if(bcf_get_format_int32(_hdr, line, "PS",&ps_array,&nps)>0)
        {// variant has samples with phase set (PS) set. we need to treat them special
            bcf1_t *new_line = bcf_dup(line);
            _line_buffer.push_back(new_line);
        }
        else
        {//no phase set. phase+flush the deque and perform standard  line-at-a-time phasing by mendelian inheritance.
            flushBuffer();
            int ret = bcf_get_genotypes(_hdr, line, &gt_array, &ngt);
            bool diploid = ret > _nsample; //are there any non-haploid genotypes??
            if (diploid)
            {
                for (int i = 0; i < _nsample; i++)
                {
                     int phase = mendelPhase(i,gt_array);
                     if(phase == -1)
                     {
                         bcf_update_info_flag(_out_hdr, line, "MENDELCONFLICT", NULL, 1);
                     }
                }
            }
            else
            {
                if (!diploid_warn)
                {
//                    cerr << "\tWARNING: found non-diploid site (" << bcf_hdr_id2name(_hdr, line->rid) << ":"
//                         << line->pos + 1 << ") was ignored. ";
//                    cerr << "You will only see this warning once." << endl;
                    diploid_warn = true;
                }
            }
            bcf_update_genotypes(_out_hdr, line, gt_array, ret);
            bcf_write1(_out_fh, _out_hdr, line);
        }
    }
    flushBuffer();

    hts_close(_out_fh);
    free(ps_array);
    free(gt_array);
    bcf_sr_destroy(_sr);

    if (nsnp > 0)
    {
        cerr << nsnp << " variants" << endl;
    }

}

int pedphase_main(int argc, char **argv)
{
    int c;
    args a;
    a.output_type = 'v';
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
    a.regions_is_file = false;
    a.targets_is_file = false;
    a.targets = a.pedigree = a.inputfile = a.include = a.regions = NULL;
    a.outfile="-";
    a.nthreads = 0;

    while ((c = getopt_long(argc, argv, "o:p:t:T:r:R:O:@:", loptions, NULL)) >= 0)
    {
        switch (c)
        {
            case 'o':
                a.outfile = optarg;
                break;
            case 'O':
                a.output_type = optarg[0];
                break;
            case 'p':
                a.pedigree = optarg;
                break;
            case 'i':
                a.include = optarg;
                break;
            case 't':
                a.targets = optarg;
                break;
            case 'T':
                a.targets = optarg;
                break;
            case 'r':
                a.regions = optarg;
                break;
            case '@':
                a.nthreads = atoi(optarg);
                break;
            case 'R':
                a.regions = optarg;
                a.regions_is_file = true;
                break;
            default:
                die("unknown argument");
        }
    }
    optind++;
    a.inputfile = argv[optind];
    if (a.inputfile == NULL)
    {
        die("no input provided");
    }

    cerr << "Output file: " << a.outfile<<endl;
    PedPhaser p(a);
    return (0);
}
