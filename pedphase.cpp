#include "pedphase.h"

using namespace std;


void usage()
{
    cerr << "usage"<<endl;
}

void swap(int & a, int & b)
{
    int tmp = a;
    a=b;
    b=tmp;
}

//performs simple trio phasing using mendelian inheritance.
//returns
//-1: mendelian inconsistent
//0:  unphaseable
//1:  phased
//int phase_trio(int kid_gt[],int dad_gt[],int mum_gt[])
int PedPhaser::phaseTrio(int idx,int *gt_arr)
{

    int mum = _ped->getMumIndex(idx);
    int dad = _ped->getDadIndex(idx);
    if(mum==-1 || dad==-1)
    {
        return(0);
    }

    int kid_gt[2],dad_gt[2],mum_gt[2];
    getGenotype(idx,gt_arr,kid_gt);
    getGenotype(dad,gt_arr,dad_gt);
    getGenotype(mum,gt_arr,mum_gt);

    bool kid_is_genotyped = is_genotyped(gt_arr, idx);
    bool dad_is_genotyped = is_genotyped(gt_arr, dad);
    bool mum_is_genotyped = is_genotyped(gt_arr, mum);
    if(!kid_is_genotyped||!mum_is_genotyped||!dad_is_genotyped)
    {
        return(0);
    }

    //these loops enumerate the 2**3 possible phase configurations and check which ones are possible
    //ie. search a binary tree with depth 2
#ifdef DEBUG
    cerr << kid_gt[0] << "/" << kid_gt[1]
         << " " << dad_gt[0] << "/" << dad_gt[1] <<  " " << mum_gt[0] << "/" << mum_gt[1] <<endl;
#endif
    vector<bool> phasetree(8,0);//indicator for the 8 possible phase configurations in a trio with two genotypes per sample.
    for (int j = 0; j < 2; j++)
    {
        if (kid_gt[0] != kid_gt[1] || j < 1)
        {
            for (int k = 0; k < 2; k++)
            {
                if (dad_gt[0] != dad_gt[1] || k < 1)
                {
                    for (int l = 0; l < 2; l++)
                    {
                        if (mum_gt[0] != mum_gt[1] || l < 1)
                        {
                            if (kid_gt[j] == dad_gt[k] and kid_gt[(j + 1) % 2] == mum_gt[l])
                            {
                                phasetree[4 * j + 2 * k + l] = 1;
                            }
                        }
                    }
                }
            }
        }
    }

    int sum = accumulate(phasetree.begin(), phasetree.end(), 0);
    if (sum > 1)//multiple solutions - cannot phase
    {
        return(0);
    }
    else if (sum == 1)
    {
        int leaf = find(phasetree.begin(), phasetree.end(), 1) - phasetree.begin();
        if ((leaf >> 2) % 2)
        {
            swap(kid_gt[0], kid_gt[1]);
        }
        if ((leaf >> 1) % 2)
        {
            swap(dad_gt[0], dad_gt[1]);
        }
        if (leaf % 2)
        {
            swap(mum_gt[0], mum_gt[1]);
        }

        gt_arr[idx * 2] = bcf_gt_phased(kid_gt[0]);
        gt_arr[idx * 2 + 1] = bcf_gt_phased(kid_gt[1]);
        gt_arr[dad * 2 ] = bcf_gt_phased(dad_gt[0]);
        gt_arr[dad * 2 + 1] = bcf_gt_phased(dad_gt[1]);
        gt_arr[mum * 2 ] = bcf_gt_phased(mum_gt[0]);
        gt_arr[mum * 2 + 1] = bcf_gt_phased(mum_gt[1]);

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
    int *gt_arr=NULL,*ps_arr=NULL,*gt_arr_dup=NULL;

    int count=0;
    vector< map< int,pair<int,int> > > flip(_nsample);
    int kid_gt[2],dad_gt[2],mum_gt[2];
    for(deque<bcf1_t *>::iterator it1=_line_buffer.begin();it1!=_line_buffer.end();it1++)
    {
        bcf1_t *line = *it1;
        assert(bcf_get_genotypes(_hdr, line, &gt_arr, &ngt)==2*_nsample);
        gt_arr_dup = (int *)realloc(gt_arr_dup,ngt*sizeof(int));
        memcpy(gt_arr_dup,gt_arr,ngt*sizeof(int));

        //wipe any existing phase information.
        for(int i=0;i<(2*_nsample);i++)
        {
            gt_arr_dup[i] = bcf_gt_unphased(bcf_gt_allele(gt_arr_dup[i]));
        }

        for(int i=0;i<_nsample;i++)
        {
            int phase=phaseTrio(i,gt_arr_dup);
            count++;
        }

        if(bcf_get_format_int32(_hdr, line, "PS", &ps_arr, &nps)>0)
        {
            for(int i=0;i<_nsample;i++)
            {
                if(ps_arr[i]!=bcf_int32_missing && bcf_gt_allele(gt_arr[2*i])!=bcf_gt_allele(gt_arr[2*i+1]) && bcf_gt_is_phased(gt_arr_dup[2*i+1]))
                {
                    if(!flip[i].count(ps_arr[i]))
                    {
                        flip[i][ps_arr[i]] = pair<int, int>(0, 0);
                    }
                    flip[i][ps_arr[i]].second++;
                    if(bcf_gt_allele(gt_arr_dup[2*i])!=bcf_gt_allele(gt_arr[2*i]))
                    {
                        flip[i][ps_arr[i]].first++;
                    }
                }
            }
        }
    }

    for(deque<bcf1_t *>::iterator it1=_line_buffer.begin();it1!=_line_buffer.end();it1++)
    {
        bcf1_t *line = *it1;

        assert(bcf_get_genotypes(_hdr, line, &gt_arr, &ngt)==2*_nsample);
        gt_arr_dup = (int *)realloc(gt_arr_dup,ngt*sizeof(int));
        memcpy(gt_arr_dup,gt_arr,ngt*sizeof(int));

        int has_ps = bcf_get_format_int32(_hdr, line, "PS", &ps_arr, &nps);
        vector<bool> flipped(_nsample,0);
        for (int i = 0; i < _nsample; i++)
        {
            int phase = phaseTrio(i,gt_arr_dup);

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
                    if ((*idx)!= -1 && has_ps>0 && ps_arr[*idx] != bcf_int32_missing )
                    {
                        if (!flipped[*idx] && flip[*idx][ps_arr[*idx]].first > flip[*idx][ps_arr[*idx]].second / 2)
                        {
                            int tmp = gt_arr[2 * (*idx)];
                            gt_arr[2 * (*idx)] = bcf_gt_phased(bcf_gt_allele(gt_arr[2 * (*idx) + 1]));
                            gt_arr[2 * (*idx) + 1] = bcf_gt_phased(bcf_gt_allele(tmp));
                            flipped[*idx] = true;
                        }
                    }
                    else if(phase==1 && (*idx)!=-1)
                    {
                        gt_arr[2*(*idx)]=gt_arr_dup[2*(*idx)];
                        gt_arr[2*(*idx)+1]=gt_arr_dup[2*(*idx)+1];
                    }
                }
            }
        }
        bcf_update_genotypes(_hdr, line, gt_arr, ngt);
    }

    while(!_line_buffer.empty())//flushes out the deque
    {
        bcf1_t *tmp_line = _line_buffer.front();
        _line_buffer.pop_front();
        bcf_write(_out_fh,_out_hdr,tmp_line);
        bcf_destroy(tmp_line);
    }

    free(gt_arr);
    free(ps_arr);
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

    bcf_hdr_append(_out_hdr,"##INFO=<ID=MENDELCONFLICT,Number=0,Type=Flag,Description=\"this variant has at least one Mendelian inconsistency\">");
    bcf_hdr_write(_out_fh, _out_hdr);
    if(a.nthreads>0)
    {
        hts_set_threads(_out_fh,a.nthreads);
    }

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
    int *gt_arr=NULL,*ps_arr=NULL;

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

        if(bcf_get_format_int32(_hdr, line, "PS",&ps_arr,&nps)>0)
        {// variant has samples with phase set (PS) set. we need to treat them special
            bcf1_t *new_line = bcf_dup(line);
            _line_buffer.push_back(new_line);
        }
        else
        {//no phase set. phase+flush the deque and perform standard  line-at-a-time phasing by mendelian inheritance.
            flushBuffer();
            int ret = bcf_get_genotypes(_hdr, line, &gt_arr, &ngt);
            bool diploid = ret == 2 * _nsample;
            if (diploid)
            {
                for (int i = 0; i < _nsample; i++)
                {
                     int phase = phaseTrio(i,gt_arr);
                     if(phase == -1)
                     {
                         bcf_update_info_flag(_out_hdr, line, "MENDELCONFLICT", NULL, 1);
                     }
//                    cerr << line->pos+1 <<" "<<phase<<endl;
                }
            }
            else
            {
                if (!diploid_warn)
                {
                    cerr << "\tWARNING: found non-diploid site (" << bcf_hdr_id2name(_hdr, line->rid) << ":"
                         << line->pos + 1 << ") was ignored. ";
                    cerr << "You will only see this warning once." << endl;
                    diploid_warn = true;
                }
            }
            bcf_update_genotypes(_out_hdr, line, gt_arr, ret);
            bcf_write1(_out_fh, _out_hdr, line);
        }
    }
    flushBuffer();

    hts_close(_out_fh);
    free(ps_arr);
    free(gt_arr);
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
