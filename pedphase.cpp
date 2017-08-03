#include "pedphase.h"

using namespace std;

typedef struct _args
{
    bool regions_is_file;
    bool targets_is_file;
    char output_type;
    char *pedigree, *inputfile, *include, *regions, *targets, *outfile;
} args;

void usage()
{
    cerr << "usage"<<endl;
}

void swap(int & a, int & b)
{
    int tmp = a;
    a = b;
    b = tmp;
}

//performs simple trio phasing using mendelian inheritance.
//returns
//-1: mendelian inconsistent
//0:  unphaseable
//1:  phased
int phase_trio(int kid_gt[],int dad_gt[],int mum_gt[])
{
    //these loops enumerate the 2**3 possible phase configurations and check which ones are possible
    //ie. search a binary tree with depth 2
    //cerr << kid_gt[0] << "/" << kid_gt[1] << " " << dad_gt[0] << "/" << dad_gt[1] <<  " " << mum_gt[0] << "/" << mum_gt[1] <<endl;
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
        return(1);
    }
    else//inconsistent with mendelian inheritance.
    {
        return(-1);
    }
}

int phaseBuffer(deque<bcf1_t*> & buf,bcf_hdr_t *hdr,sampleInfo &ped)
{
    int nsample = bcf_hdr_nsamples(hdr);
    int ngt=0,nps=0;
    int *gt_arr=NULL,*ps_arr=NULL;

    for(int i=0;i<nsample;i++)
    {
        int mum = ped.getMumIndex(i);
        int dad = ped.getDadIndex(i);
        if(dad>-1 && mum>-1)
        {
            for(deque<bcf1_t *>::iterator it1=buf.begin();it1!=buf.end();it1++)
            {
                bcf1_t *line = *it1;
                int ret = bcf_get_genotypes(hdr, line, &gt_arr, &ngt);
                if(bcf_get_format_int32(hdr, line, "PS",&ps_arr,&nps)>0)
                {
                    continue;
                }
                else
                {
                    continue;
                }
            }
        }
    }
    return(0);
}

int pedphase(args &a)
{
    //open a file.
    bcf_srs_t *sr = bcf_sr_init();

    if (a.targets != NULL)
    {
        if (bcf_sr_set_targets(sr, a.targets, a.targets_is_file, 0) < 0)
        {
            cerr << "ERROR: Failed to set targets " << a.targets << endl;
            exit(1);
        }
    }

    if (a.regions != NULL)
    {
        if (bcf_sr_set_regions(sr, a.regions, a.regions_is_file) < 0)
        {
            cerr << "ERROR: Failed to read the regions: " << a.regions << endl;;
            exit(1);
        }
    }

    if (bcf_sr_add_reader(sr, a.inputfile) != 1)
    {
        cerr << "ERROR: problem opening " << a.inputfile << endl;
        exit(1);
    }

    bcf_hdr_t *hdr = sr->readers[0].header;
    htsFile *out_fh = NULL;
    bcf_hdr_t *out_hdr = NULL;
    out_hdr = bcf_hdr_dup(hdr);
    char output_type[] = "wv";
    output_type[1] = a.output_type;
    out_fh = hts_open(a.outfile, output_type);
    bcf_hdr_append(out_hdr,"##INFO=<ID=MENDELCONFLICT,Number=0,Type=Flag,Description=\"this variant has at least one Mendelian inconsistency\">");
    bcf_hdr_write(out_fh, out_hdr);

    sampleInfo ped(a.pedigree, out_hdr);

    bcf1_t *line = bcf_init1();
    int nsample = bcf_hdr_nsamples(out_hdr);
    int ngt=0,nps=0;
    int *gt_arr=NULL,*ps_arr=NULL;

    vector<int> gt(nsample);

    bool diploid_warn = false;
    int nsnp=0;
    deque<bcf1_t *> line_buffer; //this buffers phase sets when necessary.
    vector<int32_t> phase_set(nsample,bcf_int32_missing);//stores the current phase set for each sample
    int prev_rid = -1;
    cerr << "Reading input from " << a.inputfile << endl;
    while (bcf_sr_next_line(sr))
    {
        line = bcf_sr_get_line(sr, 0);
        bcf_unpack(line, BCF_UN_ALL);

        if(bcf_get_format_int32(hdr, line, "PS",&ps_arr,&nps)>0)
        {// variant has samples with phase set (PS) set. we need to treat them special
            bcf1_t *new_line = bcf_dup(line);
            line_buffer.push_back(new_line);
        }
        else
        {//no phase set. phase+flush the deque and perform standard  line-at-a-time phasing by mendelian inheritance.
            phaseBuffer(line_buffer,out_hdr,ped);//phases all records in deque
            while(!line_buffer.empty())//flushes out the deque
            {
                bcf1_t *tmp_line = line_buffer.front();
                line_buffer.pop_front();
                bcf_write(out_fh,out_hdr,tmp_line);
                bcf_destroy(tmp_line);
            }

            int ret = bcf_get_genotypes(hdr, line, &gt_arr, &ngt);
            bool diploid = ret == 2 * nsample;
            if (diploid)
            {
                for (int i = 0; i < nsample; i++)
                {
                    int mum = ped.getMumIndex(i);
                    int dad = ped.getDadIndex(i);

                    if (ped.isPhaseableTrio(i,gt_arr) )
                    {
                        int kid_gt[2],dad_gt[2],mum_gt[2];
                        getGenotype(i,gt_arr,kid_gt);
                        getGenotype(dad,gt_arr,dad_gt);
                        getGenotype(mum,gt_arr,mum_gt);
                        int phase = phase_trio(kid_gt,dad_gt,mum_gt);
                        if(phase==1)
                        {
                            gt_arr[i * 2] = bcf_gt_phased(kid_gt[0]);
                            gt_arr[i * 2 + 1] = bcf_gt_phased(kid_gt[1]);
                            gt_arr[dad * 2 + 1] = bcf_gt_phased(dad_gt[0]);
                            gt_arr[dad * 2 + 1] = bcf_gt_phased(dad_gt[1]);
                            gt_arr[mum * 2 + 1] = bcf_gt_phased(mum_gt[0]);
                            gt_arr[mum * 2 + 1] = bcf_gt_phased(mum_gt[1]);
                        }
                        if(phase == -1)
                        {
                            bcf_update_info_flag(out_hdr, line, "MENDELCONFLICT", NULL, 1);
                        }
                    }
                    else if (dad != -1)
                    {
                        int kid_gt[2] = {bcf_gt_allele(gt_arr[i * 2]), bcf_gt_allele(gt_arr[i * 2 + 1])};
                        int dad_gt[2] = {bcf_gt_allele(gt_arr[dad * 2]), bcf_gt_allele(gt_arr[dad * 2 + 1])};
                        continue;
                    }
                    else if (mum != -1)
                    {
                        int kid_gt[2] = {bcf_gt_allele(gt_arr[i * 2]), bcf_gt_allele(gt_arr[i * 2 + 1])};
                        int mum_gt[2] = {bcf_gt_allele(gt_arr[mum * 2]), bcf_gt_allele(gt_arr[mum * 2 + 1])};
                        continue;
                    }
                }
                nsnp++;
            } else
            {
                if (!diploid_warn)
                {
                    cerr << "\tWARNING: found non-diploid site (" << bcf_hdr_id2name(hdr, line->rid) << ":"
                         << line->pos + 1 << ") was ignored. ";
                    cerr << "You will only see this warning once." << endl;
                    diploid_warn = true;
                }
            }
            bcf_update_genotypes(out_hdr, line, gt_arr, ngt);
            bcf_write1(out_fh, out_hdr, line);
        }
    }

    hts_close(out_fh);
    free(ps_arr);
    free(gt_arr);
    bcf_sr_destroy(sr);

    if (nsnp > 0)
    {
        cerr << nsnp << " variants" << endl;
    }

    return (0);
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
            {"targets",      required_argument, NULL, 't'},
            {"targets-file", required_argument, NULL, 'T'},
            {"regions-file", required_argument, NULL, 'R'},
            {"regions",      required_argument, NULL, 'r'},
            {0,          0, 0,                        0}
    };
    a.regions_is_file = false;
    a.targets_is_file = false;
    a.targets = a.pedigree = a.inputfile = a.include = a.regions = a.outfile = NULL;

    while ((c = getopt_long(argc, argv, "o:p:t:T:r:R:O:", loptions, NULL)) >= 0)
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
            case 'R':
                a.regions = optarg;
                a.regions_is_file = true;
                break;
        }
    }
    optind++;
    a.inputfile = argv[optind];
    if (a.inputfile == NULL)
    {
        die("no input provided");
    }
    if (a.pedigree == NULL)
    {
        die("the -p argument is required");
    }
    cerr << "Output file: " << a.outfile<<endl;
    pedphase(a);
    return (0);
}
