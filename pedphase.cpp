#include "pedphase.hh"

using namespace std;

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
    fprintf(stderr, "    -x, --exclude-chromosome       leave these chromosomes unphased (unphased lines will still be in in output)  eg. -x chrM,chrY\n");
    exit(1);
}

static int bcf_int32_count_missing(int32_t *array, int n)
{
    int num_missing = 0;
    for (int i = 0; i < n; i++)
    {
        if (array[i] == bcf_int32_missing)
        {
            num_missing++;
        }
    }
    return (num_missing);
}

static void bcf_int32_set_missing(int32_t *array, int n)
{
    for (int i = 0; i < n; i++)
    {
        array[i] = bcf_int32_missing;
    }
}

PedPhaser::PedPhaser(args &a)
{
    setup_io(a);
    cerr << "Reading input from " << a.inputfile << endl;    
    _num_sample = bcf_hdr_nsamples(_out_header);
    _parental_genotypes.assign(2*_num_sample,pair<int,int>(bcf_gt_missing,bcf_gt_missing));
    _gt_array = nullptr;
    _ps_array = nullptr;
    _rps_array = nullptr;
    _gt_array = nullptr;
    _gt_array_dup = nullptr;
    _num_gt=0;
    _num_ps=0;
    if (!a.exclude_chromosomes.empty())
    {
        cerr << "Will not phase chromosomes: " << a.exclude_chromosomes << " (--exclude-chromosomes)" << endl;
        vector<string> chromosomes;
        stringSplit(a.exclude_chromosomes, ',', chromosomes);
        if (chromosomes.empty())
        {
            die("Problem parsing --exclude-chromosomes: " + a.exclude_chromosomes);
        }
        for (vector<string>::iterator chromosome = chromosomes.begin(); chromosome != chromosomes.end(); chromosome++)
        {
            int rid = bcf_hdr_id2int(_in_header, BCF_DT_CTG, chromosome->c_str());
            if (rid == -1)
            {
                cerr << "WARNING: chromosome " << *chromosome << " was not in header" << endl;
            }
            else
            {
                _chromosomes_to_ignore.push_back(rid);
            }
        }
    }
    main();
}

void PedPhaser::main()
{
    bcf1_t *line;
    vector<int> gt(_num_sample);
    vector<int32_t> phase_set(_num_sample, bcf_int32_missing); //stores the current phase set for each sample
    while (bcf_sr_next_line(_bcf_reader))
    {
        line = bcf_sr_get_line(_bcf_reader, 0);
        if (chromosome_is_in_ignore_list(line))
        {
            flush_buffer(); //just in case something was sitting in buffer from previous chromosome
            bcf_write(_out_file, _out_header, line);
        }
        else
        {
            bcf_unpack(line, BCF_UN_ALL);
            if (bcf_get_format_int32(_in_header, line, "PS", &_ps_array, &_num_ps) > 0)
            { // This variant has samples with phase set (PS) set. We need to buffer these.
                bcf1_t *new_line = bcf_dup(line);
                _line_buffer.push_back(new_line);
            }
            else
            { //No phase set. We flush the deque and then perform standard phase-by-transmission on the current line.
#ifdef DEBUG
		cerr << line->pos+1<<endl;
#endif		
                flush_buffer();
                int ret = bcf_get_genotypes(_in_header, line, &_gt_array, &_num_gt);
                bool diploid = ret > _num_sample; //are there any non-haploid genotypes?
                if (diploid)
                {
		    _sample_has_been_phased.assign(_num_sample,false);		    
                    for (int i = 0; i <_num_sample ; i++)
                    {
			if (mendel_phase(i, _gt_array)==-1)
			    bcf_update_info_flag(_out_header, line, "MENDELCONFLICT", nullptr, 1);
                    }
                }
                bcf_update_genotypes(_out_header, line, _gt_array, ret);
                bcf_write(_out_file, _out_header, line);
            }
        }
    }
    flush_buffer();
}

//performs simple duo/trio phasing using mendelian inheritance.
//returns
//-1: mendelian inconsistent
//0:  unphaseable
//1:  phased
int phase_by_transmission(Genotype & kid_gt,Genotype & dad_gt,Genotype & mum_gt)
{
    const int NUM_LEAVES= 8;//maximum number of leaves on the binary tree(we are only doing duos/trios so it never gets this big)    
    int pedigree_size = 3; //we might make this dynamic later
    if ((dad_gt.isMissing() && mum_gt.isMissing()) || kid_gt.isMissing())  return (0); //unphaseable due to missingness

    //phasetree is perfect binary tree (stored as array) that enumerates every genotype configuration in the pedigree
    //the leaf of each tree is 0 if the genotype configuration is inconsistent with inheritance and 1 otherwise.
    //redundant leaves (eg. where a sample is homozygous) are also 0
    vector<bool> phasetree(pow(2, pedigree_size), 0);
    assert(phasetree.size() <= NUM_LEAVES);
    //this loop enumerates the 2**n possible phase configurations and checks which are compatible with inheritance
    for (size_t i = 0; i < phasetree.size(); i++)
    {
        bitset<NUM_LEAVES> leaf((int)i);
        bool kid_branch = leaf[0];
        bool dad_branch = leaf[1];
        bool mum_branch = leaf[2];
        if ((!kid_gt.is_phased() || !kid_branch) && (!dad_gt.is_phased() || !dad_branch) && (!mum_gt.is_phased() || !mum_branch))
        {
            bool is_inheritance_consistent = dad_gt.isMissing() || dad_gt.getGenotype(dad_branch) == kid_gt.getGenotype((kid_branch + 1) % 2);
            is_inheritance_consistent &= mum_gt.isMissing() || mum_gt.getGenotype(mum_branch) == kid_gt.getGenotype(kid_branch);
            if (is_inheritance_consistent)
            {
                phasetree[i] = 1;
            }
        }
    }
    int sum = accumulate(phasetree.begin(), phasetree.end(), 0);

    if (sum > 1)
    {
        return (0);  //multiple solutions - cannot phase
    }
    else if (sum == 1) //found a unique solution for phasing. update the genotype array.
    {
        int leaf = find(phasetree.begin(), phasetree.end(), 1) - phasetree.begin();
        if ((leaf >> 0) % 2)   kid_gt.swap();
        if ((leaf >> 1) % 2)   dad_gt.swap();
        if (leaf >> 2)         mum_gt.swap();

	kid_gt.setPhase(true);
	dad_gt.setPhase(true);
	mum_gt.setPhase(true);
        return (1);
    }
    else //inconsistent with mendelian inheritance.
    {
        return (-1);
    }
}

bool is_mendel_inconsistent(Genotype  kid,Genotype  dad,Genotype  mum)
{
    if(kid.isMissing() || !kid.is_phased()) return false;
    int k0=kid.first();
    int k1=kid.second();
    int m_transmitted = k0;
    int d_transmitted = k1;
    if(!dad.isMissing())
	d_transmitted=dad.first();
    if(!mum.isMissing())
	m_transmitted=mum.first();
    return(k0!=m_transmitted || k1!=d_transmitted);
}

bool is_mendel_inconsistent(pair<int,int> kid,pair<int,int> dad,pair<int,int> mum)
{
    if(bcf_gt_is_missing(kid.first) || bcf_gt_is_missing(kid.second))
	return false;
    int k0=kid.first;
    int k1=kid.second;
    int m_transmitted = k0;
    int d_transmitted = k1;
    if(!bcf_gt_is_missing(dad.first)&&!bcf_gt_is_missing(dad.second))
	d_transmitted=dad.first;
    if(!bcf_gt_is_missing(mum.first)&&!bcf_gt_is_missing(mum.second))
	m_transmitted=mum.first;    
    return(k0!=m_transmitted || k1!=d_transmitted);
}

int PedPhaser::mendel_phase(int kid_index, int *gt_array, int *ps_array)
{
    int dad_index = _pedigree->getDadIndex(kid_index);
    int mum_index = _pedigree->getMumIndex(kid_index);
    Genotype kid_gt(kid_index, gt_array, ps_array);
    Genotype dad_gt(dad_index, gt_array, ps_array);
    Genotype mum_gt(mum_index, gt_array, ps_array);        
    if(mum_index>=0 || dad_index>=0) _sample_has_been_phased[kid_index]=true;
    bool update_dad=false,update_mum=false;
    if(mum_index>=0)
    {
	update_mum = !_sample_has_been_phased[mum_index];        
	_sample_has_been_phased[mum_index]=true;
    }
    if(dad_index>=0)
    {
	update_dad = !_sample_has_been_phased[dad_index];	
	_sample_has_been_phased[dad_index]=true;
    }
    int ret = phase_by_transmission(kid_gt,dad_gt,mum_gt);
    if(ret)
    {
	if(dad_index>=0) _parental_genotypes[kid_index*2]   = pair<int,int>(dad_gt.first(),dad_gt.second());
	if(mum_index>=0) _parental_genotypes[kid_index*2+1] = pair<int,int>(mum_gt.first(),mum_gt.second());	
	if(update_dad) dad_gt.update_bcf_gt_array(gt_array, dad_index);
	if(update_mum) mum_gt.update_bcf_gt_array(gt_array, mum_index);
	kid_gt.update_bcf_gt_array(gt_array, kid_index);
    }    
    return(ret);
}

bool PedPhaser::chromosome_is_in_ignore_list(bcf1_t *record)
{
    for (vector<int>::iterator rid = _chromosomes_to_ignore.begin(); rid != _chromosomes_to_ignore.end(); rid++)
        if (*rid == record->rid)
            return true;

    return false;
}

int PedPhaser::flush_buffer()
{
    if (_line_buffer.empty()) return(0); 
    int _num_gt=0,_num_ps=0;
    HaplotypeBuffer hap_transmission(_num_sample,_pedigree);//stores the transmission phased haplotypes
    HaplotypeBuffer hap_phaseset(_num_sample,_pedigree);//stores the phase-set phased haplotypes
    for (deque<bcf1_t *>::iterator it1 = _line_buffer.begin(); it1 != _line_buffer.end(); it1++)
    {
	bcf1_t *line = *it1;	
	assert(bcf_get_genotypes(_out_header, line, &_gt_array, &_num_gt) == 2 * _num_sample);
	hap_transmission.push_back(_gt_array);
	int status = bcf_get_format_int32(_in_header, line,"PS", &_ps_array, &_num_ps);
	if(status==_num_sample)
	    hap_phaseset.push_back(_gt_array,_ps_array);
	else if(_num_sample<=0)
	    hap_phaseset.push_back(_gt_array,nullptr);
	else
	    die(("Invalid PS length: "+std::to_string(status)).c_str());
    }
    hap_transmission.phase();
    hap_transmission.align(hap_phaseset);
    
    int count=0;
    _rps_array = (int32_t *)realloc(_rps_array,_num_sample*sizeof(int32_t));
    //Finally, flush the buffer to the output file.
    while (!_line_buffer.empty()) 
    {
        bcf1_t *line = _line_buffer.front();
        _line_buffer.pop_front();
	hap_transmission.update_bcf1_genotypes(count,_gt_array,_ps_array,_rps_array);
	bcf_update_genotypes(_out_header, line, _gt_array, _num_gt);
	if(bcf_int32_count_missing(_ps_array,_num_sample)==_num_sample)
	    bcf_update_format_int32(_out_header, line, "PS", nullptr,0);
	else
	    bcf_update_format_int32(_out_header, line, "PS", _ps_array, _num_sample);
	if(bcf_int32_count_missing(_rps_array,_num_sample)==_num_sample)	
	    bcf_update_format_int32(_out_header, line, "RPS", nullptr,0);
	else
	    bcf_update_format_int32(_out_header, line, "RPS", _rps_array, _num_sample);

	if(!hap_transmission.is_mendel_consistent(count))
	    bcf_update_info_flag(_out_header, line, "MENDELCONFLICT", nullptr, 1);
	bcf_write(_out_file, _out_header, line);	
        bcf_destroy(line);
	count++;
    }
    return (0);
}


void PedPhaser::setup_io(args &a)
{
    //open a file.
    _bcf_reader = bcf_sr_init();

    if (a.targets != nullptr)
    {
        if (bcf_sr_set_targets(_bcf_reader, a.targets, a.targets_is_file, 0) < 0)
        {
            cerr << "ERROR: Failed to set targets " << a.targets << endl;
            exit(1);
        }
    }

    if (a.regions != nullptr)
    {
        if (bcf_sr_set_regions(_bcf_reader, a.regions, a.regions_is_file) < 0)
        {
            cerr << "ERROR: Failed to read the regions: " << a.regions << endl;
            exit(1);
        }
    }

    if (bcf_sr_add_reader(_bcf_reader, a.inputfile) != 1)
    {
        cerr << "ERROR: problem opening " << a.inputfile << endl;
        exit(1);
    }

    _in_header = _bcf_reader->readers[0].header;

    if (a.pedigree == nullptr)
    {
        _pedigree = new sampleInfo(_in_header);
    }
    else
    {
        _pedigree = new sampleInfo(a.pedigree, _in_header);
    }
    if (_pedigree->N <= 0)
    {
        die("no pedigree detected");
    }
    setup_output(a);
}

void PedPhaser::setup_output(args &a)
{
    _out_header = bcf_hdr_dup(_in_header);
    char output_type[] = "wv";
    output_type[1] = a.output_type;
    _out_file = hts_open(a.outfile, output_type);
    if (a.nthreads > 0)
    {
        bcf_sr_set_threads(_bcf_reader, a.nthreads);
        hts_set_threads(_out_file, a.nthreads);
    }

    bcf_hdr_remove(_out_header, BCF_HL_FMT, "PS"); //remove the old PS descripion
    bcf_hdr_append(_out_header, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Read-backed phase set. If missing from a phased genotype then it indicates the genotype was pedigree-phased such that children are phased as 'maternal allele | paternal allele' and parents are phased as 'allele transmitted to first child | untransmitted allele'\">");
    bcf_hdr_append(_out_header, "##FORMAT=<ID=RPS,Number=1,Type=Integer,Description=\"Read-backed phase set. The phase set (PS) value before this phased genotype was incorporated into the pedigree phase set\">");
    bcf_hdr_append(_out_header, "##INFO=<ID=MENDELCONFLICT,Number=0,Type=Flag,Description=\"This variant has at least one duo/trio with genotypes that are inconsistent with Mendelian inheritance\">");
    bcf_hdr_append(_out_header, ("##akt_pedphase_version=" + (string)AKT_VERSION).c_str());
    bcf_hdr_write(_out_file, _out_header);
}

PedPhaser::~PedPhaser()
{
    delete _pedigree;
    hts_close(_out_file);
    free(_ps_array);
    free(_gt_array);
    free(_gt_array_dup);
    bcf_sr_destroy(_bcf_reader);
    bcf_hdr_destroy(_out_header);
    if(_rps_array) free(_rps_array);
}

int pedphase_main(int argc, char **argv)
{
    int c;
    args arguments;
    arguments.output_type = 'v';
    if (argc < 3)
        usage();
    static struct option loptions[] = {
        {"out", 1, 0, 'o'},
        {"output-type", 1, 0, 'O'},
        {"pedigree", 1, 0, 'p'},
        {"threads", required_argument, nullptr, '@'},
        {"targets", required_argument, nullptr, 't'},
        {"targets-file", required_argument, nullptr, 'T'},
        {"regions-file", required_argument, nullptr, 'R'},
        {"regions", required_argument, nullptr, 'r'},
        {"exclude-chromosome", required_argument, nullptr, 'x'},
        {0, 0, 0, 0}};
    arguments.regions_is_file = false;
    arguments.targets_is_file = false;
    arguments.targets = arguments.pedigree = arguments.inputfile = arguments.include = arguments.regions = nullptr;
    arguments.outfile = "-";
    arguments.nthreads = 0;
    arguments.exclude_chromosomes = "";

    while ((c = getopt_long(argc, argv, "o:p:t:T:r:R:O:@:x:", loptions, nullptr)) >= 0)
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
        case 'x':
            arguments.exclude_chromosomes = optarg;
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
    if (arguments.inputfile == nullptr)
    {
        die("no input provided");
    }

    cerr << "Output file: " << arguments.outfile << endl;
    PedPhaser p(arguments);
    return (0);
}


