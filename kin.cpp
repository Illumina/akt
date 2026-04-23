/**
 * @file   kin.cpp
 * @Author Rudy Arthur (rudy.d.arthur@gmail.com)
 * @brief  Kinship calculator.
 *
 * A simple tool to read vcf/bcf files and calculate IBD and kinship.
 * Calculate allele frequencies from data or from input file.
 */

#include "kin.hh"

using namespace std;

/**
 * @name    read_pair
 * @brief   read input file containing sample pair names
 *
 * @param [in] in file 		handle of input list
 * @param [in] relpairs  	vector of pairs found
 * @param [in] name_to_id  	samples found in VCF
 *
 */
void read_pairs(ifstream &in, vector<pair<string, string>> &relpairs, map<string, int> &name_to_id)
{

	if (!in.is_open())
	{
		cout << "Failed to open file." << endl;
		exit(1);
	}

	string line = "";
	while (getline(in, line)) // loop through the file
	{
		stringstream is(line);
		istream_iterator<string> begin(is);
		istream_iterator<string> end;
		vector<string> tokens(begin, end);

		if (name_to_id.find(tokens[0]) == name_to_id.end())
		{
			cerr << tokens[0] << " not found in input" << endl;
			exit(1);
		}
		else if (name_to_id.find(tokens[1]) == name_to_id.end())
		{
			cerr << tokens[1] << " not found in input" << endl;
			exit(1);
		}

		pair<string, string> tmps = make_pair(tokens[0], tokens[1]);
		relpairs.push_back(tmps);
	}
}

/**
 * @name    make_pair_list
 * @brief   make a list of all possible sample pairs
 *
 * @param [in] relpairs  	vector of pairs
 * @param [in] names  		samples found in VCF
 *
 */
void make_pair_list(vector<pair<string, string>> &relpairs, vector<string> names)
{
	for (size_t j1 = 0; j1 < names.size(); ++j1)
	{
		for (size_t j2 = j1 + 1; j2 < names.size(); ++j2)
		{
			relpairs.push_back(make_pair(names[j1], names[j2]));
		}
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
	cerr << "\nAbout: Calculate kinship/IBD statistics from a multisample BCF/VCF" << endl;
	cerr << "Usage: akt kin [options] <in.bcf>" << endl;
	cerr << "Expects input.bcf to contain genotypes." << endl;

	cerr << "\nKinship calculation options:" << endl;
	cerr << "\t -k --minkin:			threshold for relatedness output (none)" << endl;
	cerr << "\t -F --freq-file:                a file containing population allele frequencies to use in kinship calculation" << endl;
	cerr << "\t -M --method:			type of estimator. 0:plink (default) 1:king-robust 2:genetic-relationship-matrix" << endl;
	umessage('a');
	umessage('@');
	cerr << "\nSite filtering options:" << endl;
	umessage('R');
	umessage('r');
	umessage('T');
	umessage('t');
	cerr << "\t    --force:			run kin without -R/-T/-F" << endl;
	cerr << "\nSample filtering options:" << endl;
	umessage('s');
	umessage('S');
	cerr << endl;
	exit(1);
}

void Kinship::addGenotypes(int *gt_arr, float p)
{
	float q = 1 - p;
	_n00 += 2 * p * p * q * q;

	_n11 += 2 * p * q;					 ///2ppq + 2qqp = 2pq(p+q) = 2pq == Hij
	_n10 += 4 * p * q * (p * p + q * q); ///4pppq + 4qqqp = 4pq(pp + qq)

	_n20 += p * p * p * p + q * q * q + 2 * p * q; ///  p^4 + q^4 + 4p^2q^2 = p^4 + q^4 + 2pq^2 + 2qp^2 = p^4 + q^4 + 2pq(p+q) = p^4 + q^4 + 2pq
	_n21 += p * p + q * q;									   ///ppp + qqq + ppq + pqq = pp(p+q) + qq(q+p) = pp + qq
	_n22 += 1;
	_af.push_back(p);

	///for each site record truth table
	///g=0  g=1  g=2  g=missing
	for (int i = 0; i < _nsample; ++i)
	{
		if (_bc == 0)
		{
			_bits[i].push_back(vector<bitset<BITSET_SIZE>>(4, bitset<BITSET_SIZE>()));
		}
		if (gt_arr[2 * i] != -1 && gt_arr[2 * i + 1] != -1)
		{
			int g = bcf_gt_allele(gt_arr[2 * i]) + bcf_gt_allele(gt_arr[2 * i + 1]);
			_bits[i].back()[g][_bc] = 1;
		}
		else
		{
			_bits[i].back()[3][_bc] = 1;
		}
	}
	///chunks of size L
	_bc = (_bc + 1) % (BITSET_SIZE);
	++_markers;
}

Kinship::Kinship(int nsample)
{
	// _lookup.resize(65536);
	// for(int i=0;i<_lookup.size();i++)
	// {
	// 	_lookup[i] = (float)((bitset<BITSET_SIZE>(i)).count());
	// }

	_nsample = nsample;
	_markers = 0;
	_bits.assign(_nsample, vector<vector<bitset<BITSET_SIZE>>>());
	_n00 = 0;
	_n10 = 0;
	_n11 = 0;
	_n20 = 0;
	_n21 = 0;
	_n22 = 0;
	_af.resize(50000); //shouldnt need more markers than this.
	_bc = 0;
}

void Kinship::estimateKinship(int j1, int j2, float &ibd0, float &ibd1, float &ibd2, float &ibd3, float &ks, int method)
{
	ibd0 = 0;
	ibd1 = 0;
	ibd2 = 0;
	ibd3 = 0;
	ks = -1;
	for (size_t i = 0; i < _bits[j1].size(); ++i)
	{
		// //opposite homozygotes. NAA,aa
		// ibd0 += _lookup[(_bits[j1][i][0] & _bits[j2][i][2]).to_ulong()] + _lookup[(_bits[j1][i][2] & _bits[j2][i][0]).to_ulong()];
		// //same genotype.
		// ibd2 += _lookup[(_bits[j1][i][0] & _bits[j2][i][0]).to_ulong()] + _lookup[(_bits[j1][i][1] & _bits[j2][i][1]).to_ulong()] + _lookup[(_bits[j1][i][2] & _bits[j2][i][2]).to_ulong()];
		// //missing in both.
		// ibd3 += _lookup[(_bits[j1][i][3] | _bits[j2][i][3]).to_ulong()];

		//opposite homozygotes. NAA,aa
		ibd0 += (_bits[j1][i][0] & _bits[j2][i][2]).count() + (_bits[j1][i][2] & _bits[j2][i][0]).count();
		//same genotype.
		ibd2 += (_bits[j1][i][0] & _bits[j2][i][0]).count() + (_bits[j1][i][1] & _bits[j2][i][1]).count() + (_bits[j1][i][2] & _bits[j2][i][2]).count();
		//missing in both.
		ibd3 += (_bits[j1][i][3] | _bits[j2][i][3]).count();
	}
	//consistent with IBD1 N - (NAA,AA + Naa,aa)
	ibd1 = _markers - ibd3 - ibd0 - ibd2;

	if (method == 0)
	{
		estimateIBD(ibd0, ibd1, ibd2, ibd3);
		ks = 0.5 * ibd2 + 0.25 * ibd1;
	}
	if (method == 1) //king
	{
		int Nhet_1 = 0, Nhet_2 = 0, Nhet_12 = 0;
		for (size_t i = 0; i < _bits[j1].size(); ++i)
		{
			auto mask = (_bits[j1][i][3] | _bits[j2][i][3]).flip();
			Nhet_1 += (mask & _bits[j1][i][1]).count();						//NAa^i
			Nhet_2 += (mask & _bits[j2][i][1]).count();						//NAa^j
			Nhet_12 += (_bits[j1][i][1] & _bits[j2][i][1]).count(); //NAa,Aa - no mask needed here
		}
		int minhet = min(Nhet_1, Nhet_2);
		ks = (Nhet_12 - 2 * ibd0) / (2 * minhet) + 0.5 - 0.25 * (Nhet_1 + Nhet_2) / minhet;
		estimateIBD(ibd0, ibd1, ibd2, ibd3);
	}
}

void Kinship::estimateIBD(float &ibd0, float &ibd1, float &ibd2, float &ibd3, bool normalise)
{
	///method of moments
	ibd0 /= _n00;
	ibd1 = (ibd1 - ibd0 * _n10) / _n11;
	ibd2 = (ibd2 - ibd0 * _n20 - ibd1 * _n21) / _n22;
	ibd3 = _n22 - ibd3;

	///_normalize i_n [0,1]
	if (normalise)
	{
		if (ibd0 > 1) //very unrelated, project to 100
		{
			ibd0 = 1;
			ibd1 = 0;
			ibd2 = 0;
		}
		if (ibd1 < 0)
		{
			ibd1 = 0;
		}
		if (ibd2 < 0)
		{
			ibd2 = 0;
		}

		float sum = ibd0 + ibd1 + ibd2;
		ibd0 /= sum;
		ibd1 /= sum;
		ibd2 /= sum;
	}
}

#define FORCE 100
int kin_main(int argc, char *argv[])
{

	int c;

	if (argc < 3)
		usage();
	static struct option loptions[] = {
		{"targets-file", 1, 0, 'T'},
		{"targets", 1, 0, 't'},
		{"regions-file", 1, 0, 'R'},
		{"regions", 1, 0, 'r'},
		{"method", 1, 0, 'M'},
		{"freq-file", 1, 0, 'F'},
		{"minkin", 1, 0, 'k'},
		{"threads", 1, 0, '@'},
		{"aftag", 1, 0, 'a'},
		{"samples", 1, 0, 's'},
		{"samples-file", 1, 0, 'S'},
		{"force", 0, 0, FORCE},
		{0, 0, 0, 0}};
	int method = 0;
	string regions = "";
	bool regions_is_file = false;
	string targets = "";
	bool targets_is_file = false;
	float min_kin = 0;
	bool tk = false;
	int thin = 1;
	int nthreads = -1;
	float min_freq = 0;
	string pairfile = "";
	string af_tag = "AF";
	sample_args sargs;

	bool force = false;
	bool used_r = false;
	bool used_R = false;
	bool used_t = false;
	bool used_T = false;

	string frq_file = "";
	while ((c = getopt_long(argc, argv, "T:t:R:r:M:F:k:h:@:m:a:s:S:f", loptions, NULL)) >= 0)
	{
		switch (c)
		{
		case 'R':
			regions = (optarg);
			used_R = true;
			regions_is_file = true;
			break;
		case 'r':
			regions = (optarg);
			used_r = true;
			break;
		case 'T':
			targets = (optarg);
			targets_is_file = true;
			used_T = true;
			break;
		case 't':
			targets = (optarg);
			used_t = true;
			break;
		case 'F':
			frq_file = optarg;
			break;
		case 'M':
			method = atoi(optarg);
			break;
		case 'k':
			tk = true;
			min_kin = atof(optarg);
			break;
		case '@':
			nthreads = atoi(optarg);
			break;
		case FORCE:
			force = true;
			break;
		case 'a':
			af_tag = string(optarg);
			break;
		case 's':
			sargs.sample_names = (optarg);
			sargs.subsample = true;
			break;
		case 'S':
			sargs.sample_names = (optarg);
			sargs.subsample = true;
			sargs.sample_is_file = 1;
			break;
		case '?':
			usage();
		default:
			cerr << "Unknown argument:" + (string)optarg + "\n"
				 << endl;
			exit(1);
		}
	}
	if (!force && targets.empty() && regions.empty() && frq_file.empty())
	{
		die("None of -R/-F/-T were provided.\n       kin does not require a dense set of markers and this can substantially increase compute time.\n       You can disable this error with --force");
	}

	if (method < 0 || method > 2)
	{
		cerr << "ERROR: method must be one of 0/1/2" << endl;
		exit(1);
	}
	if (used_r && used_R)
	{
		cerr << "-r and -R cannot be used simultaneously" << endl;
		exit(1);
	}
	if (used_t && used_T)
	{
		cerr << "-t and -T cannot be used simultaneously" << endl;
		exit(1);
	}
	assert(min_freq >= 0 && min_freq <= 1);
	if ((!targets.empty() && !regions.empty()))
	{
		cerr << "ERROR: -r/-R and -r/-R are incompatible" << endl;
		exit(1);
	}
	if (!frq_file.empty() && !regions.empty())
	{
		cerr << "ERROR: -F and -R/-r are incompatible!" << endl;
		exit(1);
	}
	if (!frq_file.empty() && regions.empty())
	{
		regions = frq_file;
		regions_is_file = true;
	}
	if (frq_file.empty())
	{
		cerr << "No frequency VCF provided (-F). Allele frequencies will be estimated from the data." << endl;
	}
	else
	{
		cerr << "Taking allele frequencies from " << frq_file << " using INFO/" << af_tag << endl;
	}
	if (nthreads < 1)
	{
		nthreads = 1;
	}
	if (!frq_file.empty() && method == 2)
	{
		die("method=2 and -F are incompatible. The GRM must estimate allele frequencies from the data.");
	}

	omp_set_num_threads(nthreads);
#pragma omp parallel
	{
		if (omp_get_thread_num() == 0)
		{
			if (omp_get_num_threads() != 1)
			{
				cerr << "Using " << omp_get_num_threads() << " threads" << endl;
				nthreads = omp_get_num_threads();
			}
		}
	}

	optind++;
	string filename = argv[optind]; ///input VCF

	int Nsamples;

	int sites = 0, num_sites = 0, num_study = 0;

	bcf_srs_t *sr = bcf_sr_init(); ///htslib synced reader.
	sr->collapse = COLLAPSE_NONE;  ///require matching ALTs
	sr->require_index = 1;		   ///require indexed VCF

	///subset regions
	if (!regions.empty())
	{
		if (bcf_sr_set_regions(sr, regions.c_str(), regions_is_file) < 0)
		{
			die("Failed to read the regions: " + regions);
		}
	}
	if (!targets.empty())
	{
		if (bcf_sr_set_targets(sr, targets.c_str(), targets_is_file, 0) < 0)
		{
			die("Failed to read the targets: " + targets);
		}
	}

	///open input VCF
	if (!(bcf_sr_add_reader(sr, filename.c_str())))
	{
		cerr << "Problem opening " << filename << endl;
		cerr << "Input file not found." << endl;
		bcf_sr_destroy(sr);
		return 0;
	}
	bcf_hdr_t *hdr = sr->readers[0].header;

	///Open file of allele freqs
	if (frq_file != "" && !(bcf_sr_add_reader(sr, frq_file.c_str())))
	{
		cerr << "Problem opening " << frq_file << endl;
		cerr << "Sites file not found." << endl;
		bcf_sr_destroy(sr);
		return 0;
	}
	///subsample input vcf
	if (sargs.subsample)
	{
		if (bcf_hdr_set_samples(hdr, sargs.sample_names, sargs.sample_is_file) != 0)
		{
			die("problem setting samples");
		}
	}
	if (bcf_hdr_nsamples(hdr) <= 0)
	{
		die("no samples!");
	}

	if (method == 2) //jump out to GRM routine.
	{
		die("method 2 (GRM) is deprecated. Try plink or GCTA.");
		//	return(grm(sr));
	}

	int N = bcf_hdr_nsamples(hdr); ///number of samples
	cerr << N << " samples" << endl;

	//     for(int i=0;i<N;i++)
	//     {
	// 	cout << i << " "<<hdr->samples[i]<<endl;
	//     }

	if (N < 50 && frq_file.empty())
	{
		cerr << "WARNING: your sample size is <50 and you have NOT provided population frequencies (-F)." << endl;
	}

	Nsamples = N;
	Kinship K(Nsamples);

	int count = 0;

	bcf1_t *line, *line2; ///bcf/vcf line structure.

	int *gt_arr = (int *)malloc(N * 2 * sizeof(int)), ngt = N * 2, ngt_arr = N * 2;
	float *af_ptr = (float *)malloc(1 * sizeof(float));
	int nval = 1;

	bool use_frq = !frq_file.empty();
	cerr << "Reading genotypes...";
	while (bcf_sr_next_line(sr)) ///read file
	{
		if (bcf_sr_has_line(sr, 0) && (!use_frq || bcf_sr_has_line(sr, 1))) ///present in the study file (and frequency file)
		{
			int nmiss = 0;
			int npres = 0;
			int sum = 0; ///AC

			line = bcf_sr_get_line(sr, 0);

			if (line->n_allele == 2 && (count++) % thin == 0) ///bi-allelic
			{
				ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
				assert(ngt == 2 * N);
				if (ngt < 0)
				{
					cerr << "Bad genotypes at " << line->pos + 1 << endl;
					exit(1);
				}
				for (int i = 0; i < 2 * N; i++) ///htslib -> int
				{
					if (bcf_gt_is_missing(gt_arr[i]) || bcf_gt_allele(gt_arr[i]) < 0 || bcf_gt_allele(gt_arr[i]) > 2)
					{
						gt_arr[i] = -1;
						++nmiss;
					}
					else
					{
						sum += bcf_gt_allele(gt_arr[i]);
						++npres;
					}
				}
				float p;
				if (frq_file.empty()) ///calculate AF from data
				{
					p = (float)sum / (float)(npres); ///allele frequency
				}
				else
				{
					assert(bcf_sr_has_line(sr, 1)); ///present in sites file.
					line2 = bcf_sr_get_line(sr, 1);
					num_sites++;
					++sites;
					int ret = bcf_get_info_float(sr->readers[1].header, line2, af_tag.c_str(), &af_ptr, &nval);
					if (ret < 0 || nval != 1)
					{
						cerr << af_tag << " read error at " << line2->rid << ":" << line->pos + 1 << endl;
						exit(1);
					}
					p = af_ptr[0];
				}
				if ((p < 0.5) ? (p > min_freq) : (1 - p > min_freq)) ///min af
				{
					K.addGenotypes(gt_arr, p);
				}
			} //thin
			++num_study;
		} //in study
	}	  //reader
	cerr << "done." << endl;
	free(gt_arr);
	free(af_ptr);

	if (frq_file.empty())
	{
		cerr << "Using " << K._markers << " markers for calculations" << endl;
	}
	else
	{
		cerr << "Kept " << K._markers << " markers out of " << sites << " in panel." << endl;
		cerr << num_study << "/" << num_sites << " of study markers were in the sites file" << endl;
	}
	cerr << "Calculating kinship values...";

//ordered lets us enforce ordereing but slows down code
//#pragma omp parallel for ordered
#pragma omp parallel for
	for (int j1 = 0; j1 < Nsamples; j1++)
	{
		for (int j2 = j1 + 1; j2 < Nsamples; j2++)
		{
			float ibd0, ibd1, ibd2, ibd3, ks;
			K.estimateKinship(j1, j2, ibd0, ibd1, ibd2, ibd3, ks, method);
			if (!tk || ks > min_kin)
			{
//#pragma omp ordered
#pragma omp critical
				{
					string id1 = hdr->samples[j1];
					string id2 = hdr->samples[j2];
					cout << id1 << " " << id2 << " " << left << " " << setprecision(5) << fixed << ibd0 << left << " " << setprecision(5) << fixed << ibd1 << left << " " << setprecision(5) << fixed << ibd2 << left << " " << setprecision(5) << fixed << ks << " " << setprecision(0) << ibd3 << "\n";
				}
			}
		}
	}

	bcf_sr_destroy(sr);
	cerr << "done." << endl;
	return 0;
}
