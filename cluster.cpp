#include <Eigen/Dense>
#include "akt.hpp"
#include "cluster.hpp"
#include "reader.hpp"

using namespace std;
using namespace Eigen;

bool cf_idx(vector<int> i, vector<int> j) { return i[0] > j[0]; } //descending order
bool cf_pair(pair<int, float> i, pair<int, float> j) { return i.first > j.first; } //descending order

//number of points less than dc away
void Cluster::localDensity(float dc_)
{
    dc = dc_;
    rho = VectorXi::Zero(N);
    float dc2 = dc * dc;
    for (int i = 0; i < N; ++i)
    {
        for (int j = i + 1; j < N; ++j)
        {
            if ((P.row(i) - P.row(j)).squaredNorm() < dc2)
            {
                ++rho[i];
                ++rho[j];
            }
        }
    }

}

//number of points less than dc away
void Cluster::minDistance()
{
    //can probably be improved by sorting
    delta = VectorXf::Ones(N) * numeric_limits<float>::infinity();    //distance of closest data point of higher density
    VectorXf maxd = VectorXf::Zero(N);    //furthest data point of lower density

    for (int i = 0; i < N; ++i)
    {
        for (int j = i + 1; j < N; ++j)
        {
            float dist = (P.row(i) - P.row(j)).squaredNorm();

            //density of i > density of j
            if (rho[i] > rho[j])
            {
                if (dist < delta[j]) { delta[j] = dist; }    //check if dist of i from j is smaller then best so far
            }
            else if (dist >  maxd[j])     //density j >= density i: check if distance from i to j is larger than max so far
            {
                maxd[j] = dist;
            }

            if (rho[j] > rho[i])
            {
                if (dist < delta[i])
                {
                    delta[i] = dist;
                }
            }
            else if (dist > maxd[i])
            {
                maxd[i] = dist;
            }
        }
    }

    for (int i = 0; i < N; ++i) {
        if (delta[i] == numeric_limits<float>::infinity()) { //there were no higher density data points
            delta[i] = maxd[i]; //assign highest density data point delta = furthest lower density point from i
        }
    }

}

//print out delta versus density graph - useful to determine parameter settings
void Cluster::densityPlot() {
    for (int i = 0; i < N; ++i) { cout << rho(i) << " " << delta(i) << "\n"; }

}

//assign clusters based on density
void Cluster::densityCluster(int min_rho, float min_delta) {

    vector<int> peak_rho;
    vector<float> peak_delta;
    vector<int> peak;
    vector<vector<int> > sorted_rho;

    //find density values at peaks
    for (int i = 0; i < N; ++i) {
        if (rho(i) > min_rho) {    //only data points at least this dense
            vector<int> tmp(2);
            tmp[0] = rho(i);
            tmp[1] = i;
            sorted_rho.push_back(tmp);
            if (delta(i) > min_delta) { //only peaks in local density at least this high
                vector<int>::iterator it = find(peak_rho.begin(), peak_rho.end(), rho(i));
                if (peak_rho.size() == 0 || it == peak_rho.end()) {
                    peak_rho.push_back(rho(i));
                    peak_delta.push_back(delta(i));
                    peak.push_back(i);
                } else {    //repeated rho value
                    int idx = distance(peak_rho.begin(), it);
                    if ((P.row(i) - P.row(peak[idx])).squaredNorm() > dc * dc) { //further than dc away - different peak
                        peak_rho.push_back(rho(i));
                        peak_delta.push_back(delta(i));
                        peak.push_back(i);
                    } else if (peak_delta[idx] < delta(i)) {    //closer than dc - overlapping peak
                        peak_delta[idx] = delta(i);
                        peak[idx] = i;
                    }
                }
            }
        }
    }
    //assign peaks arbitrary cluster numbers
    K = peak.size() + 1;
    centres.resize(K, d);
    sizes = VectorXi::Zero(K);
    for (size_t i = 0; i < peak.size(); ++i) {
        assignment[peak[i]] = i + 1;
        centres.row(i + 1) = P.row(peak[i]);
    }

    sort(sorted_rho.data(), sorted_rho.data() + sorted_rho.size(), cf_idx); //sort to avoid backtracking

    //assign peaks
    for (size_t i = 0; i < sorted_rho.size(); ++i) {

        vector<int>::iterator it = find(peak.begin(), peak.end(), sorted_rho[i][1]);
        if (it != peak.end()) { //found a peak skip it
            //0 cluster is unassigned
        } else { //not a peak
            //find nn of higher density
            float min_dist = numeric_limits<float>::infinity();
            for (size_t j = 0; j < sorted_rho.size(); ++j) {
                float dist = (P.row(sorted_rho[i][1]) - P.row(sorted_rho[j][1])).squaredNorm();
                if (assignment[sorted_rho[j][1]] != 0 && i != j && sorted_rho[j][0] >= sorted_rho[i][0] &&
                    dist < min_dist) {
                    min_dist = dist;
                    assignment[sorted_rho[i][1]] = assignment[sorted_rho[j][1]];
                }
            }
        }
        ++sizes[assignment[sorted_rho[i][1]]];
    }

}

//clusters to stdout, with labels appended
void Cluster::clustered_data_dump(vector<vector<string> > &labels) {

    for (int k = 0; k < K; k++) {
        for (int i = 0; i < N; i++) {
            if (assignment[i] == k) {
                for (int j = 0; j < d - 1; ++j) {
                    cout << P(i, j) << "\t";
                }
                cout << P(i, d - 1);
                if (silset) { cout << "\t" << sil(i); }
                cout << "\tCluster" << k << "\t";
                if (labels[i].size() > 0) {
                    for (size_t j = 0; j < labels[i].size() - 1; ++j) {
                        cout << labels[i][j] << "\t";
                    }
                    cout << labels[i].back() << "\n";
                } else {
                    cout << "\n";
                }
            }
        }
        if (k < K - 1) cout << endl;
    }

}

//assign points to clusters based on centres
void Cluster::clusterAssign()
{
    sizes = VectorXi::Zero(K);

    //assign to cluster
    for (int i = 0; i < N; ++i)
    {
        float min_dist = numeric_limits<float>::infinity();
        int c = -1;
        for (int j = 0; j < K; ++j)
        {
            float dist = (P.row(i) - centres.row(j)).squaredNorm();
            if (dist < min_dist)
            {
                min_dist = dist;
                c = j;
            }
        }
        assignment(i) = c;
        ++sizes(c);
    }
}

//k++ means initialisation
void Cluster::initialiseCentres() {    //k++ means init step

    centres.row(0) = P.row(int(N * drand48()));

    int ac = 1;
    while (ac < K) {
        vector<float> D(N); //distance of each point from the nearest centre
        float dnorm = 0.0;
        //fill in D
        for (int i = 0; i < N; ++i) {
            float min_dist = (P.row(i) - centres.row(0)).squaredNorm();
            for (int j = 1; j < ac; ++j) {
                float dist = (P.row(i) - centres.row(j)).squaredNorm();
                if (dist < min_dist) { min_dist = dist; }
            }
            D[i] = min_dist;
            dnorm += min_dist;
        }
        //choose a new centre
        float rand = drand48();    //good enough random number generator
        float sum = 0;
        for (int i = 0; i < N; ++i) {
            sum += D[i] / dnorm;
            if (rand <= sum) {
                centres.row(ac++) = P.row(i);
                break;
            }
        }
        if (rand >= sum) centres.row(ac++) = P.row(N - 1);
    }
}

//k means clustering
void Cluster::kMeans(int max_it) {

    VectorXi old_sizes(K);

    //initial cluster assignments
    clusterAssign();

    //If number of cluster members doesn't change exit.
    int it;
    for (it = 0; it < max_it; ++it) {
        //set cluster location to mean
        centres = MatrixXf::Zero(K, d);
        for (int i = 0; i < N; ++i) { centres.row(assignment(i)) += P.row(i); }
        for (int c = 0; c < K; ++c) { centres.row(c) /= sizes(c); }
        //reassign cluster
        old_sizes = sizes;
        clusterAssign();
        //check counts
        int diff = 0;
        for (int i = 0; i < K; ++i) {
            diff += abs(sizes(i) - old_sizes(i));
        }
        if (diff == 0) {
            break;
        }
    }
    cerr << "Clustering completed after " << it << " iterations" << endl;
}

//k++ means algorithm
void Cluster::kppMeans(int max_it) {

    initialiseCentres();
    kMeans(max_it);

}

// calculate weights for gaussian mixture clusters
float Cluster::EMweights() {

    //calc weights
    for (int k = 0; k < K; ++k) {
        MatrixXf MS = P.transpose();
        for (int n = 0; n < N; ++n) { MS.col(n) -= centres.row(k); }    //MS.col(i) = x_i - mu_k
        MatrixXf SMS = covariance[k].householderQr().solve(MS);    //SMS.col(i) = Sigma^-1 ( x_i - mu_k )
        float norm = 1.0 / (pow(2 * M_PI, d * 0.5) * sqrt(abs(covariance[k].determinant()))); //technically correct
        for (int n = 0; n < N; ++n) { prob(k, n) = norm * exp(-0.5 * MS.col(n).dot(SMS.col(n))); }
    }

    float newP = 0;
    for (int n = 0; n < N; ++n) {
        for (int k = 0; k < K; ++k) { weight(k, n) = alpha(k) * prob(k, n); }
        float norm = weight.col(n).sum();
        newP += log(norm);
        weight.col(n) /= norm;
    }

    return newP;
}

//assign data to clusters based on gaussian probability weights
void Cluster::EMassign() {

    for (int k = 0; k < K; ++k) { sizes(k) = 0; }
    for (int n = 0; n < N; ++n) {
        int mx;
        weight.col(n).maxCoeff(&mx);
        assignment[n] = mx;
        ++sizes[mx];
    }
}

//EM Gaussian clustering algorithm
void Cluster::EMcluster(int max_its) {

    clusterAssign();

    //init
    weight = MatrixXf::Zero(K, N);
    for (int n = 0; n < N; ++n) { weight(assignment(n), n) = 1; }

    float oldP = 0;
    for (int n = 0; n < N; ++n) {
        float norm = weight.col(n).sum();
        oldP += log(norm);
    }

    for (int it = 0; it < max_its; ++it) {

        //calc params
        centres = MatrixXf::Zero(K, d);
        for (int k = 0; k < K; ++k) {
            alpha(k) = 0;
            for (int n = 0; n < N; ++n) {
                alpha(k) += weight(k, n);
                centres.row(k) += weight(k, n) * P.row(n).transpose();
            }
            centres.row(k) *= (1.0 / alpha(k));
            covariance[k] = MatrixXf::Zero(d, d);
            for (int n = 0; n < N; ++n) {
                covariance[k].noalias() +=
                        weight(k, n) * (P.row(n) - centres.row(k)).transpose() * (P.row(n) - centres.row(k));
            }
            covariance[k] *= (1.0 / (float) alpha(k));
            alpha(k) *= (1.0 / N);
        }

        float newP = EMweights();
        if (abs(newP - oldP) / abs(oldP) < eps) {
            cerr << "Clustering completed after " << it << " iterations" << endl;
            break;
        }
        if (it + 1 == max_its) { cerr << "Max its reached: " << abs(newP - oldP) / abs(oldP) << endl; }
        oldP = newP;
    }

    EMassign();
}

//calculate silhouette scores
void Cluster::silhouette() {

    float a; //dissimilarity to assigned cluster
    float b; //dissimilarity to nearest cluster not assigned
    float s; //silhouette

    sil = VectorXf::Zero(N);
    for (int i = 0; i < N; ++i) {
        VectorXf dists = VectorXf::Zero(K);
        for (int j = 0; j < N; ++j) { dists(assignment(j)) += (P.row(i) - P.row(j)).squaredNorm(); }
        for (int k = 0; k < K; ++k) { dists(k) /= sizes(k); }
        a = dists(assignment(i));
        b = 1.0 / 0.0;
        for (int k = 0; k < K; ++k) {
            if (assignment(i) != k) {
                if (dists(k) < b) { b = dists(k); }
            }
        }
        s = (b - a) / max(a, b);
        sil(i) = s;
    }
    silset = true;
}


/**
 * @name    usage
 * @brief   print out options
 *
 * List of input options
 *
 */
static void usage() {
    cerr << "Clustering on text files" << endl;
    cerr << "Usage:   ./akt cluster input.txt" << endl;
    cerr << "\t -k --K:			number of clusters" << endl;
    cerr << "\t -i --seed:			random seed for starting values" << endl;
    cerr << "\t -a --alg:			clustering algorithm 0 = k++means, 1 = gaussian mixture, 2 = density method"
         << endl;
    umessage('c');
    cerr << "\t -C --cfile:			initial guess for cluster centres in text file" << endl;
    cerr << "\t -o --outputcfile:		assigned cluster centres" << endl;
    cerr << "\t -c --cols:			which columns to use e.g. 2-4" << endl;
    cerr << "\t -I --maxits:			max number of iterations to use: a=0,1" << endl;
    cerr << "\t -d --dc:			radius for density method: a=2" << endl;
    cerr << "\t -p --rho_min:			min density for cluster centre: a=2" << endl;
    cerr << "\t -D --delta_min:		min radius for cluster centre: a=2" << endl;
    cerr << "\t --density-plot:		plot the density and finish: a=2" << endl;
    cerr << "\t -e --silhouette:		calculate silhouette score" << endl;
    exit(1);
}


int cluster_main(int argc, char **argv) {

    int c;

    if (argc < 3) usage();
    static struct option loptions[] = {
            {"K",            1, 0, 'k'},
            {"alg",          1, 0, 'a'},
            {"cfile",        1, 0, 'C'},
            {"output_cfile", 1, 0, 'o'},
            {"cols",         1, 0, 'c'},
            {"maxits",       1, 0, 'I'},
            {"dc",           1, 0, 'd'},
            {"rho_min",      1, 0, 'p'},
            {"delta_min",    1, 0, 'D'},
            {"density-plot", 0, 0, 1},
            {"seed",         1, 0, 'i'},
            {"silhouette",   1, 0, 'e'},
            {0,              0, 0, 0}
    };
    int seed = 12345;
    int K = 2;
    int alg = 0;
    bool use_file = false;
    string cfile = "";
    string dims = "";
    float rho_min = 1;
    float delta_min = 0.01;
    bool density_plot = false;
    float dc = 0.05;
    int max_its = 100;
    //string prefix = "";
    bool dosil = false;
    string output_c = "";

    while ((c = getopt_long(argc, argv, "k:a:c:C:I:d:p:D:i:eo:", loptions, NULL)) >= 0) {
        switch (c) {
            case 'k':
                K = atoi(optarg);
                break;
            case 'a':
                alg = atoi(optarg);
                break;
            case 'I':
                max_its = atoi(optarg);
                break;
            case 'd':
                dc = atof(optarg);
                break;
            case 'p':
                rho_min = atof(optarg);
                break;
            case 'D':
                delta_min = atof(optarg);
                break;
            case 1:
                density_plot = true;
                break;
            case 'i':
                seed = atoi(optarg);
                break;
            case 'e':
                dosil = true;
                break;
            case 'C':
                use_file = true;
                cfile = (optarg);
                break;
            case 'c':
                dims = (optarg);
                break;
            case 'o':
                output_c = (optarg);
                break;
            case '?':
                usage();
            default:
                die("Unknown argument:" + (string) optarg + "\n");
        }
    }
    if (dims == "") {
        cerr << "-c argument is mandatory" << endl;
        exit(1);
    }
    switch (alg) {
        case 0:
            cerr << "K++ means" << endl;
            break;
        case 1:
            cerr << "EM Gaussian" << endl;
            break;
        case 2:
            cerr << "Density Clustering" << endl;
            break;
    }
    optind++;
    string input = argv[optind];
    cerr << "Input: " << input << endl;

    srand48(seed);

    vector<vector<float> > data;
    vector<vector<string> > labels;
    ifstream in_file(input.c_str());
    readMatrix(in_file, data, labels, dims);
    in_file.close();

    int N = data.size();
    int d = data[0].size();

    //vector to Eigen
    MatrixXf P(N, d);
    for (int i = 0; i < N; i++) { P.row(i) = VectorXf::Map(&data[i][0], d); }

    //initial centres
    MatrixXf mu;
    if (use_file) {
        ifstream in_cfile(cfile.c_str());
        vector<vector<float> > mu_data;
        vector<vector<string> > mu_lab;
        readMatrix(in_cfile, mu_data, mu_lab, "1-" + to_string(d));
        in_cfile.close();
        if (mu_data[0].size() != (size_t) d) {
            cerr << "Center init must be same dim as data" << endl;
            exit(1);
        }
        K = mu_data.size();
        mu.resize(K, d);
        for (int i = 0; i < K; i++) { mu.row(i) = VectorXf::Map(&mu_data[i][0], d); }
    }
    Cluster C(P, K);
    if (use_file) {
        C.assignCentres(mu);
    } else {
        C.initialiseCentres();
    }

    //which algorithm?
    switch (alg) {
        case 0:
            C.kMeans(max_its);
            break;
        case 1:
            C.EMcluster(max_its);
            break;
        case 2:
            C.localDensity(dc);
            C.minDistance();
            if (density_plot) {
                C.densityPlot();
                return 0;
            } else {
                C.densityCluster(rho_min, delta_min);
                cerr << "found " << C.K << " clusters" << endl;
            }
    }
    //calculate silhouette?
    if (dosil) { C.silhouette(); }
    //output
    C.clustered_data_dump(labels);

    //output centres
    if (output_c != "") {
        ofstream fk(output_c.c_str());
        int st = (alg == 2) ? 1 : 0;
        for (int k = st; k < C.K; ++k) {
            fk << C.centres.row(k) << endl;
        }
        fk.close();
    }

    return (0);
}

