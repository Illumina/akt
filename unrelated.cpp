#define __STDC_LIMIT_MACROS
#include "akt.hh"
#include "family.hh"
#include "cluster.hh"
#include "relatives.hh"

using namespace std;
using namespace Eigen;

/**
 * @name    usage
 * @brief   print out options
 *
 * List of input options
 *
 */
static void usage()
{
    cerr << "Print a list of unrelated individuals taking the output from akt kin as input." << endl;
    cerr << "Usage:" << endl;
    cerr << "./akt unrelated ibdfile" << endl;
    cerr << "\t -k --kmin:			threshold for relatedness (0.025)" << endl;
    cerr << "\t -i --its:			number of iterations to find unrelated (10)" << endl;
    exit(1);
}

///compare graphs based on number of members
struct less_than_graph
{
    inline bool operator() (const graph& struct1, const graph& struct2)
    {
        return (struct1.nv < struct2.nv);
    }
};


int unrelated_main(int argc, char* argv[])
{

    int c;

    if(argc<3) usage();
    static struct option loptions[] =    {
            {"kmin",1,0,'k'},
            {"its",1,0,'i'},
            {0,0,0,0}
    };
    float relmin = 0.025;
    int uits = 10;
    string prefix="out.";
    bool gout = false;

    while ((c = getopt_long(argc, argv, "k:i:?",loptions,NULL)) >= 0) {
        switch (c)
        {
            case 'k': relmin = atof(optarg); break;
            case 'i': uits = atoi(optarg); break;
            case '?': usage();
            default: cerr << "Unknown argument:"+(string)optarg+"\n" << endl; exit(1);
        }
    }
    optind++;
    string cfilename = argv[optind];
    cerr <<"Input: " << cfilename << endl;

    //read ibd data
    vector< vector<string> > pnames;	//sample pairs
    set<string> unames;					//sample names
    vector< vector<float> > ibd;		//ibd data
    ifstream in(cfilename.c_str());
    read_ibd1(in, ibd, pnames, unames, relmin);
    in.close();

    graph F;    //contains families only
    for (size_t i = 0; i < pnames.size(); i++)    //for all pairs
    {
        if (!F.hasvertex(pnames[i][0]))
        {
            F.add(pnames[i][0]);
        }
        if (!F.hasvertex(pnames[i][1]))
        {
            F.add(pnames[i][1]);
        }
        F.link(pnames[i][0], pnames[i][1], 2);
    }

    //deal with families
    vector<graph> DF;
    F.assign_disconnected(DF);
    sort(DF.begin(), DF.end(), less_than_graph());

    //FIND UNRELATED SET
    int uc = 0;
    //add singletons
    for (set<string>::iterator it=unames.begin(); it!=unames.end(); ++it)
    {
        if( !F.hasvertex( *it ) )
        {
            cout << *it << endl;
            ++uc;
        }
    }

//    for(size_t g=0; g<DF.size(); ++g)
//    {
//        vector<string> nm = DF[g].names();
//        vector<string> unrelated;
//        size_t ms = 0;
//        //Find a random unrelated set a few times and save the biggest one
//        //Do this for long enough, you'll find the best set...
//        for(int i=0; i<uits; ++i)
//        {
//            vector<string> ur;
//            DF[g].unrelated(ur);
//            if( ur.size() > ms )
//            {
//                ms = ur.size();
//                unrelated = ur;
//            }
//        }
//
//        for(size_t j=0; j<unrelated.size(); ++j)
//        {
//            cout << unrelated[j] << endl;
//            ++uc;
//        }
//    }

    for(size_t g=0; g<DF.size(); ++g)
    {
        vector<string> unrelated;
        DF[g].unrelatedGreedy(unrelated);
        for(size_t j=0; j<unrelated.size(); ++j)
        {
            cout << unrelated[j] << endl;
            ++uc;
        }
    }

    cerr << uc << " nominally unrelated samples." << endl;

    return 0;
}
