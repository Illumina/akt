#include "akt.hh" 
#include "family.hh"
#include "cluster.hh"

using namespace std;

/**
 * @name    read_ibd1
 * @brief   read all ibd values < cutoff
 *
 * @param [in] in  	read ibd values from here
 * @param [in] ibd	data container	- ibd values
 * @param [in] ln	data container	- sample pairs
 * @param [in] ls	data container	- unique sample names
 * @param [in] relmin	most distant relation to consider
 *
 */
void read_ibd1(ifstream &in, vector< vector<float> > &ibd, vector< vector<string> > &ln, set<string> &ls, float relmin );


/**
 * @name    read_ibd2
 * @brief   read all ibd values
 *
 * @param [in] in  	read ibd values from here
 * @param [in] ibd	data container	- ibd values
 * @param [in] ln	data container	- sample pairs
 * @return number of unique sample ids
 */
int read_ibd2(ifstream &in, vector< vector<float> >  &ibd, vector< vector<string> > &ln );