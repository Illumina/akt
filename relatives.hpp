#include "akt.hpp" 
#include "family.hpp"
#include "cluster.hpp"

using namespace std;

//read all ibd values < cutoff
void read_ibd1(ifstream &in, vector< vector<float> > &ibd, vector< vector<string> > &ln, set<string> &ls, float relmin )
{
	
	if(!in.is_open()){
		cout << "Failed to open file." << endl;
		exit(1);
	}
		
	string line = ""; 	
	while(getline(in,line)) // loop through the file
	{
		stringstream is(line);
		istream_iterator<string> begin(is);
		istream_iterator<string> end;
        vector<string> tokens(begin, end);
		
		vector<string> tmps(2); tmps[0] = tokens[0]; tmps[1] = tokens[1];
		ls.insert(tokens[0]);
		ls.insert(tokens[1]);

		vector<float> tmp(3); 
		tmp[0] = atof(tokens[2].c_str());
		tmp[1] = atof(tokens[3].c_str());
		tmp[2] = atof(tokens[5].c_str());
		if( tmp[2] > relmin ){
			ibd.push_back( tmp );
			ln.push_back(tmps);
		}
	}
	
}

int read_ibd2(ifstream &in, vector< vector<float> >  &ibd, vector< vector<string> > &ln )
{
	
	if(!in.is_open()){
		cout << "Failed to open file." << endl;
		exit(1);
	}
		
	string line = ""; 	
	set<string> unames;
	while(getline(in,line)) // loop through the file
	{
		stringstream is(line);
		istream_iterator<string> begin(is);
		istream_iterator<string> end;
        vector<string> tokens(begin, end);
		
		vector<string> tmps(2); tmps[0] = tokens[0]; tmps[1] = tokens[1];
		unames.insert(tokens[0]);
		unames.insert(tokens[1]);
		ln.push_back(tmps);
		vector<float> tmp(2); 
		tmp[0] = atof(tokens[2].c_str());
		tmp[1] = atof(tokens[3].c_str());
		ibd.push_back( tmp );
	}
	return (int)unames.size();
}

