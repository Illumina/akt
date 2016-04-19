#define __STDC_LIMIT_MACROS

#include <Eigen/Dense>
#include "akt.hpp"
#include "family.hpp"
#include "cluster.hpp"
#include "relatives.hpp"

using namespace std;
using namespace Eigen;

/**
 * @name    usage
 * @brief   print out options
 *
 * List of input options
 *
 */
static void usage(){
	cerr << "Discover relatives from IBD" << endl;	
	cerr << "Usage:" << endl;
	cerr << "./akt relatives ibdfile" << endl;
	cerr << "\t -k --kmin:			threshold for relatedness (0.05)" << endl;
	cerr << "\t -i --its:			number of iterations to find unrelated (10)" << endl;
	cerr << "\t -g --graphout:			if present output pedigree graph files" << endl;
	cerr << "\t -p --prefix:			output file prefix (out)" << endl;
	cerr << "arrow types     : solid black	= parent-child" << endl;
	cerr << "                : dotted black	= siblings" << endl;
	cerr << "                : blue 		= second order" << endl;
	cerr << "                : red		= duplicates" << endl;
	cerr << "                : directed	= from parent to child" << endl;
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


int relatives_main(int argc, char* argv[])
{
	
	int c;
    
    if(argc<3) usage();
    static struct option loptions[] =    {
        {"kmin",1,0,'k'},
        {"its",1,0,'i'},
        {"prefix",1,0,'p'},
        {"graphout",1,0,'g'},
        {0,0,0,0}
    };
    float relmin = 0.05;
    int uits = 10;
    string prefix="out.";
	bool gout = false;
	
    while ((c = getopt_long(argc, argv, "k:i:p:g",loptions,NULL)) >= 0) {  
        switch (c)
        {
        case 'k': relmin = atof(optarg); break;
        case 'i': uits = atoi(optarg); break;
        case 'g': gout = true; break;
        case 'p': prefix = (optarg); prefix += "."; break;
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
	
	int K=6;
	int d=2;
	int N = ibd.size();
	cerr << N << " ibd pairs above threshold" << endl;
	
	MatrixXf mu(K,d);	//cluster centres, known from theory
	mu(0,0) = 0; 	mu(0,1) = 1; 	//PO
	mu(1,0) = 0.25; mu(1,1) = 0.5;  //FS
	mu(2,0) = 0.5; 	mu(2,1) = 0.5;	//rel1
	mu(3,0) = 0.75; mu(3,1) = 0.25;	//rel2
	mu(4,0) = 1; 	mu(4,1) = 0;	//UN
	mu(5,0) = 0; 	mu(5,1) = 0;	//Duplicates

	//add dummy variables so default clusters have at least one element
	for(int k=0; k<K; ++k){	
		vector<float> vec(d);
		for(int i=0; i<d; ++i){ vec[i] = mu(k,i); } 
		ibd.push_back(vec);
		++N;
	}
	
	//vector to Eigen
	MatrixXf P(N,d);	
	for (int i=0; i<N; i++){ P.row(i) = VectorXf::Map(&ibd[i][0],d); }

	//assign data points to category
	Cluster C(P,K);
	C.assignCentres(mu);
	C.clusterAssign(); //don't bother to cluster, just use given centres
	
	graph F;	//contains families only
	graph Fdup; //contains duplicates only

	int ct = 0; int rels = 0; int dups = 0;
	for(size_t i=0; i<pnames.size(); i++){	//for all pairs

		int type = C.assignment[ct];
		
		if( type == 0 || type == 1 || type == 2 || type == 5 ){ //parents, sibs, 2nd order or duplicates
			++rels;
			if( !F.hasvertex(pnames[i][0]) ){ F.add(pnames[i][0]); }
			if( !F.hasvertex(pnames[i][1]) ){ F.add(pnames[i][1]); }
			F.link( pnames[i][0], pnames[i][1], type );
		} 
		if( type == 5 ){	//duplicates
			++dups;
			if( !Fdup.hasvertex(pnames[i][0]) ){ Fdup.add(pnames[i][0]); }
			if( !Fdup.hasvertex(pnames[i][1]) ){ Fdup.add(pnames[i][1]); }
			Fdup.link( pnames[i][0], pnames[i][1], type );
		}
		
		++ct;
	}
	cerr << rels << " filtered ibd pairs >= 2nd order" << endl; 
	cerr << dups << " duplicate pairs" << endl;
	
	if( gout ){	//print out the big graph
		cerr <<  prefix + "allgraph includes all relative pairs 2nd order or lower" << endl; 
		ofstream f( (prefix + "allgraph").c_str() );
		F.gviz_neato(f);
		f.close();
	}

	//deal with duplicate samples
	vector<graph> DFdup;
	Fdup.assign_disconnected(DFdup);
	sort(DFdup.begin(), DFdup.end(), less_than_graph());
	//print duplicates
	for(size_t i=0; i<DFdup.size(); ++i){
		for(viter iter=DFdup[i].vlist.begin(); iter != DFdup[i].vlist.end(); ++iter){
			cout << "Dup" << i << "\t" << (*iter).second->name << endl; 
		}
	}

	//deal with families
	vector<graph> DF;
	F.assign_disconnected(DF);
	sort(DF.begin(), DF.end(), less_than_graph());
	
	map<string, string> fam_names;
	vector< string > fam_labs;
	//print families
	for(size_t i=0; i<DF.size(); ++i){
		fam_labs.push_back("Fam" + to_string(i));
		for(viter iter=DF[i].vlist.begin(); iter != DF[i].vlist.end(); ++iter){
			fam_names[ (*iter).second->name ] = fam_labs.back();
		}
	}

    //FIND UNRELATED SET
	int uc = 0;
	//add singletons
	for (set<string>::iterator it=unames.begin(); it!=unames.end(); ++it){ 
		if( !F.hasvertex( *it ) ){
			cout << "Unrel" << uc << "\t" << *it << endl; ++uc;
		}
	}

	for(size_t g=0; g<DF.size(); ++g){
		
		vector<string> nm = DF[g].names();
		vector<string> unrelated;

		size_t ms = 0;
		//Find a random unrelated set a few times and save the biggest one 
		//Do this for long enough, you'll find the best set...
		for(int i=0; i<uits; ++i){
			vector<string> ur; 
			DF[g].unrelated(ur);	
			if( ur.size() > ms ){
				ms = ur.size();
				unrelated = ur;
			}
		}
		
		for(size_t j=0; j<unrelated.size(); ++j){
			cout << "Unrel" << uc << "\t" << unrelated[j] << endl; ++uc;
		}			
	}
	cerr << uc << " nominally unrelated samples." << endl;	
	
	cerr << "Attempting to resolve pedigrees." << endl;
	//try to read every ibd pair
	vector< vector<float> > tibd; 
	vector< vector<string> > pnamesr;
	ifstream in2(cfilename.c_str());
	int Nsamples = read_ibd2(in2, tibd, pnamesr); 
	in2.close();
	//have to have all to all data
	cerr << Nsamples << " unique samples names" << endl;
	if( Nsamples*(Nsamples-1)/2 != (int)tibd.size() ){
		cerr << "Found " << tibd.size() << " total pairs when " << Nsamples*(Nsamples-1)/2 << " expected." << endl;
		cerr << "\"akt relatives\" expects unfiltered output from \"akt kin\"." << endl;
		exit(1);
	}
			
	vector< vector< vector<string> > > all_names( fam_labs.size() );
	vector< vector< vector<float> > > ibdr( fam_labs.size() );

	//strip out singletons and non family ibd pairs.
	int sz = 0;
	for(size_t n=0; n<pnamesr.size(); ++n){
		//fam_names[ pnames[n][0] ] = which family this is
		//then find the index of this family
		int pos1 = find( fam_labs.begin(), fam_labs.end(), fam_names[ pnamesr[n][0] ] ) - fam_labs.begin();
		int pos2 = find( fam_labs.begin(), fam_labs.end(), fam_names[ pnamesr[n][1] ] ) - fam_labs.begin();
		if( (size_t)pos1 < fam_labs.size() && (size_t)pos2 < fam_labs.size() ){
			if( pos1 == pos2 ){	
				all_names[ pos1 ].push_back( pnamesr[n] );
				ibdr[ pos1 ].push_back( tibd[n] );	
				++sz;
			}
		}
	}

	N = sz;
	cerr << N << " ibd pairs to cluster" << endl;
	
	//vector to Eigen
	MatrixXf Pr(N+K,d);
	ct = 0;	
	for(size_t i=0; i<ibdr.size(); ++i){
		for(size_t j=0; j<ibdr[i].size(); ++j){
			Pr.row(ct) = VectorXf::Map(&ibdr[i][j][0],d);
			++ct;
		}
	}
	//add dummy variables so default clusters have at least one element
	for(int i=0; i<K; ++i){Pr.row(ct) = mu.row(i); ++ct; }

	Cluster Cr(Pr,K);
	Cr.assignCentres(mu);
	Cr.clusterAssign(); //don't bother to cluster, just use given centres
	
	ct = 0;
	ofstream out_file2 ( (prefix + "fam").c_str() );
	
	//big loop over all families
	for(size_t n=0; n<all_names.size(); ++n){
			
		//temporary copy of family
		graph H;
		for(size_t m=0; m<all_names[n].size(); ++m){ 
			if( ! H.hasvertex(all_names[n][m][0]) ){ H.add(all_names[n][m][0]);}
			if( ! H.hasvertex(all_names[n][m][1]) ){ H.add(all_names[n][m][1]); }
		}
		N = all_names[n].size();
		map<string, int> relationship;

		//hash table of relationships and remove duplicates from graph
		for(int i=0; i<N; i++){
			relationship[all_names[n][i][0] + all_names[n][i][1]] = Cr.assignment[ct];
			relationship[all_names[n][i][1] + all_names[n][i][0]] = Cr.assignment[ct];

			if( Cr.assignment[ct] == 5 ){
				if( H.hasvertex(all_names[n][i][1]) ){ H.remove_vertex( all_names[n][i][1] ); }	
			} 
			++ct;
		}
		
		//Add PO links to H
		for(int i=0; i<N; i++){
			if( relationship[all_names[n][i][0] + all_names[n][i][1]] == 0 &&
				H.hasvertex(all_names[n][i][0]) && H.hasvertex(all_names[n][i][1]) &&
				!H.linked(all_names[n][i][0], all_names[n][i][1]) &&
				!H.descendant(all_names[n][i][0], all_names[n][i][1]) ){
					H.link( all_names[n][i][0], all_names[n][i][1], -1 ); //1 is parent, 1 is child, don't know which is which -1
			}
		}
		
		//trio resolution
		for(viter iter=H.vlist.begin(); iter != H.vlist.end(); ++iter){
			if( (*iter).second->num_in + (*iter).second->num_out == 2){	//2 links attached to node

				string tname = (*iter).second->name;
				vector<string> tmp;	//names of nodes originating links  tmp[0] --- iter --- tmp[1]
				
				for(int i=0; i<(*iter).second->num_in; ++i){ 
					if( (*iter).second->in[i].type == -1 ){ 				
						tmp.push_back( (*iter).second->in[i].from_name ); 
					}
				}
				for(int i=0; i<(*iter).second->num_out; ++i){  		
					if( (*iter).second->out[i].type == -1 ){ 				
						tmp.push_back( (*iter).second->out[i].to_name ); 
					}
				}
				
				for(size_t i=0; i<tmp.size(); ++i){
					for(size_t j=i+1; j<tmp.size(); ++j){	
						
						if( relationship[ tmp[i] + tmp[j] ] == 1 ){ //sibs -> this is the parent of 2 sibs

							if(H.linked(tname, tmp[i])){ H.unlink(tname, tmp[i]); }
							if(H.linked(tmp[i], tname)){ H.unlink(tmp[i], tname); }
							if(H.linked(tname, tmp[j])){ H.unlink(tname, tmp[j]); }
							if(H.linked(tmp[j], tname)){ H.unlink(tmp[j], tname); }
							
								H.link( tname, tmp[i], 0 );
								H.link( tname, tmp[j], 0 );	
						} 
						//unrelated => this is the child of 2 parents
						if( relationship[ tmp[i] + tmp[j] ] == 4 || relationship[ tmp[i] + tmp[j] ] == 3 ){ 
							if(H.linked(tname, tmp[i])){ H.unlink(tname, tmp[i]); }
							if(H.linked(tmp[i], tname)){ H.unlink(tmp[i], tname); }
							if(H.linked(tname, tmp[j])){ H.unlink(tname, tmp[j]); }
							if(H.linked(tmp[j], tname)){ H.unlink(tmp[j], tname); }
								H.link( tmp[i], tname, 0 );
								H.link( tmp[j], tname, 0 );	
						}
						
					}
				}

			}
		} 
		
		//multiple link resolution
		for(viter iter=H.vlist.begin(); iter != H.vlist.end(); ++iter){
			if( (*iter).second->num_in + (*iter).second->num_out > 2){	//more than in 2 links attached to node

				string tname = (*iter).second->name;
				vector<string> tmp;	//names of nodes originating links
				
				for(int i=0; i<(*iter).second->num_in; ++i){ 
					if( (*iter).second->in[i].type == -1 ){ 				
						tmp.push_back( (*iter).second->in[i].from_name ); 
					}	
				}
				for(int i=0; i<(*iter).second->num_out; ++i){  		
					if( (*iter).second->out[i].type == -1 ){ 				
						tmp.push_back( (*iter).second->out[i].to_name ); 
					}
				}

				vector<size_t> parent;
				set<string> sib;
				int sc = 0;
				for(size_t i=0; i<tmp.size(); ++i){
					for(size_t j=i+1; j<tmp.size(); ++j){
						++sc;
						if( relationship[ tmp[i] + tmp[j] ] == 4 ){ //unrelated => parents, 
																//not 3 because 2 misclassified as 3 can happen
							parent.push_back(i);
							parent.push_back(j);
						}
						if( relationship[ tmp[i] + tmp[j] ] == 0 ){
							sib.insert( tmp[i] );
							sib.insert( tmp[j] );
						} 
						
					}
				}
				if(parent.size() == 2){
					for(size_t i=0; i<tmp.size(); ++i){
						if( i != parent[0] && i != parent[1] ){
														
							if(H.linked(tname, tmp[i])){ H.unlink(tname, tmp[i]); }
							if(H.linked(tmp[i], tname)){ H.unlink(tmp[i], tname); }
							
							H.link( (*iter).second->name, tmp[i], 0 );
						} else {
							
							if(H.linked(tname, tmp[i])){ H.unlink(tname, tmp[i]); }
							if(H.linked(tmp[i], tname)){ H.unlink(tmp[i], tname); }
							
							H.link( tmp[i], (*iter).second->name, 0 );
						}
					}
				} else if (parent.size() == 0 && (int)sib.size() == 2*sc){ //all links are siblings => vertex is parent
					for (set<string>::iterator it=sib.begin(); it!=sib.end(); ++it){
			
						if(H.linked(tname, *it)){ H.unlink(tname, *it); }
						if(H.linked(*it, tname)){ H.unlink(*it, tname); }
						
						relationship[ tname + *it ] = 1; //update relationship map;
						relationship[ *it + tname ] = 1;
						
						H.link( tname, *it, 0 );
						
					}
				} else { //one grandparent and multiple grandchildren and all other cases
					if( parent.size() > 2 ){ 
						cout << "Warning: found more than 2 unrelated 'parents' of node " << (*iter).second->name << endl;
					} 
				}
			}
		}

		//add in other relationships for graph output
		for(int i=0; i<N; i++){
			if( relationship[all_names[n][i][0] + all_names[n][i][1]] < 4 && 
			    relationship[all_names[n][i][0] + all_names[n][i][1]] > 0 &&
				H.hasvertex(all_names[n][i][0]) && H.hasvertex(all_names[n][i][1]) &&
				!H.linked(all_names[n][i][0], all_names[n][i][1]) &&
				!H.linked(all_names[n][i][1], all_names[n][i][0]) ){
					H.link( all_names[n][i][0], all_names[n][i][1], relationship[all_names[n][i][0] + all_names[n][i][1]] );
			}
		}
		
		//graph file
		if( gout ){
			ofstream out_file ( (prefix + fam_labs[n] + ".graph").c_str() );
			H.gviz_dot(out_file);
			out_file.close();
		}
		//family info
		for(viter iter=H.vlist.begin(); iter != H.vlist.end(); ++iter){
			cout << fam_labs[n] << "\t" << (*iter).second->name << endl;
		}
		int ntypes = 0;
		for(viter iter=H.vlist.begin(); iter != H.vlist.end(); ++iter){
			for(int i=0; i<(*iter).second->num_out; ++i){
				if( (*iter).second->out[i].type >= -1 && (*iter).second->out[i].type <= 3){
					cout << "Type\t" << fam_labs[n] << "\t" << (*iter).second->name << 
					"\t" << (*iter).second->out[i].to_name << "\t"; 
					if( (*iter).second->out[i].type == -1){ cout << "Parent/Child" << endl; }
					if( (*iter).second->out[i].type == 0){ cout << "Parent/Child" << endl; }
					if( (*iter).second->out[i].type == 1){ cout << "Sibling" << endl; }
					if( (*iter).second->out[i].type == 2){ cout << "Second-order" << endl; }
					if( (*iter).second->out[i].type == 3){ cout << "Higher-order" << endl; }
					++ntypes;
				}
			}	
		}
		//output fam file
		H.ped_print(out_file2, fam_labs[n]);
	}
	out_file2.close();


	return 0;
}





