#ifndef PEDIGREE_H
#define PEDIGREE_H

#include <iostream>
#include <iomanip>
#include <time.h>       
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <list>
#include <deque>
#include <iomanip>
#include <stdint.h>
#include <stdlib.h> 
#include <algorithm>
#include <sstream>
#include <ctype.h>
#include <limits> 
#include <map>
#include <algorithm>
#include <set>
#include <queue>
#include <stack>
#include "logs.hpp"

using namespace std;

//(from_name, from_id) ----type----> (to_name, to_id)
class edge{
public:
	string from_name;	
	int from_id;		
	string to_name;
	int to_id;
	//0 = parent child
	//1 = sibling
	//2 = 2nd order relative 
	//3 = 3rd order relative
	//4 = unrelated
	//5 = sample duplicate
	int type;	
	
	
	edge(){ from_name = ""; from_id = -1; to_name=""; to_id = -1;  type = 0; }
	edge(string from_name_, int from_id_, string to_name_, int to_id_, int type_){ 
		from_name = from_name_; to_name = to_name_; 
		from_id = from_id_;	to_id = to_id_;
		type = type_;
	}
};


/*
in0->-[name, id]->-out0
in1->-			->-out1
*/
class vertex{	
public:
	string name;
	int id;
	int num_in;
	int num_out;
	
	vector<edge> in;	//incoming edges
	vector<edge> out;	//outgoing edges

	int dist;			//distance marker for graph traversal

	vertex(string name_, int id_){
			name = name_;
			id = id_;
			num_in = 0;
			num_out = 0;
	}

};

//return a pointer to a copy of the vertex pointed to by in
//in -> v, out -> w, w=v
vertex* copy_vertex(vertex *in){
	
	vertex *out = new vertex(in->name, in->id);
	out->num_in = in->num_in;
	out->num_out = in->num_out;
	out->dist = in->dist;
	out->in.resize(0);
	//copy edges
	for(size_t i=0; i<in->in.size(); ++i){
			edge tmp(in->in[i].from_name, in->in[i].from_id, in->in[i].to_name, in->in[i].to_id, in->in[i].type );
			out->in.push_back(tmp);
	}
	out->out.resize(0);
	for(size_t i=0; i<in->out.size(); ++i){
			edge tmp(in->out[i].from_name, in->out[i].from_id, in->out[i].to_name, in->out[i].to_id, in->out[i].type );
			out->out.push_back(tmp);
	}	
	return out;
}

typedef map<string, vertex*>::iterator viter;
typedef map<string, vertex*>::const_iterator const_viter;

//a hash table of vertex pointers
class graph
{
	public:
		//default
		graph(){
			nv=0;
			ne=0;
		}
		//constructor
		graph(string name_){ 
			vertex *v = new vertex(name_, 0);
			vlist[name_] = v;
			nv = 1; 
			ne = 0; 
		}
		//destruction of graph
		void clear(){
			for( viter iter= vlist.begin(); iter!=vlist.end(); iter++){	
				if((*iter).second){
					delete iter->second;
				}				
			}
			vlist.clear();
		}
		//copy constructor
		graph( const graph& other ){
			nv = other.nv;
			ne = other.ne;
			
			clear();
			for( const_viter iter= other.vlist.begin(); iter!=other.vlist.end(); ++iter){	
				vertex* v = copy_vertex( (*iter).second );
				vlist[v->name] = v;
			}
			
		}
		//destructor
		~graph(){ clear(); }
		//copy operator
		//no copy swap because constructor doesn't initialise vlist
		graph& operator=(const graph& second) 
		{
			this->ne = second.ne;
			this->nv = second.nv;

			clear();
			for( const_viter iter= second.vlist.begin(); iter!=second.vlist.end(); ++iter){	
				this->vlist[ (*iter).second->name ] = copy_vertex( (*iter).second );
			}
			
			return *this;
		}

		//output functions
		void gviz_dot(ofstream& of);
		void ped_print(ofstream& of, string fam);
		void gviz_neato(ofstream& of);
		void gviz_neato_named(ofstream& of);

		//joining functions
		void add(string name_);
		void remove_vertex(string name_);
		void link(string from_, string to_, int type);
		void unlink(string from_, string to_);
		void reverse(string from_, string to_);
		void copy_links(const graph &G);

		//access functions
		int id(string name_);
		string name( int id_);
		vector<string> names();
		bool linked(string from_, string to_);
		int link_type(string from_, string to_);
		bool hasvertex(string name);
		bool can_join(string n1, string n2);

		//Graph searches 
		//Breadth first search BFS
		//depth first search DFS
		int num_disconnected(); //BFS
		void assign_disconnected(vector<graph> &D); //BFS
		bool descendant(string from, string to); //DFS
		bool relatives(string from, string to); //BFS
		void unrelated(vector<string> &I); //BFS

		//data
		map<string, vertex*> vlist; //data
		int nv;	//number of vertices
		int ne; //number of edges
};

//Gviz suitable graph output
//dot -Tpng -O offile
void graph::gviz_dot(ofstream& of){
	of << "digraph G {" << endl;
	for(viter iter=vlist.begin(); iter != vlist.end(); ++iter){
		string tmp = (*iter).second->name;
		replace( tmp.begin(), tmp.end(), '-', '_'); //replace - with _ because - will break gviz!
		of << "\t" << tmp << endl; 
			for(vector<edge>::iterator it=(*iter).second->out.begin(); it!=(*iter).second->out.end(); ++it){
				if( (*it).to_name != "" ){
					string tmp2 = (*it).to_name;
					replace( tmp2.begin(), tmp2.end(), '-', '_'); //replace - with _ because - will break gviz!
					if( (*it).type == -1 || (*it).type == 0){ of << "\t" << tmp << " -> " << tmp2 << endl; }
					if( (*it).type == -1){ of << "\t" << tmp2 << " -> " << tmp << endl; }
					if( (*it).type == 1){ 
						of << "\t" << tmp << " -> " << tmp2 << " [ style=\"dashed\" arrowhead=\"none\"]" <<endl; 
					}
					if( (*it).type == 5){ 
						of << "\t" << tmp << " -> " << tmp2 << " [ color=\"red\" arrowhead=\"none\"]" <<endl; 
					}
					if( (*it).type == 2){ 
						of << "\t" << tmp << " -> " << tmp2 << " [ color=\"blue\" arrowhead=\"none\"]" <<endl; 
					}
				}
			}	
	}
	of << "fixedsize = true\n}" << endl;
}
//Gviz suitable graph output
//neato -Tpng -O offile
//fdp -Tpng -O offile
//nodes are points
void graph::gviz_neato(ofstream& of){
	of << "graph G {" << endl;
	of << "node[label=\"\",shape=\"point\"]" << endl;
	for(viter iter=vlist.begin(); iter != vlist.end(); ++iter){
		string tmp = (*iter).second->name;
		replace( tmp.begin(), tmp.end(), '-', '_'); //replace - with _ because - will break gviz!
		of << "\t" << tmp << endl; 
			for(vector<edge>::iterator it=(*iter).second->out.begin(); it!=(*iter).second->out.end(); ++it){
				if( (*it).to_name != "" ){
					string tmp2 = (*it).to_name;
					replace( tmp2.begin(), tmp2.end(), '-', '_'); //replace - with _ because - will break gviz!
					if( (*it).type == 0 ){ of << "\t" << tmp << " -- " << tmp2 << endl; } 
					else if( (*it).type == 1 ){ of << "\t" << tmp << " -- " << tmp2 << " [style=\"dashed\"]" << endl; } 
					else if( (*it).type == 2 ){ of << "\t" << tmp << " -- " << tmp2 << " [color=\"blue\"]" << endl; } 
					else if( (*it).type == 3 ){ of << "\t" << tmp << " -- " << tmp2 << " [color=\"green\"]" << endl; } 
					else if( (*it).type == 5 ){ of << "\t" << tmp << " -- " << tmp2 << " [color=\"red\"]" << endl; } 
					else { of << "\t" << tmp << " -- " << tmp2 << " [color=\"yellow\"]" << endl; } 
				}
			}	
	}
	of << "fixedsize = true\n}" << endl;
}
//Gviz suitable graph output
//neato -Tpng -O offile
//fdp -Tpng -O offile
//nodes are sample ids
void graph::gviz_neato_named(ofstream& of){
	of << "graph G {" << endl;
	for(viter iter=vlist.begin(); iter != vlist.end(); ++iter){
		string tmp = (*iter).second->name;
		replace( tmp.begin(), tmp.end(), '-', '_'); //replace - with _ because - will break gviz!
		of << "\t" << tmp << endl; 
			for(vector<edge>::iterator it=(*iter).second->out.begin(); it!=(*iter).second->out.end(); ++it){
				if( (*it).to_name != "" ){
					string tmp2 = (*it).to_name;
					if( (*it).type == 0 || (*it).type == -1 ){ of << "\t" << tmp << " -- " << tmp2 << endl; } 
					else if( (*it).type == 1 ){ of << "\t" << tmp << " -- " << tmp2 << " [style=\"dashed\"]" << endl; } 
					else if( (*it).type == 2 ){ of << "\t" << tmp << " -- " << tmp2 << " [style=\"blue\"]" << endl; } 
					else if( (*it).type == 5 ){ of << "\t" << tmp << " -- " << tmp2 << " [color=\"red\"]" << endl; } 
					else { of << "\t" << tmp << " -- " << tmp2 << " [color=\"green\"]" << endl; } 
				}
			}	
	}
	of << "fixedsize = true\n}" << endl;
}
//output .fam file
//fam s1 p1 p2 0 2
void graph::ped_print(ofstream& of, string fam){
	for(viter iter=vlist.begin(); iter != vlist.end(); ++iter){
		of << fam << "\t" << (*iter).second->name << "\t";
		int nu = 0;
		vector<string> parents(2,"0");
		for(vector<edge>::iterator it=(*iter).second->in.begin(); it!=(*iter).second->in.end(); ++it){
			if( (*it).type == 0 || (*it).type == -1 ){ 
				if(nu < 2){ parents[nu] = (*it).from_name; }
				++nu; 
			}
		}
		if(nu < 3){
			for(int i=0; i<2; ++i){ of << parents[i] << "\t"; }
		} else {	//more than 3 parents? something fishy...
			of << "0\t0\t";
		}
		of << "0\t"; //sex
		of << nu << "\n"; //number of parents
	}
}
//samples names in graph
vector<string> graph::names(){
	vector<string> name(nv);
	for(viter iter=vlist.begin(); iter != vlist.end(); ++iter){
		name[ (*iter).second->id ] = (*iter).second->name;
	}
	return name;
}
//add a sample
void graph::add(string name_){
	vertex *v = new vertex(name_, nv++);
	vlist[name_] = v; 
}
//remove a sample
void graph::remove_vertex(string name_){
	
	viter iter = vlist.find (name_); 
	
	for(size_t i=0; i<(*iter).second->out.size(); ++i){	//unlink out
		unlink(name_, (*iter).second->out[i].to_name);	
	}
	for(size_t i=0; i<(*iter).second->in.size(); ++i){	//unlink out
		unlink((*iter).second->out[i].from_name, name_);	
	}
	
	if((*iter).second){ delete (*iter).second; }
	vlist.erase(iter);
	--nv;
}
//sample id in graph
int graph::id(string name_){
	return vlist[name_]->id;
}
//sample name from graph id
string graph::name(int id_){
	for(viter iter=vlist.begin(); iter != vlist.end(); ++iter){
		if( (*iter).second->id == id_ ) return (*iter).second->name;
	}
	return "";
}
//join samples
void graph::link(string from_, string to_, int type=0){
	edge tmp(from_, id(from_), to_, id(to_) , type);
	vlist[from_]->out.push_back(tmp); vlist[from_]->num_out++;
	vlist[to_]->in.push_back(tmp); vlist[to_]->num_in++;
	++ne;
}
//unjoin samples
void graph::unlink(string from_, string to_){

	bool waslinked = false;
	for(size_t i=0; i<vlist[from_]->out.size(); ++i){
		if(vlist[from_]->out[i].to_name == to_){
			vlist[from_]->out.erase(  vlist[from_]->out.begin() + i );
			vlist[from_]->num_out--;
			waslinked = true;
			break;
		}
	}
	for(size_t i=0; i<vlist[to_]->in.size(); ++i){
		if(vlist[to_]->in[i].from_name == from_){
			vlist[to_]->in.erase(  vlist[to_]->in.begin() + i );
			vlist[to_]->num_in--;
			waslinked = true;
			break;
		}
	}
	if( waslinked ){ --ne; }
	
}
//s1->s2 becomes s2->s1
void graph::reverse(string from_, string to_){
	unlink(from_, to_);
	link(to_,from_);
}
//check if samples are linked
bool graph::linked(string from_, string to_){
	for(size_t i=0; i<vlist[from_]->out.size(); ++i){
		if(vlist[from_]->out[i].to_name == to_) return true;
	}
	for(size_t i=0; i<vlist[from_]->in.size(); ++i){
		if(vlist[from_]->in[i].from_name == to_) return true;
	}
	return false;
}
//check if samples are linked and if so what kind of link
int graph::link_type(string from_, string to_){
	for(size_t i=0; i<vlist[from_]->out.size(); ++i){
		if(vlist[from_]->out[i].to_name == to_) return vlist[from_]->out[i].type;
	}
	for(size_t i=0; i<vlist[from_]->in.size(); ++i){
		if(vlist[from_]->in[i].from_name == to_) return vlist[from_]->in[i].type;
	}
	return -1;
}
//is this sample in the graph?
bool graph::hasvertex(string name){
	return ( vlist.find(name) != vlist.end() );
}
//add links to G if they exist in this graph
void graph::copy_links(const graph &G){
	for(const_viter iter = G.vlist.begin(); iter != G.vlist.end(); ++iter){
		if( hasvertex( (*iter).second->name ) ){
			for(vector<edge>::const_iterator it=(*iter).second->out.begin(); it!=(*iter).second->out.end(); ++it){
				if( hasvertex( (*it).to_name ) ){
					link((*iter).second->name, (*it).to_name);
				}
			}
		}
	}
}
//check if sample is ancestor 
bool graph::descendant(string from, string to){
	for(viter iter=vlist.begin(); iter != vlist.end(); ++iter){
		(*iter).second->dist = -1;
	}
	stack<vertex*> S;

	vertex* v=vlist[from];
	S.push(v);

	while( !S.empty() ){
		vertex* u = S.top();
		S.pop();
		if( u->dist == -1 ){
			u->dist = 0;
			for(size_t i=0; i<u->out.size(); ++i){
					if( u->out[i].to_name == to){
						return true;
					}
					S.push(vlist[u->out[i].to_name]);
			}
		}
	}
	return false;
}
//remove relatives until only unrelated set remains
void graph::unrelated(vector<string> &I){
	
	vector<string> ns = names();
	vertex* v= vlist[ ns[ rand() % ns.size() ] ];
	I.push_back(v->name);
	
	
	while( !ns.empty() ){
		
		ns.erase( remove( ns.begin(), ns.end(), v->name ), ns.end() ); 
		for(size_t i=0; i<v->out.size(); ++i){
			ns.erase( remove( ns.begin(), ns.end(), v->out[i].to_name ), ns.end() ); 
		}
		for(size_t i=0; i<v->in.size(); ++i){
			ns.erase( remove( ns.begin(), ns.end(), v->in[i].from_name ), ns.end() ); 
		}
		
		if(ns.size() > 0){ 
			v = vlist[ ns[ rand() % ns.size() ] ];
			I.push_back(v->name);
		}	
	}	
}
//number of components in graph
int graph::num_disconnected(){
	
	int nd = 0;
	vector<string> ns = names();
	
	while( !ns.empty() ){
		for(viter iter=vlist.begin(); iter != vlist.end(); ++iter){
			(*iter).second->dist = -1;
		}
		queue<vertex*> Q;
		
		vertex* v= vlist[ ns[0] ]; 
		v->dist = 0;
		ns.erase( remove( ns.begin(), ns.end(), v->name ), ns.end() ); 
		Q.push(v);
		
		while( !Q.empty() ){
			vertex* u = Q.front();
			Q.pop();
			for(size_t i=0; i<u->out.size(); ++i){
				if( vlist[u->out[i].to_name]->dist == -1 ){
					vlist[u->out[i].to_name]->dist = u->dist+1;
					Q.push(vlist[u->out[i].to_name]);
					ns.erase( remove( ns.begin(), ns.end(), Q.back()->name ), ns.end() ); 
				}
				
			}
			for(size_t i=0; i<u->in.size(); ++i){
				if( vlist[u->in[i].from_name]->dist == -1 ){
					vlist[u->in[i].from_name]->dist = u->dist+1;
					Q.push(vlist[u->in[i].from_name]);
					ns.erase( remove( ns.begin(), ns.end(), Q.back()->name ), ns.end() ); 
				}
			}
		}
		++nd;
	}
	return nd;
}
//each component of graph gets copied to a new graph object
void graph::assign_disconnected(vector<graph> &D){
	
	D.resize(0);
	vector<string> ns = names();
	
	while( !ns.empty() ){
		for(viter iter=vlist.begin(); iter != vlist.end(); ++iter){
			(*iter).second->dist = -1;
		}
		queue<vertex*> Q;
		graph Dt;
		
		vertex* v= vlist[ ns[0] ]; 
		
		v->dist = 0;
		ns.erase( remove( ns.begin(), ns.end(), v->name ), ns.end() ); 
		Q.push(v);
		Dt.add(v->name);
		
		while( !Q.empty() ){
			vertex* u = Q.front();
			Q.pop();
			for(size_t i=0; i<u->out.size(); ++i){
				if( vlist[u->out[i].to_name]->dist == -1 ){
					vlist[u->out[i].to_name]->dist = u->dist+1;
					Q.push(vlist[u->out[i].to_name]);
					Dt.add(Q.back()->name); 
					ns.erase( remove( ns.begin(), ns.end(), Q.back()->name ), ns.end() ); 
				}		
			}
			for(size_t i=0; i<u->in.size(); ++i){
				if( vlist[u->in[i].from_name]->dist == -1 ){
					vlist[u->in[i].from_name]->dist = u->dist+1;
					Q.push(vlist[u->in[i].from_name]);
					Dt.add(Q.back()->name); 
					ns.erase( remove( ns.begin(), ns.end(), Q.back()->name ), ns.end() ); 
				}
			}
		}

		Dt.copy_links(*this);
		D.push_back(Dt);

	}

}
//check if two samples are relatives of any kind
bool graph::relatives(string from, string to){
	

	vector<string> ns = names();
	
	for(viter iter=vlist.begin(); iter != vlist.end(); ++iter){
		(*iter).second->dist = -1;
	}
	queue<vertex*> Q;
		
	vertex* v= vlist[ from ];
	v->dist = 0;
	Q.push(v);
			
	while( !Q.empty() ){
		vertex* u = Q.front();
		Q.pop();
		for(size_t i=0; i<u->out.size(); ++i){
			if( vlist[u->out[i].to_name]->dist == -1 ){
				vlist[u->out[i].to_name]->dist = u->dist+1;
				Q.push(vlist[u->out[i].to_name]);
				if( Q.back()->name == to ){ return true; }
			}
		}
		for(size_t i=0; i<u->in.size(); ++i){
			if( vlist[u->in[i].from_name]->dist == -1 ){
				vlist[u->in[i].from_name]->dist = u->dist+1;
				Q.push(vlist[u->in[i].from_name]);
				ns.erase( remove( ns.begin(), ns.end(), Q.back()->name ), ns.end() ); 
				if( Q.back()->name == to ){ return true; }
			}
		}
	}
	
	return false;
}
//is it possible to join two samples without violating simple rules:
//can't be an ancestor of yourself
//no one can have 3 parents
bool graph::can_join(string n1, string n2){
	
	if( (n1 != n2) && (vlist[n2]->num_in < 2) && ( !linked(n1, n2) ) && ( !descendant(n2, n1) ) ){ 
		return true;
	} 	
	return false;
}


#endif
