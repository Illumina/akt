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
#include <math.h>
#include "logs.hpp"

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
	int degree() {return(num_in + num_out);}

	vertex(string name_, int id_)
	{
			name = name_;
			id = id_;
			num_in = 0;
			num_out = 0;
	}

};

vertex* copy_vertex(vertex *in);

typedef map<string, vertex*>::iterator viter;
typedef map<string, vertex*>::const_iterator const_viter;

//a hash table of vertex pointers
class graph
{
public:
	//default
	graph()
	{
		nv=0;
		ne=0;
	}
	//constructor
	graph(string name_)
	{
		vertex *v = new vertex(name_, 0);
		vlist[name_] = v;
		nv = 1;
		ne = 0;
	}
	//destruction of graph
	void clear()
	{
		for( viter iter= vlist.begin(); iter!=vlist.end(); iter++)
		{
			if((*iter).second)
			{
				delete iter->second;
			}
		}
		vlist.clear();
	}
	//copy constructor
	graph( const graph& other )
	{
		nv = other.nv;
		ne = other.ne;

		clear();
		for( const_viter iter= other.vlist.begin(); iter!=other.vlist.end(); ++iter)
		{
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
		for( const_viter iter= second.vlist.begin(); iter!=second.vlist.end(); ++iter)
		{
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
	void unrelatedGreedy(vector<string> &I); //BFS

	//data
	map<string, vertex*> vlist; //data
	int nv;	//number of vertices
	int ne; //number of edges
};


#endif
