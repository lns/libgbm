#ifndef LIBGBM_H
#define LIBGBM_H

#ifndef DEFS_H
#define DEFS_H
/* 
 * Defines NAN,INF
 * INFINITY,POS_INF,NEG_INF
 */
#ifndef NAN
	#ifdef __INTEL_COMPILER
		#define NAN ((__builtin_nanl("")))
	#else
		#include <cmath>
		#define NAN (-sqrt(-1))
	#endif
#endif

#ifndef INF
	#ifdef __INTEL_COMPILER
		#define INF ((__builtin_huge_vall()))
	#else
		#include <cmath>
		#ifdef INFINITY
			#define INF INFINITY
		#else
			#define INF ((double)1/0)
		#endif
	#endif
#endif

#ifndef INFINITY
#define INFINITY INF
#endif
#ifndef POS_INF
#define POS_INF INF
#endif
#ifndef NEG_INF
#define NEG_INF (-INF)
#endif
#endif

#include <string>

typedef float datatype;
class DataSheet
{
public:
	std::string file_name;
	unsigned p;
	unsigned n;
	datatype* y;
	datatype** x;
	std::string y_name;
	std::string* var_name;
	double* w;

	void init(unsigned N,unsigned P);
	void init();
	DataSheet(unsigned N,unsigned P);
	DataSheet(char *,char delimiter=',',bool haveY=true);
	DataSheet();
	~DataSheet();
	void shuffle(unsigned n);
	void head(unsigned n=20,int p=9);
	void read_csv(char *,char delimiter=',',bool haveY=true);
	void write_csv(char *,char delimiter=',');
	// TODO:
	void dump(char * file_name);
	void load(char * file_name);
};

class Node{
public:
	Node* left;
	Node* right;
	#ifndef _NO_NAN
	Node* nan;
	#endif
	unsigned var;
	double val;
	double er;

	double sumy;
	unsigned obs;

	Node();
};

class Parameters
{
public:
	unsigned n_iters;
	unsigned n_bags;
	double init_shrinkage;
	double var_fraction;
	double bag_fraction;
	double train_fraction;
	unsigned n_minobsinnode;
	unsigned n_depth;

	void print();
};

struct SplitResult
{
	Node* node;
	double er;
	unsigned var;
	unsigned i;
	#ifndef _NO_NAN
	double nan_obs;
	#endif
};
typedef struct SplitResult SplitResult;

class HashTable
{
public:
	unsigned p;
	unsigned n;
	unsigned * nan_index;
	unsigned** table;

	HashTable(DataSheet* ds,unsigned n);
	~HashTable();
};

class GBM
{
private:
	Node* new_node();
	void find_best_split_i(Node* node,unsigned var,SplitResult * sr);
	SplitResult * find_best_split(Node* node);
	void grow(SplitResult * sr);
	void update_node(Node* node,Node* parent,const unsigned * l,
				unsigned start, unsigned end);
	void estimate(Node* root);
	Node* CART(int depth);

public:
	DataSheet* ds;
	HashTable* ht;
	Parameters* params;

//for each in n_bag:
	datatype * y;
	double * f;
	double * z;
	Node** train_set;
	bool* var_set;

	double * train_error;
	double *  test_error;

	double init_val;
	double shrinkage;

	//cart.cpp
	Node * node_array;
	unsigned node_index;
	unsigned node_max_length;
	Node ** tree_array;
	unsigned tree_index;
	unsigned tree_max_length;

	SplitResult ** grow_candidate;
	unsigned grow_candidate_index;
	unsigned grow_candidate_max_length;
	Node ** unvisited;
	unsigned unvisited_index;
	unsigned unvisited_max_length;

	//machine.cpp
	GBM();
	GBM(DataSheet* ds, Parameters* p);
	GBM(char * file_name);
	~GBM();
	void fill(DataSheet* ds,unsigned best_iter);
	void dump(char * file_name);
	void load(char * file_name);
	unsigned find_best_iter();
	void relevance(double* array, unsigned best_iter);

	//boost.cpp
	void boost();
};
double predict(Node* root,DataSheet* ds,unsigned index);

#endif
