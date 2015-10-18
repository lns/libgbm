#ifndef GBM_H
#define GBM_H

#include "node.h"
#include "datasheet.h"
#include "hashtable.h"
#include "split_result.h"
#include "parameters.h"

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

	//gbm.cpp
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

void CrossValidate(DataSheet* ds,Parameters* params,unsigned cvfolds);
void CrossPredict(DataSheet* ds,unsigned cvfolds,unsigned best_iter);

double predict(Node* root,DataSheet* ds,unsigned index);

//boost.cpp
#endif
