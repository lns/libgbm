#include <cmath>
#include <cstdlib>
#include <cassert>
#include <cstddef>
#include  <cstdio>
#include "zlib.h"
#include "gbm.h"

GBM::GBM()
{
	tree_array = NULL;
	node_array = NULL;
	train_error = NULL;
	test_error = NULL;
}

GBM::GBM(DataSheet* datasheet, Parameters* p)
{	
	ds = datasheet;
	params = p;
	tree_max_length = p->n_iters*p->n_bags;
	tree_array = new Node*[tree_max_length];
	tree_index = 0;
	unsigned n_depth = p->n_depth;
	node_max_length = tree_max_length;
	#ifndef _NO_NAN
		#ifdef _EXTENDED
		node_max_length *= (1+3*n_depth+3*(2*n_depth+1));
		#else
		node_max_length *= (1+3*n_depth);
		#endif
	#else
		#ifdef _EXTENDED
		node_max_length *= (1+2*n_depth+2*(n_depth+1));
		#else
		node_max_length *= (1+2*n_depth);
		#endif
	#endif
	node_array = new Node[node_max_length];
	node_index = 0;
	train_error = new double[p->n_iters];
	test_error = new double[p->n_iters];
	assert(tree_array!=NULL and node_array!=NULL 
		and train_error!=NULL and test_error!=NULL);
}

GBM::~GBM()
{
	delete[] tree_array;
	delete[] node_array;
	delete[] train_error;
	delete[] test_error;
}

void GBM::fill(DataSheet* dstest,unsigned best_iter)
{
	const unsigned n = dstest->n;
	const unsigned n_bags = params->n_bags;
	#pragma omp parallel for
	for(unsigned i=0;i<n;i++)
	{
		dstest->y[i] = init_val;
		for(unsigned t=0;t<best_iter;t++)
			for(unsigned s=0;s<n_bags;s++)
				dstest->y[i] += shrinkage/n_bags*predict(tree_array[t*n_bags+s],dstest,i);
		#ifdef _BERNOULLI
		dstest->y[i] = 1/(1+exp(-dstest->y[i]));
		#endif
		if(i%10000==0)
			printf("i: %d\n",i);
	}
}

/* You should make sure that sizeof(unsigned long)==sizeof(void*) */
#define RAWADD(x,y) ((Node*)((unsigned long)(x)+(unsigned long)(y)))
#define RAWSUB(x,y) ((Node*)((unsigned long)(x)-(unsigned long)(y)))

void GBM::dump(char * file_name)
/*
 * Dump GBM to pickle file.
 */
{
	assert(sizeof(unsigned long)==sizeof(Node*));
	for(unsigned i=0;i<tree_index;i++)
		tree_array[i] = RAWSUB(tree_array[i],node_array);
	for(unsigned i=0;i<node_index;i++)
	{
		if(node_array[i].left  != NULL)
			node_array[i].left = RAWSUB(node_array[i].left,node_array);
		if(node_array[i].right != NULL)
			node_array[i].right= RAWSUB(node_array[i].right,node_array);
		#ifndef _NO_NAN
		if(node_array[i].nan   != NULL)
			node_array[i].nan  = RAWSUB(node_array[i].nan,node_array);
		#endif
	}
	Parameters * tmpparams = params;
	Node **tmptree = tree_array;
	Node * tmpnode = node_array;
	double *tmptre = train_error;
	double *tmptee = test_error;
	unsigned n_iters = params->n_iters;
	params = NULL;
	tree_array = NULL;
	node_array = NULL;
	train_error = NULL;
	test_error = NULL;
	grow_candidate = NULL;
	unvisited = NULL;
	gzFile f = gzopen(file_name,"w");
	gzbuffer(f,131072);
	gzwrite(f,this,sizeof(GBM));
	gzwrite(f,tmpparams,sizeof(Parameters));
	gzwrite(f,tmpnode,node_index*sizeof(Node));
	gzwrite(f,tmptree,tree_index*sizeof(Node*));
	gzwrite(f,tmptre,n_iters*sizeof(double));
	gzwrite(f,tmptee,n_iters*sizeof(double));
	for(int i=0;i<n_iters;i++)
		printf("train_error: %12.10lf  test_error: %12.10lf\n",tmptre[i],tmptee[i]);
	gzclose(f);
	return;
}

void GBM::load(char * file_name)
{
	assert(access(file_name,R_OK)==0);
	gzFile f = gzopen(file_name,"r");
	gzbuffer(f,131072);
	gzread(f,this,sizeof(GBM));
	printf("init_val: %lf\n",init_val);
	printf("shrinkage: %lf\n",shrinkage);
	printf("node_index: %u\n",node_index);
	printf("tree_index: %u\n",tree_index);
	printf("node_max_length: %u\n",node_max_length);
	printf("tree_max_length: %u\n",tree_max_length);
	params = new Parameters();
	gzread(f,params,sizeof(Parameters));
	params->print();
	tree_array = new Node*[tree_max_length];
	node_array = new Node[node_max_length];
	train_error = new double[params->n_iters];
	test_error = new double[params->n_iters];
	assert(tree_array!=NULL and node_array!=NULL
		and train_error!=NULL and test_error!=NULL);
	gzread(f,node_array,node_index*sizeof(Node));
	gzread(f,tree_array,tree_index*sizeof(Node*));
	for(int i=0;i<tree_index;i++)
		tree_array[i] = RAWADD(tree_array[i],node_array);
	for(int i=0;i<node_index;i++)
	{
		if(node_array[i].left!=NULL)
			node_array[i].left  = RAWADD(node_array[i].left,node_array);
		if(node_array[i].right!=NULL)
			node_array[i].right = RAWADD(node_array[i].right,node_array);
		#ifndef _NO_NAN
		if(node_array[i].nan!=NULL)
			node_array[i].nan   = RAWADD(node_array[i].nan,node_array);
		#endif
	}
	unsigned n_iters = params->n_iters;
	gzread(f,train_error,n_iters*sizeof(double));
	gzread(f,test_error,n_iters*sizeof(double));
	for(int i=0;i<params->n_iters;i++)
		printf("train_error: %12.10lf  test_error: %12.10lf\n",train_error[i],test_error[i]);
	gzclose(f);
	return;
}

GBM::GBM(char * file_name)
{
	load(file_name);
	return;
}

void update_relevance(double* array,Node* node,double shrinkage)
{
	if(node==NULL || node->er==0)//not split
		return;
	array[node->var] += shrinkage*node->er;
	update_relevance(array,node->left,shrinkage);
	update_relevance(array,node->right,shrinkage);
	#ifndef _NO_NAN
	update_relevance(array,node->nan,shrinkage);
	#endif
}

void GBM::relevance(double* array,unsigned best_iter)
{
	unsigned n_bags = params->n_bags;
	for(unsigned t=0;t<best_iter;t++)
		for(unsigned s=0;s<n_bags;s++)
			update_relevance(array,tree_array[t*n_bags+s],shrinkage);
	return;
}

unsigned GBM::find_best_iter()
{
	double best_test_error = test_error[0];
	unsigned best_iter = 0;
	for(int i=1;i<params->n_iters;i++)
		if(test_error[i]<best_test_error)
		{
			best_iter = i;
			best_test_error = test_error[i];
		}
	return best_iter;
}
