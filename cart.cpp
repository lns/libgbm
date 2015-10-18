#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <ctime>
#ifdef _OPEN_MP
#include "omp.h"
#else
#define omp_get_num_threads() 1
#define omp_get_thread_num() 0
#endif

#include "gbm.h"
/* We won't use 16+ core! */
#define MAX_OMP_THREADS 16

Node * GBM::new_node()
{
	assert(node_index<node_max_length);
	node_array[node_index].er = 0;
	return node_array+node_index++;
}

double predict(Node * root, DataSheet* ds, unsigned n)
/*
 * Make prediction for ds->y[n] using CART root.
 *
 * Please check this source and make sure this
 * prediction's logic is what you want, or you
 * should rewrite it.
 */
{
	Node * node = root;
	#ifdef _NAN_MEANS_MISSING
	while(node->er>0)
	{
		if(isnan(ds->x[node->var][n]))
		{
			if(node->nan!=NULL)
				return (node->left->obs*predict(node->left,ds,n)
					+ node->right->obs*predict(node->right,ds,n)
					+ node->nan->obs*predict(node->nan,ds,n))/node->obs;
			else
				return (node->left->obs*predict(node->left,ds,n)
					+ node->right->obs*predict(node->right,ds,n))/node->obs;
		}
		if(ds->x[node->var][n] <= node->val)
			node = node->left;
		else
			node = node->right;
	}
	#endif
	#ifdef _NAN_MEANS_UNAVAILABLE
	while(node->er>0)
	{
		if(isnan(ds->x[node->var][n]))
		{
			if(node->nan!=NULL)
				node = node->nan;
			else
				break;// need revision
		}
		if(ds->x[node->var][n] <= node->val)
			node = node->left;
		else
			node = node->right;
	}
	#endif
	#ifdef _NO_NAN
	while(node->er>0)
	{
		if(ds->x[node->var][n] <= node->val)
			node = node->left;
		else
			node = node->right;
	}
	assert(node->er==0);
	#endif
	return node->val;
}

inline double error_reduction(double  left_sumy,	double  left_obs,
							  double right_sumy,	double right_obs,
							  double   nan_sumy, 	double   nan_obs,
							  double coeff)
/*
 * Calculate the Error Reduction corresponding to a split.
 *
 * For a production release, this function should be
 * cross optimized with "CART_Factory::find_best_split_i".
 * (eg. rewrite er /= (left_obs+right_obs+nan_obs) 
 * 		as er *= coeff; (const double coeff = 1/node->obs) )
 *
 */
{
	double left_ybar = left_sumy/left_obs;
	double right_ybar = right_sumy/right_obs;
	double er = left_obs*right_obs
				*(left_ybar-right_ybar)
				*(left_ybar-right_ybar);
	#ifdef _NAN_MEANS_MISSING
	#endif
	#ifdef _NAN_MEANS_UNAVAILABLE
	if(nan_obs>0)
	{
		double nan_ybar = nan_sumy/nan_obs;
		er += left_obs*nan_obs
			*(left_ybar-nan_ybar)
			*(left_ybar-nan_ybar);
		er += right_obs*nan_obs
			*(right_ybar-nan_ybar)
			*(right_ybar-nan_ybar);
	}
	#endif
	#ifdef _NO_NAN
	#endif
	er *= coeff;
	return er;
}

inline void GBM::find_best_split_i(Node* node, unsigned var, SplitResult * sr)
/* 
 * Find best split position(split_val) for a specified split_var.
 * Result returned to sr.
 *
 * Note:	if best_er of this split_var is not better than sr->er,
 * 			then nothing in sr would be changed.
 *
 * Makesure these vars are already calculated:
 * 		node->sumy
 * 		node->obs
 * Makesure train_set is correctly set.
 */
{
	const datatype* x = ds->x[var];
	const unsigned* l = ht->table[var];
	#ifndef _NO_NAN
	/* 1. Count nan_sum and nan_obs */
	double    nan_sum = 0;
	unsigned  nan_obs = 0;
	for(int i=ht->nan_index[var];i<ht->n;i++)
		if(train_set[l[i]]==node)
		{
			nan_sum += z[l[i]];
			nan_obs++;
		}
	double nan_ybar = nan_sum/nan_obs;
	#else
	#define nan_sum 0
	#define nan_obs 0
	#endif
	/* 2. Prepare */
	const double both_sum = node->sumy - nan_sum;
	const double both_obs = node->obs - nan_obs;
	if(both_obs<2*params->n_minobsinnode)
		return;
	double right_sum,left_sum = 0;
	unsigned right_obs,left_obs = 0;
	unsigned next_i = 0;
	while(train_set[l[next_i]]!=node)
		next_i++;
	/* 3. Find Best Split */
	double best_er = sr->er;
	unsigned best_split_i;
	const double coeff = 1.0/node->obs;
	for(;;)
	{
		unsigned i = next_i++;
		left_sum += z[l[i]];
		left_obs++;
		while(left_obs+params->n_minobsinnode <= both_obs)
		{
			while(train_set[l[next_i]]!=node)
				next_i++;
			if(x[l[next_i]]==x[l[i]])
			{
				i = next_i++;
				left_sum += z[l[i]];
				left_obs++;
			}else
				break;
		}
		right_sum = both_sum - left_sum;
		right_obs = both_obs - left_obs;
		if(left_obs < params->n_minobsinnode)
			continue;
		if(right_obs < params->n_minobsinnode)
			break;
		// You may want to try a different split measurement.
		// then redefine "error_reduction"
		double er = error_reduction(left_sum,left_obs,
									right_sum,right_obs,
									nan_sum,nan_obs,
									coeff);
		if(er>best_er)
		{
			best_er = er;
			best_split_i = i;
		}
	}
	if(best_er>sr->er)
	{
		sr->er = best_er;
		sr->var = var;
		sr->i = best_split_i;
		#ifndef _NO_NAN
		sr->nan_obs  = nan_obs;
		#endif
	}
}

SplitResult * GBM::find_best_split(Node* node)
{
	SplitResult * sr = new SplitResult;
	assert(sr!=NULL);
	sr->node = node;
	sr->er = 0;
	if(node->obs<2*params->n_minobsinnode)
		return sr;
#pragma omp parallel
{
	SplitResult * threadsr[MAX_OMP_THREADS];
	unsigned n_tid;
	#pragma omp single
	{
		n_tid = omp_get_num_threads();
		{
			for(int tid=0;tid<n_tid;tid++)
			{
				threadsr[tid] = new SplitResult;
				threadsr[tid]->node = node;
				threadsr[tid]->er = 0;
			}
		}
	}
	#pragma omp for
	for(int j=0;j<ds->p;j++)
	{
		if(var_set[j]==false)//Skip
			continue;
		unsigned tid = omp_get_thread_num();
		find_best_split_i(node,j,threadsr[tid]);
	}
	#pragma omp single
	{
		for(int tid=0;tid<n_tid;tid++)
		{
			if(threadsr[tid]->er > sr->er)
				*sr = *threadsr[tid];
			delete threadsr[tid];
		}
	}
}
	return sr;
}

void GBM::update_node(Node* t,Node* parent,const unsigned * l,
				unsigned start,unsigned end)
/*
 * update train_set and t->sumy t->obs
 */
{
	assert(t->sumy==0);
	assert(t->obs==0);
	for(int i=start;i<end;i++)
	{
		unsigned k = l[i];
		if(train_set[k]==parent)
		{
			train_set[k] = t;
			t->sumy += z[k];
			t->obs++;
		}
	}
}

/*
 * Estimate the predict value of terminal node.
 */
#ifdef _GAUSSIAN
void GBM::estimate(Node* root)
{
	if(root->left)
		estimate(root->left);
	if(root->right)
		estimate(root->right);
	#ifndef _NO_NAN
	if(root->nan)
		estimate(root->nan);
	#endif
	if(root->er==0)
	{
		assert(root->var==~0);
		assert(isnan(root->val));
		root->val = root->sumy/root->obs;
	}
}
#endif

#ifdef _BERNOULLI
void GBM::estimate(Node* root)
{
	Node * node = root;
	unsigned nnd = params->n_depth;
	nnd = 1+3*nnd+3*(2*nnd+1);
	double * tmp = new double[2*nnd];
	assert(tmp!=NULL);
	for(int i=0;i<2*nnd;i++)
		tmp[i] = 0;
	for(int i=0;i<ht->n;i++)
	{
		if(train_set[i]==NULL)
			continue;
		double p = 1/(1+exp(-f[i]));
		unsigned offset = train_set[i]-root;
		tmp[2*offset] += y[i] - p;
		tmp[2*offset+1] += p*(1-p);
	}
	for(int i=0;i<nnd;i++)
		if(tmp[2*i+1]!=0)
		{
			assert(node[i].er==0);
			assert(node[i].var==~0);
			assert(isnan(node[i].val));
			node[i].val = tmp[2*i]/tmp[2*i+1];
		}
	delete[] tmp;
}
#endif

#ifdef _RSQUARE
/*
void GBM::estimate(Node* root)
{
	Node * node = root;
	unsigned nnd = params->n_depth;
	nnd = 1+3*nnd+3*(2*nnd+1); 
	double * tmp = new double[2*nnd];
	assert(tmp!=NULL);
	for(int i=0;i<2*nnd;i++)
		tmp[i] = 0;
	for(int i=0;i<ht->n;i++)
	{
		if(train_set[i]==NULL)
			continue;
		unsigned offset = train_set[i]-root;
		tmp[2*offset] += y[i];
		tmp[2*offset+1] += f[i];
	}
	for(int i=0;i<nnd;i++)
		if(tmp[2*i+1]!=0)
		{
			if(node[i].er!=0)
				*(int*)NULL = 0;
			assert(node[i].er==0);
			assert(node[i].var==~0);
			assert(isnan(node[i].val));	
			node[i].val = (fbar*nobs-tmp[2*i+1])/(nobs-node[i].obs)
						 + (tmp[2*i]/node[i].obs-ybar)/ratio*nobs/(nobs-node[i].obs);
		}
	delete[] tmp;
}
*/

void GBM::estimate(Node* root)
{
	if(root->left)
		estimate(root->left);
	if(root->right)
		estimate(root->right);
	#ifndef _NO_NAN
	if(root->nan)
		estimate(root->nan);
	#endif
	if(root->er==0)
	{
		assert(root->var==~0);
		assert(isnan(root->val));
		root->val = root->sumy/root->obs;
	}
}

#endif

void GBM::grow(SplitResult * sr)
/*
 * Make new child node according to SplitResult.
 */
{
	Node* node = sr->node;
	assert(node->var==~0);
	assert(isnan(node->val));
	assert(node->er == 0);
	node->var = sr->var;
	node->er = sr->er;
	const unsigned * l = ht->table[sr->var];
	const unsigned nan_index = ht->nan_index[sr->var];
	unsigned next_i = sr->i+1;
	// split_val = (this_val + next_val) / 2
	while(train_set[l[next_i]]!=node)
		next_i++;
	const datatype* x = ds->x[node->var];
	node->val = (x[l[sr->i]]+x[l[next_i]])/2;
	node->left = new_node();
	update_node(node->left,node,l,0,sr->i+1);
	node->right = new_node();
	update_node(node->right,node,l,sr->i+1,nan_index);
	#ifndef _NO_NAN
	if(sr->nan_obs>0)// Should this be n_minobsinnode?
	{
		node->nan = new_node();
		update_node(node->nan,node,l,nan_index,ht->n);
	}
	#endif
}

Node* GBM::CART(int depth)
{
	grow_candidate_index = 0;
	grow_candidate_max_length = 2*depth;
	grow_candidate = new SplitResult*[2*depth];
	unvisited_index = 0;
	unvisited_max_length = 3;
	unvisited = new Node*[3];
	assert(grow_candidate!=NULL and unvisited !=NULL);
	Node* root = new_node();
	/* 2. Prepare root node */
	root->sumy = 0;
	root->obs = 0;
	for(int i=0;i<ht->n;i++)
		if(train_set[i]==(Node*)0x1)
		{
			train_set[i] = root;
			root->sumy += z[i];
			root->obs ++;
		}
	unvisited[unvisited_index++] = root;
	for(;depth>0;depth--)
	{
		/* 1.Calculate Split_Result for each node in "unvisited". */
		while(unvisited_index>0)
			grow_candidate[grow_candidate_index++] = 
				find_best_split(unvisited[--unvisited_index]);
		assert(unvisited_index==0);
		/* 2.Grow best Split_Result. */
		double best_er = 0;
		unsigned best_split_index;
		for(int i=0;i<grow_candidate_index;i++)
			if(grow_candidate[i]->er>best_er)
			{
				best_split_index = i;
				best_er = grow_candidate[i]->er;
			}
		if(best_er==0)
			break;
		grow(grow_candidate[best_split_index]);
		SplitResult * best_split = grow_candidate[best_split_index];
		//SplitResult * sr = best_split;
		/* 3. Put new nodes in "unvisited". */
		unvisited[unvisited_index++] = best_split->node->left;
		unvisited[unvisited_index++] = best_split->node->right;
		#ifndef _NO_NAN
		if(best_split->node->nan!=NULL)
			unvisited[unvisited_index++] = best_split->node->nan;
		#endif
		delete grow_candidate[best_split_index];
		grow_candidate[best_split_index] = grow_candidate[--grow_candidate_index];
	}
	#ifdef _EXTENDED
	for(int i=0;i<grow_candidate_index;i++)
		grow(grow_candidate[i]);
	#endif
	estimate(root);
	for(int i=0;i<grow_candidate_index;i++)
		delete grow_candidate[i];
	delete[] grow_candidate;
	delete[] unvisited;
	return root;
}

