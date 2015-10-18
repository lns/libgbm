#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <iostream>
#include <fstream>
#include <ctime>

#include "gbm.h"

template <typename T>
inline T * subsample(unsigned m,unsigned n,
						T a, T b)
{
	T * set = new T[n];
	assert(set!=NULL);
	const unsigned N = n;
	for(int i=0;i<N;i++)
	{
		bool result = n*drand48()<m;
		if(result)
		{
			set[i] = a;
			m--;
		}else{
			set[i] = b;
		}
		n--;
	}
	return set;
}

inline void update_f(double* f, DataSheet* ds,Node** tree,unsigned n_bags,
						unsigned iter,double shrinkage)
{
	for(int i=0;i<ds->n;i++)
		for(int k=0;k<n_bags;k++)
		{
			double pred = predict(tree[n_bags*iter+k],ds,i);
			assert(!isnan(pred));
			f[i] += shrinkage*pred;
		}
}
inline void update_z(datatype* y,double* f,double* z,unsigned start,unsigned end,GBM* pgbm)
{
	for(int i=start;i<end;i++)
	{
		#ifdef _GAUSSIAN
		z[i] = y[i]-f[i];
		#endif
		#ifdef _BERNOULLI
		z[i] = y[i]-1/(1+exp(-f[i]));
		#endif
		#ifdef _RSQUARE
		#define LAMBDA_1 0.1
		z[i] = y[i]-pgbm->ybar-pgbm->ratio*(f[i]-pgbm->fbar)+LAMBDA_1/pgbm->nobs*(pgbm->ybar-pgbm->fbar);
		#endif
	}
}

#if (defined _BERNOULLI) || (defined _GAUSSIAN)
inline double loss(double& variance,datatype* y,double* f,unsigned start,unsigned end)
{
	double res = 0;
	for(int i=start;i<end;i++)
	{
		#ifdef _GAUSSIAN
		double d = (y[i]-f[i])*(y[i]-f[i]);
		res += d;
		variance += d*d;
		#endif
		#ifdef _BERNOULLI
		double d = -2*(y[i]*f[i]-log(1+exp(f[i])));
		res += d;
		variance += d*d;
		#endif
	}
	variance = (variance - res*res/(end-start))/(end-start);
	return res/(end-start);
}
#endif
#ifdef _RSQUARE
inline double loss(double& variance,datatype* y,double* f,unsigned start,unsigned end)
{
	double dyy=0,dyf=0,dff=0,sy=0,sf=0;
	for(int i=start;i<end;i++)
	{
		sy += y[i];
		sf += f[i];
	}
	sy /= end-start;
	sf /= end-start;
	for(int i=start;i<end;i++)
	{
		dyy += (y[i]-sy)*(y[i]-sy);
		dyf += (y[i]-sy)*(f[i]-sf);
		dff += (f[i]-sf)*(f[i]-sf);
	}
	double res = dyf/sqrt(dyy*dff);
	variance = NAN;
	return res;
}
#endif

inline double error_reduction(DataSheet* ds,unsigned start,unsigned end,
							Node * root,datatype*y,double*f,double*z)
{
	double res = 0;
	for(int i=start;i<end;i++)
	{
		double pred = predict(root,ds,i);
		res += pred*z[i];
	}
	return res/(end-start);
}

inline double error_reduction(DataSheet* ds,Node** set,unsigned n,
							Node * root,datatype*y,double*f,double*z)
{
	//only count i: set[i]==0
	double res = 0;
	unsigned count = 0;
	for(int i=0;i<n;i++)
	{
		if(set[i]==0)
		{
			double pred = predict(root,ds,i);
			res += pred*z[i];
			count++;
		}
	}
	return res/count;
}

inline double initial_value(unsigned m,datatype* y)
{
	double avg_y = 0;
	for(int i=0;i<m;i++)
		avg_y += y[i];
	avg_y /= m;
	#ifdef _GAUSSIAN
	return avg_y;
	#endif
	#ifdef _BERNOULLI
	assert(avg_y!=1);
	return log(avg_y/(1-avg_y));
	#endif
	#ifdef _RSQUARE
	return avg_y;
	#endif
}

void GBM::boost()
{
	/* 0.Generate Hash Table for training data. */
	const unsigned n = ds->n;
	const unsigned m = n*params->train_fraction;
	// GBM won't use the test_set in any way.
	// It's safe to use the test_set for a later validation.
	y = ds->y;
	z = new double[n];
	f = new double[n];
	assert(z!=NULL and f!=NULL);
	ht = new HashTable(ds,m);
	printf("Boost: HashTable generated.\n");
	/* 1.Fill with Initial value */
	params->print();
	/* 1.1 Compute Initial value */
	init_val = initial_value(m,y);
	printf("Init_val: %lf\n",init_val);
	for(int i=0;i<n;i++)
	{
		f[i] = init_val;
		#ifdef _GAUSSIAN
		z[i] = y[i] - init_val;
		#endif
		#ifdef _BERNOULLI
		z[i] = y[i] - 1/(1+exp(-init_val));
		#endif
		#ifdef _RSQUARE
		z[i] = y[i] - init_val;
		#endif
	}
	/* 2.Main Loop: fit and update. */
	shrinkage = params->init_shrinkage;
	printf("  Iter\t  TrainLoss\t  TestLoss\t");
	#if (defined _BERNOULLI) || (defined _GAUSSIAN)
	printf("  StdDeviation");
	#endif
	#ifdef _RSQUARE
	printf("    Ratio     ");
	#endif
	printf("\t  Shrinkage\t  Improve\n");
	#ifdef _RSQUARE
	ybar = 0;
	for(int i=0;i<m;i++)
		ybar += y[i];
	nobs = m;
	ybar /= nobs;
	fbar = ybar;
	ratio = 1.0;
	#endif
	for(int iter=0;iter<params->n_iters;iter++)
	{
		#pragma omp parallel for
		for(int k=0;k<params->n_bags;k++)
		{
			/* 2.1 Set variables and samples for BasicLearner_k */
			//label t out of m as 1, others 0
			unsigned t = m*params->bag_fraction;
			train_set = subsample(t,m,(Node*)0x1,(Node*)0x0);
			//label q out of p as true, others false
			unsigned q = ds->p*params->var_fraction;
			var_set = subsample(q,ds->p,true,false);
			/* 2.2 Fit BasicLearner_k to z, return results to *_er */
			Node* root = CART(params->n_depth);
			tree_array[tree_index] = root;
			tree_index++;
			delete[] train_set;
			delete[] var_set;
		}
		/* 3.1 Update f and z, Calculate Test Loss */
		double train_error_var, test_error_var;
		update_f(f,ds,tree_array,params->n_bags,iter,shrinkage/params->n_bags);
		#ifdef _RSQUARE
		double dyf=0,dff=0;
		fbar = 0;
		for(int i=0;i<m;i++)
			fbar += f[i];
		nobs = m;
		fbar /= nobs;
		for(int i=0;i<m;i++)
		{
			dyf += (y[i]-ybar)*(f[i]-fbar);
			dff += (f[i]-fbar)*(f[i]-fbar);
		}
		printf("nobs: %u \t ybar: %lf \t fbar: %lf \t dyf: %lf \t dff: %lf \n",
				nobs,ybar,fbar,dyf,dff);
		ratio = dyf/dff;
		#endif
		update_z(y,f,z,0,m,this);
		train_error[iter] = loss(train_error_var,y,f,0,m);
		test_error[iter]  = loss(test_error_var, y,f,m,n);
		/* 3.2 Print results */
		if(iter<=100 or iter%10==0)
		{
			printf("%5d",iter);
			printf("\t%12.8lf",train_error[iter]);
			printf("\t%12.8lf",test_error[iter]);
			#if (defined _BERNOULLI) || (defined _GAUSSIAN)
			printf("\t%12.8lf",sqrt(test_error_var/(n-m)/params->n_bags));
			#endif
			#ifdef _RSQUARE
			printf("\t%12.8le",ratio);
			#endif
			printf("\t%12.8lf",shrinkage);
			printf("\t%12.8lf\n",iter>0?
				test_error[iter-1]-test_error[iter] : NAN);
		}
	}
	printf("node_max_length: %u, node_index: %u\n",
			node_max_length,node_index);
	/*
	fill(ds,50);
	for(int i=0;i<ds->n;i++)
		if((y[i]-f[i])*(y[i]-f[i])>1e-10)
		{
			printf("%u %le %le\n",i,y[i],f[i]);
			exit(-1);
		}
	*/
	delete[] z;
	delete[] f;
}

void CrossValidate(DataSheet* ds,Parameters* params,unsigned cvfolds)
/*
 * Cross Validate
 */
{
	DataSheet* cvds = new DataSheet(ds->n,ds->p);
	double * train_error = new double[params->n_iters*cvfolds];
	double * test_error = new double[params->n_iters*cvfolds];
	for(int fold=0;fold<cvfolds;fold++)
	{
		printf("===== FOLD %02d =====\n",fold);
		unsigned idx = 0;
		const unsigned mark1 = ds->n*fold/cvfolds;
		const unsigned mark2 = ds->n*(fold+1)/cvfolds;
		/* 1. Making DataSheet */
		idx = cvds->overwrite(idx,ds,0,mark1);
		idx = cvds->overwrite(idx,ds,mark2,ds->n);
		idx = cvds->overwrite(idx,ds,mark1,mark2);
		assert(idx==ds->n);
		params->train_fraction = (double)(ds->n+mark1-mark2)/ds->n;
		/* 2. Start Training */
		char name[200];
		sprintf(name,"gbmcv_fold%02d.gz",fold);
		GBM* gbmcv = new GBM(cvds,params);
		gbmcv->boost();
		/* 3. Save Results */
		for(int i=0;i<params->n_iters;i++)
		{
			train_error[params->n_iters*fold+i] = gbmcv->train_error[i];
			test_error [params->n_iters*fold+i] = gbmcv->test_error[i];
		}
		gbmcv->dump(name);
	}
	delete cvds;
	/* 4. Output CV Error Result */
	FILE* f = fopen("cverror.csv","w");
	fprintf(f,"iter");
	for(int fold=0;fold<cvfolds;fold++)
		fprintf(f,",train_fold%02d,test_fold%02d",fold,fold);
	fprintf(f,"\n");
	for(int i=0;i<params->n_iters;i++)
	{
		fprintf(f,"%d",i);
		for(int fold=0;fold<cvfolds;fold++)
			fprintf(f,",%12.8lf,%12.8lf",
					train_error[params->n_iters*fold+i],
					test_error [params->n_iters*fold+i]);
		fprintf(f,"\n");
	}
	fclose(f);
	delete[] train_error;
	delete[] test_error;
}

void CrossPredict(DataSheet* ds,unsigned cvfolds,unsigned best_iter)
/*
 * Make Cross Prediction with CV Models
 */
{
	DataSheet* cvds = new DataSheet(ds->n,ds->p);
	FILE* f = fopen("cvresult.csv","w");
	double** rele = new double*[cvfolds];
	for(int fold=0;fold<cvfolds;fold++)
		rele[fold] = new double[ds->p];
	fprintf(f,"prediction\n");
	for(int fold=0;fold<cvfolds;fold++)
	{
		printf("===== PREDICT FOLD %02d =====\n",fold);
		unsigned mark1 = ds->n*fold/cvfolds;
		unsigned mark2 = ds->n*(fold+1)/cvfolds;
		cvds->overwrite(0,ds,mark1,mark2);
		cvds->n = mark2-mark1;
		char name[200];
		sprintf(name,"gbmcv_fold%02d.gz",fold);
		GBM* gbmcv = new GBM(name);
		gbmcv->fill(cvds,best_iter);
		gbmcv->relevance(rele[fold],best_iter);
		for(int i=0;i<cvds->n;i++)
			fprintf(f,"%12.8lf\n",cvds->y[i]);
		cvds->n = ds->n;
	}
	fclose(f);
	FILE* g = fopen("cvrele.csv","w");
	fprintf(g,"Var,");
	for(int fold=0;fold<cvfolds;fold++)
	{
		fprintf(g,"VarImportance%02d",fold);
		if(fold==cvfolds-1)
			fprintf(g,"\n");
		else
			fprintf(g,",");
	}
	for(int i=0;i<ds->p;i++)
	{
		fprintf(g,"%s,",ds->var_name[i].c_str());
		for(int fold=0;fold<cvfolds;fold++)
		{
			fprintf(g,"%12.8lf",rele[fold][i]);
			if(fold==cvfolds-1)
				fprintf(g,"\n");
			else
				fprintf(g,",");
		}
	}
	fclose(g);
	delete cvds;
	for(int fold=0;fold<cvfolds;fold++)
		delete[] rele[fold];
	delete[] rele;
}
