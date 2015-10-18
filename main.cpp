#include <cstdio>
#include <cstdlib>
#include "gbm.h"

int main()
{
/*
	DataSheet* origin = new DataSheet("/home/linus/playground/SAS/features/features01");
	DataSheet* ds = new DataSheet(origin->n,origin->p-1);
	ds->y = origin->x[0];
	for(int i=0;i<ds->p;i++)
		ds->x[i] = origin->x[i+1];
	ds->head();
	ds->shuffle(ds->n);
*/
/*
	FILE* f = fopen("target.csv","w");
	fprintf(f,"target\n");
	for(int i=0;i<ds->n;i++)
		fprintf(f,"%12.8lf\n",ds->y[i]);
	fclose(f);
	CrossPredict(ds,10,1846);
	exit(0);
*/
/*
	DataSheet* ds = new DataSheet("/media/Main/kaggle/MMAC/Data/TrainingSet/ACT13_competition_training_ordered.csv");
	ds->head();
	delete[] ds->y;
	ds->y = ds->x[0];
	for(int i=0;i<ds->p-2;i++)
		ds->x[i] = ds->x[i+1];
	ds->p = ds->p-1;
	//for(int i=ds->p/2;i<ds->p;i++)
	//	delete[] ds->x[i];
	ds->p = ds->p/2;
*/	
/*	FILE* f = fopen("predict.gbm6","r");
	for(int i=0;i<ds->n;i++)
	{
		double temp;
		fscanf(f,"%lf\n",&temp);
		ds->y[i] -= temp;
	}
*/
/*
	ds->shuffle(ds->n);
	Parameters * params = new Parameters;
	params->n_iters = 500;
	params->n_bags = 1;
	params->train_fraction = 0.9;
	params->bag_fraction = 0.5;
	params->var_fraction = 0.5;
	params->init_shrinkage = 0.01;
	params->n_minobsinnode = 30;
	params->n_depth = 30;
	GBM* gbm = new GBM(ds,params);
	gbm->boost();
//	GBM* gbm = new GBM("mmac_13a.gz");
//	gbm->params->print();
	float * target = new float[ds->n];
	for(int i=0;i<ds->n;i++)
		target[i] = ds->y[i];
	unsigned best_iter = gbm->find_best_iter();
	printf("Best_iter: %u, score: %lf\n",best_iter,gbm->test_error[best_iter]);
	gbm->fill(ds,best_iter+1);
	double dxx,dxy,dyy,xbar,ybar;
	dxx=0;dxy=0;dyy=0;xbar=0;ybar=0;
	unsigned n = 0;
	for(int i=ds->n*params->train_fraction+1;i<ds->n;i++)
	{
		xbar += target[i];
		ybar += ds->y[i];
		n++;
	}
	xbar /= n;
	ybar /= n;
	for(int i=ds->n*params->train_fraction+1;i<ds->n;i++)
	{
		dxx += (target[i]-xbar)*(target[i]-xbar);
		dxy += (target[i]-xbar)*(ds->y[i]-ybar);
		dyy += (ds->y[i]-ybar)*(ds->y[i]-ybar);
	}
	double r2 = dxy*dxy/dxx/dyy;
	printf("Rsquare: %lf\n",r2);
	delete[] target;
	gbm->dump("mmac_13c.gz");
*/
///*
	DataSheet* ds = new DataSheet("train.csv",',',true);
	ds->head();
	Parameters * params = new Parameters;
	params->n_iters = 2000;//800;//1500;
	params->n_bags = 1;
	params->train_fraction = 0.8;
	params->bag_fraction = 0.5;
	params->var_fraction = 1.0;
	params->init_shrinkage = 0.1;
	params->n_minobsinnode = 30;//50;
	params->n_depth = 6;//18;
	GBM* gbm = new GBM(ds,params);

	ds->shuffle(ds->n);
	srand48(7);
	ds->shuffle(ds->n);

/*
	DataSheet* dstest = new DataSheet(ds->n-params->train_fraction*ds->n,ds->p);
	unsigned offset = params->train_fraction*ds->n;
	for(int i=offset;i<ds->n;i++)
	{
		//dstest->y[i-offset] = ds->y[i];
		for(int j=0;j<ds->p;j++)
			dstest->x[j][i-offset] = ds->x[j][i];
	}
*/
	gbm->boost();
	gbm->dump("gbm_ipinyou.gz");
	exit(0);
/*
	gbm->fill(dstest,params->n_iters);
	unsigned err = 0;
	for(int i=offset;i<ds->n;i++)
		if( abs((dstest->y[i-offset]>0.5?1:0)-ds->y[i])>0.1 )
			err += 1;
	printf("%u/%u: %lf\n",err,dstest->n,(double)err/dstest->n);
*/
//*/
//CrossValidate(ds,params,10);
//	GBM* gbm = new GBM("gbm_test.gz");
//	gbm->fill(ds,100);
//	ds->head();
/*
	delete gbm;
	gbm = new GBM("gbm_heart_before.gz");
	double * rele = new double[42];
	unsigned best_iter;
	best_iter = gbm->find_best_iter();
	printf("best_iter: %u\n",best_iter);
	gbm->relevance(rele,best_iter);
	for(int i=0;i<42;i++)
		printf("%s : %lf\n",gbm->ds->var_name[i].c_str(),rele[i]);
*/
	return 0;
}
