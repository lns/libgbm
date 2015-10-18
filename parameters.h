#ifndef PARAMETERS_H
#define PARAMETERS_H

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

#endif
