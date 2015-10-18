#include <cstdio>
#include "parameters.h"

void Parameters::print()
{
    printf("Parameters:\n");
    printf("         n_iters: %6u\n",n_iters);
    printf("         n_bags : %6u\n",n_bags);
    printf("    var_fraction: %6.3lf\n",var_fraction);
    printf("    bag_fraction: %6.3lf\n",bag_fraction);
    printf("  train_fraction: %6.3lf\n",train_fraction);
    printf("  init_shrinkage: %8.4lf\n",init_shrinkage);
    printf("  n_minobsinnode: %6u\n",n_minobsinnode);
    printf("         n_depth: %6u\n",n_depth);
}

