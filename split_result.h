#ifndef SPLIT_RESULT_H
#define SPLIT_RESULT_H
#include "node.h"

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

#endif
