#include "node.h"

Node::Node()
{
	left  = NULL;
	right = NULL;
	#ifndef _NO_NAN
	nan   = NULL;
	#endif
	
	var   = ~0;
	val   = NAN;
	er    = -INF;
	sumy  = 0;
	obs   = 0;
}

#include <cstdio>
#include <cassert>
void Node::print(int depth)
{
	for(int i=0;i<depth;i++)
		printf("    ");
	if(var!=~0)
	{
		printf("Var[%5u]<>%8le ",var,val);
		printf("ER:%8le ",er);
	}else{
		printf("Pred:%8le ",val);
		assert(er==0);
	}
	printf("sumz: %le obs: %u\n",sumy,obs);
	if(left!=NULL)
		left->print(depth+1);
	if(right!=NULL)
		right->print(depth+1);
	#ifndef _NO_NAN
	if(nan!=NULL)
		nan->print(depth+1);
	#endif
}
