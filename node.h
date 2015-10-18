#ifndef NODE_H
#define NODE_H
#include "defs.h"
#include "datasheet.h"

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
	void print(int depth=0);
};

#endif
