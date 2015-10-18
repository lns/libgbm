#ifndef HASHTABLE_H
#define HASHTABLE_H
#include "datasheet.h"

class HashTable
{
public:
	unsigned p;
	unsigned n;
	unsigned * nan_index;
	unsigned** table;

	HashTable(DataSheet* ds,unsigned n);
	~HashTable();
};

#endif
