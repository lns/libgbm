#include <cstdio>
#include <cmath>
#include "datasheet.h"
#include "hashtable.h"

void quicksort(datatype* x,unsigned* l,unsigned n)
/*
 * Sort l such that x[i] is ascending for i=1:len(l)
 *
 * x is ds->x[j] , the jth variable
 * l is ds->table[j], the jth variable's hash index.
 */
{
	if(n<2)
		return;
	const unsigned m = n/2;
	unsigned tmp = l[m];
	datatype pivot = x[tmp];
	l[m] = l[n-1];
	unsigned *p=l,*q=l+n-1;
	for(;;)
	{
		while(p<q and x[*p]<pivot)
			p++;
		if(p==q)
			break;
		*q = *p;
		q--;
		while(p<q and x[*q]>pivot)
			q--;
		if(p==q)
			break;
		*p = *q;
		p++;
	}
	*p = tmp;
	unsigned *p1=p-1,*p2=p+1;
	while(p1>=l and x[*p1]==x[*p])
		p1--;
	while(p2<=l+n-1 and x[*p2]==x[*p])
		p2++;
	if(p1>l)
		quicksort(x,l,p1+1-l);
	if(p2<l+n-1)
		quicksort(x,p2,l+n-p2);
}

#ifndef NO_NAN
unsigned precipitation(datatype* x,unsigned* l,unsigned n)
/*
 * Push all 'nan' to the end of array
 * Return index of the first nan x
 *
 * x is indeed ds->x[j]
 * l is indeed ds->table[j]
 */
{
	unsigned tmp=l[n-1];
	unsigned *p=l,*q=l+n-1;
	for(;;)
	{
		while(p<q and !isnan(x[*p]))
			p++;
		if(p==q)
			break;
		*q = *p;
		q--;
		while(p<q and isnan(x[*q]))
			q--;
		if(p==q)
			break;
		*p = *q;
		p++;
	}
	*p = tmp;
	return (p-l)+(isnan(x[tmp])?0:1);
}
#endif

HashTable::HashTable(DataSheet* ds,unsigned N)
{
	n = N;
	p = ds->p;
	nan_index = new unsigned[p];
	table = new unsigned*[p];
	printf("Sorting..\n");
	#pragma omp parallel for
	for(int j=0;j<p;j++)
	{
		table[j] = new unsigned[n];
		for(int i=0;i<n;i++)
			table[j][i] = i;
		#ifndef _NO_NAN
		nan_index[j] = precipitation(ds->x[j],table[j],n);
		quicksort(ds->x[j],table[j],nan_index[j]);
		#else
		nan_index[j] = n;
		quicksort(ds->x[j],table[j],n);
		#endif
	}
}

HashTable::~HashTable()
{
	delete [] nan_index;
	for(int j=0;j<p;j++)
		delete [] table[j];
	delete [] table;
}
