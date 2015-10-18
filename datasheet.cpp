#include "datasheet.h"
#include "zlib.h"
#include <cstdio>
#include <cmath>
#include <cassert>
#include <cstdlib>

void DataSheet::init(unsigned N,unsigned P)
{
	n = N;
	p = P;
	y = new datatype[N];
	assert(y!=NULL);
	x = new datatype*[P];
	assert(x!=NULL);
	file_name = "";
	y_name = "Y";
	var_name = new std::string[P];
	assert(var_name!=NULL);
	for(int j=0;j<P;j++)
	{
		x[j] = new datatype[N];
		assert(x[j]!=NULL);
	}
	w = new double[N];
	assert(w!=NULL);
}

void DataSheet::init()
{
	n = 0;
	p = 0;
	y = NULL;
	x = NULL;
	w = NULL;
	file_name = "";
	y_name = "Y";
	var_name = NULL;
}

DataSheet::DataSheet(unsigned N,unsigned P)
{
	init(N,P);
}

DataSheet::DataSheet()
{
	init();
}

DataSheet::~DataSheet()
{
	if(y!=NULL)
		delete[] y;
	if(x!=NULL)
	{
		for(int j=0;j<p;j++)
		{
			if(x[j]!=NULL)
			delete[] x[j];
		}
		delete[] x;
	}
	if(w!=NULL)
		delete[] w;
	if(var_name!=NULL)
		delete[] var_name;
}

void DataSheet::shuffle(unsigned n)
{
	for(unsigned i=0;i<n;i++)
	{
		unsigned m = floor(drand48()*(n-1-i));
		if(m==0)
			continue;
		datatype tmp;
		double tmp_w;
		/* Swap y */
		tmp    = y[i];
		y[i]   = y[i+m];
		y[i+m] = tmp;
		/* Swap w */
		tmp_w    = w[i];
		w[i]   = w[i+m];
		w[i+m] = tmp_w;
		/* Swap x */
		for(int j=0;j<p;j++)
		{
			tmp       = x[j][i];
			x[j][i]   = x[j][i+m];
			x[j][i+m] = tmp;
		}
	}
}

void DataSheet::head(unsigned N,int P)
{
	printf("%s\n",file_name.c_str());
	if(0>P or P>p)
		P = p;
	if(N>n)
		N = n;
	printf("%8s",y_name.c_str());
	for(int j=0;j<P;j++)
	{
		if(var_name!=NULL)
			printf(" \t%8s",var_name[j].c_str());
		else
			printf(" \t  X[%d]  ",j);
	}
	printf("\n");
	for(int i=0;i<N;i++)
	{
		printf("%8le",y[i]);
		for(int j=0;j<P;j++)
			printf(" \t%8le",x[j][i]);
		printf("\n");
	}
}

unsigned DataSheet::overwrite(unsigned index,
						DataSheet* from_ds,unsigned start,unsigned end)
/*
 * Overwrite self[index:index+(end-start)] with data from_ds[start:end]
 * Return index+(end-start)
 */
{
	assert(index+end-start<=n and end<=from_ds->n);
	unsigned idx = index;
	for(int i=start;i<end;i++)
	{
		y[idx] = from_ds->y[i];
		for(int j=0;j<p;j++)
			x[j][idx] = from_ds->x[j][i];
		idx++;
	}
	return idx;
}

void DataSheet::read_csv(char * input_file_name,char delimiter,bool haveY)
{
	assert(y==NULL and x==NULL and w==NULL and var_name==NULL);
	assert(access(input_file_name,R_OK)==0);
	FILE* f = fopen(input_file_name,"r");
	file_name = std::string(input_file_name);
	size_t buffer_size = 32768;
	char* buffer = (char*)malloc(sizeof(char)*buffer_size);
	assert(buffer!=NULL);
	unsigned n_line = 0;
	unsigned n_col = 1;
	while(!feof(f))
	{
		getline(&buffer,&buffer_size,f);
		n_line++;
		//Now we are at file lineno == i
		if(n_line==1)
		{
			char* q = buffer;
			while(*q!='\n')
			{
				if(*q==delimiter)
					n_col++;
				q++;
			}
			if(haveY)
				p = n_col-1;
			else
				p = n_col;
			printf("File has %u columns, p=%u\n",n_col,p);
		}
	}
	n_line--;
	printf("File has %u lines, n=%u\n",n_line,n_line-1);
	init(n_line-1,p);
	assert(var_name!=NULL);
	printf("Assuming file has data head\n");
	fseeko(f,0,SEEK_SET);
	int i,j;
	//read var_name
	getline(&buffer,&buffer_size,f);
	char* q = buffer;
	if(haveY)
	{
		while(*q!=delimiter)
			q++;
		*q = '\0';
		y_name = std::string(buffer);
		q++;
	}
	for(j=0;j<p-1;j++)
	{
		char* start = q;
		while(*q!=delimiter)
			q++;
		*q = '\0';
		var_name[j] = std::string(start);
		q++;
	}
	char* start = q;
	while(*q!='\n')
		q++;
	*q = '\0';
	var_name[p-1] = std::string(start);
	//read data
	for(i=0;i<n_line-1;i++)
	{
		getline(&buffer,&buffer_size,f);
		char* q = buffer;
		if(haveY)
		{
			while(*q!=delimiter)
				q++;
			*q = '\0';
			y[i] = atof(buffer);
			q++;
		}
		for(j=0;j<p-1;j++)
		{
			char* start = q;
			while(*q!=delimiter)
				q++;
			*q = '\0';
			if(start==q)
				x[j][i] = NAN;
			else
				x[j][i] = atof(start);
			q++;
		}
		char* start = q;
		while(*q!='\n')
			q++;
		*q = '\0';
		x[p-1][i] = atof(start);
	}
	free(buffer);
}

DataSheet::DataSheet(char * input_file_name,char delimiter,bool haveY)
{
	init();
	read_csv(input_file_name,delimiter,haveY);
}

void DataSheet::write_csv(char * output_file_name,char delimiter)
{
	FILE* f = fopen(output_file_name,"w");
	fprintf(f,y_name.c_str());
	for(int j=0;j<p;j++)
	{
		fprintf(f,"%c",delimiter);
		fprintf(f,"%s",var_name[j].c_str());
	}
	fprintf(f,"\n");
	for(int i=0;i<n;i++)
	{
		fprintf(f,"%le",y[i]);
		for(int j=0;j<p;j++)
		{
			fprintf(f,"%c",delimiter);
			fprintf(f,"%le",x[j][i]);
		}
		fprintf(f,"\n");
	}
	fclose(f);
}

void DataSheet::dump(char * dump_file_name)
/*
 * Dump DataSheet to a Gzipped file
 */
{
	gzFile f = gzopen(dump_file_name,"w");
	gzbuffer(f,131072);
	//ds.file_name
	unsigned uint = file_name.length();
	gzwrite(f,&uint,sizeof(unsigned));
	gzwrite(f,&(file_name[0]),uint*sizeof(char));
	//ds.n, ds.p
	gzwrite(f,&p,sizeof(unsigned));
	gzwrite(f,&n,sizeof(unsigned));
	//ds.y
	if(y==NULL)
	{
		uint = 0;
		gzwrite(f,&uint,sizeof(unsigned));
	}else{
		uint = 1;
		gzwrite(f,&uint,sizeof(unsigned));
		gzwrite(f,y,n*sizeof(datatype));
	}
	//ds.x
	assert(x!=NULL);
	for(int i=0;i<p;i++)
	{
		if(x[i]==NULL)
		{
			uint = 0;
			gzwrite(f,&uint,sizeof(unsigned));
		}else{
			uint = 1;
			gzwrite(f,&uint,sizeof(unsigned));
			gzwrite(f,x[i],n*sizeof(datatype));
		}
	}
	//ds.y_name
	uint = y_name.length();
	gzwrite(f,&uint,sizeof(unsigned));
	gzwrite(f,&(y_name[0]),uint*sizeof(char));
	//ds.var_name
	for(int i=0;i<p;i++)
	{
		uint = var_name[i].length();
		gzwrite(f,&uint,sizeof(unsigned));
		gzwrite(f,&(var_name[i][0]),uint*sizeof(char));
	}
	//ds.w
	if(w==NULL)
	{
		uint = 0;
		gzwrite(f,&uint,sizeof(unsigned));
	}else{
		uint = 1;
		gzwrite(f,&uint,sizeof(unsigned));
		gzwrite(f,w,n*sizeof(w[0]));
	}
	//
	gzclose(f);
	return;
}

void DataSheet::load(char * load_file_name)
/*
 * Load Gzipped DataSheet
 */
{
	assert(access(load_file_name,R_OK)==0);
	gzFile f = gzopen(load_file_name,"r");
	gzbuffer(f,131072);
	unsigned uint;
	char buffer[512];
	//file_name
	gzread(f,&uint,sizeof(uint));
	gzread(f,&buffer,uint*sizeof(char));
	file_name = std::string(buffer,uint);
	// n,p
	gzread(f,&p,sizeof(p));
	gzread(f,&n,sizeof(n));
	// Initialize
	init(n,p);
	// y
	gzread(f,&uint,sizeof(unsigned));
	if(uint==1)
		gzread(f,y,n*sizeof(datatype));
	else
	{
		delete[] y;
		y = NULL;
	}
	// x
	for(int i=0;i<p;i++)
	{
		gzread(f,&uint,sizeof(unsigned));
		if(uint==1)
			gzread(f,x[i],n*sizeof(datatype));
		else
		{
			delete[] x[i];
			x[i] = NULL;
		}
	}
	// y_name
	gzread(f,&uint,sizeof(unsigned));
	gzread(f,&buffer,uint*sizeof(char));
	y_name = std::string(buffer,uint);
	// var_name
	for(int i=0;i<p;i++)
	{
		gzread(f,&uint,sizeof(unsigned));
		gzread(f,&buffer,uint*sizeof(char));
		var_name[i] = std::string(buffer,uint);
	}
	// w
	gzread(f,&uint,sizeof(unsigned));
	if(uint==1)
		gzread(f,w,n*sizeof(w[0]));
	else
	{
		delete[] w;
		w = NULL;
	}
	//
	gzclose(f);
	return;
}

