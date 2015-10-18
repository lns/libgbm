#ifndef DATASHEET_H
#define DATASHEET_H
#include <string>

typedef float datatype;
class DataSheet
{
public:
	std::string file_name;
	unsigned p;
	unsigned n;
	datatype* y;
	datatype** x;
	std::string y_name;
	std::string* var_name;
	double* w;

	void init(unsigned N,unsigned P);
	void init();
	DataSheet(unsigned N,unsigned P);
	DataSheet(char *,char delimiter=',',bool haveY=true);
	DataSheet();
	~DataSheet();
	void shuffle(unsigned n);
	void head(unsigned n=20,int p=9);
	unsigned overwrite(unsigned index,
						DataSheet* from_ds,unsigned start,unsigned end);
	void read_csv(char *,char delimiter=',',bool haveY=true);
	void write_csv(char *,char delimiter=',');
	void dump(char * file_name);
	void load(char * file_name);
};
#endif
