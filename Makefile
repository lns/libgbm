VERSION = 3

CXX = g++
CXXFLAGS = -Wall -O3 -fpic -g -D_BERNOULLI -D_NO_NAN

GBMOBJS = node.o datasheet.o hashtable.o parameters.o cart.o gbm.o boost.o

main: libgbm.so main.o
	$(CXX) $(CXXFLAGS) main.o libgbm.so -o main
libgbm.so: $(GBMOBJS)
	$(CXX) -shared -o libgbm.so $(GBMOBJS) -lz

clean: 
	-rm -f *.o
