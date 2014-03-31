#include "my_test_base.h"
#include "my_heap.h"
#include <iostream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <vector>
#include <cstdlib>

class TestPlot : public my_lib::StatPlot<TestPlot&>{
	typedef my_lib::StatPlot<TestPlot&> Super;
	unsigned int start = 100;
	unsigned int d = 3;
public:
	TestPlot() : Super(1, *this){
		srand(time(0));
		maxNumPoints = 10000;
	}
	bool operator()(double *values){
		my_lib::NaryHeap<unsigned int> heap(d);
		clock_t begin = clock();
		for(unsigned int i = 0;i < 10;++ i){
			for(unsigned int i = 0;i < start;++ i)
				heap.add(rand());
			for(unsigned int i = 0;i < start / 2;++ i)
				heap.pop();
			for(unsigned int i = 0;i < start;++ i)
				heap.add(rand());
			for(unsigned int i = 0;i < start;++ i)
				heap.pop();
		}
		values[0] = start;
		values[1] = (clock() - begin);
		start += 10;
		return true;
	}
};
int main(int argc, char *argv[]){
	TestPlot test;
	test.run(argc, argv);
	return 0;
}
