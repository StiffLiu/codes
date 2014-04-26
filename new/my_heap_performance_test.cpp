#include "my_test_base.h"
#include "my_heap.h"
#include "my_access_counted_object.h"
#include <iostream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <vector>
#include <cstdlib>

#define ARSIZE(a) (sizeof(a) / (sizeof *(a)))
namespace{
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
struct NotUseFloyd : public my_lib::NaryHeapAlgorithmTraits{
	static const bool useFloydMethod = false;
};
class DiffNaryHeapCompare : public my_lib::StatPlot<DiffNaryHeapCompare&>{
	typedef my_lib::StatPlot<DiffNaryHeapCompare&> Super;
	unsigned int start = 100;
	unsigned int dimensions[4] = {2, 3, 5, 7};
public:
	DiffNaryHeapCompare(): Super(2 * ARSIZE(dimensions), *this){
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
			0.0, 0.0, 1.0,
			1.0, 1.0, 0.0,
			
			0.5, 0.0, 0.0,
			0.0, 0.5, 0.0,
			0.0, 0.0, 0.5,
			0.5, 0.5, 0.0,
		};
		colors.assign(c, c + ARSIZE(c));
	}
	bool operator()(double *values){
		using namespace my_lib;
		std::vector<unsigned int> vals(start);
		for(decltype(start) i = 0;i < start;++ i)
			vals[i] = i;
		std::random_shuffle(vals.begin(), vals.end());

		typedef AccessCountedFunctor<unsigned int, 
			std::less<unsigned int> > LessThan;
		typedef NaryHeap<unsigned int, std::vector<unsigned int>,
		       LessThan> Heap;	
		auto n = ARSIZE(dimensions);
		for(decltype(n) i = 0;i < n;++ i){
			Heap heap(dimensions[i], vals.begin(), vals.end()); 
			while(!heap.empty())
				heap.pop();

			auto index = 2 * i;
			values[index] = start;
			values[index + 1] = heap.getComparator().getCounter();
		}	
		typedef NaryHeap<unsigned int, std::vector<unsigned int>,
		       LessThan, NotUseFloyd> HeapNotFloyd;	
		for(decltype(n) i = 0;i < n;++ i){
			HeapNotFloyd heap(dimensions[i], vals.begin(), 
				vals.end()); 
			while(!heap.empty())
				heap.pop();

			auto index = 2 * i + 8;
			values[index] = start;
			values[index + 1] = heap.getComparator().getCounter();
		}
		start += 10;
		return true;
	}
};
class ConstructionTimePercentage : public my_lib::StatPlot<ConstructionTimePercentage&>{
	typedef my_lib::StatPlot<ConstructionTimePercentage&> Super;
	unsigned int start = 1000;
	unsigned int dimensions[4] = {2, 3, 5, 7};
public:
	ConstructionTimePercentage(): Super(ARSIZE(dimensions), *this){
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
			0.0, 0.0, 1.0,
			1.0, 1.0, 0.0,
		};
		colors.assign(c, c + ARSIZE(c));
	}
	bool operator()(double *values){
		using namespace my_lib;
		std::vector<unsigned int> vals(start);
		for(decltype(start) i = 0;i < start;++ i)
			vals[i] = i;
		std::random_shuffle(vals.begin(), vals.end());

		typedef AccessCountedFunctor<unsigned int, 
			std::less<unsigned int> > LessThan;
		typedef NaryHeap<unsigned int, std::vector<unsigned int>,
		       LessThan> Heap;	
		auto n = ARSIZE(dimensions);
		for(decltype(n) i = 0;i < n;++ i){
			
			Heap heap(dimensions[i], vals.begin(), vals.end()); 
			double construction = heap.getComparator().getCounter(); 
			while(!heap.empty())
				heap.pop();
			auto index = 2 * i;
			values[index] = start;
			values[index + 1] = construction / heap.getComparator().getCounter();
		}	
		start += 10;
		return true;
	}
};
class PathologicalCases : public my_lib::StatPlot<PathologicalCases&>{
	typedef my_lib::StatPlot<PathologicalCases&> Super;
	unsigned int start = 1000;
	unsigned int dimension = 3;
public:
	PathologicalCases(): Super(4, *this){
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
			0.0, 0.0, 1.0,
			1.0, 1.0, 0.0,
		};
		colors.assign(c, c + ARSIZE(c));
	}
	bool operator()(double *values){
		using namespace my_lib;

		typedef AccessCountedFunctor<unsigned int, 
			std::less<unsigned int> > LessThan;
		typedef NaryHeap<unsigned int, std::vector<unsigned int>,
		       LessThan, NotUseFloyd> Heap;	

		Heap heap1(dimension), heap2(dimension), heap3(dimension), heap4(dimension);
		
		for(decltype(start) i = 0;i < start;++ i)
			heap1.add(i);
		for(decltype(start) i = start;i > 0;-- i)
			heap2.add(i);
		for(decltype(start) i = 0;i < start;++ i)
			heap3.add(0);
		for(decltype(start) i = 0;i < start;++ i)
			heap4.add(i > start / 2 ? 0 : 1);

		values[0] = start;
		values[1] = heap1.getComparator().getCounter();
		values[2] = start;
		values[3] = heap2.getComparator().getCounter();
		values[4] = start;
		values[5] = heap3.getComparator().getCounter();
		values[6] = start;
		values[7] = heap4.getComparator().getCounter();
		
		start += 10;
		return true;
	}

};
}
int main(int argc, char *argv[]){
	PathologicalCases test;
	test.run(argc, argv);
	return 0;
}
