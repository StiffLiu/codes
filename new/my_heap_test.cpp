#include "my_heap.h"
#include "my_min_max_heap.h"
#include "my_median_heap.h"
#include "my_math.h"
#include <cassert>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <random>

#define ARSIZE(ar) (sizeof(ar) / sizeof (*(ar)))

static int validateHeap(int argc, char *argv[]){
	using namespace my_lib;
	std::random_device rd;
	typedef int Type;
	std::default_random_engine dre(rd());
	std::uniform_int_distribution<Type> urd(0, 100);
	const auto count = 100;
	Type values[count];
	Type sorted[count];

	for(auto i = 0;i < count;++ i) 
		values[i] = urd(dre);

	std::copy(values, values + count, sorted);
	std::sort(sorted, sorted + count);
	
	auto d = 3;
	NaryHeap<Type> heap(d, values, values + count);
	NaryMinMaxHeap<Type> minMaxHeap(d, values, values + count);
	NaryMedianHeap<Type> medianHeap(d, values, values + count);

	assert(heap.size() == count);
	assert(minMaxHeap.size() == count);
	assert(medianHeap.size() == count);
	assert(heap.top() == sorted[0]);

	for(auto i = 0;i < count;++ i){
		assert(heap.top() == sorted[i]);
		heap.pop();
	}
	assert(heap.empty());

	assert(minMaxHeap.min() == sorted[0]);
	assert(minMaxHeap.max() == sorted[count - 1]);
	for(auto i = 0, j = count - 1;i <= j;++ i, -- j){
		assert(minMaxHeap.min() == sorted[i]);
		assert(minMaxHeap.max() == sorted[j]);
		minMaxHeap.popMin();
		minMaxHeap.popMax();
		assert(minMaxHeap.isValid());
	}
	assert(minMaxHeap.empty());
	
	auto start1 = (count - 1) / 2, start2 = count / 2;
	if(start1 == start2){
		assert(medianHeap.top() == sorted[start1]);
		medianHeap.pop();
		--start1, ++start2;
	}
	//copy(sorted, sorted + count, std::ostream_iterator<Type>(std::cout, " "));
	//std::cout << std::endl;
	assert(medianHeap.isValid());
	for(;start2 < count;-- start1, ++ start2){
		//std::cout << medianHeap.size() << std::endl;
		//std::cout << medianHeap.top() << ',' << sorted[start2] << std::endl;
		assert(medianHeap.top() == sorted[start2]);
		medianHeap.pop();
		assert(medianHeap.top() == sorted[start1]);
		medianHeap.pop();
		assert(medianHeap.isValid());	
	}
	assert(medianHeap.empty());
	return 0;
}
static int testDiscreteDist(int argc, char *argv[]){
	double probabilities[] = {0.1, 0.2, 0.3, 0.4, 0.5};
	auto count = ARSIZE(probabilities);
	my_lib::DiscreteDistribution dd(probabilities, probabilities + count);
	int counts[count];
	auto iteration = 1000000;
	double sum = std::accumulate(probabilities, probabilities + count, 0.0);
	for(decltype(count) i = 0;i < count;++ i){
		counts[i] = 0;
	}
	for(auto i = 0;i < iteration;++ i){
		++ counts[dd.next()];
	}
	for(decltype(count) i = 0;i < count;++ i){
		std::cout << counts[i] << "\texpected : " 
			<< (int)(iteration * probabilities[i] / sum) << std::endl;
	}
	return 0;
}
static int outputHeap(int argc, char *argv[]){
	using namespace my_lib;
	using namespace std;
	std::initializer_list<double> l = {10.0, 4.0, 5.0, 8.0, 3.0, 4.0};
	NaryHeap<double > heap(2, l);
	NaryMinMaxHeap<double> minMaxHeap(3, l);
	NaryMedianHeap<double> medianHeap(3, l);

	assert(minMaxHeap.isValid());
	assert(medianHeap.isValid());
	cout << "all keys : ";	
	copy(l.begin(), l.end(), ostream_iterator<double>(cout, " "));
	cout << endl;
	while(!minMaxHeap.empty()){
		cout << "min, max : " <<  minMaxHeap.min() << ", " << minMaxHeap.max() << endl;
		minMaxHeap.popMin();
		minMaxHeap.popMax();
		assert(minMaxHeap.isValid());
	}
	while(!medianHeap.empty()){
		cout << "median : " << medianHeap.top() << endl;
		medianHeap.pop();
		assert(medianHeap.isValid());
	}
	return 0;
}
int heapTest(int argc, char *argv[]){
	for(auto i = 0;i < 100;++ i)
	validateHeap(argc, argv);
	outputHeap(argc, argv);
	testDiscreteDist(argc, argv);
	return 0;
}
