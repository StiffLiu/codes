#include "my_rand_string_generator.h"
#include "my_list_symbol_table.h"
#include "my_bin_search_symbol_table.h"
#include "my_access_counted_object.h"
#include "my_test_base.h"
#include <set>
#include <ctime>
#include <vector>
#include <iostream>
#include <iterator>
#include <vector>
#include <fstream>

namespace{
using namespace std;
using namespace my_lib;
#define ARSIZE(ar) (sizeof (ar) / sizeof *(ar))
class TestPlot : public my_lib::StatPlot<TestPlot&>{
	typedef my_lib::StatPlot<TestPlot&> Super;
	unsigned int start = 100;
	double totalTime1 = 0;
	double totalOp1 = 0;
	double totalTime2 = 0;
	double totalOp2 = 0;
public:
	TestPlot() : Super(4, *this){
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
			0.0, 0.0, 1.0,
			1.0, 1.0, 0.0,
		};
		colors.assign(c, c + ARSIZE(c));
	}	
	bool operator()(double *values){
		vector<string> strs;
		auto func = [&strs](size_t i1, size_t i2){
			return strs[i1] < strs[i2];
		};
		set<unsigned int, decltype(func)> strSet(func);
		RandStringGenerator gen(2, 50);
		
		while(strs.size() < start){
			strs.push_back(string());	
			gen.randStringAlphaNum(strs.back());
			if(!strSet.insert(strs.size() - 1).second)
				strs.pop_back();
			
		}

		clock_t begin = 0;
		clock_t end = 0;
		ListSymbolTable<string, bool> listST;
		BinSearchSymbolTable<string, bool> binSearchST;
		unsigned int half = start / 2; 
		double ops = 0;
		double duration = 0;

		begin = clock();
		for(unsigned int i = 0;i < half;++ i){
			listST.put(strs[i], false);	
		}
		for(auto i = 0;i < 10;++ i){
			for(unsigned int j = 0;j < half;++ j){
				listST.get(strs[j]);
				listST.get(strs[j + half]);
			}
		}
		end = clock();
		
		duration = end - begin;
		ops = 21 * half; 
		totalTime1 += duration;
		totalOp1 += ops;
		values[0] = start;
		values[1] = duration / ops;
		values[2] = start;
		values[3] = totalTime1 / totalOp1;

		begin = clock();
		for(unsigned int i = 0;i < half;++ i){
			binSearchST.put(strs[i], false);	
		}
		for(auto i = 0;i < 10;++ i){
			for(unsigned int j = 0;j < half;++ j){
				binSearchST.get(strs[j]);
				binSearchST.get(strs[j + half]);
			}
		}
		end = clock();
		
		duration = end - begin;
		ops = 21 * half; 
		totalTime2 += duration;
		totalOp2 += ops;
		values[4] = start;
		values[5] = duration / ops;
		values[6] = start;
		values[7] = totalTime2 / totalOp2;

		start += 2;
		return true;
	}

};
template<template<class U, class V> class T>
class TestOne  : public my_lib::StatPlot<TestOne<T>&> {
	typedef my_lib::StatPlot<TestOne<T>&> Super;
	unsigned int start = 100;
	unsigned int totalIteration = 10;
	unsigned int currentIteration = 0;
	double totalTime = 0;
	double totalOp = 0;
public:
	TestOne() : Super(2, *this){
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
		};
		Super::colors.assign(c, c + ARSIZE(c));
		Super::maxNumPoints = 0;
		Super::interval = 1;
	}
	bool operator()(double *values){
		vector<string> strs;
		auto func = [&strs](size_t i1, size_t i2){
			return strs[i1] < strs[i2];
		};
		set<unsigned int, decltype(func)> strSet(func);
		RandStringGenerator gen(2, 50);
		
		while(strs.size() < start){
			strs.push_back(string());	
			gen.randStringAlphaNum(strs.back());
			if(!strSet.insert(strs.size() - 1).second)
				strs.pop_back();
			
		}

		clock_t begin = 0;
		clock_t end = 0;
		T<string, bool> st;
		unsigned int half = start / 2; 
		double ops = 0;
		double duration = 0;

		begin = clock();
		for(unsigned int i = 0;i < half;++ i){
			st.put(strs[i], false);	
		}
		for(unsigned int j = 0;j < half;++ j){
			st.get(strs[j]);
			st.get(strs[j + half]);
		}
		end = clock();

		duration = end - begin;
		++currentIteration;
		if(currentIteration >= totalIteration){
			currentIteration = 0;
			totalTime = 0;
			totalOp = 0;
			++ start;
		}else{
			if(!Super::graphs[1].empty()){
				Super::graphs[1].pop_back();
				Super::graphs[1].pop_back();
			}
		}

		ops = 3 * half;
		values[0] = start;
		values[1] = duration / ops;
		totalTime += duration;
		totalOp += ops;
		values[2] = start;
		values[3] = totalTime / totalOp;
		return true;
	}
};
template<template<class U, class V> class T>
class AverageCost : public my_lib::StatPlot<AverageCost<T>&> {
	typedef my_lib::StatPlot<AverageCost<T>&> Super;
	unsigned int start = 100;
	unsigned int totalIteration = 10;
	unsigned int currentIteration = 0;
	double totalTime = 0;
	double totalOp = 0;
	size_t current = 0;
	std::vector<clock_t> times;
public:
	AverageCost() : Super(2, *this){
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
		};
		Super::colors.assign(c, c + ARSIZE(c));

		if(false){
			std::vector<std::string> vec;
			if(false){
				std::ifstream in("/development/documents/books/algos4/algs4.cs.princeton.edu/31elementary/tale.txt");
				std::copy(std::istream_iterator<std::string>(in), 
					std::istream_iterator<std::string>(),
					std::back_inserter(vec));
			}else{
				auto count = (1 << 16);
				RandStringGenerator gen(5, 10);
				cout << "generating " << count << " strings with length between 5 and 10 " << endl;
				for(auto i = 0;i < count;++ i){
					vec.push_back(string());
					gen.randStringPrintable(vec.back());
					//cout << vec.back() << endl;
				}
				sort(vec.begin(), vec.end());
				reverse(vec.begin(), vec.end());
			}

			T<std::string, bool> st;
			times.resize(vec.size());
			//clock_t start = clock();
			for(size_t i = 0;i < vec.size();++ i){
				st.put(vec[i], false);
				times[i] = st.getComparator().getCounter();//clock() - start;
			}
			std::cout << "total operations : " << vec.size() << std::endl;
			std::cout << "total keys : " << st.size() << std::endl;
			std::cout << "total time consumed : " << times.back() << std::endl;
			Super::maxNumPoints = vec.size();
		}else{
			std::vector<unsigned int> vec;
			auto count = (1 << 16);
			T<unsigned int, bool> st;
			for(auto i = 0;i < count;++ i)
				vec.push_back(i);
			//random_shuffle(vec.begin(), vec.end());
			reverse(vec.begin(), vec.end());
			times.resize(vec.size());
			//clock_t start = clock();
			for(size_t i = 0;i < vec.size();++ i){
				st.put(vec[i], false);
				times[i] = st.getComparator().getCounter();//clock() - start;
			}
			std::cout << "total operations : " << vec.size() << std::endl;
			std::cout << "total keys : " << st.size() << std::endl;
			std::cout << "total time consumed : " << times.back() << std::endl;
			Super::maxNumPoints = vec.size();

		}
		Super::interval = 1;
	}
	bool operator()(double *values){
		auto i = current;
		if(i >= times.size())
			return false;
		values[0] = values[2] = i;
		values[1] = (i == 0 ? times[i] : times[i] - times[i - 1]);
		values[3] = times[i] / (double)(i + 1);
		++ current;
		return true;
	}
};
	
template<class K, class V>
class DefaultListST : public ListSymbolTable<K, V, 
	std::list<std::pair<K, V> >, AccessCountedFunctor<unsigned long, std::less<K> > >{
};
template<class K, class V>
class DefaultBinST : public BinSearchSymbolTable<K, V>{
};
struct Node{
	Node *l = nullptr, *r = nullptr;
};
struct Func{
	bool left(Node *n, Node*& l){
		if(n->l == nullptr)
			return false;
		l = n->l;
		return true;
	}
	bool right(Node *n, Node*& r){
		if(n->r == nullptr)
			return false;
		r = n->r;
		return true;
	}
};

Node *randTree(unsigned int n){
	if(n == 0)
		return nullptr;
	if(n == 1)
		return new Node();
	n -= 1;
	unsigned int t = rand() % n;
	Node *node = new Node();
	node->l = randTree(t);
	node->r = randTree(n - t);
	return node;
}
}
int main(int argc, char *argv[]){
	//TestOne<DefaultBinST> testPlot;
	//AverageCost<DefaultListST> testPlot;	
	TreePlot testPlot;
	testPlot.calculate(randTree(100), Func());
	testPlot.run(argc, argv);
	return 0;
}
