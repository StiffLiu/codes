#include <random>
#include <thread>
#include <vector>
#include <iostream>
#include <algorithm>

template<class T, class U = T, class Comparator>
void merge(T ar1, unsigned int n1, T ar2, unsigned int n2, U output, Comparator comparator){
	decltype(n1) i = 0;
	decltype(n2) j = 0;
	auto count = n1 + n2;
	decltype(count) k = 0;
	while(i < n1 && j < n2){
		if(comparator(ar2[j], ar1[i])){
			output[k] = ar2[j];
			++ j;
		}else{
			output[k] = ar1[i];
			++ i;
		}
		++ k;
	}
	while(i < n1)
		output[k ++] = ar1[i ++ ];
	while(j < n2)
		output[k ++] = ar2[j ++ ];

}

template<class T, bool optimize, class U = T, class Comparator>
void mergeSort(T ar, unsigned int n, U output, Comparator comparator){
	if(n <= 1){
		output[0] = ar[0];
		return;
	}

	decltype(n) n1 = n / 2;
	decltype(n) n2 = n - n1;
	T ar2 = ar + n1;
	mergeSort(ar, n1, output, comparator);
	mergeSort(ar2, n2, output, comparator);
	if(optimize && comparator(ar[n1], ar[n1 - 1])){
		for(decltype(n) i = 0;i < n;++ i)
			output[i] = ar[i];
		merge(output, n1, output, n2, ar, comparator);
	}
}
template<class T, class U = T, class Comparator>
void bottomUpMergeSort(T ar, unsigned int n, U output, Comparator comparator){
	bool isArSorted = true;
	for(unsigned int sz = 1;sz < n;sz += sz){
		unsigned int cnt = sz + sz;
		unsigned int j = sz;
		while(j < n){
			merge(ar + (j - cnt), sz, ar + (j - sz), sz, output, comparator);
			j += cnt;
		}
		merge(ar + (j - cnt), sz, ar + (j - sz), n - j + sz, output, comparator);
		std::swap(ar, output);
		isArSorted = (!isArSorted);
	}
	if(!isArSorted){
		for(decltype(n) i = 0;i < n;++ i)
			ar[i] = output[i];
	}
}
template<class T, class Comparator>
unsigned int firstIndexOfDisorder(T ar, unsigned int n, Comparator comparator){
	for(decltype(n) i = 1;i < n;++ i){
		if(comparator(ar[i], ar[i - 1]))
			return i;
	}
	return 0;
}

int testMergeSort(int argc, char *argv[]){
	unsigned int count;
	return 0;
}

