#include <cstdlib>
#include <iostream>
#include <iostream>
#include <iterator>
#include <ctime>
#include <cassert>
#include <algorithm>
using namespace std;

template<class T, class Comparator>
unsigned int insertionSort(T *values, unsigned int n, Comparator comparator) {

	unsigned int arrayAccesses = 0;
	//loop invariant:
	//	the range [0, i) is sorted.
	for (unsigned int i = 1; i < n; ++i) {
		T val = values[i];
		unsigned int j = i;
		//loop invariant:
		//	the values in the range (j, i] is greater than "val"
		//  the ranges [0, j ) and (j, i] is sorted.
		for (; j > 0 && (arrayAccesses++, values[j - 1] > val); --j) {
			values[j] = values[j - 1];
		}
		arrayAccesses += ((i - j) * 2 + 1);

		//When the for loop is terminated, there are two cases:
		//         j == 0,  we have reached beginning of the array.
		//         values[j - 1] <= val, since the range [0, i) is sorted we could
		//                  conclude that all the element in the range [0, j - 1] is no greater than "val".
		if (j != i) {
			values[j] = val;
			arrayAccesses++;
		}
	}
	return arrayAccesses;
}
template<class T>
unsigned int insertionSort(T *values, unsigned int n) {
	return insertionSort(values, n, less<T>());
}
template<class T, class Comparator>
unsigned int medianNSelect(T *values, unsigned int n, Comparator comparator) {
	const unsigned int num = 5;
	const unsigned int median = (num / 2);
	unsigned int arrayAccesses = 0;

	while (n > 1) {
		unsigned int g = n / num;
		unsigned int startIndex = 0;
		for (unsigned int i = 0; i < g; ++i, startIndex += num) {
			arrayAccesses += insertionSort(values + startIndex, num,
					comparator);
			values[i] = values[startIndex + median];
		}
		arrayAccesses += 2 * g;

		unsigned int remaining = n - startIndex;
		if (remaining != 0) {
			arrayAccesses += insertionSort(values, remaining, comparator);
			values[g] = values[startIndex + remaining / 2];
			arrayAccesses += 2;
			++g;
		}

		n = g;
	}
	return arrayAccesses;
}
template<class T>
unsigned int medianNSelect(T *values, unsigned int n) {
	return medianNSelect(values, n, less<T>());
}

/**
 * partition the elements in values, into two parts.
 * The partition index "p" with pivot "pivot" is defined as:
 *      the elements in the range [0, p) is no greater than "pivot" and
 *      the elements in the range [p, n) is greater than "pivot".
 * If every elements in the range [0, n) is no greater than "pivot", then "n" is returned,
 *     else the partition index is returned.
 * Note that these implementations are just for experimentations, thus some extra values
 *     are returned, such as the number of array accesses, used by this algorithm.
 */
template<class T, class Comparator>
unsigned int partition(T *values, T pivot, unsigned int n, unsigned int& arrayAccesses,
		Comparator comparator) {
	unsigned int i = 0;
	unsigned int j = 0;
	//loop invariant:
	//      the elements in the range [0, i) is no greater than "pivot"
	//      and the elements in the range [i, j) is greater than "pivot"
	while(j < n){
		arrayAccesses ++;
		if(!comparator(pivot, values[j])){
			if(i != j){
				T tmp = values[i];
				values[i] = values[j];
				values[j] = tmp;
				arrayAccesses += 3;
			}
			++ i;
		}
		++ j;
	}
	return i;
}
/**
 * partition using the default comparator.
 */
template<class T>
unsigned int partition(T *values, T pivot, unsigned int n, unsigned int& arrayAccesses) {
	return partition(values, pivot, n, arrayAccesses, less<T>());
}
/**
 * partition the input array into three ranges using the last element as the pivot element.
 * Assume that the returned index of this procedure is "p".
 * The elements in the first range [0, p) are no greater than the "pivot".
 * The element in the second range [p, p] is equal to the "pivot".
 * The element in the last range (p, n) are greater than the "pivot".
 */
template<class T, class Comparator>
unsigned int partition(T *values, unsigned int n, unsigned int& arrayAccesses,
		Comparator comparator){
	if(n == 1)
		return 0;
	unsigned int p = partition(values, values[n - 1], n - 1, arrayAccesses, comparator);
	if(p != n - 1){
		T tmp = values[p];
		values[p] = values[n - 1];
		values[n - 1] = tmp;
		arrayAccesses += 3;
	}
	return p;
}
/**
 * partition using the default comparator.
 */
template<class T>
unsigned int partition(T *values, unsigned int n, unsigned int& arrayAccesses){
	return partition(values, n, arrayAccesses, less<T>());
}
/**
 * partition using a random value in the input array as the "pivot"
 */
template<class T, class Comparator>
unsigned int randPartition(T *values, unsigned int n, unsigned int& arrayAccesses, Comparator comparator){
	unsigned int index = rand() % n;
	if(index != n -1){
		T tmp = values[index];
		values[index] = values[n - 1];
		values[n - 1] = tmp;
		arrayAccesses += 3;
	}
	return partition(values, n, arrayAccesses, comparator);
}
/**
 * partition using the default comparator.
 */
template<class T>
unsigned int randPartition(T *values, unsigned int n, unsigned int& arrayAccesses){
	return randPartition(values, n, arrayAccesses, less<T>());
}
template<class T, class Comparator>
unsigned int medianNPartiton(T *values, unsigned int n, unsigned int& arrayAccesses, Comparator comparator){
	if(n <= 1)
		return 0;
	arrayAccesses += medianNSelect(values, n, comparator);

	T tmp = values[0];
	values[0] = values[n - 1];
	values[n - 1] = tmp;
	arrayAccesses += 3;
	return partition(values, n, arrayAccesses, comparator);
}
/**
 * partition using the default comparator.
 */
template<class T>
unsigned int medianNPartiton(T *values, unsigned int n, unsigned int& arrayAccesses){
	return medianNPartiton(values, n, arrayAccesses, less<T>());
}
void randInts(unsigned int *values, unsigned int n, unsigned int maxValue);
extern bool isSorted(unsigned int *values, unsigned int n);
template<class T, class Comparator>
bool isPartitionedByIndex(T *values, unsigned int n, unsigned int p, Comparator comparator){
	for(unsigned int i = 0;i < p;++ i)
		if(comparator(values[p], values[i]))
			return false;
	for(unsigned int i = p + 1;i < n;++ i)
		if(!comparator(values[p], values[i]))
			return false;
	return true;
}
template<class T>
bool isPartitionedByIndex(T *values, unsigned int n, unsigned int p){
	return isPartitionedByIndex(values, n, p, less<T>());
}
int testSelectMedian(int argc, char *argv[]) {
	const unsigned int n = 10000;
	unsigned int values[n];
	srand(time(0));
	for (int i = 0; i < 100; ++i) {
		randInts(values, n, n);
		insertionSort(values, n);
		assert(isSorted(values, n));
	}
	for (unsigned int i = 0; i < 100; ++i) {
		randInts(values, n, n * n);
		unsigned int comparisions = medianNSelect(values, n);
		unsigned int median = values[0];
		unsigned int num = (n <= 20 ? 0 : (n * 3 / 10 - 6));
		sort(values, values + n);
		assert(median >= values[num - 1] && median <= values[n - num]);
		cout << "ratio : " << (comparisions / (double) n) << endl;
	}
	return 0;
}
int testPartitionRandom(int argc, char *argv[]){
	const unsigned int n = 1000;
	unsigned int values[n];
	unsigned int array1[n];
	unsigned int array2[n];
	unsigned int array3[n];

	unsigned int maxValue = 1000;

	for(unsigned int i = 0;i < 10;++ i){

		randInts(values, n, maxValue);
		//sort(values, values + n);
		copy(values, values + n, array1);
		copy(values, values + n, array2);
		copy(values, values + n, array3);

		unsigned int arrayAccesses1 = 0, arrayAccesses2 = 0, arrayAccesses3 = 0;
		unsigned int p1 = partition(array1, n, arrayAccesses1);
		unsigned int p2 = randPartition(array2, n, arrayAccesses2);
		unsigned int p3 = medianNPartiton(array3, n, arrayAccesses3);

		assert(isPartitionedByIndex(array1, n, p1));
		assert(isPartitionedByIndex(array2, n, p2));
		assert(isPartitionedByIndex(array3, n, p3));

		cout << "parition :  index - " << p1 << ", number of array accesses - " <<  arrayAccesses1 << endl;
		cout << "random parition :  index - " << p2 << ", number of array accesses - " <<  arrayAccesses2 << endl;
		cout << "median n parition :  index - " << p3 << ", number of array accesses - " <<  arrayAccesses3 << endl;

		cout << endl;
	}
	return 0;
}
int testPartition(int argc, char *argv[]){
	const unsigned int n = 20;
	unsigned int values[n];
	unsigned int arrayAccesses = 0;
	unsigned int p;
	randInts(values, n, 100);
	cout << "before partition :";
	copy(values, values +n, ostream_iterator<int>(cout, " "));
	cout << endl;
	p = partition(values, n, arrayAccesses);
	cout << "partition index : " << p << endl;
	cout << "after partition : ";
	copy(values, values +n, ostream_iterator<int>(cout, " "));
	cout << endl;
	return testPartitionRandom(argc, argv);
}
