/**
 * This file is for the test, visualization, comparison of some basic sorting algorithms,
 * that is selection sort, insertion sort, shell sort and merge sort.
 *
 * The algorithms in this file comes from the book "algorithms 4th edition", section 2.1.
 * Many of the codes are solutions to the exercise of section 2.1.
 *
 * Features of c++11 is applied in this source file,
 * to compile it, you have to set some compiler options.
 * For g++, you have to use "-std=c++11"
 */
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <random>
#include <algorithm>
#include <iostream>
#include <random>
#include <thread>
#include <cmath>
#include <GL/glut.h>
#include <mutex>
#include <cassert>
#include <condition_variable>
#include <unordered_map>
#include <map>
#include <chrono>
#include <tuple>
#include <queue>
#include <iterator>
#include <cstring>
#include "common.h"
using namespace std;
/**
 * This template class is for performance test purpose.
 * It's a wrapper for arrays and it counts the number
 * of array accesses. This number is used to measure
 * the performance of an algorithm.
 */
template<class T, class Counter = unsigned int *>
class AccessCountedArray {
	T *val;
	Counter counter;
public:
	AccessCountedArray(T *val = nullptr, Counter counter = Counter()) :
			val(val), counter(counter) {
	}
	T& operator[](unsigned int index) {
		if (counter != NULL)
			++*counter;
		return val[index];
	}
	AccessCountedArray& operator++() {
		++val;
		return *this;
	}
	T& operator*() {
		return *val;
	}
	AccessCountedArray operator++(int) {
		AccessCountedArray old = *this;
		++val;
		return old;
	}
	bool operator==(AccessCountedArray array) {
		return val == array.val;
	}
	bool operator!=(AccessCountedArray array) {
		return val != array.val;
	}
	int operator-(AccessCountedArray array) {
		return val - array.val;
	}
	operator T*() {
		return val;
	}
	operator const T*() {
		return val;
	}
	AccessCountedArray operator+(unsigned int n) {
		return {val + n, counter};
	}
	AccessCountedArray operator-(unsigned int n) {
		return {val - n, counter};
	}
	const T& operator[](unsigned int index) const {
		if (counter != NULL)
			++*counter;
		return val[index];

	}
};
template<class T>
struct CountedLessThan {
	unsigned long long *counter;
	CountedLessThan(unsigned long long * counter = nullptr) :
			counter(counter) {
	}
	bool operator()(T val1, T val2) {
		if (counter != nullptr)
			++*counter;
		return val1 < val2;
	}
};
/**
 * selection sort
 * C++ version of Algorithm 2.1 in book "Algorithms 4th edition".
 * @param array The array to sort.
 * @param n Number of elements in the array.
 * @param comparator The "less" relation between two elements in the array.
 *
 * The best-case, average-case and worse-case cost of this algorithm is O(n * n).
 */
template<class T, class Comparator>
void selectionSort(T array, unsigned int n, Comparator comparator) {
	for (unsigned int i = 0; i < n; i++) {
		unsigned int min = i;
		for (unsigned int j = i + 1; j < n; j++)
			if (comparator(array[j], array[min]))
				min = j;
		if (min != i) {
			swap(array[i], array[min]);
		}
	}
}
/**
 * insertion sort
 * C++ version of Algorithm 2.2 in book "Algorithms 4th edition".
 * @param array The array to sort.
 * @param n Number of elements in the array.
 * @param comparator The "less" relation between two elements in the array.
 *
 * best-case cost of this algorithm is O(n), such as when the array is already sorted.
 * average-case cost of this algorithm is O(n * n) if the elements in the input array
 *     aren't all identical and are randomly ordered.
 * worst-case of this algorithm is O(n * n), such as when the array is in reversed order.
 */
template<class T, class Comparator>
void insertionSort(T a, unsigned int n, Comparator comparator) {
	for (decltype(n) i = 1; i < n; i++) {
		for (auto j = i; j > 0 && comparator(a[j], a[j - 1]); j--) {
			swap(a[j], a[j - 1]);
		}
	}
}
/**
 * shell sort
 * C++ version of Algorithm 2.3 in book "Algorithms 4th edition".
 * @param array The array to sort.
 * @param n Number of elements in the array.
 * @param comparator The "less" relation between two elements in the array.
 *
 */
template<class T, class Comparator>
void shellSort(T a, unsigned int n, Comparator comparator = Comparator()) {
	decltype(n) h = 1;
	while (h < n / 3)
		h = 3 * h + 1;
	while (h >= 1) {
		for (auto i = h; i < n; i++) {
			for (auto j = i; j >= h && comparator(a[j], a[j - h]); j -= h) {
				swap(a[j], a[j - h]);
			}
		}
		h = h / 3;
	}
}
/**
 * shell sort with increments stored in an array, rather than computed it.
 * @param a The sort to sort
 * @param n Number of elements in the array
 * @param increments The increment used by shell sort, its expected to
 *        be in decreasing order. And At least one element is "1" .
 * @param m Number of increments
 * @param comparator The "less" relation between two elements in the array.
 */
template<class T, class U, class Comparator>
void shellSort(T a, unsigned int n, U increments, unsigned int m,
		Comparator comparator = Comparator()) {
	for (auto k = 0; k < m; ++k) {
		auto h = increments[k];
		for (auto i = h; i < n; i++) {
			for (auto j = i; j >= h && comparator(a[j], a[j - h]); j -= h) {
				swap(a[j], a[j - h]);
			}
		}
	}
}

/**
 * Merge two sorted array into one sorted array.
 * If {@code ar1} and {@code ar2} contains the some "equal" elements,
 *  the elements in {@code ar1} will be put before the elements in {@code ar2}
 *  in the merged sorted array.
 *  This is important for the stability of algorithms like merge sort.
 *
 * @param ar1 One of the arrays to be merged
 * @param n1 Number of elements in array {@code ar1}
 * @param ar2 The other array to be merged
 * @param n2 Number of elements in array {@code ar2}
 * @param output The merged sorted array.
 * @param comparator The "less than" relation between every two element.
 *
 * The expected number of elements remaining in the other array when the first
 * array is exhausted is about 2(1- 1/ceil(n/2)).
 */
template<class T, class U, class Comparator>
void merge(T ar1, unsigned int n1, T ar2, unsigned int n2, U output,
		Comparator comparator) {
	decltype(n1) i = 0;
	decltype(n2) j = 0;
	auto count = n1 + n2;
	decltype(count) k = 0;
	while (i < n1 && j < n2) {
		if (comparator(ar2[j], ar1[i])) {
			output[k] = ar2[j];
			++j;
		} else {
			output[k] = ar1[i];
			++i;
		}
		++k;
	}

	while (i < n1)
		output[k++] = ar1[i++];
	while (j < n2)
		output[k++] = ar2[j++];
}
/**
 * Top-down merge sort.
 * @param ar The array to be sorted.
 * @param n Number of elements in the array.
 * @param output An auxiliary array that will be used by this sorting procedure,
 * 		the length of it should be at list the same as {@code ar}.
 * 	@comparator The "less than" relation between any two element in the array
 *
 * After returned from this procedure, the elements in {@code ar} will be sorted.
 *
 * Note the template argument {@code optimize} indicates whether to used the
 * optimized version of the algorithm. {@code cutoff} is a threshold value,
 * when {@code n} is not greater than this value insertion sort will be used.
 */
template<class T, class U, class Comparator, bool optimize = true,
		unsigned int cutoff = 10>
void mergeSort(T ar, unsigned int n, U output, Comparator comparator) {
	//if only one element no need to sort
	if (n <= 1)
		return;
	//if use the optimized version and the elements to be sorted is less than the "cutoff" value
	//insertion sort will be used.
	if (optimize && n <= cutoff) {
		insertionSort(ar, n, comparator);
		return;
	}
	decltype(n) n1 = n / 2;
	decltype(n) n2 = n - n1;
	T ar2 = ar + n1;

	//sort the left half
	mergeSort(ar, n1, output, comparator);
	//sort the right half
	mergeSort(ar2, n2, output + n1, comparator);
	//if use the optimized version and the smallest element in the right half
	//is greater than the largest element in the left half, then no need to merge.
	//The possibility that this optimization succeeds depend on the "cutoff" value,
	//For a cutoff value of "10", the possibility is less than 0.002
	//	(about 1.0/C(cutoff, cutoff / 2), where C(m,n) is the combination of "m" and "n"
	//Several tests have validated this assertion.
	if (optimize && !comparator(ar[n1], ar[n1 - 1])) {
		return;
	}

	//merge the left and left sorted subarrays.
	merge(ar, n1, ar2, n2, output, comparator);
	for (decltype(n) i = 0; i < n; ++i)
		ar[i] = output[i];
}
/**
 * Merge two sorted array and into one sorted array and count the number of inversions in the array
 * formed by concanating {@code ar2} to {@code ar1}.
 *
 * @param ar1 The first sorted array.
 * @param n1 Number of elements in the first array.
 * @param ar2 The second sorted array.
 * @param n2 Number of elements in the second array.
 * @param output An auxiliary array that will be used by this algorithm.
 * @param comparator The "less than" relation between any two elements.
 * @return The number of inversions.
 */
template<class T, class U, class Comparator>
unsigned int mergeWithInversionCount(T ar1, unsigned int n1, T ar2,
		unsigned int n2, U output, Comparator comparator) {
	decltype(n1) i = 0;
	decltype(n2) j = 0;
	auto count = n1 + n2;
	decltype(count) k = 0;
	unsigned int counter = 0;
	while (i < n1 && j < n2) {
		if (comparator(ar2[j], ar1[i])) {
			output[k] = ar2[j];
			counter += (n1 - i);
			++j;
		} else {
			output[k] = ar1[i];
			++i;
		}
		++k;
	}

	while (i < n1)
		output[k++] = ar1[i++];
	while (j < n2)
		output[k++] = ar2[j++];
	return counter;
}

/**
 * Count the number of inversions in the input array. The cost of this algorithm is linearithmic.
 * @param ar The array in which the number of inversions will be counted.
 * @param n Number of elements in the array.
 * @param output An auxiliary array that will be used by this algorithm.
 * @param comparator The "less than" relation between any two elements.
 * @return The number of inversions.
 */
template<class T, class U, class Comparator>
unsigned int countInversions(T ar, unsigned int n, U output,
		Comparator comparator) {
	if (n <= 1)
		return 0;
	decltype(n) n1 = n / 2;
	decltype(n) n2 = n - n1;
	T ar2 = ar + n1;
	unsigned int counter = countInversions(ar, n1, output, comparator);
	counter += countInversions(ar2, n2, output + n1, comparator);

	counter += mergeWithInversionCount(ar, n1, ar2, n2, output, comparator);
	for (decltype(n) i = 0; i < n; ++i)
		ar[i] = output[i];
	return counter;
}
/**
 * Find the first index, say {@var i},  in an array of length n, starting from index {@var offset},
 * such that the range [offset, i) is in increasing order, but the range [offset, i] is not in increasing order.
 * If the range [offset, n) is in increasing order(we define that an empty range is always in increasing order),
 * then offset will be returned.
 * We define this index as : first disorder index.
 *
 * The cost of this algorithm is linear.
 *
 * @param values The array bo be searched.
 * @param n Number of elements in the array.
 * @param comparator The "less than" relation between two elements in the array
 * @param offset The index from which to start the search.
 * @return The first disorder index.
 */
template<class T, class Comparator>
unsigned int firstDisorderIndex(T values, unsigned int n, Comparator comparator,
		unsigned int offset = 0) {
	for (unsigned int i = offset + 1; i < n; ++i)
		if (comparator(values[i], values[i - 1]))
			return i;
	return offset;
}
/**
 * Find the number of ordered sub-arrays.
 *
 * An ordered sub-array of an array {@var A}, starting from index {@var i} and ending at index {@var j},
 * is the range [i, j) and A[i] < A[i - 1] if i > 0, A[j] < A[j - 1].
 *
 * The cost of this algorithm is linear.
 *
 * @param values The array to be searched.
 * @param n Number of elements in the array.
 * @param comparator The "less than" relation between two elements.
 */
template<class T, class Comparator>
unsigned int numberOfOrderedSubarrays(T values, unsigned int n,
		Comparator comparator) {
	if (n == 0)
		return 0;
	unsigned int count = 0;
	unsigned int i = 0;
	unsigned int tmp = 0;
	while ((tmp = firstDisorderIndex(values, n, comparator, i)) != i) {
		++count;
		i = tmp;
	}
	++count;
	return count;
}
template<class T, class U, class Comparator>
void bottomUpMergeSort(T ar, unsigned int n, U output, Comparator comparator) {
	bool isArSorted = true;
	for (unsigned int sz = 1; sz < n; sz += sz) {
		unsigned int cnt = sz + sz;
		unsigned int j = 0;
		while (j + cnt < n) {
			merge(ar + j, sz, ar + (j + sz), sz, output + j, comparator);
			j += cnt;
		}
		if (j + sz < n) {
			merge(ar + j, sz, ar + (j + sz), n - sz - j, output + j,
					comparator);
		} else {
			for (unsigned int i = j; i < n; ++i)
				output[i] = ar[i];
		}
		std::swap(ar, output);
		isArSorted = (!isArSorted);
	}
	if (!isArSorted) {
		for (decltype(n) i = 0; i < n; ++i)
			output[i] = ar[i];
	}
}
template<class T, class U, class Comparator>
void naturalMergeSort(T ar, unsigned int n, U output, Comparator comparator) {
	unsigned int last = 0;
	unsigned int cur = 0;
	unsigned int offset = 0;
	unsigned int passes = 1;
	bool isArSorted = true;
	while ((cur = firstDisorderIndex(ar, n, comparator, last + offset))
			!= offset) {
		if (cur == last + offset) {
			while (last < n) {
				output[last] = ar[last];
				++last;
			}
			std::swap(ar, output);
			isArSorted = (!isArSorted);
			last = 0;
			offset = 2 * offset + 1;
			++passes;
			continue;
		}

		unsigned int next = firstDisorderIndex(ar, n, comparator, cur + offset);
		if (next == cur + offset)
			next = n;
		merge(ar + last, cur - last, ar + cur, next - cur, output + last,
				comparator);
		last = next;
	}
	if (!isArSorted) {
		for (decltype(n) i = 0; i < n; ++i)
			output[i] = ar[i];
	}
}
template<class T, class Comparator>
bool isSorted(T values, unsigned int n, Comparator comparator) {
	return firstDisorderIndex(values, n, comparator) == 0;
}

/**
 * This algorithm is the C++ version of the algorithm 2.5 from the book "Algorithms 4th Edition"
 * This function partitions the input array by the first element in it.
 * Say the first element in the array is {@var a}.
 * This function returns an partition index {@var i}, such that :
 * 		1. the elements in the range [0, i] is not greater than than {@var a}.
 * 		2. the i-th element equals {@var a}.
 * 		3. the elements in the range [i+1, n) is not less than {@var a}.
 * @param a The input array, which will be partitioned.
 * @param n Number of elements in the input array.
 * @param comparator The "less than" relation between two elements in the array.
 * @return The partition index.
 */
template<class T, class Comparator>
unsigned int partition(T a, unsigned int n, Comparator comparator) {
	if (n == 0)
		return 0;
	unsigned int i = 0, j = n;
	while (true) {
		while (comparator(a[++i], a[(unsigned int) 0]))
			if (i == n - 1)
				break;
		while (comparator(a[(unsigned int) 0], a[--j]))
			if (j == 0)
				break;
		if (i >= j)
			break;
		std::swap(a[i], a[j]);
	}
	std::swap(a[(unsigned int) 0], a[j]);
	return j;
}
/**
 * This algorithm is the C++ version of the algorithm 2.5 from the book "Algorithms 4th Edition"
 * The very famous quick sort.
 *
 * @param a The array to be sorted.
 * @param n Number of elements in the array.
 * @param comparator The "less than" relation between two elements in the array.
 */
//static unsigned int smallArrayCount;
template<class T, class Comparator, class PartitionMethod,
		bool optimize = false, unsigned int cutoff = 5,
		unsigned int toBeIgnored = 1>
unsigned int quickSortBase1(T a, unsigned int n, Comparator comparator,
		PartitionMethod partition) {
	if (n <= toBeIgnored)
		return 1;
	if (optimize && n <= cutoff) {
		insertionSort(a, n, comparator);
		return 1;
	}
	unsigned int j = partition(a, n, comparator);
	//if(j <= 2)
	//	++ smallArrayCount;
	//if(j + 1 + 2 >= n)
	//	++smallArrayCount;
	unsigned int l = 1, r = 1;
	if (j > 1)
		l += quickSortBase1(a, j, comparator, partition);
	if (j + 2 < n)
		r += quickSortBase1(a + (j + 1), n - j - 1, comparator, partition);
	return std::max(l, r);
}
unsigned int recursionDepth = 0;
template<class T, class Comparator, class PartitionMethod,
		bool optimize = false, unsigned int cutoff = 5,
		unsigned int toBeIgnored = 1>
void quickSortBase(T a, unsigned int n, Comparator comparator,
		PartitionMethod partition) {
	recursionDepth = quickSortBase1<T, Comparator, PartitionMethod, optimize,
			cutoff, toBeIgnored>(a, n, comparator, partition);
}
template<class T, class Comparator, bool optimize = false, unsigned int cutoff =
		5>
void quickSort(T a, unsigned int n, Comparator comparator) {
	typedef unsigned int (*Func)(T, unsigned int, Comparator);
	Func func = partition;
	quickSortBase<T, Comparator, Func, optimize, cutoff>(a, n, comparator,
			func);
}
template<class T, class Comparator>
void nonRecursiveQuickSort(T a, unsigned int n, Comparator comparator) {
	if (n <= 1)
		return;
	unsigned int count = ceilLog2(n);
	unsigned int stacks[count * 2];
	int current = 0;
	stacks[0] = 0;
	stacks[1] = n;
	while (current >= 0) {
		assert(current < 2 * count);
		unsigned int s = stacks[current], e = stacks[current + 1];
		n = e - s;
		unsigned int j = partition(a + s, e - s, comparator);
		current -= 2;
		unsigned int n1 = 0;
		unsigned int n2 = 0;
		if (j > 1)
			n1 = j;
		if (j + 2 < n)
			n2 = n - j - 1;

		// This is important to make the stack size less than 2* log(n)
		if (n1 >= n2) {
			if (n1 > 0) {
				current += 2;
				stacks[current] = s;
				stacks[current + 1] = s + j;
			}
			if (n2 > 0) {
				current += 2;
				stacks[current] = s + j + 1;
				stacks[current + 1] = e;
			}
		} else {
			current += 2;
			stacks[current] = s + j + 1;
			stacks[current + 1] = e;
			if (n1 > 0) {
				current += 2;
				stacks[current] = s;
				stacks[current + 1] = s + j;
			}
		}
	}
}

template<class T, class Comparator>
void specialQuickSort(T a, unsigned int n, Comparator comparator) {
	typedef unsigned int (*Func)(T, unsigned int, Comparator);
	Func func = partition;
	quickSortBase<T, Comparator, Func, false, 0, 5>(a, n, comparator, func);
	insertionSort(a, n, comparator);
}
/**
 * the difference between this partition method and {@see partition} is 
 * this partition method assumes a sentinel element at the end of the input array, that is
 * there is an element at position {@code a[n]} or  {@code a[n - 1]} that is greater than {@code a[0]}.
 * With a sentinel, bounds checking could be eliminated,
 */
template<class T, class Comparator>
unsigned int partitionWithSentinel(T a, unsigned int n, Comparator comparator) {
	if (n == 0)
		return 0;
	unsigned int i = 0, j = n;
	while (true) {
		while (comparator(a[++i], a[(unsigned int) 0]))
			;
		while (comparator(a[(unsigned int) 0], a[--j]))
			;
		if (i >= j)
			break;
		std::swap(a[i], a[j]);
	}
	std::swap(a[(unsigned int) 0], a[j]);
	return j;
}
/**
 * This version of quick sort puts an element that is not less than {@code a[0]} at
 * the end of the array and use the the {@code partitionWithSentinel}
 * partition method to sort the array.
 */
template<class T, class Comparator>
void quickSortWithSentinal(T a, unsigned int n, Comparator comparator) {
	unsigned int (*func)(T, unsigned int, Comparator) = partitionWithSentinel;
	unsigned int index = 0;
	for (unsigned int i = 1; i < n; ++i)
		if (comparator(a[index], a[i])) {
			index = i;
			break;
		}
	swap(a[index], a[n - 1]);
	quickSortBase(a, n, comparator, func);
}
/****
 * the difference between this partition method and {@see partition} is that
 * instead of choosing the first element, the median of the first two and the last elements
 *  is chosen. The largest of the three is put at the end of the array as a sentinel so that bounds
 *  checking could be eliminated.
 */
template<class T, class Comparator>
void median3Partition(T a, unsigned int n, Comparator comparator) {
	if (n >= 3) {
		if (comparator(a[1], a[0]))
			swap(a[1], a[0]);
		if (comparator(a[n - 1], a[1]))
			swap(a[n - 1], a[1]);
		if (comparator(a[1], a[0]))
			swap(a[1], a[0]);
	}
	partitionWithSentinel(a, n, comparator);
}
template<class T, class Comparator>
void median5Partition(T a, unsigned int, Comparator comparator) {

}

/**
 * partition the array into three ranges, the first range contains the elements
 * that are less than the elements in the second range, the second range contains elements
 * that are all the same, the third range contains elements that are greater than the elements
 * in the second range.
 * @param a The array to be partitioned
 * @param n Number of elements in the array
 * @param comparator "less than" relation
 */
template<class T, class Comparator>
std::pair<unsigned int, unsigned int> fast3wayPartition(T a, unsigned int n,
		Comparator comparator) {
	if (n > 1) {
		unsigned int lo = 0, p = 1, i = 1, j = n - 1, q = n, hi = n;

		//loop invariant
		//[lo, p) = a[lo]
		//[p, i) < a[lo]
		//[q, hi) = a[lo]
		//(j, q) > a[lo]
		//(i, j) unchecked
		while (i <= j) {
			//loop invariant
			//[lo, p) = a[lo]
			//[p, i) < a[lo]
			//(i, j) unchecked
			while (i <= j && !comparator(a[lo], a[i])) {
				if (!comparator(a[i], a[lo])) {
					if (i != p)
						swap(a[i], a[p]);
					++p;
				}
				++i;
			}
			//if the while is finished, then "i > j or a[lo] < a[i]
			//loop invariant
			//[q, hi) = a[lo]
			//(j, q) > a[lo]
			//(i, j) unchecked
			while (i <= j && !comparator(a[j], a[lo])) {
				if (!comparator(a[lo], a[j])) {
					--q;
					if (j != q)
						swap(a[j], a[q]);
				}
				--j;
			}
			//if the while is finished, then "i > j or a[j] < a[lo]
			if (i > j)
				break;
			swap(a[i], a[j]);
			++i;
			--j;
		}
		while (lo < p) {
			--p;
			--i;
			swap(a[p], a[i]);
		}
		while (q < hi) {
			++j;
			swap(a[j], a[q]);
			++q;
		}
		return {i, j};
	}
	return {0, 0};
}
template<class T, class Comparator>
void fast3wayQuickSort(T a, unsigned int n, Comparator comparator) {
	if (n <= 1)
		return;
	std::pair<unsigned int, unsigned int> tmp = fast3wayPartition(a, n,
			comparator);
	if (tmp.first > 1) {
		fast3wayQuickSort(a, tmp.first, comparator);
	}

	unsigned int start = tmp.second + 1;
	if (start + 1 < n) {
		fast3wayQuickSort(a + start, n - start, comparator);
	}
}
/**
 * This class is for calculate the array accesses in each iteration of one increment
 */
class ShellSortIncrements {
	unsigned int *indices;
	vector<unsigned int> *arrayAccesses;
	unsigned int *counter;
public:
	ShellSortIncrements(unsigned int *indices,
			vector<unsigned int> *arrayAccesses, unsigned int *counter) :
			indices(indices), arrayAccesses(arrayAccesses), counter(counter) {
	}
	ShellSortIncrements& operator=(unsigned int *indices) {
		this->indices = indices;
		return *this;
	}
	unsigned int operator[](unsigned int index) {
		arrayAccesses->push_back(*counter);
		return indices[index];
	}
	auto done()->void {
		arrayAccesses->push_back(*counter);
	}
};

class InsertionSortAndSelectionSortAnimation {
	vector<unsigned int> values;
	unsigned int maxValue;
	std::thread sortingThread;
	struct SignalData {
		std::mutex m;
		std::condition_variable cv;
		bool done;
		bool isReady;
		SignalData() :
				done(false), isReady(false) {
		}
	} signalData;
	class UIntsArrayForAnimation {
	public:
		class UInt;
		class UInt {
			unsigned int *value;
			SignalData *signalData;
		public:
			UInt(unsigned int *value, SignalData *signalData) :
					value(value), signalData(signalData) {
			}
			bool operator<(const UInt& val) const {
				return *value < *val.value;
			}
			operator unsigned int() {
				return *value;
			}

			friend void swap(UInt va, UInt vb) {
				if (&va != &vb) {
					if (!va.signalData->done) {
						std::unique_lock < std::mutex > lk(va.signalData->m);
						va.signalData->cv.wait(lk,
								[va] {return va.signalData->isReady;});
						va.signalData->isReady = false;
					}
					unsigned int tmp = *va.value;
					*va.value = *vb.value;
					*vb.value = tmp;
				}
			}
		};
		UIntsArrayForAnimation(unsigned int *values, SignalData *signalData) :
				uints(values), signalData(signalData) {
		}
		UInt operator[](unsigned int index) const {
			return {uints + index, signalData};
		}
		UInt operator[](unsigned int index) {
			return {uints + index, signalData};
		}
	private:
		unsigned *uints;
		SignalData *signalData;
	};
	static void display() {
		instance->show();
	}
	static void timer(int value) {
		instance->alarm(value);
	}
	static void keyboard(unsigned char ch, int x, int y) {
		instance->kbd(ch, x, y);
	}
	void kbd(unsigned char ch, int x, int y) {
		int method = -1;
		switch (ch) {
		case '0':
			method = 0;
			break;
		case '1':
			method = 1;
			break;
		case '2':
			method = 2;
			break;
		}
		if (method >= 0) {
			done();
			sortingThread.join();
			generateData();
			sortingThread = std::thread(sort, this, method);
			signalData.done = false;
			signalData.isReady = false;
			glutTimerFunc(interval, timer, 1000);
		}
	}
	void show() {
		double xStart = -0.95;
		double yStart = -0.95;
		double xEnd = 0.95;
		double yEnd = 0.95;
		glClear(GL_COLOR_BUFFER_BIT);
		glColor3f(0.314, 0.314, 0.000); // JMU Purple
		glLineWidth(2.0);
		PerformancePloter::drawBars(xStart, xEnd, yStart, yEnd, values,
				maxValue);
		glFlush();
	}
	static void sort(InsertionSortAndSelectionSortAnimation *instance,
			int method) {
		UIntsArrayForAnimation ar(&instance->values[0], &instance->signalData);
		switch (method) {
		case 0:
			shellSort(ar, (unsigned int) (instance->values.size()),
					std::less<UIntsArrayForAnimation::UInt>());
			break;
		case 1:
			selectionSort(ar, (unsigned int) (instance->values.size()),
					std::less<UIntsArrayForAnimation::UInt>());
			break;
		case 2:
			insertionSort(ar, (unsigned int) (instance->values.size()),
					std::less<UIntsArrayForAnimation::UInt>());
			break;
		}
		instance->signalData.done = true;
	}
	void notify() {
		std::lock_guard < std::mutex > lk(signalData.m);
		signalData.cv.notify_one();
		signalData.isReady = true;
	}
	void done() {
		std::lock_guard < std::mutex > lk(signalData.m);
		signalData.cv.notify_one();
		signalData.isReady = true;
		signalData.done = true;
	}
	void alarm(int value) {
		if (!signalData.done) {
			glutTimerFunc(interval, timer, value);
			glutPostRedisplay();
			notify();
		}
	}
	void generateData() {
		unsigned int count = 100;
		std::random_device rd;
		maxValue = (count * 2);
		randUInts(values, count, maxValue);
	}
public:
	InsertionSortAndSelectionSortAnimation() {
		instance = this;
		generateData();
		sortingThread = std::thread(sort, this, 2);
	}
	int run(int argc, char *argv[]) {
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize(640, 480);
		glutInitWindowPosition(0, 0);
		glutCreateWindow("Test");
		glutDisplayFunc(display);
		glutTimerFunc(interval, timer, 1000);
		glutKeyboardFunc(keyboard);
		PerformancePloter::openGLInit();
		glutMainLoop();
		done();
		sortingThread.join();
		return 0;
	}
	static const int interval = 200;
	static thread_local InsertionSortAndSelectionSortAnimation *instance;
};
thread_local InsertionSortAndSelectionSortAnimation *InsertionSortAndSelectionSortAnimation::instance =
NULL;
class ShellSortTraces {
	vector<vector<unsigned int>> traces;
	unsigned int maxValue;
	static void display() {
		instance->show();
	}
	static void keyboard(unsigned char ch, int x, int y) {
		instance->kbd(ch, x, y);
	}
public:
	ShellSortTraces() {
		instance = this;
		reset();
	}
	void show() {
		double yStart = -0.95;
		double yEnd = 0.95;
		double yStep = (yEnd - yStart) / traces.size();
		glClear(GL_COLOR_BUFFER_BIT);
		glColor3f(0.314, 0.314, 0.000); // JMU Purple
		glLineWidth(2.0);
		for (size_t i = 0; i < traces.size(); ++i) {
			double y = yStart + i * yStep;
			double y1 = y + yStep * 0.95;
			PerformancePloter::drawBars(-0.95, 0.95, y, y1, traces[i],
					maxValue);
		}
		glFlush();
	}
	void reset() {
		traces.clear();
		vector<unsigned int> values;
		unsigned int count = 100;
		maxValue = (count * 2);
		randUInts(values, count, maxValue);
		unsigned int h = 1;
		traces.push_back(values);
		while (h < count / 3)
			h = 3 * h + 1;
		while (h >= 1) {
			shellSort(&values[0], count, &h, 1, std::less<unsigned int>());
			traces.push_back(values);
			h /= 3;
		}
	}
	void kbd(unsigned char ch, int x, int y) {
		if (ch == 'r') {
			reset();
			glutPostRedisplay();
		}
	}
	int run(int argc, char *argv[]) {
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize(640, 480);
		glutInitWindowPosition(0, 0);
		glutCreateWindow("Test");
		glutDisplayFunc(display);
		glutKeyboardFunc(keyboard);
		PerformancePloter::openGLInit();
		glutMainLoop();
		return 0;
	}
	static thread_local ShellSortTraces *instance;
};
thread_local ShellSortTraces *ShellSortTraces::instance;

int testShellSort(int argc, char *argv[]) {
	unsigned int counter = 0;
	decltype(counter) maxValue = (1 << 28);
	for (decltype(counter) start = 10; start <= 1000000; start *= 10) {
		vector<decltype(counter)> incrementValues;
		{
			decltype(start) h = 1;
			while (h < start / 3)
				h = 3 * h + 1;
			while (h >= 1) {
				incrementValues.push_back(h);
				h /= 3;
			}
		}
		decltype(counter) *values = new decltype(counter)[start];

		vector<decltype(counter)> countOfEachIteration;
		ShellSortIncrements increments(&incrementValues[0],
				&countOfEachIteration, &counter);
		randUInts(values, values + start, maxValue);
		shellSort(AccessCountedArray<decltype(counter)>(values, &counter),
				start, increments, incrementValues.size(),
				std::less<unsigned int>());
		increments.done();
		cout << "number of elemens to sort : " << start << endl;
		for (decltype(countOfEachIteration.size()) i = 1;
				i < countOfEachIteration.size(); ++i) {
			cout
					<< (countOfEachIteration[i] - countOfEachIteration[i - 1])
							/ (double) start << "\t";
		}
		cout << endl;
		cout << "average : "
				<< countOfEachIteration.back()
						/ ((double) start * incrementValues.size()) << endl;
		assert(isSorted(values, start, std::less<unsigned int>()));
		delete[] values;
	}
	return 0;
}
class DoublingTestForSort {
	int method = 0;
	clock_t last = -1;
	vector<vector<unsigned int>> values;
public:
	DoublingTestForSort(int method) :
			method(method) {
	}
	void problemSize(unsigned int size, unsigned int maxValue) {
		values.clear();
		unsigned int trials = 10;
		for (unsigned int i = 0; i < trials; ++i) {
			values.push_back(vector<unsigned int>());
			randUInts(values.back(), size, maxValue);
		}

	}
	void problemSize(unsigned int size) {
		problemSize(size, size * 2);
	}
	void operator()() {
		for (size_t i = 0; i < values.size(); ++i)
			switch (method) {
			case 0:
				shellSort(&values[i][0], (unsigned int) (values[i].size()),
						std::less<unsigned int>());
				break;
			case 1:
				selectionSort(&values[i][0], (unsigned int) (values[i].size()),
						std::less<unsigned int>());
				break;
			case 2:
				insertionSort(&values[i][0], (unsigned int) (values[i].size()),
						std::less<unsigned int>());
				break;
			case 3: {
				unsigned int tmpValues[values[i].size()];
				bottomUpMergeSort(&values[i][0],
						(unsigned int) (values[i].size()), tmpValues,
						std::less<unsigned int>());
				break;
			}
			case 4: {
				unsigned int tmpValues[values[i].size()];
				mergeSort(&values[i][0], (unsigned int) (values[i].size()),
						tmpValues, std::less<unsigned int>());
				break;
			}
			}
	}
	void add(unsigned int size, clock_t duration) {
		cout << "number of elements to sort : " << size << " duration : "
				<< duration << ", average : "
				<< (duration / (double) values.size());
		if (last != -1) {
			cout << ", ratio : " << duration / (double) last;
		}
		cout << endl;
		last = duration;
	}
};
int validateSortAlgorithms(int argc, char *argv[]) {
	const unsigned int count = 30;
	unsigned int values[count];
	unsigned int toBeSorted[count];
	unsigned int output[count];
	unsigned int maxValue = count / 5;
	randUInts(values, values + count, maxValue);

	copy(values, values + count, toBeSorted);
	selectionSort(toBeSorted, count, std::less<unsigned int>());
	assert(isSorted(toBeSorted, count, std::less<unsigned int>()));

	copy(values, values + count, toBeSorted);
	insertionSort(toBeSorted, count, std::less<unsigned int>());
	assert(isSorted(toBeSorted, count, std::less<unsigned int>()));

	copy(values, values + count, toBeSorted);
	shellSort(toBeSorted, count, std::less<unsigned int>());
	assert(isSorted(toBeSorted, count, std::less<unsigned int>()));

	copy(values, values + count, toBeSorted);
	bottomUpMergeSort(toBeSorted, count, output, std::less<unsigned int>());
	assert(isSorted(toBeSorted, count, std::less<unsigned int>()));

	copy(values, values + count, toBeSorted);
	mergeSort(toBeSorted, count, output, std::less<unsigned int>());
	assert(isSorted(toBeSorted, count, std::less<unsigned int>()));

	copy(values, values + count, toBeSorted);
	naturalMergeSort(toBeSorted, count, output, std::less<unsigned int>());
	assert(isSorted(toBeSorted, count, std::less<unsigned int>()));

	copy(values, values + count, toBeSorted);
	quickSort(toBeSorted, count, std::less<unsigned int>());
	assert(isSorted(toBeSorted, count, std::less<unsigned int>()));

	copy(values, values + count, toBeSorted);
	quickSortWithSentinal(toBeSorted, count, std::less<unsigned int>());
	assert(isSorted(toBeSorted, count, std::less<unsigned int>()));

	copy(values, values + count, toBeSorted);
	nonRecursiveQuickSort(toBeSorted, count, std::less<unsigned int>());
	assert(isSorted(toBeSorted, count, std::less<unsigned int>()));

	copy(values, values + count, toBeSorted);
	fast3wayQuickSort(toBeSorted, count, std::less<unsigned int>());
	assert(isSorted(toBeSorted, count, std::less<unsigned int>()));
	return 0;
}
int doublingTestForSorts(int argc, char *argv[]) {

	DoublingTestForSort buMerge(3);
	cout << "doubling test for bottom up merge sort:" << endl;
	doublingTest((unsigned int) 128, (unsigned int) (1 << 21), buMerge);

	DoublingTestForSort merge(4);
	cout << "doubling test for merge sort:" << endl;
	doublingTest((unsigned int) 128, (unsigned int) (1 << 21), merge);

	DoublingTestForSort shell(0);
	cout << "doubling test for shell sort:" << endl;
	doublingTest((unsigned int) 128, (unsigned int) (1 << 21), shell);

	DoublingTestForSort insertion(2);
	cout << "doubling test for insertion sort:" << endl;
	doublingTest((unsigned int) 128, (unsigned int) (1 << 21), insertion);

	DoublingTestForSort selection(1);
	cout << "doubling test for selection sort:" << endl;
	doublingTest((unsigned int) 128, (unsigned int) (1 << 21), selection);

	return 0;
}
/**
 * The total number of inversions over all binary strings of length n is :
 * 		(n * n - n) * (2 ^ (n - 3)).
 * So the expected cost of insertion sort of a randomly generated binary string of length n is :
 * 		(n * n - n) / 8;
 * This test is to validate this hypothesis.
 */
int testForSortingZeroOneArray(int argc, char *argv[]) {
	//The minimum length of binary string we begin the test
	unsigned int start = 20;
	//The maximum length of binary string we end the test
	//This value should not be too big since the number of
	//possible binary strings of length "n" is "2^n".
	//If we have too large a value, we have to do a very
	//large number of trials to make the result accurate.
	unsigned int end = 80;
	//vector to store the binary string as unsigned integer zero and one.
	vector<unsigned int> values;
	while (start <= end) {
		values.resize(start);

		//The counter to calculate the total array accesses.
		unsigned long long counter = 0;
		AccessCountedArray<unsigned int, unsigned long long*> ar(&values[0],
				&counter);
		//Number of trials. The total number of binary string of length "start" is "2 ^ start".
		//However this number is too large. So we will only have "start ^ 3" randomly generated
		//binary strings out of the "2 ^ start" binary strings. If "start" is too large, the result
		//won't be accurate, because of inadequate sampling.
		unsigned int trials = start * start * start;
		for (unsigned int i = 0; i < trials; ++i) {
			//generate random binary string.
			randUInts(values, start, 1);
			//use insertion sort to sort the binary string.
			insertionSort(ar, (unsigned int) (values.size()),
					std::less<unsigned int>());
		}
		//"average ratio" should not vary too much from a constant.
		cout << "number of elements : " << start << ", average ratio : "
				<< counter / (double) (trials * start * (start - 1) / 8)
				<< endl;
		++start;
	}
	return 0;
}
template<int N, class T = double>
struct TTuple {
	static_assert(N > 0, "N should be greater than zero");
	T values[N];
	static const int size = N;
	typedef T type;
	template<typename ... U>
	TTuple(U ... vals) :
			values { vals... } {

	}
};
template<class T, int N = 2>
class RunningTimePlotter {
	std::thread workingThread;
	std::mutex m;
	bool isDone = false;
	T val;
	vector<vector<double>> values;
	vector<unsigned int> counts;
	TTuple<3 * N> colors;
	static void display() {
		instance->show();
	}
	static void timer(int value) {
		glutTimerFunc(interval, timer, 1000);
		glutPostRedisplay();
	}
	static void func(RunningTimePlotter * const i) {
		i->work();
	}
	void work() {
		while (!isDone) {
			TTuple<N> tmp = val();
			{
				std::lock_guard < std::mutex > lk(m);
				for (size_t i = 0; i < values.size(); ++i) {
					vector<double>& v = values[i];
					unsigned int& count = counts[i];
					if (count == v.size()) {
						v.erase(v.begin());
						v.push_back(tmp.values[i]);
					} else {
						v[count] = tmp.values[i];
						++count;
					}

				}
			}
			std::chrono::milliseconds dura(50);
			std::this_thread::sleep_for(dura);
		}
	}
	double getMaxValue() {
		double maxValue = 0.0;
		for (size_t i = 0; i < values.size(); ++i) {
			vector<double>& v = values[i];
			maxValue = std::max(maxValue,
					*std::max_element(v.begin(), v.end()));
		}
		if (maxValue == 0.0)
			maxValue = 1.0;
		return maxValue;
	}
public:
	RunningTimePlotter(const T& t, const TTuple<3 * N>& colors) :
			val(t), colors(colors) {
		instance = this;
		values.resize(N);
		for (size_t i = 0; i < N; ++i) {
			values[i].resize(1000);
			counts.push_back(0);
		}
	}

	void show() {
		glClear(GL_COLOR_BUFFER_BIT);
		{
			std::lock_guard < std::mutex > lk(m);
			double maxValue = getMaxValue();
			glPointSize(3.0);
			for (size_t i = 0; i < values.size(); ++i) {
				vector<double>& v = values[i];
				glColor3dv(colors.values + 3 * i);
				PerformancePloter::plotPoints(NULL, v, maxValue, counts[i]);

			}
		}
		glFlush();

	}
	int run(int argc, char *argv[]) {
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
		glutInitWindowSize(640, 480);
		glutInitWindowPosition(0, 0);
		glutCreateWindow("Test");
		glutDisplayFunc(display);
		glutTimerFunc(interval, timer, 1000);
		PerformancePloter::openGLInit();
		for (size_t i = 0; i < counts.size(); ++i)
			counts[i] = 0;
		isDone = false;
		workingThread = std::thread(func, this);
		glutMainLoop();
		workingThread.join();
		return 0;
	}
	static const int interval = 200;
	static thread_local RunningTimePlotter *instance;
};
class SortingTime {
	unsigned int start = 100;
	unsigned int end = 10000;
	unsigned int step = 10;
	unsigned int method = 0;
public:
	SortingTime(unsigned int method = 0) :
			method(method) {
	}
	TTuple<2> operator()() {
		vector<unsigned int> values;
		randUInts(values, start, start * 2);
		double ratio = 200000.0;
		clock_t begin = clock();
		typedef AccessCountedArray<unsigned int, unsigned long long *> Array;
		void (*func)(Array, unsigned int, Array,
				std::less<unsigned int>) = nullptr;
		switch (method) {
		case 0:
			insertionSort(values, (unsigned int) (values.size()),
					std::less<unsigned int>());
			start += step;
			return {(double)(clock() - begin), (start - step) * (start - step) / ratio};
		case 1:
			selectionSort(values, (unsigned int) (values.size()),
					std::less<unsigned int>());
			start += step;
			return {(double)(clock() - begin), (start - step) * (start - step) / ratio};
		case 2:
			shellSort(values, (unsigned int) (values.size()),
					std::less<unsigned int>());
			start += step;
			return {(double)(clock() - begin), (start - step) * std::pow(start - step, 1 / 2.0) / ratio};
		case 3:
			func = bottomUpMergeSort;
		case 4:
			if (func == nullptr)
				func = mergeSort;
		case 5: {
			if (func == nullptr)
				func = naturalMergeSort;

			unsigned int tmpValues[values.size()];
			unsigned long long counter = 0;
			double val = 6 * start * std::log(start);
			Array ar1(&values[0], &counter), ar2(tmpValues, &counter);
			func(ar1, start, ar2, std::less<unsigned int>());
			start += step;
			return {(double)counter, val};
			//return {(double)(clock() - begin), val / ratio};
		}
		case 6: {
			unsigned int comp = 0;
			unsigned int count = values.size();
			//double expected = count * 5.0 / 6.0;
			//smallArrayCount = 0;
			quickSort(&values[0], count,
					[&comp](unsigned int v1, unsigned int v2) {++comp; return v1 < v2;});
			//cout << "expected : " << expected << ", actual : " << smallArrayCount << ", ratio : " << smallArrayCount / expected <<  endl;
			start += step;
			return {(double)comp, 2 * count * log(count)};
			break;
		}
		}
		return {0.0, 0.0};
	}
};
class ArmortizedCost {
	unsigned int problemSize = 1000;
	unsigned int counts = 0;
	unsigned long long totalTime = 0;
	unsigned int method;
public:
	ArmortizedCost(unsigned int method = 0, unsigned int problemSize = 1000) :
			problemSize(problemSize), method(method) {

	}
	TTuple<2> operator()() {
		vector<unsigned int> values;
		randUInts(values, problemSize, problemSize * 2);
		clock_t begin = clock();
		switch (method) {
		case 0:
			insertionSort(values, (unsigned int) (values.size()),
					std::less<unsigned int>());
			break;
		case 1:
			selectionSort(values, (unsigned int) (values.size()),
					std::less<unsigned int>());
			break;
		case 2:
			shellSort(values, (unsigned int) (values.size()),
					std::less<unsigned int>());
			break;
		}
		clock_t duration = clock() - begin;
		totalTime += duration;
		++counts;
		return {(double)duration, totalTime / (double)counts};
	}
};
class CompareDistributions {
	unsigned int start = 100;
	unsigned int method = 0;
	unsigned int distribution = 0;
public:
	CompareDistributions(unsigned int method) :
			method(method) {
	}
	template<class Func>
	TTuple<4> testOne(Func func) {
		TTuple<4> result;
		if (func != nullptr) {
			unsigned int values[start];
			unsigned long long counter = 0;
			AccessCountedArray<unsigned int, unsigned long long *> ar(values,
					&counter);
			static std::default_random_engine rd;
			static std::uniform_int_distribution<unsigned int> discrete(0,
					2 * start);
			static std::poisson_distribution<unsigned int> poisson(start);
			static std::geometric_distribution<unsigned int> geometric(
					1 / log(start));
			static std::normal_distribution<double> normal(start, start * 0.01);

			std::generate(values, values + start, [&] {return discrete(rd);});
			func(ar, start, std::less<unsigned int>());
			result.values[0] = counter;

			counter = 0;
			std::generate(values, values + start, [&] {return poisson(rd);});
			func(ar, start, std::less<unsigned int>());
			result.values[1] = counter;

			counter = 0;
			std::generate(values, values + start, [&] {return geometric(rd);});
			func(ar, start, std::less<unsigned int>());
			result.values[2] = counter;

			counter = 0;
			std::generate(values, values + start,
					[&] {return (unsigned int)normal(rd);});
			func(ar, start, std::less<unsigned int>());
			result.values[3] = counter;
		}
		return result;
	}
public:

	TTuple<4> operator()() {
		void (*func)(AccessCountedArray<unsigned int, unsigned long long *>,
				unsigned int, std::less<unsigned int>) = nullptr;
		switch (method) {
		case 0:
			func = selectionSort;
			break;
		case 1:
			func = insertionSort;
			break;
		case 2:
			func = shellSort;
			break;
		}
		start += 10;
		return testOne(func);
	}
};
class SpecialDistributions {
	unsigned int start = 100;
	unsigned int method = 0;
	unsigned int distribution = 0;
public:
	SpecialDistributions(unsigned int method = 0) :
			method(method) {
	}
	static void generate0(unsigned int *values, unsigned int count) {
		for (unsigned int i = 0; i < count; ++i)
			values[i] = (i % 2);
		std::random_shuffle(values, values + count);
	}
	static void generate1(unsigned int *values, unsigned int count) {
		for (unsigned int i = 0; i < count / 2; ++i)
			values[i] = 0;
		randUInts(values + count / 2, values + count, count);
		std::random_shuffle(values, values + count);
	}
	static void generate2(unsigned int *values, unsigned int count) {
		unsigned int i = 0;
		unsigned int *p = values;
		unsigned int n = count;
		while (count != 1) {
			unsigned int tmp = count / 2;
			for (unsigned int j = 0; j < tmp; ++j)
				values[j] = i;
			values += tmp;
			count -= tmp;
			++i;
		}
		*values = i;
		std::random_shuffle(p, p + n);
	}
	template<class Func>
	TTuple<3> testOne(Func func) {
		TTuple<3> result;
		if (func != nullptr) {
			unsigned int values[start];
			unsigned long long counter = 0;
			AccessCountedArray<unsigned int, unsigned long long *> ar(values,
					&counter);

			generate0(values, start);
			func(ar, start, std::less<unsigned int>());
			result.values[0] = counter;

			counter = 0;
			generate1(values, start);
			func(ar, start, std::less<unsigned int>());
			result.values[1] = counter;

			counter = 0;
			generate2(values, start);
			func(ar, start, std::less<unsigned int>());
			result.values[2] = counter;
		}
		return result;
	}
public:

	TTuple<3> operator()() {
		void (*func)(AccessCountedArray<unsigned int, unsigned long long *>,
				unsigned int, std::less<unsigned int>) = nullptr;
		switch (method) {
		case 0:
			func = selectionSort;
			break;
		case 1:
			func = insertionSort;
			break;
		case 2:
			func = shellSort;
			break;
		}
		start += 10;
		return testOne(func);
	}
};
template<int N>
class ShellSortIncrementCompare {
	unsigned int start = 100;
public:
	typedef std::vector<vector<unsigned int>> Increments;
	ShellSortIncrementCompare(const Increments& increments) :
			increments(increments) {
		for (size_t i = 0; i < increments.size(); ++i) {
			vector<unsigned int>& vals = this->increments[i];
			std::sort(vals.begin(), vals.end(), std::greater<unsigned int>());
			for (int j = 0; j < vals.size(); ++j)
				cout << vals[j] << ' ';
			cout << endl;
		}
	}
	TTuple<N> operator()() {
		unsigned int origin[start];
		unsigned int toBeSorted[start];
		clock_t duration = 0;
		randUInts(origin, origin + start, start * 2);
		unsigned long long counter = 0;
		AccessCountedArray<unsigned int, unsigned long long *> ar(toBeSorted,
				&counter);
		TTuple<N> result;
		for (size_t i = 0; i < increments.size(); ++i) {
			std::copy(origin, origin + start, toBeSorted);
			counter = 0;
			duration = clock();
			shellSort(ar, start, &increments[i][0], increments[i].size(),
					less<unsigned int>());
			duration = clock() - duration;
			result.values[i] = counter;
		}
		start += 10;
		return result;
	}
private:
	Increments increments;
};
class ShellSortIncrementCompare1: public ShellSortIncrementCompare<4> {
public:
	typedef vector<unsigned int> Irts;
	ShellSortIncrementCompare1() :
			ShellSortIncrementCompare<4>(
					{
					//formed by merging together the se	quences "9* (4 ^ k) - 9 * (2 ^ k) + 1" and "(4 ^ k - 3 * (2 ^ k) + 1"
							Irts { 1, 5, 19, 41, 109, 209, 505, 929, 2161, 3905,
									8929, },
							//formed from the sequence "(3 ^ k - 1) / 2"
							Irts { 1, 4, 13, 40, 121, 364, 1093, 3280, 9841, },
							//formed from the sequence "(2^(k - 1))"
							Irts { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, },
							//prime numbers of the form "(2^k - 1)"
							Irts { 1, 3, 7, 31, 127, 8191, }, }) {
	}
};
template<int N>
class MergeSortCompare {
	unsigned int start = 100;
	bool measureTime = false;
public:
	static void generate(unsigned int *values, unsigned int count) {
		randUInts(values, values + count, 2 * count);
		//SpecialDistributions::generate2(origin, start);
		//std::sort(values, values + count);

		//unsigned int tmp = count * 0.1;
		//for(unsigned int i = 0;i < tmp;++ i)
		//	swap(values[rand() % count], values[rand() % count]);
	}
	typedef AccessCountedArray<unsigned int, unsigned long long *> Array;
	typedef void (*Func)(Array, unsigned int, Array, std::less<unsigned int>);
	TTuple<N> operator()() {
		unsigned int origin[start];
		unsigned int toBeSorted[start];
		unsigned int tmpValues[start];
		clock_t duration = 0;
		unsigned long long counter = 0;
		Func *func = funcs.values;
		Array ar1(toBeSorted, &counter), ar2(tmpValues, &counter);
		TTuple<N> result;

		generate(origin, start);
		for (size_t i = 0; i < N; ++i) {
			std::copy(origin, origin + start, toBeSorted);
			counter = 0;
			duration = clock();
			func[i](ar1, start, ar2, std::less<unsigned int>());
			duration = clock() - duration;
			if (!measureTime)
				duration = counter;
			result.values[i] = duration;
		}
		start *= 1.01;
		//start += 10;
		return result;
	}
	void set(Func fs[N]) {
		for (int i = 0; i < N; ++i) {
			funcs.values[i] = fs[i];
		}
	}
	//MergeSortCompare(const TTuple<N, Func>& funcs) : funcs(funcs){
	//}
protected:
	TTuple<N, Func> funcs;
};
class DifferentMergeCompare: public MergeSortCompare<3> {
public:
	DifferentMergeCompare() {
		Func tmps[3] = { bottomUpMergeSort, mergeSort, naturalMergeSort };
		set(tmps);
	}

};
class DifferentCutOffCompare: public MergeSortCompare<7> {
public:
	DifferentCutOffCompare() {
		Func tmps[] = {
				mergeSort<Array, Array, std::less<unsigned int>, true, 5>,
				mergeSort<Array, Array, std::less<unsigned int>, true, 8>,
				mergeSort<Array, Array, std::less<unsigned int>, true, 11>,
				mergeSort<Array, Array, std::less<unsigned int>, true, 14>,
				mergeSort<Array, Array, std::less<unsigned int>, true, 17>,
				mergeSort<Array, Array, std::less<unsigned int>, true, 20>,
				mergeSort<Array, Array, std::less<unsigned int>, true, 23> , };
		set(tmps);
	}

};
class ShellSortIncrementCompare2: public ShellSortIncrementCompare<12> {
public:
	typedef vector<unsigned int> Irts;
	ShellSortIncrementCompare2() :
			ShellSortIncrementCompare<12>( {
			//formed from the sequence (2)^(k - 1))
					Irts { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048,
							4096, },
					//formed from the sequence (3)^(k - 1))
					Irts { 1, 3, 9, 27, 81, 243, 729, 2187, 6561, },
					//formed from the sequence (4)^(k - 1))
					Irts { 1, 4, 16, 64, 256, 1024, 4096, },
					//formed from the sequence (5)^(k - 1))
					Irts { 1, 5, 25, 125, 625, 3125, },
					//formed from the sequence (6)^(k - 1))
					Irts { 1, 6, 36, 216, 1296, 7776, },
					//formed from the sequence (7)^(k - 1))
					Irts { 1, 7, 49, 343, 2401, },
					//formed from the sequence (8)^(k - 1))
					Irts { 1, 8, 64, 512, 4096, },
					//formed from the sequence (9)^(k - 1))
					Irts { 1, 9, 81, 729, 6561, },
					//formed from the sequence (10)^(k - 1))
					Irts { 1, 10, 100, 1000, },
					//formed from the sequence (11)^(k - 1))
					Irts { 1, 11, 121, 1331, },
					//formed from the sequence (12)^(k - 1))
					Irts { 1, 12, 144, 1728, },
					//prime numbers of the form "(2^k - 1)"
					Irts { 1, 3, 7, 31, 127, 8191, }, }) {
	}
	static void output() {
		for (unsigned int i = 2; i <= 12; ++i) {
			cout << "//formed from the sequence (" << i << ")^(k - 1))" << endl;
			cout << "Irts{";
			unsigned int t = 1;
			while (t <= 8000) {
				cout << t << ", ";
				t *= i;
			}
			cout << "}," << endl;
		}

	}
};
template<unsigned int N>
class QuickSortCompare {
protected:
	unsigned int start = 100;
	unsigned int method = 0;
public:
	typedef AccessCountedArray<unsigned int, unsigned long long *> Array;
	typedef CountedLessThan<unsigned int> LessThan;
	typedef void (*Func)(Array, unsigned int, LessThan);
	TTuple<N> operator()() {
		//dangerous!!! platform dependent, if "start" is larger than the stacck size limit
		//of the platform, program will crash.
		//unsigned int origin[start];
		//unsigned int toBeSorted[start];
		vector<unsigned int> originVec(start, 0);
		vector<unsigned int> toBeSortedVec(start, 0);
		unsigned int *origin = &originVec[0];
		unsigned int *toBeSorted = &toBeSortedVec[0];
		clock_t duration = 0;
		unsigned long long accessCounter = 0, compareCounter = 0;
		Func *func = funcs.values;
		Array ar(toBeSorted, &accessCounter);
		TTuple<N> result;
		LessThan lessThan(&compareCounter);

		//generate(origin, start);
		for (size_t i = 0; i < N; ++i) {
			std::copy(origin, origin + start, toBeSorted);
			accessCounter = 0;
			compareCounter = 0;
			duration = clock();
			func[i](ar, start, lessThan);
			duration = clock() - duration;
			switch (method % 3) {
			case 0:
				result.values[i] = duration;	// + 0.01;
				break;
			case 1:
				result.values[i] = accessCounter;	// + 0.01;
				break;
			default:
				result.values[i] = compareCounter;	// + 0.01;
				break;
			}
			//assert(isSorted(toBeSorted, start, lessThan));
		}
		//for(size_t i = 1;i < N;++ i){
		//	result.values[i]  = (result.values[i] + 1e-8) / (result.values[0] + 1e-8);
		//}
		//result.values[0] = 1;
		start *= 1.01;
		//start += 10;
		return result;
	}
protected:
	void set(Func fs[N]) {
		for (int i = 0; i < N; ++i) {
			funcs.values[i] = fs[i];
		}
	}
	TTuple<N, Func> funcs;
};
class QuickSortCompare1: public QuickSortCompare<4> {
public:
	QuickSortCompare1(unsigned int method = 0) {
		Func tmps[4] = { quickSort, quickSortWithSentinal,
				nonRecursiveQuickSort, specialQuickSort };
		set(tmps);
		this->method = method;
	}
};
class QuickSortCompare2: public QuickSortCompare<4> {
public:
	QuickSortCompare2(unsigned int method = 0) {
		Func tmps[4] = { quickSort<Array, LessThan, true, 5>, quickSort<Array,
				LessThan, true, 8>, quickSort<Array, LessThan, true, 10>,
				quickSort<Array, LessThan, true, 15> };
		set(tmps);
		this->method = method;
	}
};
template<class T, class Comparator>
unsigned int randPartition(T a, unsigned int n, Comparator comparator) {
	unsigned int index = rand() % n;
	std::swap(a[0], a[index]);
	return partition(a, n, comparator);
}
template<class T, class Comparator>
void randShuffleSort(T a, unsigned int n, Comparator comparator) {
	random_shuffle(a, a + n);
	return quickSort(a, n, comparator);
}
template<class T, class Comparator>
void randPartitionSort(T a, unsigned int n, Comparator comparator) {
	typedef unsigned int (*Func)(T, unsigned int, Comparator);
	Func func = randPartition;
	quickSortBase<T, Comparator, Func, false, 5>(a, n, comparator, func);

}
class QuickSortCompare3: public QuickSortCompare<2> {
public:
	QuickSortCompare3(unsigned int method = 0) {
		Func tmps[2] = { randShuffleSort, randPartitionSort, };
		set(tmps);
		this->method = method;
	}
};
template<class T, int N>
thread_local RunningTimePlotter<T, N> *RunningTimePlotter<T, N>::instance =
		nullptr;
class RecursionDepth {
	unsigned int start = 100;
public:
	TTuple<2> operator()() {
		vector<unsigned int> values;
		unsigned int tmpDepth = 0;
		unsigned int count = start;
		for (unsigned int i = 0; i < 10; ++i) {
			randUInts(values, count, count * 2);

			::recursionDepth = 0;
			quickSort(&values[0], count, std::less<unsigned int>());
			tmpDepth += recursionDepth;
		}
		start *= 1.01;
		auto expected = 4 * log(count), average = tmpDepth / 10.0;
		cout << "expected : " << expected << ", actual : " << average
				<< ", ratio : " << average / expected << endl;
		return {average, expected};
	}
};
class DebugRunningTime {
	unsigned int start = 0;
	unsigned int N = 10000;
public:
	TTuple<2> operator()() {
		unsigned int counts[N];
		static std::default_random_engine rd;
		static std::uniform_int_distribution<unsigned int> discrete(0, N);
		memset(counts, 0, sizeof(counts));
		start += 10;
		for (unsigned int i = 0; i < start; ++i)
			++counts[discrete(rd)];
		double expected = (1 - 1 / exp(start * 1.0 / N));
		double count = 0;
		for (unsigned int i = 0; i < N; ++i)
			if (counts[i] != 0)
				count += 1;
		return {expected, count / N};
	}
};
int testPlotting(int argc, char *argv[]) {
	//InsertionSortAndSelectionSortAnimation test;
	//ShellSortTraces test;
//	ShellSortIncrementCompare1 tmp;
//	RunningTimePlotter<ShellSortIncrementCompare1, 4> test(tmp,
//				{
//						1.0, 0.0, 0.0,
//						0.0, 1.0, 0.0,
//						1.0, 1.0, 0.0,
//						1.0, 1.0, 1.0,
//				});
//	ShellSortIncrementCompare2 tmp;
//	RunningTimePlotter<ShellSortIncrementCompare2, 12> test(tmp,
//			{
//					0.0, 0.0, 0.0,
//					0.1, 0.1, 0.0,
//					0.2, 0.2, 0.0,
//					0.3, 0.3, 0.0,
//					0.4, 0.4, 0.0,
//					0.5, 0.5, 0.0,
//					0.6, 0.6, 0.0,
//					0.7, 0.7, 0.0,
//					0.8, 0.8, 0.0,
//					0.9, 0.9, 0.0,
//					1.0, 1.0, 0.0,
//					1.0, 0.0, 0.0,
//			});
	SpecialDistributions tmp(1);
	RunningTimePlotter<SpecialDistributions, 3> test(tmp, { 1.0, 0.0, 0.0, 0.0,
			1.0, 0.0, 1.0, 1.0, 0.0 });
//	DebugRunningTime tmp;
//	RunningTimePlotter<DebugRunningTime, 2> test(tmp, { 1.0, 0.0, 0.0, 0.0, 1.0,
//			0.0 });
//	return test.run(argc, argv);
	//SortingTime tmp(6);
	//RunningTimePlotter<SortingTime, 2> test(tmp,
	//		{ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 });
//	DifferentCutOffCompare tmp;
//	RunningTimePlotter<DifferentCutOffCompare, 7> test(tmp,
//			{
//										0.1, 0.1, 0.1,
//										0.2, 0.2, 0.0,
//										0.3, 0.3, 0.0,
//										0.4, 0.4, 0.0,
//										0.5, 0.5, 0.0,
//										0.6, 0.6, 0.0,
//										0.7, 0.7, 0.7,
//			});
//	QuickSortCompare3 tmp(1);
//	RunningTimePlotter<QuickSortCompare3, 2> test(tmp, { 1.0, 0.0, 0.0, 0.0,
//			1.0, 0.0, });

	//RecursionDepth tmp;
	//RunningTimePlotter<RecursionDepth, 2> test(tmp, { 1.0, 0.0, 0.0, 0.0, 1.0,
	//		0.0, });
	return test.run(argc, argv);
}

struct NaturalOrderAnalysis {
	/**
	 * Count the number of cases, denoted by C(n, i),  that there are i ordered subarrays in an array of length n.
	 *
	 * @param n Length of an array.
	 * @param counts {@code counts[i]} denote the number of cases that that there
	 * 		are (i+1) ordered subarray in an array of length {@code n}
	 *
	 * 	Note :
	 * 		The following recurrence relation is used :
	 * 			C(n + 1, i) = i * C(n, i) + (n + 2 - i) * C(n, i - 1), for 2 <= i <= n
	 * 			C(n + 1, 1) = C(n + 1, n + 1) = 1
	 */
	static void count(unsigned int n, unsigned long long *counts) {
		if (n == 0)
			return;
		counts[0] = 1;
		if (n == 1)
			return;

		counts[1] = 1;
		for (unsigned int j = 3; j <= n; ++j) {

			unsigned long long last = counts[0];
			for (unsigned int i = 1; i < j - 1; ++i) {
				unsigned long long tmp = counts[i];
				counts[i] = i * tmp + (j + 1 - i) * last;
				last = tmp;
			}
			counts[j - 1] = 1;
		}

	}
	static void validateAlgorithm() {
		unsigned int values1[] = { 0, 1, 2, 3, 4 };
		unsigned int values2[] = { 4, 3, 2, 1, 0 };
		unsigned int values3[] = { 0, 3, 1, 2, 5, 4 };
		unsigned int values4[] = { 1, 1, 1, 0, 0, 0, 1, 1, 0 };
		assert(
				numberOfOrderedSubarrays(values1,
						(sizeof values1 / sizeof(*values1)),
						std::less<unsigned int>()) == 1);
		assert(
				numberOfOrderedSubarrays(values2,
						(sizeof values2 / sizeof(*values2)),
						std::less<unsigned int>()) == 5);
		assert(
				numberOfOrderedSubarrays(values3,
						(sizeof values3 / sizeof(*values3)),
						std::less<unsigned int>()) == 3);
		assert(
				numberOfOrderedSubarrays(values4,
						(sizeof values4 / sizeof(*values4)),
						std::less<unsigned int>()) == 3);
	}
	/**
	 * Given an randomly generated array of length n, what's the expected number of ordered subarrays in it?
	 * Denote the expected number of ordered subarrays in an array of length n as E(n).
	 * Then we will have the following recurrence relation:
	 * 		E(n+1) = (n / (n+1))E(n) + 1;
	 * 		The base case is E(1) = 1.
	 * 	After solving this recurrence relation, we will have
	 * 		E(n) = (n + 1) / 2
	 *
	 */
	static int test(unsigned int argc, char *argv[]) {
		validateAlgorithm();

		unsigned int start = 3;
		while (start <= 10000) {
			unsigned int values[start];
			unsigned int iteration = 1;
			unsigned int count = 0;
			for (unsigned int i = 0; i < start; ++i)
				values[i] = i;
			for (unsigned int i = 0; i < iteration; ++i) {
				std::random_shuffle(values, values + start);
				count += numberOfOrderedSubarrays(values, start,
						std::less<unsigned int>());
			}

			double average = count / (double) iteration;
			cout << "expected : " << (start + 1) / 2 << ", actual : " << average
					<< ", ratio : " << average * 2 / (start + 1) << endl;
			start += 10;
		}
		return 0;
	}
};
struct ExpectedInversions {
	static void validateAlgorithm() {
		unsigned int values1[] = { 1, 2, 3, 4, 5, 6 };
		unsigned int values2[] = { 6, 5, 4, 3, 2, 1 };
		unsigned int values3[] = { 3, 1, 2, 5, 4, 1 };
		unsigned int values4[] = { 1, 1, 1, 0, 0, 0 };
		unsigned int aux[20];
		assert(
				countInversions(values1, (sizeof values1 / sizeof(*values1)),
						aux, std::less<unsigned int>()) == 0);
		assert(
				countInversions(values2, (sizeof values2 / sizeof(*values2)),
						aux, std::less<unsigned int>()) == 15);
		assert(
				countInversions(values3, (sizeof values3 / sizeof(*values3)),
						aux, std::less<unsigned int>()) == 7);
		assert(
				countInversions(values4, (sizeof values4 / sizeof(*values4)),
						aux, std::less<unsigned int>()) == 9);
	}
	/**
	 * For an randomly generated array of length {@var n}, the expected number of
	 * inversions in it is {@code n * (n - 1) / 2}
	 */
	static int test(unsigned int argc, char *argv[]) {
		validateAlgorithm();
		unsigned int start = 3;
		while (start <= 10000) {
			unsigned int values[start];
			unsigned int aux[start];
			unsigned int iteration = 100;
			unsigned long long count = 0;
			for (unsigned int i = 0; i < start; ++i)
				values[i] = i;
			for (unsigned int i = 0; i < iteration; ++i) {
				std::random_shuffle(values, values + start);
				count += countInversions(values, start, aux,
						std::less<unsigned int>());
			}

			double average = count / (double) iteration;
			double expected = start * (start - 1) / 4;
			cout << "expected : " << expected << ", actual : " << average
					<< ", ratio : " << average / expected << endl;
			start += 10;
		}
		return 0;
	}
};
template<class T>
struct LinkedList {
	struct Link {
		T val;
		Link *next;
		Link * prev;
		Link(const T&val, Link *n = nullptr, Link *p = nullptr) :
				val(val), next(n), prev(p) {

		}
	};
	template<class Comparator = std::less<T> >
	void sort(Comparator comparator = Comparator()) {
		if (count <= 1)
			return;
		LinkedList newList;
		split(count / 2, newList);
		sort(comparator);
		newList.sort(comparator);
		merge(newList, comparator);
	}
	template<class Comparator = std::less<T> >
	bool isSorted(Comparator comparator = Comparator()) {
		if (count <= 1)
			return true;
		Link *i = head;
		Link *j = i->next;
		while (j != nullptr) {
			if (comparator(j->val, i->val))
				return false;
			i = j;
			j = i->next;
		}
		return true;
	}
	template<class Random>
	void shuffle(Random random) {
		if (count <= 1)
			return;
		LinkedList newList;
		split(count / 2, newList);
		shuffle(random);
		newList.shuffle(random);
		randomMerge(newList, random);
	}
	bool empty() {
		return count == 0;
	}
	void swap(LinkedList& l) {
		std::swap(head, l.head);
		std::swap(tail, l.tail);
		std::swap(count, l.count);
	}
	template<class Comparator = std::less<T> >
	void merge(LinkedList& l, Comparator comparator = Comparator()) {
		LinkedList newList;
		merge(*this, l, newList, comparator);
		swap(newList);
	}

	template<class Random>
	void randomMerge(LinkedList& l, Random random) {
		LinkedList newList;
		randomMerge(*this, l, newList, random);
		swap(newList);
	}
	void split(unsigned int count, LinkedList& l) {
		Link *link = head;
		if (count == 0) {
			swap(l);
			return;
		}
		for (unsigned int i = 0; i < count && link != nullptr;
				++i, link = link->next)
			;
		if (link != nullptr) {
			l.head = link;
			l.tail = tail;
			tail = link->prev;
			tail->next = nullptr;
			l.head->prev = nullptr;
			l.count = this->count - count;
			this->count = count;
		}
	}
	template<class Comparator = std::less<T> >
	LinkedList() :
			head(nullptr), tail(nullptr), count(0) {
	}
	~LinkedList() {
		clear();
	}
	void clear() {
		Link *link = head;
		while (link != nullptr) {
			Link *next = link->next;
			delete link;
			link = next;
		}
		head = tail = nullptr;
		count = 0;
	}
	void push_back(const T& val) {
		add(new Link(val, nullptr, nullptr));
	}
	const Link* getHead() const {
		return head;
	}
	LinkedList& operator=(const LinkedList&) = delete;
	LinkedList(const LinkedList&) = delete;
	template<class Comparator = std::less<T> >
	static void merge(LinkedList& l1, LinkedList& l2, LinkedList& l,
			Comparator comparator = Comparator()) {
		while (!l1.empty() && !l2.empty()) {
			Link *i1 = l1.head;
			Link *i2 = l2.head;
			if (comparator(i2->val, i1->val)) {
				l2.remove(i2);
				l.add(i2);
			} else {
				l1.remove(i1);
				l.add(i1);
			}
		}
		l.addBack(l1);
		l.addBack(l2);
	}
	template<class Random>
	static void randomMerge(LinkedList& l1, LinkedList& l2, LinkedList& l,
			Random random) {
		while (!l1.empty() && !l2.empty()) {
			unsigned int tmp = random() % (l2.count + 1);
			Link *i1 = l1.head;
			LinkedList newList;
			l2.split(tmp, newList);
			l.addBack(l2);
			newList.swap(l2);
			l1.remove(i1);
			l.add(i1);
		}
		l.addBack(l1);
		l.addBack(l2);
	}
private:
	void add(Link *link) {
		if (head == nullptr) {
			assert(count == 0);
			assert(tail == nullptr);
			head = tail = link;
		} else {
			assert(count != 0);
			assert(tail != nullptr);
			link->next = tail->next;
			link->prev = tail;
			tail->next = link;
			tail = link;
		}
		++count;
	}
	void remove(Link *link) {
		if (link == head) {
			head = link->next;
			if (head != nullptr) {
				assert(count > 1);
				assert(tail != link);
				head->prev = nullptr;
			} else {
				assert(count == 1);
				assert(tail == link);
			}
		} else if (link == tail) {
			assert(tail != head);
			tail = link->prev;
			tail->next = nullptr;
			assert(count > 1);
		} else {
			Link *prev = link->prev;
			Link *next = link->next;
			prev->next = next;
			next->prev = prev;
		}
		--count;
		return;
	}
	void addBack(LinkedList& l) {
		if (head == nullptr) {
			swap(l);
		} else if (l.head != nullptr) {
			l.head->prev = tail;
			tail->next = l.head;
			tail = l.tail;
			count += l.count;
			l.head = l.tail = nullptr;
			l.count = 0;
		}
	}
	Link *head, *tail;
	unsigned int count;
};
struct LinkedListTest {
	static void generate(unsigned int count, LinkedList<unsigned int>& l) {
		std::random_device rd;
		std::uniform_int_distribution<unsigned int> generator(0, 2 * count);
		l.clear();
		for (unsigned int i = 0; i < count; ++i)
			l.push_back(generator(rd));
	}
	static void output(const LinkedList<unsigned int>& l, ostream& os) {
		auto begin = l.getHead();
		while (begin != nullptr) {
			os << begin->val << ' ';
			begin = begin->next;
		}
		os << endl;
	}
	static int test(int argc, char *argv[]) {
		for (unsigned int start = 2; start <= 1000; start += 10) {
			LinkedList<unsigned int> l;
			generate(start, l);
			l.sort();
			if (start <= 40)
				output(l, cout);
			assert(l.isSorted());
		}
		unsigned int iteration = 10000;
		unsigned int count = 10;
		unsigned int counts[count][count];
		LinkedList<unsigned int> l;
		for (auto& ar : counts)
			for (auto& val : ar)
				val = 0;

		for (unsigned int i = 0; i < count; ++i)
			l.push_back(i);
		output(l, cout);
		cout << "shuffle test : " << endl;
		for (unsigned int i = 0; i < iteration; ++i) {
			l.shuffle(rand);
			output(l, cout);
			unsigned int j = 0;
			auto begin = l.getHead();
			while (begin != nullptr) {
				++counts[begin->val][j];
				++j;
				begin = begin->next;
			}
		}

		for (auto& ar : counts) {
			for (auto& val : ar) {
				cout << val / (double) iteration << '\t';
			}
			cout << endl;
		}
		return 0;
	}
};
template<class Algorithm>
void sortRank(unsigned int n, unsigned int *r, Algorithm algorithm) {
	std::vector<unsigned int> indices(n, 0);
	for (unsigned int i = 0; i < n; ++i)
		indices[i] = i;
	algorithm(&indices[0], n);
	for (unsigned int i = 0; i < n; ++i)
		r[indices[i]] = i;
}
template<typename T>
class SortRank {
	T items;
	T aux;
	unsigned int method;
public:
	SortRank(T items, T aux, unsigned int method) :
			items(items), aux(aux), method(method) {
	}
	void operator()(unsigned int *values, unsigned int n) {
		auto func =
				[&](unsigned int i, unsigned int j) {return this->items[i] < this->items[j];};
		switch (method % 9) {
		case 0:
			selectionSort(values, n, func);
			return;
		case 1:
			insertionSort(values, n, func);
			return;
		case 2:
			shellSort(values, n, func);
			return;
		case 3:
			quickSort(values, n, func);
			return;
		case 4:
			fast3wayQuickSort(values, n, func);
			return;
		case 5:
			mergeSort(&values[0], n, &aux[0], func);
			return;
		case 6:
			bottomUpMergeSort(&values[0], n, &aux[0], func);
			return;
		case 7:
			naturalMergeSort(&values[0], n, &aux[0], func);
			return;
		case 8: {
			auto func1 = [&](unsigned int i, unsigned int j) {return this->items[j] < this->items[i];};
			std::priority_queue<unsigned int, vector<unsigned int>, decltype(func1)> pq(&values[0], &values[0] + n, func1);
			unsigned int i = 0;
			while(!pq.empty()){
				values[i] = pq.top();
				pq.pop();
				++ i;
			}
			return ;
		}
		}
	}
};

int testRank(int argc, char *argv[]) {
	typedef vector<unsigned int> UVec;
	unsigned int count = 8;
	UVec values(count, 0);
	UVec tmp(count, 0), result(count, 0);
	SortRank<UVec> ranks[] = { { values, tmp, 0 }, { values, tmp, 1 }, { values,
			tmp, 2 }, { values, tmp, 3 }, { values, tmp, 4 },
			{ values, tmp, 5 }, { values, tmp, 6 }, { values, tmp, 7 }, { values, tmp, 8 },};
	const char *descs[] = { "selectionSort", "insertionSort", "shellSort",
			"quickSort", "fast3wayQuickSort", "mergeSort", "bottomUpMergeSort",
			"naturalMergeSort", "heap sort"};
	for (unsigned int i = 0; i < (sizeof(ranks) / sizeof(*ranks)); ++i) {
		sortRank(count, &result[0], ranks[i]);
		cout << "sort algorithm : " << descs[i] << endl;
		std::copy(result.begin(), result.end(),
				std::ostream_iterator<unsigned int>(cout, " "));

		cout << endl;
	}
	return 0;
}
template<class T, class Comparator>
int kendallTauDistance(T items1, T items2, unsigned int n,
		Comparator comparator) {
	typedef vector<unsigned int> UVec;
	UVec orders(n, 0);
	UVec ranks(n, 0);
	UVec indices(n, 0);
	UVec aux(n, 0);
	for (unsigned int i = 0; i < n; ++i) {
		orders[i] = i;
		indices[i] = i;
		aux[i] = i;
	}
	fast3wayQuickSort(&orders[0], n,
			[&](unsigned int i, unsigned int j) {return items1[i] < items1[j];});
	fast3wayQuickSort(&aux[0], n,
			[&](unsigned int i, unsigned int j) {return items2[i] < items2[j];});
	for (unsigned int i = 0; i < n; ++i)
		ranks[aux[i]] = i;
	unsigned int inversions =
			countInversions(&indices[0], n, &aux[0],
					[&](unsigned int i, unsigned int j) {
						if(!comparator(items2[i], items2[j]) && !comparator(items2[j], items2[i]))
						return false;
						return orders[ranks[i]] < orders[ranks[j]];
					});
	////////assertions///////////////
	for (unsigned int i = 0; i < n; ++i) {
		assert(
				!comparator(items1[i], items2[aux[i]])
						&& !comparator(items2[aux[i]], items1[i]));
	}
	////////////////////////////////
	return inversions;
}
void testKendallTauDistance(unsigned int *values1, unsigned int *values2,
		unsigned int count) {
	std::random_shuffle(values1, values1 + count);
	std::random_shuffle(values2, values2 + count);
	cout << "sequence 1 :" << endl;
	std::copy(values1, values1 + count,
			std::ostream_iterator<unsigned int>(cout, " "));
	cout << "\nsequence 2 :" << endl;
	std::copy(values2, values2 + count,
			std::ostream_iterator<unsigned int>(cout, " "));
	cout << "\ninversions : "
			<< kendallTauDistance(values1, values2, count,
					std::less<unsigned int>()) << endl;
}
int testKendallTauDistance(int argc, char *argv[]) {
	unsigned int count = 8;
	unsigned int values1[count];
	unsigned int values2[count];
	for (unsigned int i = 0; i < count; ++i) {
		values1[i] = i;
		values2[i] = i;
	}
	testKendallTauDistance(values1, values2, count);
	for (unsigned int i = 0; i < count; ++i) {
		values1[i] = i % 3;
		values2[i] = i % 3;
	}
	testKendallTauDistance(values1, values2, count);
	return 0;
}

int ext21Main(int argc, char *argv[]) {
	validateSortAlgorithms(argc, argv);
	//return testRank(argc, argv);
	return testPlotting(argc, argv);
}
