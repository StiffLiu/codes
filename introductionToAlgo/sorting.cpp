#include <cassert>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <iterator>
#include <iostream>
using namespace std;
/**
 * The counting sort.
 * This algorithm comes from the book "Introduction to algorithms", second edition, chapter 8.
 * This algorithm sorts an input array by the calculated "key", where the range of "key" is between
 * 0 and "keyCount".
 * For an input array with n elements, the time complexity of this algorithm is O(n + keyCount).
 * The space complexity of this algorithm is also O(n + keyCount).
 * This algorithm is stable, the stability is very import for some other algorithms, such as "radix sort".
 */
template<class KeyCalculator>
class StableCountingSort {
private:
	//avoid copy constructor
	StableCountingSort(const StableCountingSort&);
	//avoid assignment
	StableCountingSort & operator=(const StableCountingSort&);
protected:
	unsigned int *counts;
	unsigned int keyCount;
	KeyCalculator keyCalculator;
public:
	StableCountingSort(unsigned int keyCount,
			const KeyCalculator& keyCaculator = KeyCalculator()) {
		counts = 0;
		this->keyCount = keyCount;
		this->keyCalculator = keyCalculator;
	}
	~StableCountingSort() {
		delete[] counts;
	}
	KeyCalculator& getKeyCalculator() {
		return keyCalculator;
	}
	const KeyCalculator& getKeyCalculator() const {
		return keyCalculator;
	}
	void setKeyCalculator(const KeyCalculator& keyCaculator) {
		this->keyCalculator = keyCalculator;
	}
	template<class T>
	void sort(T *values, unsigned int count) {
		if (keyCount <= 1)
			return;
		if (counts == 0)
			counts = new unsigned int[keyCount];

		//the number of appearance of each key is initialized to zero.
		for (unsigned int i = 0; i < keyCount; ++i)
			counts[i] = 0;

		//count the number of appearance of each key.
		T *tmpValues = new T[count];
		for (unsigned int i = 0; i < count; ++i) {
			unsigned int index = keyCalculator(values[i]);
			tmpValues[i] = values[i];
			assert(index < keyCount);
			++counts[index];
		}

		//calculate the last index of each key in the sorted array.
		for (unsigned int i = 1; i < keyCount; ++i)
			counts[i] += counts[i - 1];

		//put every value  in right position.
		//The back-to-start processing is essential to maintain the stability of this algorithm
		for (unsigned int i = count; i > 0; --i) {
			unsigned int index = keyCalculator(tmpValues[i - 1]);
			unsigned int& j = counts[index];
			assert(j > 0);
			--j;
			values[j] = tmpValues[i - 1];
		}
	}
};
template<class KeyCalculator>
class UnStableCountingSort: private StableCountingSort<KeyCalculator> {
public:
	UnStableCountingSort(unsigned int keyCount,
			const KeyCalculator& keyCaculator = KeyCalculator()) :
			StableCountingSort<KeyCalculator>(keyCount, keyCaculator) {

	}
	using StableCountingSort<KeyCalculator>::getKeyCalculator;
	using StableCountingSort<KeyCalculator>::setKeyCalculator;
	template<class T>
	int sort(T *values, unsigned int count) {
		unsigned int& keyCount = StableCountingSort<KeyCalculator>::keyCount;
		unsigned int *&counts = StableCountingSort<KeyCalculator>::counts;
		KeyCalculator& keyCalculator =
				StableCountingSort<KeyCalculator>::keyCalculator;
		if (keyCount <= 1)
			return 0;
		if (counts == 0)
			counts = new unsigned int[keyCount * 2];

		unsigned int *lowBound = counts + keyCount;
		//the number of appearance of each key is initialized to zero.
		for (unsigned int i = 0; i < keyCount; ++i)
			counts[i] = 0;

		//count the number of appearance of each key.
		for (unsigned int i = 0; i < count; ++i) {
			unsigned int index = keyCalculator(values[i]);
			assert(index < keyCount);
			++counts[index];
		}

		//calculate the last index of each key in the sorted array.
		lowBound[0] = counts[0];
		for (unsigned int i = 1; i < keyCount; ++i) {
			counts[i] += counts[i - 1];
			lowBound[i] = counts[i];
		}

		//put every value in right position.
		//The back-to-start processing is essential to maintain the stability of this algorithm
		int k = 0;
		for (unsigned int i = count; i > 0; ++k) {
			unsigned int j = i - 1;
			T& v1 = values[j];
			unsigned int index = keyCalculator(v1);
			if (j >= lowBound[index] && j < counts[index]) {
				--i;
				continue;
			}
			--lowBound[index];
			swap(values[lowBound[index]], values[j]);
		}

		return 2 * keyCount + count + k;
	}
};
struct RadixKeyCalculator {
	unsigned int radix;
	unsigned int modulee;
	unsigned int operator()(unsigned int val) {
		return val / modulee % radix;
	}
};
struct IntKeyCalculator {
	unsigned int operator()(unsigned int val) {
		return val;
	}
};
/*
 * radix sort.
 */
void radixSort(unsigned int *values, unsigned int n) {

	const unsigned int radix = 10;
	StableCountingSort<RadixKeyCalculator> countingSort(radix);
	RadixKeyCalculator& keyCalculator = countingSort.getKeyCalculator();
	keyCalculator.radix = 10;
	keyCalculator.modulee = 1;

	unsigned int maxValue = *std::max_element(values, values + n);
	while (maxValue != 0) {
		countingSort.sort(values, n);
		keyCalculator.modulee *= keyCalculator.radix;
		maxValue /= keyCalculator.radix;
	}
}
//void sortString(const char **strs, unsigned int n, unsigned int startIndex = 0){
//	unsigned int minLen = 0;
//	for(unsigned int i = 0;i < n;++ i){
//		unsigned j = 0
//	}
//}
void randInts(unsigned int *values, unsigned int n, unsigned int maxValue) {
	for (unsigned int i = 0; i < n; ++i)
		values[i] = rand() % maxValue;
}
bool isSorted(unsigned int *values, unsigned int n) {
	for (unsigned int i = 1; i < n; ++i)
		if (values[i] < values[i - 1])
			return false;
	return true;
}
int testSort(int argc, char *argv[]) {
	const unsigned int n = 100;
	unsigned int values[n];
	unsigned int maxValue = 1000;
	srand(time(0));
	randInts(values, n, maxValue);
	radixSort(values, n);
	assert(isSorted(values, n));

	maxValue = 10;
	for (int i = 0; i < 100; ++i) {
		randInts(values, n, maxValue);
		UnStableCountingSort<IntKeyCalculator> countingSort(maxValue);
		int iterations = countingSort.sort(values, n);
		assert(isSorted(values, n));
		copy(values, values + n, ostream_iterator<unsigned int>(cout, " "));
		cout << endl;
		cout << iterations / (double) n << endl;
	}
	return 0;
}

