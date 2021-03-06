#include <cstdlib>
#include <iostream>
#include <iostream>
#include <iterator>
#include <ctime>
#include <cassert>
#include <algorithm>
#include <cmath>
using namespace std;
/*Every procedure
 *        in this file is intended as experiments. The costs of the procedures will be measured
 *        and will be compared with the problem size. The measurement of almost every
 *        procedure in this file is the number of array access. That's why many of the procedures
 *        in this file counts the number of array accesses.
 *
 *Many of the cost assertions of algorithms in this file is based on the hypothesis that the elements
 *        of an array are distinct or almost distinct.Otherwise, say every(or almost every )element in an array is the same,
 *        If not stated specifically,  then the cost of worst case should apply if the cost is stated.
 */

/**
 *
 * insertion sort.
 * The cost of worst case insertion sort is O(n * n), for example when the array is ordered in "decreasing" order.
 * The cost of best case insertion sort is O(n), when the array is(or almost is) sorted in "non-decreasing" order.
 * However, the cost of insertion sort on average is O(n * n).
 * @param values The array to be sorted.
 * @param n Number of arrays in the input array.
 * @param comparator A functor that defines the "smaller" relation between two elements.
 *               If the first element is "smaller" than the second, then "true" should be returned.
 * @return The number of array accesses is returned.
 */
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
		for (; j > 0 && comparator(val, (arrayAccesses++, values[j - 1]));
				--j) {
			values[j] = values[j - 1];
		}
		arrayAccesses += ((i - j) * 2 + 1);

		//there are two cases, when the for loop is terminated :
		//        case 1 :  j == 0,  we have reached the beginning of the array.
		//        case 2 :  values[j - 1] <= val, since the range [0, i) is sorted we could
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

			T tmp = values[i];
			values[i] = values[startIndex + median];
			values[startIndex + median] = tmp;
		}
		arrayAccesses += 3 * g;

		unsigned int remaining = n - startIndex;
		if (remaining != 0) {
			arrayAccesses += insertionSort(values, remaining, comparator);

			T tmp = values[g];
			values[g] = values[startIndex + remaining / 2];
			values[startIndex + remaining / 2] = tmp;
			arrayAccesses += 3;
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
 * partition the elements in "values", into two parts.
 * The partition index "p" with pivot "pivot" is defined as:
 *      the elements in the range [0, p) is no greater than "pivot" and
 *      the elements in the range [p, n) is greater than "pivot".
 * If every elements in the range [0, n) is no greater than "pivot", then "n" is returned,
 *     else the partition index is returned.
 * Note that these implementations are just for experimentations, thus some extra values
 *     are returned, such as the number of array accesses used by this procedure.
 */
template<class T, class Comparator>
unsigned int partition(T *values, T pivot, unsigned int n,
		unsigned int& arrayAccesses, Comparator comparator) {
	unsigned int i = 0;
	unsigned int j = 0;
	//loop invariant:
	//      the elements in the range [0, i) is no greater than "pivot"
	//      and the elements in the range [i, j) is greater than "pivot"
	while (j < n) {
		arrayAccesses++;
		if (!comparator(pivot, values[j])) {
			if (i != j) {
				T tmp = values[i];
				values[i] = values[j];
				values[j] = tmp;
				arrayAccesses += 3;
			}
			++i;
		}
		++j;
	}
	return i;
}
/**
 * partition using the default comparator.
 */
template<class T>
unsigned int partition(T *values, T pivot, unsigned int n,
		unsigned int& arrayAccesses) {
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
		Comparator comparator) {
	if (n == 1)
		return 0;
	unsigned int p = partition(values, values[n - 1], n - 1, arrayAccesses,
			comparator);
	if (p != n - 1) {
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
unsigned int partition(T *values, unsigned int n, unsigned int& arrayAccesses) {
	return partition(values, n, arrayAccesses, less<T>());
}
/**
 * partition using a random value in the input array as the "pivot"
 */
template<class T, class Comparator>
unsigned int randPartition(T *values, unsigned int n,
		unsigned int& arrayAccesses, Comparator comparator) {
	unsigned int index = rand() % n;
	if (index != n - 1) {
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
unsigned int randPartition(T *values, unsigned int n,
		unsigned int& arrayAccesses) {
	return randPartition(values, n, arrayAccesses, less<T>());
}
template<class T, class Comparator>
unsigned int medianNPartition(T *values, unsigned int n,
		unsigned int& arrayAccesses, Comparator comparator) {
	if (n <= 1)
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
unsigned int medianNPartition(T *values, unsigned int n,
		unsigned int& arrayAccesses) {
	return medianNPartition(values, n, arrayAccesses, less<T>());
}

template<class T, class Comparator>
bool isPartitionedByIndex(T *values, unsigned int n, unsigned int p,
		Comparator comparator) {
	for (unsigned int i = 0; i < p; ++i)
		if (comparator(values[p], values[i]))
			return false;
	for (unsigned int i = p + 1; i < n; ++i)
		if (!comparator(values[p], values[i]))
			return false;
	return true;
}
template<class T>
bool isPartitionedByIndex(T *values, unsigned int n, unsigned int p) {
	return isPartitionedByIndex(values, n, p, less<T>());
}
template<class T, class Comparator>
bool isPartitionedByIndex0(T *values, unsigned int n, unsigned int p,
		Comparator comparator) {
	for (unsigned int i = 0; i < p; ++i)
		if (comparator(values[p], values[i]))
			return false;
	for (unsigned int i = p + 1; i < n; ++i)
		if (comparator(values[i], values[p]))
			return false;
	return true;
}
template<class T>
bool isPartitionedByIndex0(T *values, unsigned int n, unsigned int p) {
	return isPartitionedByIndex0(values, n, p, less<T>());
}
template<class T, class U>
void calculateSum(T *values, unsigned int start, unsigned int end, U& val) {
	for (unsigned int i = start; i < end; ++i)
		val += values[i].second;
}
/*
 * This procedure assumes that type "T" contains a member variable, which is named "second"
 *  and is of type "U".This member variable represents the weight.
 *  "weightMedian" is a value that will be used by this procedure.
 *  "weightTotal" must be the sum of weights of all the elements in the input array.
 *  This procedure partitions the input array, and returns the partition index.
 *  Assume that the returned partition index is "p", then the following properties hold after the partition is done :
 *  	Every element in the range [0, p) is no greater than "values[p]", and the sum of weights of elements in this range is less than "weightMedian".
 *  	Every element in the range (p + 1, n) is no less than "values[p]", and the sum of weights of elements in this range is no greater than "weightMedian".
 */
template<class T, class U, class Comparator>
unsigned int weightedMedian(T *values, U weightMedian, U weightTotal,
		unsigned int n, unsigned int& arrayAccesses, Comparator comparator) {
	if (n <= 1)
		return 0;;
	unsigned int startIndex = 0;
	unsigned int endIndex = n;
	U lowerHalfSum = U();
	while (startIndex < endIndex - 1) {
		unsigned int p = startIndex
				+ /*randPartition*/medianNPartition(values + startIndex,
						endIndex - startIndex, arrayAccesses, comparator);
		U tmpSum = lowerHalfSum;
		calculateSum(values, startIndex, p, tmpSum);
		arrayAccesses += (p - startIndex);
		if (tmpSum < weightMedian) {
			if (!(weightMedian + tmpSum + (arrayAccesses++, values[p]).second
					< weightTotal))
				return p;
			startIndex = p;
			lowerHalfSum = tmpSum;
		} else {
			endIndex = p;
		}
	}
	return startIndex;
}
template<class T>
unsigned int weightedMedian(pair<T, double> *values, unsigned int n,
		unsigned int& arrayAccesses) {
	typedef pair<T, double> Type;
	class Func {
	public:
		static bool func(const Type& v1, const Type& v2) {
			return v1.first < v2.first;
		}
	};
	return weightedMedian(values, 0.5, 1.0, n, arrayAccesses, Func::func);
}
/**
 * The kth order statistic of an array is the kth smallest element in the array
 * @param values The input array.
 * @param n Number of elements in the array.
 * @param k The index at which the kth order statitistic will be placed.
 * @param comparator The functor used to compare element.
 *        If the first element is less than the second element, this functor
 *        should return true.
 *
 * @return The number of array accesses used by this procedure.
 *
 * The cost of this algorithms is decided by the choosing of pivot element when we partition the input array.
 * If at every step a pivot is chosen such that the partition produces a very bad split,
 * for example when the smallest or the greatest element is selected every time,
 * the cost of would be O(n*n).
 * However if every time(or almost every time) a pivot element is selected that partition the input array
 * into two parts such that the the number of elements in the two partitions are "a*n" and "(1-a)*n".
 * The expected cost of this algorithm is O(n).
 *
 *
 */
template<class T, class Comparator>
unsigned int kthOrderStatistic(T *values, unsigned int n, unsigned int k,
		Comparator comparator) {
	if (k >= n)
		return 0;
	unsigned int arrayAccesses = 0;
	while (n > 0) {
		unsigned int p = randPartition(values, n, arrayAccesses, comparator);
		if (p == k)
			break;
		if (p < k) {
			values += (p + 1);
			n -= (p + 1);
			k -= (p + 1);
			continue;
		}
		n = p;
	}
	return arrayAccesses;
}
template<class T>
unsigned int kthOrderStatistic(T *values, unsigned int n, unsigned int k) {
	return kthOrderStatistic(values, n, k, less<T>());
}
/*
 * Calculate a, b, and m, such that :
 *           n = a * m + b * (m - 1),
 * 			 k = a + b
 */
inline void calculate(unsigned int n, unsigned int k, unsigned int& a,
		unsigned int& b, unsigned int& m) {
	m = n / k;
	a = k;
	b = 0;
	if (n % k != 0) {
		++m;
		b = k * m - n;
		a -= b;
	}
}

template<class T, class Comparator>
unsigned int kQuantiles(T *values, unsigned int n, unsigned int k,
		Comparator comparator) {
	if (n <= 1 || k <= 1)
		return 0;
	assert(k <= n);
	unsigned int a, b, m;
	calculate(n, k, a, b, m);
	unsigned int i = 0;
	unsigned int startIndex = 0;
	if (2 * a * m > n) {
		for (; i < a && 2 * i * m < n; ++i)
			;

		startIndex = i * m;
	} else {
		i = a;
		startIndex = i * m;
		for (; i < k && 2 * startIndex < n; ++i, startIndex += (m - 1))
			;
	}
	//select the "startIndex"-th order statistics.
	unsigned int arrayAccesses = kthOrderStatistic(values, n, startIndex,
			comparator);
	return arrayAccesses + kQuantiles(values, startIndex, i, comparator)
			+ kQuantiles(values + startIndex, n - startIndex, k - i, comparator);
}
template<class T>
unsigned int kQuantiles(T *values, unsigned int n, unsigned int k) {
	return kQuantiles(values, n, k, less<T>());
}
/*
 * This procedure calculates the k numbers in an array that are closest to the median of the array.
 * These k elements are placed in the range [n / 2 - k / 2,  n / 2 - k / 2  + k).
 * "comparator" is a functor that defines the "smaller" relation between two elements.
 * The expected running time of this algorithm is O(n) if the elements in the input array are (or almost) distinct.
 * To measure the performance of the algorithm, the number of array accesses is returned.
 */
template<class T, class Comparator>
unsigned int closestToMedian(T *values, unsigned int n, unsigned int k,
		Comparator comparator) {
	if (k >= n || k <= 0)
		return 0;
	unsigned int start = n / 2 - k / 2;
	unsigned int arrayAccesses = kthOrderStatistic(values, n, start,
			comparator);
	if (k > 1)
		arrayAccesses += kthOrderStatistic(values + start, n, k - 1,
				comparator);
	return arrayAccesses;
}
template<class T>
unsigned int closestToMedian(T *values, unsigned int n, unsigned int k) {
	return closestToMedian(values, n, k, less<T>());
}
/*
 * Calculate the kth order statistic of two sorted arrays.
 * This cost of this algorithm is O(lg(max(n1, n2))) in worst case.
 * Where "n1" and "n2" are the number of elements in the two sorted arrays.
 */
template<class T, class Comparator>
std::pair<T*, unsigned int> kthOrderStatisticOf2OrderedArrays(T *ordered1,
		unsigned int n1, T *ordered2, unsigned int n2, unsigned int k,
		unsigned int& arrayAccesses, Comparator comparator) {
	if (k > n1 + n2 && k != 0)
		k = n1 + n2 - 1;
	unsigned int s1 = 0, e1 = n1, s2 = 0, e2 = n2;
	while (k > 0 && e1 > s1 && e2 > s2) {
		int m1 = (e1 + s1) / 2;
		int m2 = (e2 + s2) / 2;
		int r1 = m1 - s1;
		int r2 = m2 - s2;
		int r = r1 + r2;
		if ((arrayAccesses += 2, comparator(ordered2[m2], ordered1[m1]))) {
			if (k > r) {
				s2 = m2 + 1;
				k -= (r2 + 1);
			} else {
				e1 = m1;
			}
		} else if ((arrayAccesses += 2, comparator(ordered1[m1], ordered2[m2]))) {
			if (k > r) {
				s1 = m1 + 1;
				k -= (r1 + 1);
			} else {
				e2 = m2;
			}
		} else {
			if (r == k || r + 1 == k) {
				return make_pair(ordered1, m1);
			}
			if (k < r) {
				e1 = m1;
				e2 = m2;
			} else {
				s1 = m1 + 1;
				s2 = m2 + 1;
				k -= (r + 2);
			}
		}
	}

	if (e1 <= s1)
		return make_pair(ordered2, s2 + k);
	if (e2 <= s2)
		return make_pair(ordered1, s1 + k);
	if (comparator(ordered1[s1], ordered2[s2]))
		return make_pair(ordered1, s1);
	return make_pair(ordered2, s2);
}
template<class T>
std::pair<T*, unsigned int> kthOrderStatisticOf2OrderedArrays(T *ordered1,
		unsigned int n1, T *ordered2, unsigned int n2, unsigned int k,
		unsigned int& arrayAccesses) {
	return kthOrderStatisticOf2OrderedArrays(ordered1, n1, ordered2, n2, k,
			arrayAccesses, less<T>());
}
template<class T, class Comparator>
bool iskQuantiles(T *values, unsigned int n, unsigned int k,
		Comparator comparator) {
	if (n <= 1 || k <= 1)
		return true;
	assert(k <= n);
	unsigned int a, b, m;
	calculate(n, k, a, b, m);

	unsigned int startIndex = m;
	while (startIndex < a * m) {
		if (!isPartitionedByIndex0(values + startIndex - m, 2 * m, m,
				comparator))
			return false;
		startIndex += m;
	}
	if (startIndex < n) {
		if (!isPartitionedByIndex0(values + startIndex - m, 2 * m - 1, m,
				comparator))
			return false;
		startIndex += m - 1;
		while (startIndex < n) {
			if (!isPartitionedByIndex0(values + startIndex - m + 1, 2 * m - 2,
					m - 1, comparator))
				return false;
			startIndex += (m - 1);
		}
	}
	return true;
}
template<class T>
bool iskQuantiles(T *values, unsigned int n, unsigned int k) {
	return iskQuantiles(values, n, k, less<T>());
}
extern void randInts(unsigned int *values, unsigned int n,
		unsigned int maxValue);
extern bool isSorted(unsigned int *values, unsigned int n);
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
int testPartitionRandom(int argc, char *argv[]) {
	const unsigned int n = 1000;
	unsigned int values[n];
	unsigned int array1[n];
	unsigned int array2[n];
	unsigned int array3[n];

	unsigned int maxValue = 1000;

	for (unsigned int i = 0; i < 10; ++i) {

		randInts(values, n, maxValue);
		//sort(values, values + n);
		copy(values, values + n, array1);
		copy(values, values + n, array2);
		copy(values, values + n, array3);

		unsigned int arrayAccesses1 = 0, arrayAccesses2 = 0, arrayAccesses3 = 0;
		unsigned int p1 = partition(array1, n, arrayAccesses1);
		unsigned int p2 = randPartition(array2, n, arrayAccesses2);
		unsigned int p3 = medianNPartition(array3, n, arrayAccesses3);

		assert(isPartitionedByIndex(array1, n, p1));
		assert(isPartitionedByIndex(array2, n, p2));
		assert(isPartitionedByIndex(array3, n, p3));

		cout << "parition :  index - " << p1 << ", number of array accesses - "
				<< arrayAccesses1 << endl;
		cout << "random parition :  index - " << p2
				<< ", number of array accesses - " << arrayAccesses2 << endl;
		cout << "median n parition :  index - " << p3
				<< ", number of array accesses - " << arrayAccesses3 << endl;

		cout << endl;
	}
	return 0;
}
int testPartition(int argc, char *argv[]) {
	const unsigned int n = 20;
	unsigned int values[n];
	unsigned int arrayAccesses = 0;
	unsigned int p;
	randInts(values, n, 100);
	cout << "before partition :";
	copy(values, values + n, ostream_iterator<int>(cout, " "));
	cout << endl;
	p = partition(values, n, arrayAccesses);
	cout << "partition index : " << p << endl;
	cout << "after partition : ";
	copy(values, values + n, ostream_iterator<int>(cout, " "));
	cout << endl;
	return testPartitionRandom(argc, argv);
}
template<class T>
void generateRandom(pair<T, double> *values, unsigned int n,
		unsigned int maxValue) {
	double sum = 0;
	for (unsigned int i = 0; i < n; ++i) {
		values[i].first = rand() % maxValue;
		values[i].second = rand() % maxValue;
		sum += values[i].second;
	}

	double current = 0.0;
	for (unsigned int i = 0; i < n; ++i) {
		if (current >= 1.0)
			values[i].second = 0.0;
		else
			values[i].second /= sum;
		current += values[i].second;
	}
}
int testWeightedMedian(int argc, char *argv[]) {
	typedef pair<unsigned int, double> Type;
	const unsigned int n = 1000;
	const unsigned int maxValue = 1000;
	Type values[n];
	double maxRatio = 0.0;
	double minRatio = n;
	double averageRatio = 0.0;
	const unsigned iteration = 100;
	for (unsigned int i = 0; i < iteration; ++i) {
		double tSum = 0.0;
		generateRandom(values, n, maxValue);
		unsigned int arrayAccesses = 0;
		unsigned int p = weightedMedian(values, n, arrayAccesses);
		double lSum = 0.0, rSum = 0.0;
		double ratio = (arrayAccesses / (double) n);
		calculateSum(values, 0, p, lSum);
		calculateSum(values, 0, n, tSum);
		calculateSum(values, p + 1, n, rSum);
		cout << "array accesses : " << arrayAccesses << ", ratio : " << ratio
				<< ", " << lSum << ", " << rSum << endl;
		assert(lSum < 0.5 && rSum <= 0.5);
		maxRatio = max(maxRatio, ratio);
		minRatio = min(minRatio, ratio);
		averageRatio += ratio;
	}
	cout << "maxRatio : " << maxRatio << ", minRatio : " << minRatio
			<< ", average ratio : " << (averageRatio / iteration) << endl;
	return 0;
}
int testKthOrderStatistic(int argc, char *argv[]) {
	const unsigned int n = 10000;
	unsigned int values[n];
	unsigned int maxValue = 10000;
	double maxRatio = 0;
	double minRatio = n;
	double averageRatio = 0;
	double iteration = n;
	srand(time(0));
	for (unsigned int i = 0; i < iteration; ++i) {
		double ratio = 0;
		randInts(values, n, maxValue);
		ratio = kthOrderStatistic(values, n, i) / (double) n;
		cout << i << "-th order statistic ratio : " << ratio << endl;
		assert(isPartitionedByIndex0(values, n, i));
		maxRatio = max(maxRatio, ratio);
		minRatio = min(minRatio, ratio);
		averageRatio += ratio;
	}
	cout << "maxRatio : " << maxRatio << ", minRatio : " << minRatio
			<< ", average ratio : " << (averageRatio / iteration) << endl;
	return 0;
}
int testKQuantilesRandom(int argc, char *argv[]) {
	const unsigned int n = 10000;
	unsigned int values[n];
	unsigned int maxValue = 10000;
	double maxRatio = 0;
	double minRatio = n;
	double averageRatio = 0;
	double iteration = n;
	srand(time(0));
	for (unsigned int i = 0; i < iteration; ++i) {
		double ratio = 0;
		randInts(values, n, maxValue);
		ratio = kQuantiles(values, n, i + 1) / (double) n;
		cout << i << " quantiles ratio : " << ratio << endl;
		assert(iskQuantiles(values, n, i + 1));
		maxRatio = max(maxRatio, ratio);
		minRatio = min(minRatio, ratio);
		averageRatio += ratio;
	}
	cout << "maxRatio : " << maxRatio << ", minRatio : " << minRatio
			<< ", average ratio : " << (averageRatio / iteration) << endl;
	return 0;
}
int testKQuantiles(int argc, char *argv[]) {
	const unsigned int n = 30;
	unsigned int values[n];
	unsigned int maxValue = 100;
	randInts(values, n, maxValue);
	kQuantiles(values, n, n);
	int ret = testKQuantilesRandom(argc, argv);
	copy(values, values + n, ostream_iterator<unsigned int>(cout, " "));
	cout << endl;
	return ret;
}
int testStaticsticOf2OrderedArrays(int argc, char *argv[]) {
	const unsigned int n = 1000;
	unsigned int input1[n];
	unsigned int input2[n];
	unsigned int input3[2 * n];
	unsigned int maxInt = n * n;
	double maxRatio = 0;
	double minRatio = n;
	double averageRatio = 0;
	double iteration = n;
	double log2n = log(n) / log(2);
	srand(time(0));
	for (unsigned int j = 0; j < iteration; ++j) {
		for (unsigned int i = 0; i < n; ++i) {
			input1[i] = rand() % maxInt;
			input2[i] = rand() % maxInt;
			input3[2 * i] = input1[i];
			input3[2 * i + 1] = input2[i];
		}

		sort(input1, input1 + n);
		sort(input2, input2 + n);
		sort(input3, input3 + 2 * n);

		unsigned int arrayAccesses = 0;
		unsigned int k = rand() % (2 * n);
		pair<unsigned int*, int> result = kthOrderStatisticOf2OrderedArrays(
				input1, n, input2, n, k, arrayAccesses);
		if (result.first[result.second] != input3[k])
			result = kthOrderStatisticOf2OrderedArrays(input1, n, input2, n, k,
					arrayAccesses);
		assert(result.first[result.second] == input3[k]);
		double ratio = arrayAccesses / log2n;

		maxRatio = max(maxRatio, ratio);
		minRatio = min(minRatio, ratio);
		averageRatio += ratio;
		cout << j << " : " << ratio << endl;
	}
	cout << "maxRatio : " << maxRatio << ", minRatio : " << minRatio
			<< ", average ratio : " << (averageRatio / iteration) << endl;
	return 0;
}
