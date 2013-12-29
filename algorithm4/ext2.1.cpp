/**
 * This file is for the test, visualization, comparison of some basic sorting algorithms,
 * that is selection sort, insertion sort and shell sort.
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
#include <chrono>
#include <tuple>
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
	operator T*() {
		return val;
	}
	operator const T*() {
		return val;
	}
	AccessCountedArray operator+(unsigned int n) {
		return {val + n, counter};
	}
	const T& operator[](unsigned int index) const {
		if (counter != NULL)
			++*counter;
		return val[index];

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
 * merge two sorted array into one sorted array
 *
 * @param ar1 One of the arrays to be merged
 * @param n1 Number of elements in array {@code ar1}
 * @param ar2 The other array to be merged
 * @param n2 Number of elements in array {@code ar2}
 * @param output The merged sorted array.
 * @param comparator The "less than" relation between every two element.
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

template<class T, class U, class Comparator, bool optimize = true,
		unsigned int cutoff = 10>
void mergeSort(T ar, unsigned int n, U output, Comparator comparator) {
	if (n <= 1)
		return;
	if (optimize && n <= cutoff) {
		insertionSort(ar, n, comparator);
		return;
	}
	decltype(n) n1 = n / 2;
	decltype(n) n2 = n - n1;
	T ar2 = ar + n1;
	mergeSort(ar, n1, output, comparator);
	mergeSort(ar2, n2, output + n1, comparator);
	if (optimize && !comparator(ar[n1], ar[n1 - 1])) {
		return;
	}

	merge(ar, n1, ar2, n2, output, comparator);
	for (decltype(n) i = 0; i < n; ++i)
		ar[i] = output[i];
}
template<class T, class Comparator>
unsigned int firstDisorderIndex(T values, unsigned int n, Comparator comparator,
		unsigned int offset = 0) {
	for (unsigned int i = offset + 1; i < n; ++i)
		if (comparator(values[i], values[i - 1]))
			return i;
	return offset;
}
template<class T, class U, class Comparator>
void bottomUpMergeSort(T ar, unsigned int n, U output, Comparator comparator) {
	bool isArSorted = true;
	unsigned int passes = 0;
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
		}
		std::swap(ar, output);
		isArSorted = (!isArSorted);
		++passes;
	}
	if (!isArSorted) {
		for (decltype(n) i = 0; i < n; ++i)
			output[i] = ar[i];
	}
	double expected = ceil(std::log(n) / std::log(2));
	cout << "bottom up actual passes : " << passes << ", expected passes : " << expected << ", ratio : " << passes / expected << endl;
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
			++ passes;
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
	double expected = ceil(std::log(n) / std::log(2));
	cout << "natural actual passes : " << passes << ", expected passes : " << expected << ", ratio : " << passes / expected << endl;
}
template<class T, class Comparator>
static bool isSorted(T values, unsigned int n, Comparator comparator) {
	return firstDisorderIndex(values, n, comparator) == 0;
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
						std::unique_lock<std::mutex> lk(va.signalData->m);
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
		std::lock_guard<std::mutex> lk(signalData.m);
		signalData.cv.notify_one();
		signalData.isReady = true;
	}
	void done() {
		std::lock_guard<std::mutex> lk(signalData.m);
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
	const unsigned int count = 1000;
	unsigned int values[count];
	unsigned int toBeSorted[count];
	unsigned int output[count];
	unsigned int maxValue = 2 * count;
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
				std::lock_guard<std::mutex> lk(m);
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
			std::lock_guard<std::mutex> lk(m);
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
	bool measureTime = true;
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
			if(!measureTime)
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
			mergeSort<Array, Array,
			std::less<unsigned int>, true, 5>,
			mergeSort<Array, Array,
			std::less<unsigned int>, true, 8>,
			mergeSort<Array, Array,
			std::less<unsigned int>, true, 11>,
			mergeSort<Array, Array,
			std::less<unsigned int>, true, 14>,
			mergeSort<Array, Array,
			std::less<unsigned int>, true, 17>,
			mergeSort<Array, Array,
			std::less<unsigned int>, true, 20>,
			mergeSort<Array, Array,
			std::less<unsigned int>, true, 23>,
		};
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
template<class T, int N>
thread_local RunningTimePlotter<T, N> *RunningTimePlotter<T, N>::instance =
		nullptr;

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
//	SpecialDistributions tmp(1);
//	RunningTimePlotter<SpecialDistributions, 3> test(tmp, { 1.0, 0.0, 0.0, 0.0,
//			1.0, 0.0, 1.0, 1.0, 0.0 });
//	return test.run(argc, argv);
//	SortingTime tmp(3);
//	RunningTimePlotter<SortingTime, 2> test(tmp,
//			{ 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 });
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
	DifferentMergeCompare tmp;
	RunningTimePlotter<DifferentMergeCompare, 3> test(tmp,
			{
										1.0, 0.0, 0.0,
										0.0, 1.0, 0.0,
										1.0, 1.0, 0.0,
			});
	return test.run(argc, argv);
}
template<int N>
struct NumberOfOrderedSubarray{
	unsigned int counts[N + 1][N + 1][N + 1];
	NumberOfOrderedSubarray(){
		for(unsigned int i = 0;i <= N;++ i)
			for(unsigned int j = 0;j <= N;++ j)
				for(unsigned int k = 0;k <= N;++ k)
					counts[i][j][k] = 0;
		counts[0][0][0] = 1;

		for(unsigned int i = 1;i <= N;++ i)
			for(unsigned int j = 0;j <= N;++ j){
				for(unsigned int k = 1;k <=N;++ k){
					//for(unsigned int s = 1;);
				}
			}
	}
};
int main(int argc, char *argv[]) {
	return 0;
}
