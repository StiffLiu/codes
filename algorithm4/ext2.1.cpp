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
	AccessCountedArray(T *val = NULL, Counter counter = Counter()) :
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
template<class T>
static bool isSorted(T values, unsigned int n) {
	for (unsigned int i = 1; i < n; ++i)
		if (values[i] < values[i - 1])
			return false;
	return true;
}
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
		assert(isSorted(values, start));
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
int doublingTestForSorts(int argc, char *argv[]) {
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
struct TTuple{
	static_assert(N > 0, "N should be greater than zero");
	T values[N];
	static const int size = N;
	typedef T type;
	template <typename... U>
	TTuple(U ... vals) : values{vals...}{

	}
};
template<class T, int N = 2>
class RunningTimePlotter {
	std::thread workingThread;
	std::mutex m;
	bool isDone = false;
	T val;
	vector<vector<double>> values;
	vector<unsigned int > counts;
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
				for(size_t i = 0;i < values.size();++ i){
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
	double getMaxValue(){
		double maxValue = 0.0;
		for(size_t i = 0;i < values.size();++ i){
			vector<double>& v = values[i];
			maxValue = std::max(maxValue, *std::max_element(v.begin(),
										v.end()));
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
		for(size_t i = 0;i < N;++ i){
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
			for(size_t i = 0;i < values.size();++ i){
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
		for(size_t i = 0;i < counts.size();++ i)
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
		clock_t begin = clock();
		switch (method) {
		case 0:
			insertionSort(values, (unsigned int) (values.size()),
					std::less<unsigned int>());
			start += step;
			return {clock() - begin, (start - step) * (start - step) / 120.0};
			break;
		case 1:
			selectionSort(values, (unsigned int) (values.size()),
					std::less<unsigned int>());
			start += step;
			return {clock() - begin, (start - step) * (start - step) / 120.0};
			break;
		case 2:
			shellSort(values, (unsigned int) (values.size()),
					std::less<unsigned int>());
			start += step;
			return {clock() - begin, (start - step) * std::pow(start - step, 1 / 2.0) / 120.0};
			break;
		}
		return {0, 0};
	}
};
class ArmortizedCost{
	unsigned int problemSize = 1000;
	unsigned int  counts = 0;
	unsigned long long totalTime = 0;
	unsigned int method;
public:
	ArmortizedCost( unsigned int method = 0, unsigned int problemSize = 1000) : problemSize(problemSize), method(method){

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
		return {duration, totalTime / (double)counts};
	}
};

class ShellSortIncrementCompare{
	//formed by merging together the se	quences "9* (4 ^ k) - 9 * (2 ^ k) + 1" and "(4 ^ k - 3 * (2 ^ k) + 1"
	TTuple<11, unsigned int> i1{1, 5, 19, 41, 109, 209, 505, 929, 2161, 3905, 8929,};
	//formed from the sequence "(3 ^ k - 1) / 2"
	TTuple<11, unsigned int> i2{1, 4, 13, 40, 121, 364, 1093, 3280, 9841, };
	//formed from the sequence "(2^(k - 1))"
	TTuple<11, unsigned int> i3{1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024,};
	unsigned int start = 100;
	template<class U>
	void reverse(U& u){
		for(unsigned int i = 0, j = U::size - 1;i < j;++ i, --j)
			swap(u.values[i], u.values[j]);
	}
public:
	ShellSortIncrementCompare(){
		reverse(i1);
		reverse(i2);
		reverse(i3);
	}
	TTuple<3> operator()() {
		unsigned int origin[start];
		unsigned int toBeSorted[start];
		clock_t duration1 = 0, duration2 = 0, duration3 = 0;
		randUInts(origin, origin + start, start * 2);
		unsigned long long counter = 0;
		double counter1 = 0, counter2 = 0, counter3 = 0;
		AccessCountedArray<unsigned int, unsigned long long *> ar(toBeSorted, &counter);

		std::copy(origin, origin + start, toBeSorted);
		counter = 0;
		duration1 = clock();
		shellSort(ar, start, i1.values, i1.size, less<unsigned int>());
		duration1 = clock() - duration1;
		std::copy(origin, origin + start, toBeSorted);
		counter1 = counter;
		counter = 0;
		duration2 = clock();
		shellSort(ar, start, i2.values, i2.size, less<unsigned int>());
		duration2 = clock() - duration2;
		std::copy(origin, origin + start, toBeSorted);
		counter2 = counter;
		counter = 0;
		duration3 = clock();
		shellSort(ar, start, i3.values, i3.size, less<unsigned int>());
		duration3 = clock() - duration3;
		counter3 = counter;
		start += 10;
		return {counter1, counter2, counter3};
	}
};
template<class T, int N>
thread_local RunningTimePlotter<T, N> *RunningTimePlotter<T, N>::instance = nullptr;
int main(int argc, char *argv[]) {
	//InsertionSortAndSelectionSortAnimation test;
	//ShellSortTraces test;
	ShellSortIncrementCompare tmp;
	RunningTimePlotter<ShellSortIncrementCompare, 3> test(tmp, {1, 0, 0, 0, 1, 0, 1, 1, 0});
	return test.run(argc, argv);
	//return testForSortingZeroOneArrays(argc, argv);
}
