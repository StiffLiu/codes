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
#include "common.h"
using namespace std;
/**
 * This template class is for performance test purpose.
 * It's a wrapper for arrays and it counts the number
 * of array accesses. This number is used to measure
 * the performance of an algorithm.
 */
template<class T>
class AccessCountedArray {
	T *val;
	unsigned int *counter;
public:
	AccessCountedArray(T *val = NULL, unsigned int *counter = NULL) :
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
		std::random_device rd;
		vector<decltype(counter)> countOfEachIteration;
		ShellSortIncrements increments(&incrementValues[0],
				&countOfEachIteration, &counter);
		std::uniform_int_distribution<decltype(counter)> generator(0, maxValue);
		std::generate(values, values + start, [&] {return generator(rd);});
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
	void problemSize(unsigned int size) {
		values.clear();
		for (int i = 0; i < 10; ++i) {
			values.push_back(vector<unsigned int>());
			randUInts(values.back(), size, 2 * size);
		}
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
		cout << "duration : " << duration << ", average : "
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
}
int main(int argc, char *argv[]) {
	//InsertionSortAndSelectionSortAnimation test;
	//ShellSortTraces test;
	//return test.run(argc, argv);
	return doublingTestForSorts(argc, argv);
}
