#include <vector>
#include <functional>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <cassert>
using namespace std;

template<class T, class Seq = std::vector<T>, class Less = std::less<T> >
class NaryHeap {
protected:
	Seq data;
	size_t s;
	Less comparator;
	unsigned int d;
	void promote(size_t index) {
		if (index <= 0)
			return;
		size_t parent = (index - 1) / d;
		while (parent > 0 && comparator(data[index], data[parent])) {
			swap(data[index], data[parent]);
			index = parent;
			parent = (index - 1) / d;
		}
		if (parent == 0 && comparator(data[index], data[parent]))
			swap(data[index], data[parent]);
	}
public:
	typedef T DataType;
	NaryHeap(unsigned int d) {
		this->d = ((d < 2) ? 2 : d);
		s = 0;
	}
	template<class ForwardIterator>
	NaryHeap(unsigned int d, ForwardIterator begin, ForwardIterator end) {
		this->d = ((d < 2) ? 2 : d);
		s = 0;
		makeHeap(begin, end);
	}
	template<class ForwardIterator>
	void makeHeap(ForwardIterator begin, ForwardIterator end) {
		data.clear();
		s = 0;
		while (begin != end) {
			data.push_back(*begin);
			++s;
			++begin;
		}
		if (s > 1) {
			size_t maxParent = (s - 2) / d;
			size_t i = maxParent;
			for (; i > 0; --i)
				heapify(i);
			heapify(i);
		}
	}
	virtual void heapify(size_t index) {
		while (index < s) {
			size_t smallest = index;
			size_t start = d * index + 1;
			size_t end = min(s, start + d);
			for (size_t i = start; i < end; ++i)
				if (comparator(data[i], data[smallest]))
					smallest = i;

			if (smallest == index)
				break;
			swap(data[index], data[smallest]);
			index = smallest;
		}
	}
	bool extract(T& ret) {
		if (empty())
			return false;
		ret = data[0];
		--s;
		data[0] = data[s];
		data.pop_back();
		heapify(0);
		return true;
	}
	const T& top() const {
		return data[0];
	}
	bool empty() const {
		return s == 0;
	}
	const T& get(size_t index) const {
		return data[index];
	}
	void set(size_t index, const T& item) {
		if (comparator(item, data[index])) {
			data[index] = item;
			promote(index);
		} else if (comparator(data[index], item)) {
			data[index] = item;
			heapify(index);
		} else {
			data[index] = item;
		}
	}
	void add(const T& item) {
		data.push_back(item);
		++s;
		promote(s - 1);
	}
	void remove(size_t index) {
		if (index >= s)
			return;
		if (s <= 0)
			return;
		--s;
		set(index, data[s]);
		data.pop_back();
	}
	size_t size() const {
		return s;
	}
	template<class U>
	void sort(U& u){
		while(!empty()){
			DataType val;
			extract(val);
			u.push_back(val);
		}
	}
};
template<class T, class Seq = std::vector<T>, class Less = std::less<T> >
class BinaryHeap: public NaryHeap<T, Seq, Less> {
public:
	typedef NaryHeap<T, Seq, Less> SuperClass;
	BinaryHeap() :
			SuperClass(2) {
	}
	template<class ForwardIterator>
	BinaryHeap(ForwardIterator begin, ForwardIterator end) :
			SuperClass(2, begin, end) {
	}
	void heapify(size_t index) {
		size_t s = SuperClass::size();
		Seq& data = SuperClass::data;
		while (index < s) {
			size_t smallest = index;
			size_t l = 2 * index + 1;
			size_t r = 2 * index + 2;
			if (l < s && SuperClass::comparator(data[l], data[smallest]))
				smallest = l;
			if (r < s && SuperClass::comparator(data[r], data[smallest]))
				smallest = r;
			if (smallest == index)
				break;
			swap(data[index], data[smallest]);
			index = smallest;
		}
	}
};
template<class T>
struct Range{
	T start;
	T end;
	bool operator<(const Range& another) const{
		if(start == end)
			return another.start != another.end;
		if(another.start == another.end)
			return false;
		return *start < *another.start;
	}
};
template<class T,class OutputIterator>
void kWayMerge(Range<T> *ranges, int k, OutputIterator output){
	typedef BinaryHeap<Range<T> > Heap;
	Heap heap(ranges, ranges + k);
	while(!heap.empty()){
		typename Heap::DataType val;
		heap.extract(val);
		if(val.start != val.end){
			*output = *val.start;
			++output;
			++ val.start;
		}
		if(val.start != val.end)
			heap.add(val);
	}
}
template<class Heap>
void testHeap(const vector<int>& numbers, const vector<int>& toAdd,
		const vector<int>& toRemove, Heap& heap) {
	vector<int> sorted(numbers);
	vector<int> heapSorted;
	cout << "original: ";
	std::sort(sorted.begin(), sorted.end());
	copy(numbers.begin(), numbers.end(), ostream_iterator<int>(cout, " "));
	cout << endl;
	cout << "original sorted: ";
	copy(sorted.begin(), sorted.end(), ostream_iterator<int>(cout, " "));
	cout << endl;
	cout << "after make heap: ";
	heap.makeHeap(numbers.begin(), numbers.end());
	heap.sort(heapSorted);
	copy(heapSorted.begin(), heapSorted.end(), ostream_iterator<int>(cout, " "));
	cout << endl;

	assert(sorted == heapSorted);
	cout << "after new elements added: ";
	heap.makeHeap(numbers.begin(), numbers.end());
	for(size_t i = 0;i < toAdd.size();++ i)
		heap.add(toAdd[i]);
	heapSorted.clear();
	heap.sort(heapSorted);
	copy(heapSorted.begin(), heapSorted.end(), ostream_iterator<int>(cout, " "));
	cout << endl;

	cout << "after elements removed: ";
	heap.makeHeap(numbers.begin(), numbers.end());
	for(size_t i = 0;i < toAdd.size();++ i)
			heap.remove(toRemove[i]);
	heapSorted.clear();
	heap.sort(heapSorted);
	copy(heapSorted.begin(), heapSorted.end(), ostream_iterator<int>(cout, " "));
	cout << endl;
}
int test21(int argc, char *argv[]) {
	const int numNumbers = 20;
	BinaryHeap<int> binHeap;
	NaryHeap<int> ternaryHeap(3);
	NaryHeap<int> naryHeap(5);
	vector<int> numbers;
	vector<int> toAdd;
	vector<int> toRemove;
	srand(time(0));
	for (int i = 0; i < numNumbers; ++i)
		numbers.push_back(rand() % 100);
	toAdd.push_back(rand() % 100);
	toRemove.push_back(rand() % numNumbers);
	cout << "binary heap:" << endl;
	testHeap(numbers, toAdd, toRemove, binHeap);
	cout << endl << "ternary heap:" << endl;
	testHeap(numbers, toAdd, toRemove, ternaryHeap);
	cout << endl << "5-ary heap:" << endl;
	testHeap(numbers, toAdd, toRemove, naryHeap);
	return 0;
}
int test22(int argc, char *argv[]){
	const int k = 10;
	const int least = 10;
	vector<vector<int> > numbers;
	vector<Range<vector<int>::iterator> > ranges;
	int totalNum = 0;
	numbers.resize(k);
	ranges.resize(k);
	for(int i = 0;i < k;++ i){
		vector<int>& one = numbers[i];
		int num = least + rand() % 30;
		for(int j = 0;j < num;++ j)
			one.push_back(rand() % 1000);
		sort(one.begin(), one.end());
		totalNum += num;
		ranges[i].start = one.begin();
		ranges[i].end = one.end();
		cout << "array " << (i + 1) << " \t:\t";
		copy(one.begin(), one.end(), ostream_iterator<int>(cout, " "));
		cout << endl;
	}

	vector<int> output;
	output.resize(totalNum);
	kWayMerge(&ranges[0], k, output.begin());

	cout << "merged :\t";
	copy(output.begin(), output.end(), ostream_iterator<int>(cout, " "));
	cout << endl;
	return 0;
}
