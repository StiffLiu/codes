#ifndef MY_LIB_NARY_MIN_MAX_HEAP_H
#define MY_LIB_NARY_MIN_MAX_HEAP_H
#include <vector>
#include <cstring>
#include <functional>
#include <algorithm>
#include <initializer_list>

namespace my_lib{
/**
 *  A minmax heap is a data structure that supports at least the following five operations: 
 *  <li>Insert a key into the heap.</li> 
 *  <li>Extract the minimum key from the heap.</li>
 *  <li>Remove the minimum key from the heap.</li> 
 *  <li>Extract the maximum key from the heap.</li>
 *  <li>Remove the maximum key from the heap.</li>
 * 
 * this class implements the minimum heap data struct and has the following two properties:
 * 	1. it uses a balanced full n-ary tree, where a key is stored in each node.
 * 	2. the key in the parent node is not greater than the key in any of its children. 
 * it utilizes a container that has <i>random-access</i> ability.
 * the indices of the children of the node with index {@var j} are :
 * 	<i>n * j + c, where 1 <= c <= d</i>
 *
 * the index of the parent of the node with index {@var j} is :
 * 	<i>(j - 1) / n, if j > 0</i>
 *
 * the node with index "0" is the root.
 * 
 * the index of a key is the index of the node which stores the key.
 *
 * the time complexity of the three operations mentioned above is o(lg(n)), where n
 * 	is the number of nodes in the heap.
 *
 * @param T Type of key stored in nodes of the heap.
 * @param Seq Type of random-access container used by this implementation.
 * @param Less The "less than" relation between any two keys.
 */
template<class T, class Seq = std::vector<T>, class Less = std::less<T> >
class NaryMinMaxHeap {
protected:
	/**
	 * The underlying random-access container used by the heap.
	 */
	Seq data;
	std::vector<size_t> indexToMinHeap;
	std::vector<size_t> indexToMaxHeap;
	std::vector<size_t> minHeap;
	std::vector<size_t> maxHeap;	
	/**
	 * Number of nodes in the heap.
	 */
	size_t s;
	/**
	 * The "less than" relation between two keys.
	 */
	Less comparator;
	/**
	 * Dimension of the heap, which is the maximum number of children each node could has.
	 */
	unsigned int d;
	/*
	 *  Assume that every node, except that the key in the node at index {@var index} may be smaller than
	 *  	that the key in its parent, statisfies the 2nd property mentioned above.
	 *  Starting from the node at index {@var index}, this function iteratively promoting the key by 
	 *  	swapping it with the key in the parent node and stops when it's not smaller than the key 
	 *  	in the parent node or when the root is reached.
	 *  @param index Index of the node, from which promoting starts.
	 *  @return Index of the node, at which promoting ends.
	 */
	size_t promoteMin(size_t index) {
		if (index <= 0)
			return index;
		size_t parent = (index - 1) / d;
		while (parent > 0 && comparator(data[minHeap[index]], data[minHeap[parent]])) {
			std::swap(indexToMinHeap[minHeap[index]], indexToMinHeap[minHeap[parent]]);
			std::swap(minHeap[index], minHeap[parent]);
			index = parent;
			parent = (index - 1) / d;
		}
		if (parent == 0 && comparator(data[minHeap[index]], data[minHeap[parent]])){
			std::swap(indexToMinHeap[minHeap[index]], indexToMinHeap[minHeap[parent]]);
			std::swap(minHeap[index], minHeap[parent]);
			index = parent;
		}
		return index;
	}
	size_t promoteMax(size_t index) {
		if (index <= 0)
			return index;
		size_t parent = (index - 1) / d;
		while (parent > 0 && comparator(data[maxHeap[parent]], data[maxHeap[index]])) {
			std::swap(indexToMaxHeap[maxHeap[index]], indexToMaxHeap[maxHeap[parent]]);
			std::swap(maxHeap[index], maxHeap[parent]);
			index = parent;
			parent = (index - 1) / d;
		}
		if (parent == 0 && comparator(data[maxHeap[parent]], data[maxHeap[index]])){
			std::swap(indexToMaxHeap[maxHeap[index]], indexToMaxHeap[maxHeap[parent]]);
			std::swap(maxHeap[index], maxHeap[parent]);
			index = parent;
		}
		return index;
	}
	/*
	 *  Assume that every node, except that the key in the node at index {@var index} may be greater than
	 *  	that the keys in its children, statisfies the 2nd property mentioned above.
	 *  Starting from the node at index {@var index}, this function iteratively sinking a key by 
	 *  	swapping it with the smallest key in its children and stops when it's not greater than
	 *  	any key in its children or when a leaf is reached.
	 *  @param index Index of the node, from which sinking starts.
	 *  @return Index of the node, at which sinking ends.
	 */
	virtual size_t heapifyMin(size_t index) {
		while (index < s) {
			size_t smallest = index;
			size_t start = d * index + 1;
			size_t end = std::min(s, start + d);
			for (size_t i = start; i < end; ++i)
				if (comparator(data[minHeap[i]], data[minHeap[smallest]]))
					smallest = i;

			if (smallest == index)
				break;
			std::swap(indexToMinHeap[minHeap[index]], indexToMinHeap[minHeap[smallest]]);
			std::swap(minHeap[index], minHeap[smallest]);
			index = smallest;
		}
		return index;
	}
	virtual size_t heapifyMax(size_t index) {
		while (index < s) {
			size_t smallest = index;
			size_t start = d * index + 1;
			size_t end = std::min(s, start + d);
			for (size_t i = start; i < end; ++i)
				if (comparator(data[maxHeap[smallest]], data[maxHeap[i]]))
					smallest = i;

			if (smallest == index)
				break;
			std::swap(indexToMaxHeap[maxHeap[index]], indexToMaxHeap[maxHeap[smallest]]);
			std::swap(maxHeap[index], maxHeap[smallest]);
			index = smallest;
		}
		return index;
	}
public:
	bool isValid(){
		for(size_t i = 0;i < s;++ i){
			if(indexToMinHeap[minHeap[i]] != i || minHeap[indexToMinHeap[i]] != i ||
				indexToMaxHeap[maxHeap[i]] != i || maxHeap[indexToMaxHeap[i]] != i)
				return false;
		}
		return true;
	}
	/**
	 * Type definition, type of the key
	 */
	typedef T KeyType;
	/**
	 * Consturct a d-ary heap.
	 * @param d Dimension of the heap, must be greater than 2.
	 * 	If {@var d} is less than 2, 2 will be used.
	 */
	NaryMinMaxHeap(unsigned int d, const Less& comparator = Less()) : comparator(comparator) {
		this->d = ((d < 2) ? 2 : d);
		s = 0;
	}
	/**
	 * Construct a d-ary heap with keys from the range {@code [begin, end)}.
	 * @param d Dimension of the heap, must be greater than 2.
	 * 	If {@var d} is less than 2, 2 will be used.
	 * @param begin A forward iterator that represents the beginning of the range, inclusive.
	 * @param end A forward iterator that represents the ending of the range, exclusive.
	 */
	template<class ForwardIterator>
	NaryMinMaxHeap(unsigned int d, ForwardIterator begin, ForwardIterator end, 
		const Less& comparator = Less()) : comparator(comparator) {
		this->d = ((d < 2) ? 2 : d);
		s = 0;
		makeHeap(begin, end);
	}
	
	NaryMinMaxHeap(unsigned int d, std::initializer_list<T> l, const Less& comparator = Less()) 
		: NaryMinMaxHeap(d, l.begin(), l.end(), comparator){
	}

	virtual ~NaryMinMaxHeap(){
	}
	/**
	 * Make the keys in the range {@code [begin, end)} into a min-max heap.
	 * @param begin A forward iterator that represents the beginning of the range, inclusive.
	 * @param end A forward iterator that represents the ending of the range, exclusive.
	 */
	template<class ForwardIterator>
	void makeHeap(ForwardIterator begin, ForwardIterator end) {
		data.clear();
		minHeap.clear();
		maxHeap.clear();
		indexToMinHeap.clear();
		indexToMaxHeap.clear();
		s = 0;
		while (begin != end) {
			data.push_back(*begin);
			++s;
			++begin;
		}
		for(size_t i = 0;i < s;++ i){
			minHeap.push_back(i);
			indexToMinHeap.push_back(i);
			maxHeap.push_back(i);
			indexToMaxHeap.push_back(i);
		}
		if (s > 1) {
			size_t maxParent = (s - 2) / d;
			for (size_t i = maxParent;i > 0; --i){
				heapifyMin(i);
				heapifyMax(i);
			}
			heapifyMin(0);
			heapifyMax(0);
		}
	}
	
	/**
	 * Remove the smallest key from the heap.
	 */
	 void popMin() {
		if (empty())
			return;
		--s;
		if(s > 0){
			size_t index = minHeap[0];
			minHeap[0] = minHeap[s];
			indexToMinHeap[minHeap[0]] = 0;
			
			size_t indexToMax = indexToMaxHeap[index];
			maxHeap[indexToMax] = maxHeap[s];
			indexToMaxHeap[maxHeap[indexToMax]] = indexToMax;
			if(promoteMax(indexToMax) == indexToMax)
				heapifyMax(indexToMax);
			heapifyMin(0);
			if(index != s){
				size_t lastKeyIndex;
				lastKeyIndex = indexToMinHeap[s];
				minHeap[lastKeyIndex] = index;
				indexToMinHeap[index] = lastKeyIndex;

				lastKeyIndex = indexToMaxHeap[s];
				maxHeap[lastKeyIndex] = index;
				indexToMaxHeap[index] = lastKeyIndex;
				std::swap(data[index], data[s]);
			}
		}
		data.pop_back();
		indexToMinHeap.pop_back();
		indexToMaxHeap.pop_back();
		minHeap.pop_back();
		maxHeap.pop_back();
	}

	 void popMax() {
		if (empty())
			return;
		--s;
		if(s > 0){
			size_t index = maxHeap[0];
			maxHeap[0] = maxHeap[s];
			indexToMaxHeap[maxHeap[0]] = 0;
			
			size_t indexToMin = indexToMinHeap[index];
			minHeap[indexToMin] = minHeap[s];
			indexToMinHeap[minHeap[indexToMin]] = indexToMin;
			if(promoteMin(indexToMin) == indexToMin)
				heapifyMin(indexToMin);
			heapifyMax(0);
			if(index != s){
				size_t lastKeyIndex;
				lastKeyIndex = indexToMinHeap[s];
				minHeap[lastKeyIndex] = index;
				indexToMinHeap[index] = lastKeyIndex;

				lastKeyIndex = indexToMaxHeap[s];
				maxHeap[lastKeyIndex] = index;
				indexToMaxHeap[index] = lastKeyIndex;
				std::swap(data[index], data[s]);
			}
		}
		data.pop_back();
		indexToMinHeap.pop_back();
		indexToMaxHeap.pop_back();
		minHeap.pop_back();
		maxHeap.pop_back();
	}
	 
	/**
	 * @return The smallest key in the heap.
	 * @note If the heap is empty, the behaviour is undefined.
	 */
	const T& min() const {
		return data[minHeap[0]];
	}
	const T& max() const{
		return data[maxHeap[0]];
	}
	/**
	 *@return {@code true} if the heap is empty, false otherwise.
	 */
	bool empty() const {
		return s == 0;
	}
	/*
	 * Insert the key {@var item} into the heap.
	 * @param item The key.
	 * @return Index of the key inserted.
	 */
	void add(const T& item) {
		data.push_back(item);
		++s;
		promoteMin(s - 1);
		promoteMax(s - 1);
	}

	/**
	 * @return The number of nodes in the heap.
	 */
	size_t size() const {
		return s;
	}
};
/*template<class T, class Seq = std::vector<T>, class Less = std::less<T> >
class BinaryHeap: public NaryMinMaxHeap<T, Seq, Less> {
protected:
	virtual size_t heapify(size_t index) {
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
		return index;
	}
	typedef NaryMinMaxHeap<T, Seq, Less> SuperClass;
public:
	BinaryHeap() :
			SuperClass(2) {
	}
	template<class ForwardIterator>
	BinaryHeap(ForwardIterator begin, ForwardIterator end) :
			SuperClass(2, begin, end) {
	}
};
*/
}
#endif //MY_LIB_NARY_MIN_MAX_HEAP_H
