#ifndef MY_LIB_NARY_HEAP_H
#define MY_LIB_NARY_HEAP_H
#include <vector>
#include <cstring>
#include <functional>
#include <algorithm>
#include <initializer_list>
#include <cassert>
#include "my_math.h"

namespace my_lib{
/**
 *  A minimum heap is a data structure that supports at least the following three operations: 
 *  <li>Insert a key into the heap.</li> 
 *  <li>Extract the minimum key from the heap.</li>
 *  <li>Remove the minimum key from the heap.</li> 
 * 
 * This class implements the minimum heap data struct and has the following two properties:
 * 	1. It uses a balanced full n-ary tree, where a key is stored in each node.
 * 	2. The key in the parent node is not greater than the key in any of its children. 
 * It utilizes a container that has <i>random-access</i> ability.
 * The indices of the children of the node with index {@var j} are :
 * 	<i>n * j + c, where 1 <= c <= d</i>
 *
 * The index of the parent of the node with index {@var j} is :
 * 	<i>(j - 1) / n, if j > 0</i>
 *
 * The node with index "0" is the root.
 * 
 * The index of a key is the index of the node which stores the key.
 *
 * The time complexity of the three operations mentioned above is O(lg(n)), where n
 * 	is the number of nodes in the heap.
 *
 * @param T Type of key stored in nodes of the heap.
 * @param Seq Type of random-access container used by this implementation.
 * @param Less The "less than" relation between any two keys.
 */
struct NaryHeapAlgorithmTraits{
	enum{ORDINARY_PROMOTE, LESS_COMPARE_PROMOTE};
	static const int promoteMethod = ORDINARY_PROMOTE;
	enum{ORDINARY_HEAPIFY, NO_SWAP_HEAPIFY};
	static const int heapifyMethod = ORDINARY_HEAPIFY; 
	static const bool useFloydMethod = true;
	static constexpr bool isPromoteMethodValid(int method){
		return method == ORDINARY_PROMOTE || method == LESS_COMPARE_PROMOTE;
	}
	static constexpr bool isHeapifyMethodValid(int method){
		return method == ORDINARY_HEAPIFY || method == NO_SWAP_HEAPIFY;
	}
};
template<class T, class Seq = std::vector<T>, class Less = std::less<T>, 
	class Traits = NaryHeapAlgorithmTraits >
class NaryHeap {
protected:
	/**
	 * The underlying random-access container used by the heap.
	 */
	Seq data;
	/**
	 * Number of nodes in the heap.
	 */
	size_t s;
	/**
	 * The "less than" relation between two keys.
	 */
	Less comparator;
	/**
	 * Dimension of the heap, which is the maximum number of children 
	 * 	each node could has.
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
	size_t promote(size_t index){
		static_assert(NaryHeapAlgorithmTraits::isPromoteMethodValid(Traits::promoteMethod),
			"unsupported promoting method");
		if(Traits::promoteMethod == NaryHeapAlgorithmTraits::ORDINARY_PROMOTE)
			return promoteOrdinary(index);
		else if(Traits::promoteMethod == NaryHeapAlgorithmTraits::LESS_COMPARE_PROMOTE)
			return promoteWithLessCompare(index);
	}
	size_t promoteOrdinary(size_t index) {
		if (index <= 0)
			return index;
		size_t parent = (index - 1) / d;
		T item = data[index];
		while (parent > 0 && comparator(item, data[parent])) {
			data[index] = data[parent];
			index = parent;
			parent = (index - 1) / d;
		}
		if (parent == 0 && comparator(item, data[parent])){
			data[index] = data[parent];
			index = parent;
		}
		data[index] = item;		
		return index;
	}
	/**
	 * an optimization that uses less than {@code lg(lg(n))} compares, where 
	 * 	{@var n} is the number of keys in the heap.
	 */
	size_t promoteWithLessCompare(size_t index){
		auto n = ceilLog(d, index + 1);
		size_t ancestors[n + 1];
		if(n <= 0)
			return index;
		int i = n;
		ancestors[i] = index;
		while(ancestors[i] > 0){
			-- i;
			ancestors[i] = (ancestors[i + 1] - 1) / d;
		}
		ancestors[i] = 0;

		auto end = ancestors + n;
		auto pos = std::upper_bound(ancestors + i, end, index, 
			[&](size_t i, size_t j){ return this->comparator(data[i], data[j]);});
		if(pos != end){
			T item = data[*end];
			while(end > pos){ 
				-- end;
				data[*(end + 1)] = data[*end];
			}
			data[*end] = item;
			index = *pos;
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
	virtual size_t heapify(size_t index){
		assert(index < s);
		static_assert(NaryHeapAlgorithmTraits::isHeapifyMethodValid(Traits::heapifyMethod),
			"unsupported heapify method");
		if(Traits::heapifyMethod == NaryHeapAlgorithmTraits::NO_SWAP_HEAPIFY)
			return heapifyNoSwap(index);
		else if(Traits::heapifyMethod == NaryHeapAlgorithmTraits::ORDINARY_HEAPIFY)
			return heapifyOrdinary(index);
	}
	size_t heapifyOrdinary(size_t index){
		while(index < s){
			size_t start = d * index + 1;
			size_t end = std::min(s, start + d);
			size_t smallest = index;
			for(size_t i = start;i < end;++ i)
				if(comparator(data[i], data[smallest]))
					smallest = i;
			if(index == smallest)
				break;
			std::swap(data[index], data[smallest]);
			index = smallest;
		}
		return index;
	}
	/**
	 * Heapify without swapping keys.
	 */
	size_t heapifyNoSwap(size_t index) {
		T item = data[index];
		while (index < s) {
			size_t start = d * index + 1;
			if(start >= s){
				data[index] = item;
				break;
			}
			
			size_t end = std::min(s, start + d);
			size_t smallest = start;
			for (size_t i = start + 1; i < end; ++i)
				if (comparator(data[i], data[smallest]))
					smallest = i;
			if (!comparator(data[smallest], item)){
				data[index] = item;
				break;
			}
			data[index] = data[smallest];
			index = smallest;
		}
		assert(index < s);
		return index;
	}
	/**
	 * Heapfiy without swapping keys and use Floyd's method to save the number of compares.
	 * This method could only be used when the heap is already constructed, that is it should
	 * not be called directly or indirectly by {@code makeHeap}.
	 */
	size_t heapifyNoSwapFloyd(size_t index) {
		T item = data[index];
		while (index < s) {
			size_t start = d * index + 1;
			if(start >= s){
				data[index] = item;
				index = promote(index);
				break;
			}
			
			size_t end = std::min(s, start + d);
			size_t smallest = start;
			for (size_t i = start + 1; i < end; ++i)
				if (comparator(data[i], data[smallest]))
					smallest = i;
			data[index] = data[smallest];
			index = smallest;
		}
		assert(index < s);
		return index;
	}
public:
	/**
	 * Type definition, type of the key
	 */
	typedef T KeyType;
	/**
	 * Consturct a d-ary heap.
	 * @param d Dimension of the heap, must be greater than 2.
	 * 	If {@var d} is less than 2, 2 will be used.
	 */
	NaryHeap(unsigned int d, Less comparator = Less()) : comparator(comparator) {
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
	NaryHeap(unsigned int d, ForwardIterator begin, ForwardIterator end, 
		Less comparator = Less()) : comparator(comparator) {
		this->d = ((d < 2) ? 2 : d);
		s = 0;
		makeHeap(begin, end);
	}
	
	NaryHeap(unsigned int d, std::initializer_list<T> l, Less comparator = Less()) 
		: NaryHeap(d, l.begin(), l.end(), comparator){
	}

	virtual ~NaryHeap(){
	}

	/**
	 * Make the keys in the range {@code [begin, end)} into a heap.
	 * @param begin A forward iterator that represents the beginning of the range, inclusive.
	 * @param end A forward iterator that represents the ending of the range, exclusive.
	 */
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
	
	/**
	 * Remove the smallest key from the heap.
	 */
	 void pop() {
		if (empty())
			return;
		--s;
		if(s > 0){
			data[0] = data[s];
			if(Traits::useFloydMethod)
				heapifyNoSwapFloyd(0);
			else
				heapify(0);
		}
		data.pop_back();
	}

	/**
	 * @return The smallest key in the heap.
	 * @note If the heap is empty, the behaviour is undefined.
	 */
	const T& top() const {
		return data[0];
	}

	/**
	 *@return {@code true} if the heap is empty, false otherwise.
	 */
	bool empty() const {
		return s == 0;
	}

	/**
	 * @param index The index of node whose key will be returned.
	 * @return The key in the node with index {@var index}.
	 * @note {@var index} must be in the range {@code [0, size())}, else 
	 * 	the behaviour is undefined.
	 */
	const T& get(size_t index) const {
		return data[index];
	}

	/**
	 * Set the key in the node with the given index.
	 * @param index Index of the node.
	 * @param item New key for the node.
	 * @return After the key of the node is set, the 2nd property of 
	 * 	may be violated, and the key will be moved to revalidate
	 * 	the 2nd property. The returned value is the final index of 
	 * 	the key.
	 * @note {@var index} must be in the range {@code [0, size())}.
	 */
	size_t set(size_t index, const T& item) {
		data[index] = item;

		size_t ret = promote(index);
		if(ret == index)
			return Traits::useFloydMethod ? heapifyNoSwapFloyd(index) 
					: heapify(index);
		return ret;
	}

	/*
	 * Insert the key {@var item} into the heap.
	 * @param item The key.
	 * @return Index of the key inserted.
	 */
	size_t add(const T& item) {
		data.push_back(item);
		++s;
		return promote(s - 1);
	}

	/*
	 * Remove the key at the index {@var index} from the heap.
	 * @param index Index of the key to be removed.
	 * @note {@var index} must be in the range {@code [0, size())}.
	 */
	void remove(size_t index) {
		if (index >= s)
			return;
		if (s <= 0)
			return;
		--s;
		set(index, data[s]);
		data.pop_back();
	}
	/**
	 * @return The number of keys in the heap.
	 */
	size_t size() const {
		return s;
	}
	/**
	 * @return The less than relation used by this heap.
	 */
	const Less& getComparator()const{
		return comparator;
	}

	/*
	 * get the underlying sequence used ty this heap.
	 */
	const Seq& getSeq(){
		return seq;
	}
};
template<class T, class Seq = std::vector<T>, class Less = std::less<T> >
class BinaryHeap: public NaryHeap<T, Seq, Less> {
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
public:
	typedef NaryHeap<T, Seq, Less> SuperClass;
	BinaryHeap() :
			SuperClass(2) {
	}
	template<class ForwardIterator>
	BinaryHeap(ForwardIterator begin, ForwardIterator end) :
			SuperClass(2, begin, end) {
	}
};

}
#endif //MY_LIB_NARY_HEAP_H
