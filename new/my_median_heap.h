#ifndef MY_LIB_NARY_MEDIAN_HEAP_H
#define MY_LIB_NARY_MEDIAN_HEAP_H
#include "my_heap.h"

namespace my_lib{
template<class T, class Seq = std::vector<T>, class Less = std::less<T> >
class NaryMedianHeap{
protected:
	Less comparator;
	struct LeftHalfLess{
		NaryMedianHeap* heap;
		LeftHalfLess(NaryMedianHeap* heap = nullptr) : heap(heap){
		}
		bool operator()(const T& val1, const T& val2){
			return heap->comparator(val2, val1);
		}
	};
	struct RightHalfLess{
		NaryMedianHeap* heap;
		RightHalfLess(NaryMedianHeap* heap = nullptr) : heap(heap){
		}
		bool operator()(const T& val1, const T& val2){
			return heap->comparator(val1, val2);
		}
	};
	NaryHeap<T, Seq, LeftHalfLess> leftHalf;
	NaryHeap<T, Seq, RightHalfLess> rightHalf;
public:
	NaryMedianHeap(unsigned int d, const Less& comparator = Less()) 
		:comparator(comparator), leftHalf(d, LeftHalfLess(this)), 
			rightHalf(d, RightHalfLess(this)){
	}
	NaryMedianHeap(unsigned int d, std::initializer_list<T> l, const Less& comparator = Less())
		: NaryMedianHeap(d, l.begin(), l.end(), comparator){
	}
	template<class ForwardIterator>
	NaryMedianHeap(unsigned int d, ForwardIterator begin, ForwardIterator end, 
		const Less& comparator = Less()) :comparator(comparator), 
			leftHalf(d, LeftHalfLess(this)), rightHalf(d, RightHalfLess(this)){ 
		while(begin != end){
			add(*begin);
			++ begin;
		}
	}
	bool isValid(){
		return leftHalf.size() == rightHalf.size() || 
			leftHalf.size() + 1 == rightHalf.size() ||
			leftHalf.size() == rightHalf.size() + 1;
	}
	void pop(){
		if(rightHalf.size() >= leftHalf.size()){
			rightHalf.pop();
			return;
		}
		leftHalf.pop();

	}
	const T& top(){
		if(rightHalf.size() >= leftHalf.size()){
			return rightHalf.top();
		}
		return leftHalf.top();
	}
	void add(const T& item){
		if(!leftHalf.empty() && comparator(item, leftHalf.top())){
			rightHalf.add(leftHalf.top());
			leftHalf.pop();
			leftHalf.add(item);
		}else{
			rightHalf.add(item);
		}
		if(rightHalf.size() > leftHalf.size() + 1){
			leftHalf.add(rightHalf.top());
			rightHalf.pop();	
		}
	}
	
	bool empty(){
		return leftHalf.empty() && rightHalf.empty();
	}
	size_t size(){
		return leftHalf.size() + rightHalf.size();
	}
	
};
}
#endif //MY_LIB_NARY_MEDIAN_HEAP_H
