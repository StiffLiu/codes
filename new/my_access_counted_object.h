#ifndef MY_LIB_MY_ACCESS_COUNTED_OBJECT_H
#define MY_LIB_MY_ACCESS_COUNTED_OBJECT_H
#include <cassert>
namespace my_lib{

/**
 * A class that contains a counter.
 * The counter could be used to count the number of array accessed, 
 * 	the number that a function(or functor) is called and 
 * 	the number a iterator is dereferenced etc.
 */
template<class Counter>
class AccessCountedObject{
protected:
	mutable Counter counter;
public:
	AccessCountedObject(const Counter counter = Counter()) : counter(counter){
	}
	const Counter& getCounter() const{
		return counter;
	}
	void setCounter(const Counter& c) const{
		this->counter = c;
	}
};

template<class T>
class Incrementor{
	T *ptr;
public:
	Incrementor(T *ptr = nullptr) : ptr(ptr){
	}
	operator const T&() const{
		assert(ptr != nullptr);
		return *ptr;
	}
	operator T&(){
		assert(ptr != nullptr);
		return *ptr;
	}
	Incrementor& operator++(){
		assert(ptr != nullptr);
		++ *ptr;
		return *this;
	}
	Incrementor& operator++(int){
		assert(ptr != nullptr);
		++ *ptr;
		return *this;
	}
};


template<class Counter, class T>
class AccessCountedArray : public AccessCountedObject<Counter>{
	typedef AccessCountedObject<Counter> Super;
protected:
	T array;
public:
	AccessCountedArray(const Counter& counter = Counter(), const T& array = T()) 
		: Super(counter), array(array){
	}
	template<class Index>
	auto operator[](Index index) ->decltype((array[index])){
		++ Super::counter;
		return array[index];
	}
	template<class Index>
	auto operator[](Index index) const ->decltype((array[index])){
		++ Super::counter;
		return array[index];
	}
};
template<class Counter, class Functor>
class AccessCountedFunctor : public AccessCountedObject<Counter>{
	typedef AccessCountedObject<Counter> Super;
protected:
	Functor functor;
public:
	AccessCountedFunctor(const Counter& counter = Counter(), const Functor& functor = Functor())
	       : Super(counter), functor(functor){
	}
	/*auto operator()(){
		++ Super::counter;
		return functor();
	}*/
	template<class T>
	auto operator()(const T& val)->decltype(functor(val)){
		++ Super::counter;
		return functor(val);
	}
	template<class T1, class T2>
	auto operator()(const T1& val1, const T2& val2)->decltype(functor(val1, val2)){
		++ Super::counter;
		return functor(val1, val2);
	}
	template<class T1, class T2, class T3>
	auto operator()(const T1& val1, const T2& val2, const T3& val3)
		->decltype(functor(val1, val2, val3)){
		++ Super::counter;
		return functor(val1, val2);
	}
	template<class T1, class T2, class T3, class T4>
	auto operator()(const T1& val1, const T2& val2, const T3& val3, const T4& val4)
		->decltype(functor(val1, val2, val3, val4)){
		++ Super::counter;
		return functor(val1, val2, val3, val4);
	}
	/*auto operator()()const{
		++ Super::counter;
		return functor();
	}*/
	template<class T>
	auto operator()(const T& val)const ->decltype(functor(val)){
		++ Super::counter;
		return functor(val);
	}
	template<class T1, class T2>
	auto operator()(const T1& val1, const T2& val2)const ->decltype(functor(val1, val2)){
		++ Super::counter;
		return functor(val1, val2);
	}
	template<class T1, class T2, class T3>
	auto operator()(const T1& val1, const T2& val2, const T3& val3)const
		->decltype(functor(val1, val2, val3)){
		++ Super::counter;
		return functor(val1, val2);
	}
	template<class T1, class T2, class T3, class T4>
	auto operator()(const T1& val1, const T2& val2, const T3& val3, const T4& val4)const
		->decltype(functor(val1, val2, val3, val4)){
		++ Super::counter;
		return functor(val1, val2, val3, val4);
	}
};
}
#endif //MY_LIB_MY_ACCESS_COUNTED_OBJECT_H
