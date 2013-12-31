/*
 * SpecialDoublyLinkedList.h
 *
 *  Created on: 2013Äê11ÔÂ26ÈÕ
 *      Author: tmlg
 */
#include <cassert>
#ifndef SPECIALDOUBLYLINKEDLIST_H_
#define SPECIALDOUBLYLINKEDLIST_H_
template<class T, typename Alloc =  std::allocator<T> >
class SpecialDoublyLinkedList {
	struct Link{
		T data;
		uintptr_t np;
		Link(const T& val, uintptr_t prev, uintptr_t next) : data(val), np(prev ^ next){
		}
		uintptr_t getAnother(uintptr_t p){
			return np ^ p;
		}
	};
	typedef typename Alloc::template rebind<Link >::other LinkAlloc;
	struct LinkImpl : public LinkAlloc{
	public:
	Link *head;
	Link *tail;
	size_t size;
	LinkImpl() : head(NULL), tail(NULL), size(0){
	}
	} impl;
public:
	template<class U>
	class _iterator{
	protected:
		U prev;
		U current;
		_iterator() : prev(NULL), current(NULL){
		}
	public:
		typedef ptrdiff_t                          difference_type;
		typedef std::bidirectional_iterator_tag    iterator_category;
    	typedef T                                value_type;
    	typedef const T*                         pointer;
    	typedef const T&                         reference;

		bool operator==(const _iterator& next){
			return prev == next.prev && current == next.current;
		}
		bool operator!=(const _iterator& next){
			return !this->operator==(next);
		}
		_iterator& operator ++ (){
			U next = (U)current->getAnother((uintptr_t)prev);
			prev = current;
			current = next;
			return *this;
		}
		_iterator& operator ++ (int){
			_iterator tmp = *this;
			this->operator++();
			return tmp;
		}
	};
	class iterator : public _iterator<Link*>{
		friend class SpecialDoublyLinkedList;
		typedef  _iterator<Link*> _base;
		iterator(){}
		iterator(Link* prev, Link* current){
			_base::prev = prev;
			_base::current = current;
		}
	public:
		T& operator*(){
			return _base::current->data;
		}
	};
	class const_iterator : public _iterator<const Link*>{
		friend class SpecialDoublyLinkedList;
		typedef  _iterator<Link*> _base;
		const_iterator(){}
		const_iterator(const Link* prev, const Link* current){
			_base::prev = prev;
			_base::current = current;
		}
	public:
		const T& operator*(){
			return _base::current->data;
		}
	};
	SpecialDoublyLinkedList(){
	}
	void push_back(const T& val){
		Link *next = impl.allocate(1);
		impl.construct(next, Link(val, (uintptr_t)impl.tail, (uintptr_t)NULL));
		if(impl.tail != NULL){
			assert(impl.head != NULL);
			Link* prev = (Link*)impl.tail->getAnother((uintptr_t)NULL);
			impl.tail->np = ((uintptr_t)next)^((uintptr_t)prev);
			impl.tail = next;
		}else{
			assert(impl.head == NULL);
			impl.tail = impl.head = next;
		}
		++impl.size;
	}
	size_t size(){
		return impl.size;
	}
	T& back(){
		return impl.tail->data;
	}
	bool empty(){
		assert(impl.tail == NULL ? (impl.head == NULL) : (impl.head != NULL));
		return impl.tail == NULL;
	}
	void pop_back(){
		if(impl.tail != NULL){
			assert(impl.head != NULL);
			Link* next = impl.tail;
			Link* prev = (Link*)impl.tail->getAnother((uintptr_t)NULL);
			impl.tail = prev;
			if(impl.tail != NULL){
				prev = (Link*)impl.tail->getAnother((uintptr_t)next);
				impl.tail->np = ((uintptr_t)NULL)^((uintptr_t)prev);
			}else{
				assert(impl.size == 1);
				impl.head = NULL;
			}
			impl.destroy(next);
			impl.deallocate(next, 0);
			--impl.size;
		}
	}
	iterator begin(){
		return iterator(NULL, impl.head);
	}
	iterator end(){
		return iterator(impl.tail, NULL);
	}
	iterator rbegin(){
			return iterator(NULL, impl.tail);
		}
	iterator rend(){
			return iterator(impl.head, NULL);
	}
	const_iterator begin() const{
		return const_iterator(NULL, impl.head);
	}
	const_iterator end() const{
		return const_iterator(impl.tail, NULL);
	}
	const_iterator rbegin() const{
		return const_iterator(NULL, impl.tail);
	}
	const_iterator rend() const{
		return const_iterator(impl.head, NULL);
	}
	void erase(iterator pos){
		if(pos.current == NULL)
			return;

		if(pos.current == impl.tail){
			impl.tail = pos.prev;
			if(impl.tail == NULL){
				impl.head = NULL;
			}else{
				Link *prev = (Link*)pos.prev->getAnother((uintptr_t)pos.current);
				impl.tail->np = ((uintptr_t)NULL) ^ ((uintptr_t)prev);
			}
		}else if(pos.current == impl.head){
			/*if we reach here, then "head" isn't equal to "tail"*/
			assert(pos.prev == NULL);
			impl.head =  (Link*)pos.current->getAnother((uintptr_t)pos.prev);
			Link* next =(Link*) impl.head->getAnother((uintptr_t)pos.current);
			impl.head->np = ((uintptr_t)NULL) ^ ((uintptr_t)next);
		}else{
			/*we are in the middle of the list*/
			Link *prev = pos.prev;
			Link *next = (Link*)pos.current->getAnother((uintptr_t)pos.prev);
			Link *prevPrev = (Link*)prev->getAnother((uintptr_t)pos.current);
			Link *nextNext = (Link*)next->getAnother((uintptr_t)pos.current);
			prev->np = ((uintptr_t)prevPrev) ^ ((uintptr_t)next);
			next->np = ((uintptr_t)prev) ^ ((uintptr_t)nextNext);
		}
		impl.destroy(pos.current);
		impl.deallocate(pos.current, 0);
		--impl.size;
	}
	void reverse(){
		std::swap(impl.head, impl.tail);
	}
};

#endif /* SPECIALDOUBLYLINKEDLIST_H_ */
