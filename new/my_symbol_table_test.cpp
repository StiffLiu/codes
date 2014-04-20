#include <set>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <cstdio>
namespace{
template<class T>
struct MyAllocator : public std::allocator<T>{
private:
	typedef std::allocator<T> Super;
public:
	template<class T1>
	struct rebind{
		typedef MyAllocator<T1> other;
	};
	T* allocate(typename Super::size_type n, const void* =0){
		return static_cast<T*>(malloc(n * sizeof (T)));
	}
	void deallocate(T* p, typename Super::size_type){
		free(p);
	}
};
struct MemAllocator{
	std::set<void*, std::less<void*>, MyAllocator<void*> > allocated;
	~MemAllocator(){
		assert(allocated.empty());
	}
}memAllocator;
}
void* operator new(unsigned int size){
	auto p = malloc(size);
	std::cout << "allocate memory : " << p << std::endl;
	memAllocator.allocated.insert(p);
	return p;
}
void operator delete(void *p){
	assert(memAllocator.allocated.find(p) != memAllocator.allocated.end());
	memAllocator.allocated.erase(p);
	std::cout << "deallocate memory : " << p << std::endl;
	free(p);
}

#include "my_list_symbol_table.h"
#include "my_bin_search_symbol_table.h"

using namespace std;
using namespace my_lib;
namespace{
#define ARSIZE(ar) (sizeof(ar) / sizeof(*(ar)))
void fillSymbolTable(SymbolTable<unsigned int, unsigned int>& st){
	unsigned int values[]={4, 3, 3, 2, 9, 10, 3, 6, 8, 7, };
	auto count = ARSIZE(values);
	for(decltype(count) i = 0;i < count;i += 2){
		st.put(values[i], values[i + 1]);
	}
}
int validateSymbolTable(SymbolTable<unsigned int, unsigned int>& st){
	std::cout << "test symbol table" << std::endl;
	assert(*st.get(4) == 3);
	assert(*st.get(3) == 6);
	assert(*st.get(8) == 7);
	assert(*st.get(9) == 10);
	assert(st.size() == 4);
	assert(!st.isEmpty());
	assert(st.contains(8));
	assert(!st.contains(5));
	st.remove(3);
	assert(st.size() == 3);
	st.remove(4);
	st.remove(8);
	st.remove(9);
	assert(st.isEmpty());
	assert(st.size() == 0);		
	return 0;
}

int validateOrderedSymbolTable(OrderedSymbolTable<unsigned int, unsigned int>& st){
	assert(*st.min() == 3);
	assert(*st.max() == 9);
	assert(*st.floor(3) == 3);
	assert(*st.ceil(4) == 4);
	assert(*st.floor(5) == 4);
	assert(*st.ceil(5) == 8);
	assert(st.rank(*st.min()) == 0);
	assert(st.rank(*st.max()) == st.size() - 1);
	assert(st.rank(*st.max() + 1) == st.size());
	assert(*st.select(0) == *st.min());
	assert(*st.select(st.size() - 1) == *st.max());
	assert(st.select(st.size()) == nullptr);
	assert(st.size(*st.min(), *st.max()) == st.size());
	assert(st.size(4, 8) == 2);
	assert(st.size(5, 10) == 2);
	assert(st.size(2, 4) == 2);
	assert(st.size(1, 2) == 0);

	unsigned int size = st.size();
	for(unsigned int  i = 0;i < size;++ i){
		assert(st.rank(*st.select(i)) == i);
	}
	st.removeMin();
	assert(*st.min() == 4);
	st.removeMax();
	assert(*st.max() == 8);
	assert(st.size() == 2);
	st.put(3, 6);
	st.put(9, 10);
	return 0;
}

int validateIterator(const SymbolTable<unsigned int, unsigned int>& st){
	cout << "validate iterator" << endl;
	for(auto i : st){
		cout << i.key() << " : " << i.value() << endl;
	}
	return 0;
}

void doublingTestOfSymbolTable(){
	
}

int testListSymbolTable(int argc, char *argv[]){
	ListSymbolTable<unsigned int, unsigned int> listST;
	fillSymbolTable(listST);
	validateIterator(listST);
	validateSymbolTable(listST);
	cout << "---------------------------" << endl;
	BinSearchSymbolTable<unsigned int, unsigned int> binST;
	fillSymbolTable(binST);
	validateIterator(binST);
	validateOrderedSymbolTable(binST);
	validateSymbolTable(binST);
	doublingTestOfSymbolTable();
	return 0;
}

}
int testSymbolTable(int argc, char *argv[]){
	return testListSymbolTable(argc, argv);
}

