#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <iterator>
#include <set>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <ctime>
#include <limits>

#define MY_MEM_TEST
#ifdef MY_MEM_TEST 
namespace{

template<typename T>
class MyAllocator {
public : 
	//    typedefs
	typedef T value_type;
	typedef value_type* pointer;
	typedef const value_type* const_pointer;
	typedef value_type& reference;
	typedef const value_type& const_reference;
	typedef std::size_t size_type;
	typedef std::ptrdiff_t difference_type;

public : 
	//    convert an allocator<T> to allocator<U>
	template<typename U>
	struct rebind {
		typedef MyAllocator<U> other;
	};

public : 
	inline explicit MyAllocator() {}
	inline ~MyAllocator() {}
	inline explicit MyAllocator(MyAllocator const&) {}
	template<typename U>
	inline explicit MyAllocator(MyAllocator<U> const&) {}

	//    address
	inline pointer address(reference r) { return &r; }
	inline const_pointer address(const_reference r) { return &r; }

	//    memory allocation
	inline pointer allocate(size_type cnt, 
			typename std::allocator<void>::const_pointer = 0) { 
		return reinterpret_cast<pointer>(malloc(cnt * sizeof (T))); 
	}
	inline void deallocate(pointer p, size_type) { 
		free(p); 
	}

	//    size
	inline size_type max_size() const { 
		return std::numeric_limits<size_type>::max() / sizeof(T);
	}

	//    construction/destruction
	inline void construct(pointer p, const T& t) {
	       	new(p) T(t);
       	}
	inline void destroy(pointer p) { p->~T(); }

	inline bool operator==(MyAllocator const&) { return true; }
	inline bool operator!=(MyAllocator const& a) { return !operator==(a); }
};    //    end of class Allocator 

struct MemAllocator{
	std::set<void*, std::less<void*>, MyAllocator<void*> > allocated;
	~MemAllocator(){
		assert(allocated.empty());
	}
	static MemAllocator& get(){
		static MemAllocator instance;
		return instance;
	}
	static void* allocate(unsigned int size){
		auto p = malloc(size);
		//printf("allocate memory : %p\n", p);
		get().allocated.insert(p);
		return p;
	}
	
	static void deallocate(void *p){
		MemAllocator& m = get();
		assert(m.allocated.find(p) != m.allocated.end());
		m.allocated.erase(p);
		//printf("deallocate memory : %p\n", p);
		free(p);
	}
};
}
void* operator new(unsigned int size){
	return MemAllocator::allocate(size);
}
void operator delete(void *p){
	return MemAllocator::deallocate(p);
}
#endif

#include "my_list_symbol_table.h"
#include "my_bin_search_symbol_table.h"
#include "my_rand_string_generator.h"
#include "my_binary_search_tree.h"

using namespace std;
using namespace my_lib;
namespace{
#define ARSIZE(ar) (sizeof(ar) / sizeof(*(ar)))
typedef double KeyType;
typedef double ValueType;
void fillSymbolTable(SymbolTable<KeyType, ValueType>& st){
	unsigned int values[]={4, 3, 3, 2, 9, 10, 3, 6, 8, 7, };
	auto count = ARSIZE(values);
	for(decltype(count) i = 0;i < count;i += 2){
		st.put(values[i], values[i + 1]);
	}
}
int validateSymbolTable(SymbolTable<KeyType, ValueType>& st){
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
	assert(!st.contains(3));
	assert(st.size() == 3);
	st.remove(4);
	st.remove(8);
	st.remove(9);
	assert(st.isEmpty());
	assert(st.size() == 0);		
	return 0;
}

int validateOrderedSymbolTable(OrderedSymbolTable<KeyType, ValueType>& st){
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

void moreOrderedTest(OrderedSymbolTable<KeyType, ValueType>& ost){
	const unsigned int count = 100;
	double nums[count];
	unsigned int indices[count];
	for(unsigned int i = 0;i < count;++ i){
		nums[i] = i;
		indices[i] = i;
	}

	std::random_shuffle(indices, indices + count);
	ost.clear();
	assert(ost.isEmpty());
	assert(ost.size() == 0);
	for(unsigned int i = 0;i < count;++ i)
		ost.put(nums[indices[i]], nums[indices[i]]);

	for(unsigned int i = 0;i < count;++ i){
		assert(ost.contains(nums[indices[i]]));
		assert(nums[indices[i]] == *ost.get(nums[indices[i]]));
	}

	unsigned int index = 0;
	for(auto i : ost){
		assert(nums[index] == i.key());
		assert(nums[index] == i.value());
		++ index;
	}

		
	cout << "test : floor, ceil, select, rank and size(Key, Key)" << endl;
	for(unsigned int i = 0;i < count;++ i){
		assert(ost.rank(*ost.select(nums[i])) == i);
		assert(*ost.select(i) == nums[i]);
		assert(*ost.floor(nums[i]) == nums[i]);
		assert(*ost.ceil(nums[i]) == nums[i]);
		assert(*ost.floor(nums[i] + 0.5) == nums[i]);
		assert(*ost.ceil(nums[i] - 0.5) == nums[i]);
		if(i > 0)
			assert(i > 0 && *ost.floor(nums[i] - 0.5) == nums[i - 1]);
		else{
			assert(ost.floor(nums[i] - 0.5) == nullptr);
		}
		if(i < count - 1)
			assert(i < count - 1 && *ost.ceil(nums[i] + 0.5) == nums[i + 1]);
		else
			assert(ost.ceil(nums[i] + 0.5) == nullptr);
		assert(ost.size(nums[i], nums[i]) == 1);
		assert(ost.size(nums[i] + 0.5, nums[i] - 0.5) == 0);
		assert(ost.size(nums[i] + 0.5, nums[i] + 0.5) == 0);
		if(i + 5 < count){
			assert(ost.size(nums[i], nums[i + 5]) == 6);
			assert(ost.size(nums[i], nums[i + 5] + 0.5) == 6);
			assert(ost.size(nums[i] - 0.5, nums[i + 5] + 0.5) == 6);
			assert(ost.size(nums[i] - 0.5, nums[i + 5]) == 6);
		}
	}

	cout << "test : min, removeMin" << endl;
	const int toRemove = 5;
	for(unsigned int i = 0;i < toRemove;++ i){
		assert(nums[i] == *ost.min());
		ost.removeMin();
	}
	

	cout << "test : contains" << endl;
	for(unsigned int i = 0;i < toRemove;++ i)
		assert(!ost.contains(nums[i]));

	cout << "test : max, removeMax" << endl;
	for(unsigned int i = count - 1, j = 0;j < toRemove;--i, ++ j){
		assert(nums[i] == *ost.max());
		ost.removeMax();
		assert(!ost.contains(nums[i]));
	}

	assert(ost.size() == (count - 2 * toRemove));

	cout << "test : iterator" << endl;
	index = toRemove;
	for(auto i : ost){
		assert(nums[index] == i.key());
		assert(nums[index] == i.value());
		++ index;
	}

	BinarySearchTree<KeyType, ValueType>* bst = 
		dynamic_cast<BinarySearchTree<KeyType, ValueType>*>(&ost);
	if(bst != nullptr)
		assert(bst->isValid());

	cout << "test : remove" << endl;
	index = ost.size();
	for(unsigned int i = toRemove + 4;i < count - toRemove;i += 2){
		ost.remove(nums[i]);
		if(bst != nullptr)
			assert(bst->isValid());
		--index;
	}
	
	cout << "test : get, put" << endl;
	assert(index == ost.size());
	
	for(auto i : ost){
		assert(i.value() == i.key());
		ost.put(i.key(), 2 * i.key());
		assert(*ost.get(i.key()) == i.key() * 2);
	}
}
int validateIterator(const SymbolTable<KeyType, ValueType>& st){
	cout << "validate iterator" << endl;
	for(auto i : st){
		cout << i.key() << " : " << i.value() << endl;
	}
	return 0;
}

clock_t testOne(SymbolTable<std::string, bool>& st, const std::vector<string>& strs, unsigned int n){
	clock_t start = clock();
	for(unsigned int i = 0;i < n;++ i)
		st.put(strs[i], true);
	return clock() - start;
}

template<class T>
void doublingTestOfSymbolTable(const std::vector<std::string>& vec){
	clock_t last = 0;
	for(unsigned int i = 1;i < vec.size(); i *= 2){
		T st;
		clock_t current = testOne(st, vec, i);
		std::cout << "count = " << i << ", ticks = " << current;
		if(i != 1){
			if(last == 0){
				cout << ", ratio = -";
			}else{
				cout << ", ratio = " << current / (double)last;
			}
		}
		cout << endl;
		last = current;
	}
}
void doublingTestOfSymbolTable(){
	std::vector<std::string> vec;
	if(false){
		std::ifstream in("/development/documents/books/algos4/algs4.cs.princeton.edu/31elementary/tale.txt");
		std::copy(std::istream_iterator<std::string>(in), 
			std::istream_iterator<std::string>(),
			std::back_inserter(vec));
	}else{
		auto count = (1 << 17);
		RandStringGenerator gen(5, 10);
		cout << "generating " << count << " strings with length between 5 and 10 " << endl;
		for(auto i = 0;i < count;++ i){
			vec.push_back(string());
			gen.randStringPrintable(vec.back());
			//cout << vec.back() << endl;
		}
	}
	cout << "-------------doubling test for binary search symbol table-------------" << endl;
	doublingTestOfSymbolTable<BinSearchSymbolTable<std::string, bool> >(vec);
	cout << "-------------doubling test for list symbol table-------------" << endl;
	doublingTestOfSymbolTable<ListSymbolTable<std::string, bool> >(vec);

}
int testSymbolTables(int argc, char *argv[]){
	ListSymbolTable<KeyType, ValueType> listST;
	cout << "-----------symbol table with sequential search implementation---------------" << endl;
	fillSymbolTable(listST);
	validateIterator(listST);
	validateSymbolTable(listST);
	cout << "------------symbol table with binary search implementation----------------" << endl;
	BinSearchSymbolTable<KeyType, ValueType> binST;
	fillSymbolTable(binST);
	validateIterator(binST);
	validateOrderedSymbolTable(binST);
	validateSymbolTable(binST);
	moreOrderedTest(binST);
	cout << "------------binary search tree---------------" << endl;
	BinarySearchTree<KeyType, ValueType> bst;
	assert(bst.isValid());
	fillSymbolTable(bst);
	assert(bst.isValid());
	validateIterator(bst);
	assert(bst.isValid());
	validateOrderedSymbolTable(bst);
	assert(bst.isValid());
	validateSymbolTable(bst);
	assert(bst.isValid());
	moreOrderedTest(bst);
	//doublingTestOfSymbolTable();
	return 0;
}

}
int testSymbolTable(int argc, char *argv[]){
	return testSymbolTables(argc, argv);
}

