#ifndef MY_LIST_SYMBOL_TABLE_H
#define MY_LIST_SYMBOL_TABLE_H
#include "my_symbol_table.h"
#include <list>
#include <utility>
#include <algorithm>
#include <functional>
#include <cassert>

namespace my_lib{
/**
 * A symbol table implementation that uses list as its
 * 	internal representation
 */
template<class K, class V, class Table = std::list<std::pair<K, V> >, class C=std::equal_to<K> >
class ListSymbolTable : public SymbolTable<K, V>{
	typedef SymbolTable<K, V> Super;
protected:
	typedef typename Table::value_type Pair;
	Table table;	
	Pair* find(const K& k){
		for(auto s = table.begin(), e = table.end();s != e;++ s)
			if(comparator(s->first, k))
				return &*s;
		return nullptr;
	}
	const Pair* find(const K& k) const {
		for(auto s = table.begin(), e = table.end();s != e;++ s)
			if(comparator(s->first, k))
				return &*s;
		return nullptr;
	}
	C comparator;	
	struct ListIteratorImpl : public Super::IteratorImpl{
		typedef typename Table::const_iterator Itor;
		Itor itor;
		
		ListIteratorImpl(Itor itor) : itor(itor){
		}

		void next() override {
			++itor;
		}
		
		bool equals(const typename Super::IteratorImpl& i) const {
			const ListIteratorImpl *p = dynamic_cast<const ListIteratorImpl*>(&i);
			assert(p != nullptr);
			return p->itor == itor;
		}
		
		void assign(const typename Super::IteratorImpl& i){
			const ListIteratorImpl *p = dynamic_cast<const ListIteratorImpl*>(&i);
			assert(p != nullptr);
			itor = p->itor;
		}

		const K& key() const override {
			return itor->first; 
		}

		const V& value() const override {
			return itor->second; 
		}

		ListIteratorImpl* copy() const override {
			return new ListIteratorImpl(itor); 
		}
	};	

	ListIteratorImpl* implBegin() const override {
		return new ListIteratorImpl(table.begin());
	}

	ListIteratorImpl* implEnd() const override {
		return new ListIteratorImpl(table.end());
	}
public:
	typedef C Comparator;
	ListSymbolTable(const Comparator& comparator = Comparator()) 
		: comparator(comparator){
	}
	void put(const K& k, const V& v) override {
		Pair* p = find(k);
		if(p != nullptr){
			p->second = v;
			return;
		}
		table.push_back(Pair{k, v});
	}
	
	const V* get(const K& k) const override {
		auto p = find(k);
		if(p == nullptr)
			return nullptr;
		return &p->second;
	}

	void remove(const K& k) override {
		for(auto s = table.begin(), e = table.end();s != e;++ s)
			if(s->first == k){
				table.erase(s);
				return;
			}
	}

	bool isEmpty() const override {
		return table.empty();
	}

	unsigned int size() const override {
		return table.size();
	}
	
	const Comparator& getComparator() const {
		return comparator;
	}
};
}
#endif //MY_LIST_SYMBOL_TABLE_H
