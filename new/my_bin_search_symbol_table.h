#ifndef MY_BIN_SEARCH_SYMBOL_TABLE_H
#define MY_BIN_SEARCH_SYMBOL_TABLE_H
#include "my_ordered_symbol_table.h"
#include <algorithm>
#include <vector>
#include <utility>

namespace my_lib{
/**
 */
template<class K, class V, class Table = std::vector<std::pair<K, V> >, class C = std::less<K> >
class BinSearchSymbolTable : public OrderedSymbolTable<K, V>{
	typedef OrderedSymbolTable<K, V> Super;
protected:
	typedef typename Table::value_type Pair;
	struct CompareKey{
		C comparator;
		CompareKey(const C& comparator) : comparator(comparator){
		}
		bool operator()(const Pair& p1, const Pair& p2) const {
			return comparator(p1.first, p2.first);
		}
		bool operator()(const K& k1, const K& k2) const {
			return comparator(k1, k2);
		}
		bool operator()(const Pair& p, const K& k) const {
			return comparator(p.first, k);
		}
		bool operator()(const K& k, const Pair& p) const {
			return comparator(k, p.first);
		}
	};
	Table table;
	CompareKey compareKey;
	
	struct BinSearchIteratorImpl: public Super::IteratorImpl{
		typedef typename Table::const_iterator Itor;
		Itor itor;
		
		BinSearchIteratorImpl(Itor itor) : itor(itor){
		}

		void next() override {
			++itor;
		}
	
		bool equals(const typename Super::IteratorImpl& i) const {
			const BinSearchIteratorImpl *p = dynamic_cast<const BinSearchIteratorImpl*>(&i);
			assert(p != nullptr);
			return p->itor == itor;
		}
		
		void assign(const typename Super::IteratorImpl& i){
			const BinSearchIteratorImpl *p = dynamic_cast<const BinSearchIteratorImpl*>(&i);
			assert(p != nullptr);
			itor = p->itor;
		}

		const K& key() const override {
			return itor->first; 
		}

		const V& value() const override {
			return itor->second; 
		}

		BinSearchIteratorImpl* copy() const override {
			return new BinSearchIteratorImpl(itor); 
		}
	};	

	BinSearchIteratorImpl* implBegin() const override {
		return new BinSearchIteratorImpl(table.begin());
	}

	BinSearchIteratorImpl* implEnd() const override {
		return new BinSearchIteratorImpl(table.end());
	}
public:
	BinSearchSymbolTable(const C& comparator = C()) : compareKey(comparator){
	}

	void put(const K& k, const V& v) override {
		auto pos = std::lower_bound(table.begin(), table.end(), k, compareKey);
		if(pos != table.end() && !compareKey.comparator(k, pos->first) &&
			!compareKey.comparator(k, pos->first)){
			pos->second = v;
			return;
		}
		table.insert(pos, Pair(k, v));
	}

	const V* get(const K& k) const override{
		auto pos = std::lower_bound(table.begin(), table.end(), k, compareKey);
		if(pos != table.end() && !compareKey.comparator(k, pos->first) &&
			!compareKey.comparator(k, pos->first)){
			return &pos->second;
		}
		return nullptr;
	}
	
	void remove(const K& k) override {
		auto pos = std::lower_bound(table.begin(), table.end(), k, compareKey);
		if(pos != table.end() && !compareKey.comparator(k, pos->first) &&
			!compareKey.comparator(k, pos->first)){
			table.erase(pos);
		}
	}	

	bool isEmpty() const override {
		return table.empty();
	}

	unsigned int size() const override {
		return table.size();
	}
	
	const K* min() const override {
		size_t s = table.size();
		if(s <= 0)
			return nullptr;
		return &table[0].first;
	}

	const K* max() const override {
		size_t s = table.size();
		if(s <= 0)
			return nullptr;
		return &table[s - 1].first;
	}	

	const K* floor(const K& k) const override {
		auto pos = std::upper_bound(table.begin(), table.end(), k, compareKey);
		if(pos != table.begin()){
			-- pos;
			return &pos->first;
		}
		return nullptr;
	}

	const K* ceil(const K& k) const override {
		auto pos = std::lower_bound(table.begin(), table.end(), k, compareKey);
		if(pos != table.end())
			return &pos->first;
		return nullptr;				
	}

	unsigned int rank(const K& k) const override {
		auto pos = std::lower_bound(table.begin(), table.end(), k, compareKey);
		return pos - table.begin();
	}

	const K* select(unsigned int index) const override {
		if(index >= table.size())
			return nullptr;
		return &table[index].first;
	}

	void removeMin() override {
		table.erase(table.begin());
	}

	void removeMax() override {
		table.erase(-- table.end());
	}

	unsigned int size(const K& l, const K& h) const override {
		auto p1 = std::lower_bound(table.begin(), table.end(), l, compareKey);
		auto p2 = std::upper_bound(table.begin(), table.end(), h, compareKey);
		if(p1 < p2)
			return p2 - p1;
		return 0;
	}

	bool isValid() const override {
		auto s = table.begin(), e = table.end();
		auto s0 = s;
		++ s;
		while(s != e){
			if(compareKey(*s, *s0))
return false;
		++ s;
		++ s0;
}
return true;
}

void clear() override {
table.clear();
}

const C& getComparator() const {
	return compareKey.comparator;
}

};
}
#endif //MY_BIN_SEARCH_SYMBOL_TABLE_H
