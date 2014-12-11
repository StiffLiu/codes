#include "my_symbol_table.h"
#include "my_list_symbol_table.h"
#include <vector>
#include <cassert>
#include <iostream>

namespace my_lib{
template<class K>
struct TypeConverter{
	template<class U>
	K operator()(const U& u) const{
		return K(u);
	}
};
template<class K, class V, class C=std::equal_to<K>, class Hash=TypeConverter<unsigned int>,
       class Chain = ListSymbolTable<K, V, std::vector<std::pair<K, V> >, C > >
class SeparateChainingHashTable : public SymbolTable<K, V>{
	typedef SymbolTable<K, V> Super;
protected:
	typedef std::vector<Chain> Table;
	Table table;
	float loadFactor = 2;
	Hash hash;
	C comparator;
	unsigned int n = 0;
	unsigned int resize(unsigned int numChain){
		Table oldTable;
		oldTable.swap(table);

		table.resize(numChain, Chain(comparator));
		n = 0;
		for(auto& chain : oldTable)
			for(auto kv : chain)
				put(kv.key(), kv.value());
		return numChain;
	}
	struct HashTableIteratorImpl : public Super::IteratorImpl{
		typedef typename Table::const_iterator Itor;
		typedef typename Chain::Iterator ChainItor;
		Itor itor;
		ChainItor chainItor;
		unsigned int remainingInChain = 0;
		const Table *table = nullptr;
		
		HashTableIteratorImpl(Itor itor, const Table *table) : itor(itor), table(table){
			moveToElement();
		}

		HashTableIteratorImpl(){
		}

		void next() override {
			//avoid memory allocation in "itor->end()"
			//if(chainItor != itor->end()){
			if(remainingInChain > 0){
				++ chainItor;
				--remainingInChain;
				if(remainingInChain > 0)
					return;
			}
			++itor;
			moveToElement();
		}
		
		bool equals(const typename Super::IteratorImpl& i) const {
			const HashTableIteratorImpl *p = dynamic_cast<const HashTableIteratorImpl*>(&i);
			assert(p != nullptr);
			assert(table != nullptr);
			assert(p->table == table);
			if(table == nullptr)
				return table == p->table;
			if(p->table != table)
				return false;
			if(p->itor != itor)
				return false;
			if(itor == table->end())
				return true;
			return chainItor == p->chainItor;
		}
		
		void assign(const typename Super::IteratorImpl& i){
			const HashTableIteratorImpl *p = dynamic_cast<const HashTableIteratorImpl*>(&i);
			assert(p != nullptr);
			assert(table != nullptr);
			assert(p->table == table);
			itor = p->itor;
			chainItor = p->chainItor;
			table = p->table;
		}

		const K& key() const override {
			return chainItor.key(); 
		}

		const V& value() const override {
			return chainItor.value();
		}

		HashTableIteratorImpl* copy() const override {
			return new HashTableIteratorImpl(itor, chainItor, table); 
		}
	private:
		HashTableIteratorImpl(Itor itor, ChainItor chainItor, const Table *table) 
			: itor(itor), chainItor(chainItor), table(table){
		}
		inline void moveToElement(){
			assert(table != nullptr);
			if(table == nullptr)
				return;
			Itor end = table->end();
			while(itor != end && itor->isEmpty())
				++ itor;
			if(itor != end){
				remainingInChain = itor->size();
				chainItor = itor->begin();
			}
		}

	};	
	HashTableIteratorImpl* implBegin() const override{
		return new HashTableIteratorImpl(table.begin(), &table);
	}
	HashTableIteratorImpl* implEnd() const override{
		return new HashTableIteratorImpl(table.end(), &table);
	}
	inline unsigned int assurePut(){
		auto numChain = table.size();
		if(numChain == 0){
			table.resize(1, Chain(comparator));
			numChain = 1;
		}
		if(n / (float)numChain >= loadFactor){
			numChain = 2 * numChain + 1;
			numChain = resize(numChain);		
		}
		return numChain;
	}
	inline void assureRemove(){
		auto numChain = table.size();
		if(n == 0){
			table.clear();
			return;
		}
		if(n / (float)numChain <= loadFactor / 2){
			resize(numChain / 2);
		}
	}
public:
	SeparateChainingHashTable(float loadFactor = 2, const Hash& hash = Hash(), 
		const C& comparator = C()) : loadFactor(loadFactor), hash(hash), comparator(comparator){
	}	
	void put(const K& k, const V& v) override{
		auto index = hash(k);
		auto numChain = assurePut();
		auto& chain = table[index % numChain];
		auto chainOldSize = chain.size();
		chain.put(k, v);
		n += (chain.size() - chainOldSize);
	}

	const V* get(const K& k) const override{
		if(n == 0)
			return nullptr;
		auto& chain = table[hash(k) % table.size()];
		return chain.get(k);
	}

	bool remove(const K& k) override{
		auto numChain = table.size();
		bool isRemoved = false;
		if(numChain <= 0)
			return isRemoved;
		auto& chain = table[hash(k) % numChain];
		auto chainOldSize = chain.size();
		isRemoved = chain.remove(k);
		n -= (chainOldSize - chain.size());
		assureRemove();
		return isRemoved;
	}

	unsigned int size() const override{
		return n;
	}

	void clear() override{
		table.clear();
		n = 0;
	}

	unsigned int getChainCount() const{
		return table.size();
	}

	unsigned int getChainSize(unsigned int i) const{
		return table[i].size();
	}

	float getLoadFactor() const{
		return loadFactor;
	}

	const C& getComparator() const{
		return comparator;
	}
};

template<class K, class V, class Hash1, class C=std::equal_to<K>, class Hash2=TypeConverter<unsigned int>,
       class Chain = ListSymbolTable<K, V, std::vector<std::pair<K, V> >, C > >
class DoubleHashSeparateChaining: public SeparateChainingHashTable<K, V, C, Hash2, Chain>{
	typedef SeparateChainingHashTable<K, V, C, Hash2, Chain> Super;
protected:
	Hash1 hash;
public:
	DoubleHashSeparateChaining(float loadFactor = 2, const Hash1& hash1 = Hash1(), const Hash2& hash2 = Hash2(), 
		const C& comparator = C()) : Super(loadFactor, hash2, comparator), hash(hash1){
	}
	void put(const K& k, const V& v) override{
		auto numChain = Super::assurePut();
		auto index1 = hash(k) % numChain;
		auto chain = &Super::table[index1];
		if(chain->get(k) != nullptr){
			chain->put(k, v);
			return;
		}
		auto index2 = Super::hash(k) % numChain;
		if(index1 != index2){
			auto chain2 = &Super::table[index2];
			if(chain2->get(k) != nullptr){
				chain2->put(k, v);
				return;
			}
			if(chain2->size() < chain->size())
				chain = chain2;
		}
		auto chainOldSize = chain->size();
		chain->put(k, v);
		Super::n += (chain->size() - chainOldSize);
	}

	const V* get(const K& k) const override{
		if(Super::n == 0)
			return nullptr;
		auto index1 = hash(k) % Super::table.size();
		auto value = Super::table[index1].get(k);
		if(value != nullptr)
			return value;
		auto index2 = Super::hash(k) % Super::table.size();
		if(index1 == index2)
			return nullptr;
		return Super::table[index2].get(k);
	}

	bool remove(const K& k) override{
		auto numChain = Super::table.size();
		if(numChain <= 0)
			return false;
		auto index1 = hash(k) % numChain;
		if(Super::table[index1].remove(k)){
			-- Super::n;
			Super::assureRemove();
			return true;
		}
		auto index2 = Super::hash(k) % numChain; 
		if(index1 != index2 && Super::table[index2].remove(k)){
			-- Super::n;
			Super::assureRemove();
			return true;
		}
		return false;
	}
	
};

template<class K, class V, class C=std::equal_to<K>, class Hash=TypeConverter<unsigned int> >
class LinearProbingHashTable : public SymbolTable<K, V>{
	typedef SymbolTable<K, V> Super;
protected:
	typedef std::vector<std::pair<K, V> > Table;
	typedef std::vector<bool> Usage;
	Table table;
	Usage usage;
	float loadFactor = 0.5;
	Hash hash;
	C comparator;
	unsigned int n = 0;
	unsigned int resize(unsigned int numSlots){
		Table oldTable;
		Usage oldUsage;
		
		oldTable.swap(table);
		oldUsage.swap(usage);
		table.resize(numSlots);
		usage.resize(numSlots);
		n = 0;

		for(auto i = 0;i < oldUsage.size();++ i)
			if(oldUsage[i]){
				auto& kv = oldTable[i];
				put(kv.first, kv.second);
			}
		return numSlots;
	}
	void assertValid(){
		assert(loadFactor > 0);
	}
	struct HashTableIteratorImpl;
	friend struct HashTableIteratorImpl;
	struct HashTableIteratorImpl : public Super::IteratorImpl{
		const LinearProbingHashTable* table;
		unsigned int n;
		
		HashTableIteratorImpl(const LinearProbingHashTable* table, unsigned int n) 
			: table(table), n(n){
		}

		void next() override {
			++ n;
			auto count = table->usage.size();
			while(n < count && !table->usage[n])
				++n;
		}
		
		bool equals(const typename Super::IteratorImpl& i) const {
			const HashTableIteratorImpl *p = dynamic_cast<const HashTableIteratorImpl*>(&i);
			assert(p != nullptr);
			assert(table == p->table);
			if(table != p->table)
				return false;
			if(n >= table->usage.size())
				return p->n >= table->usage.size();
			return p->n == n;
		}
		
		void assign(const typename Super::IteratorImpl& i){
			const HashTableIteratorImpl *p = dynamic_cast<const HashTableIteratorImpl*>(&i);
			assert(p != nullptr);
			assert(table == p->table);
			table = p->table;
			n = p->n;
		}

		const K& key() const override {
			return table->table[n].first; 
		}

		const V& value() const override {
			return table->table[n].second; 
		}

		HashTableIteratorImpl* copy() const override {
			return new HashTableIteratorImpl(table, n); 
		}
	};	

	HashTableIteratorImpl* implBegin() const override {
		return new HashTableIteratorImpl(this, 0);
	}

	HashTableIteratorImpl* implEnd() const override {
		return new HashTableIteratorImpl(this, usage.size());
	}
public:
	LinearProbingHashTable(float loadFactor = 0.5, const Hash& hash = Hash(), const C& comparator = C()) 
	       : loadFactor(loadFactor > 1.0 ? 1.0 : loadFactor), hash(hash), comparator(comparator){
	}	
	void put(const K& k, const V& v) override{
		auto numSlots = usage.size();
		if(n / loadFactor + 1 > numSlots){
			numSlots = (numSlots == 0 ? 2 : 2 * numSlots + 1);
			numSlots = resize(numSlots);
		}
		unsigned int h = hash(k) % numSlots;
		for(;usage[h];h = (h + 1) % numSlots)
			if(comparator(table[h].first, k)){
				//Do we need to use copy constructor,
				//instead of assignment opeator?
				table[h].second = v;
				return;
			}
		//Do we need to use copy constructor,
		//instead of assignment opeator?
		table[h] = {k, v};
		usage[h] = true;
		++ n;
	}

	const V* get(const K& k) const override{
		if(n > 0){
			auto numSlots = usage.size();
			for(auto h = hash(k) % numSlots;usage[h];h = (h + 1) % numSlots)
				if(comparator(table[h].first, k))
					return &table[h].second;
		}
		return nullptr;
	}

	bool remove(const K& k) override{
		auto numSlots = usage.size();
		for(auto h = hash(k) % numSlots;usage[h];h = (h + 1) % numSlots)
			if(comparator(table[h].first, k)){
				if(n == 1){
					clear();
					return true;
				}
				usage[h] = false;
				if(2 * n / loadFactor + 1 < usage.size()){
					resize((usage.size() - 1) / 2);
				}else{
					--n;
					h = (h + 1) % numSlots;
					while(usage[h]){
						usage[h] = false;
						-- n;
						put(table[h].first, table[h].second);
						h = (h + 1) % numSlots;
					}
				}
				return true;
			}
		return false;
	}

	unsigned int size() const override{
		return n;
	}

	unsigned int getSlotSize() const{
		return table.size();
	}

	const Usage& getSlotUsage() const{
		return usage;
	}

	const K* getKeyAtSlot(unsigned int slotIndex) const{
		if(slotIndex >= usage.size() || !usage[slotIndex])
			return nullptr;
		return &table[slotIndex].first;
	}

	float getLoadFactor() const{
		return loadFactor;
	}

	const C& getComparator() const{
		return comparator;
	}

	void clear() override{
		table.clear();
		usage.clear();
		n = 0;
	}
};

template<class K, class V, class C=std::equal_to<K>, class Hash=TypeConverter<unsigned int> >
class LazyDeleteLinearProbingHashTable : public SymbolTable<K, V>{
	typedef SymbolTable<K, V> Super;
protected:
	typedef std::vector<std::pair<K, V> > Table;
	typedef std::vector<bool> Usage;
	typedef std::vector<bool> HasElement;
	Table table;
	Usage usage;
	HasElement hasElement;
	float loadFactor = 0.5;
	Hash hash;
	C comparator;
	unsigned int n = 0;
	unsigned int resize(unsigned int numSlots){
		Table oldTable;
		HasElement oldHasElement;
		
		oldTable.swap(table);
		oldHasElement.swap(hasElement);
		table.resize(numSlots);

		hasElement.resize(numSlots);
		usage.clear();
		usage.resize(numSlots);
		n = 0;

		for(auto i = 0;i < oldHasElement.size();++ i)
			if(oldHasElement[i]){
				auto& kv = oldTable[i];
				put(kv.first, kv.second);
			}
		return numSlots;
	}
	void assertValid(){
		assert(loadFactor > 0);
	}
	struct HashTableIteratorImpl;
	friend struct HashTableIteratorImpl;
	struct HashTableIteratorImpl : public Super::IteratorImpl{
		const LazyDeleteLinearProbingHashTable* table;
		unsigned int n;
		
		HashTableIteratorImpl(const LazyDeleteLinearProbingHashTable* table, unsigned int n) 
			: table(table), n(n){
		}

		void next() override {
			++ n;
			auto count = table->hasElement.size();
			while(n < count && !table->hasElement[n])
				++n;
		}
		
		bool equals(const typename Super::IteratorImpl& i) const {
			const HashTableIteratorImpl *p = dynamic_cast<const HashTableIteratorImpl*>(&i);
			assert(p != nullptr);
			assert(table == p->table);
			if(table != p->table)
				return false;
			if(n >= table->hasElement.size())
				return p->n >= table->hasElement.size();
			return p->n == n;
		}
		
		void assign(const typename Super::IteratorImpl& i){
			const HashTableIteratorImpl *p = dynamic_cast<const HashTableIteratorImpl*>(&i);
			assert(p != nullptr);
			assert(table == p->table);
			table = p->table;
			n = p->n;
		}

		const K& key() const override {
			return table->table[n].first; 
		}

		const V& value() const override {
			return table->table[n].second; 
		}

		HashTableIteratorImpl* copy() const override {
			return new HashTableIteratorImpl(table, n); 
		}
	};	

	HashTableIteratorImpl* implBegin() const override {
		return new HashTableIteratorImpl(this, 0);
	}

	HashTableIteratorImpl* implEnd() const override {
		return new HashTableIteratorImpl(this, hasElement.size());
	}
public:
	LazyDeleteLinearProbingHashTable(float loadFactor = 0.5, const Hash& hash = Hash(), const C& comparator = C()) 
	       : loadFactor(loadFactor > 1.0 ? 0.9 : loadFactor), hash(hash), comparator(comparator){
	}	
	void put(const K& k, const V& v) override{
		auto numSlots = usage.size();
		if(n / loadFactor + 1 > numSlots){
			numSlots = (numSlots == 0 ? 2 : 2 * numSlots + 1);
			numSlots = resize(numSlots);
		}
		unsigned int h = hash(k) % numSlots;
		for(;usage[h];h = (h + 1) % numSlots)
			if(!hasElement[h]){
				hasElement[h] = true;
				table[h] = {k, v};
				++ n;
				return;

			}else if(comparator(table[h].first, k)){
				//Do we need to use copy constructor,
				//instead of assignment opeator?
				table[h].second = v;
				return;
			}
		//Do we need to use copy constructor,
		//instead of assignment opeator?
		table[h] = {k, v};
		usage[h] = true;
		hasElement[h] = true;
		++ n;
	}

	const V* get(const K& k) const override{
		if(n > 0){
			auto numSlots = usage.size();
			for(auto h = hash(k) % numSlots;usage[h];h = (h + 1) % numSlots)
				if(hasElement[h] && comparator(table[h].first, k))
					return &table[h].second;
		}
		return nullptr;
	}

	bool remove(const K& k) override{
		auto numSlots = usage.size();
		for(auto h = hash(k) % numSlots;usage[h];h = (h + 1) % numSlots)
			if(hasElement[h] && comparator(table[h].first, k)){
				if(n == 1){
					clear();
					return true;
				}
				-- n;
				hasElement[h] = false;
				if(loadFactor * n + 1 < usage.size()){
					resize((usage.size() - 1) / 2);
				}
				return true;
			}
		return false;
	}

	unsigned int size() const override{
		return n;
	}

	unsigned int getSlotSize() const{
		return table.size();
	}

	const Usage& getSlotUsage() const{
		return usage;
	}

	const K* getKeyAtSlot(unsigned int slotIndex) const{
		if(slotIndex >= usage.size() || !usage[slotIndex])
			return nullptr;
		return &table[slotIndex].first;
	}

	void clear() override{
		table.clear();
		usage.clear();
		hasElement.clear();
		n = 0;
	}
};
}
