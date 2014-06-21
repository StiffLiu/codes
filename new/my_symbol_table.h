#ifndef MY_SYMBOL_TABLE_H
#define MY_SYMBOL_TABLE_H 
namespace my_lib{
/**
 * A symbol table is a data structure that must support the following interfaces:
 * 1. put a key-value pair into the symbol table
 * 2. get the value associated with a key.
 * 3. delete a key-value pair from the symbol table.
 * 3. test whether a key is contained in the symbol table.
 *
 * A key could only be associated with one value. If the symbol table already contains 
 * a key, and a new key-value pair in put in the symbol table with that key, the old 
 * value will replace the new value. 
 * 
 *To make the symbol table more convinient to use, the following interfaces is  
 *also supported:
 * 1. get the number of key value pairs in the symbol table.
 * 2. test whether the symbol table is empty.
 * 3. an iterator used to iterate through the symbol table.
 */
template<class K, class V>
class SymbolTable{
protected:
	struct IteratorImpl{
		virtual void next() = 0;
		virtual bool equals(const IteratorImpl& i) const = 0;
		virtual void assign(const IteratorImpl& i) = 0;
		virtual const K& key() const = 0;
		virtual const V& value() const = 0;
		virtual IteratorImpl* copy() const = 0;
		virtual ~IteratorImpl(){
		}
	};
	virtual IteratorImpl* implBegin() const = 0;
	virtual IteratorImpl* implEnd() const = 0;
public:
	typedef K KeyType;
	typedef V ValueType;
	class Iterator;
	struct kv{
		const K& key(){
			return impl->key();
		}
		const V& value(){
			return impl->value();
		}
	private:
		friend class Iterator;
		kv(IteratorImpl *impl) : impl(impl){
		}
		IteratorImpl *impl;
	};
	class Iterator{
	public:
		Iterator(const Iterator& i) : impl(i.impl->copy()){
		}
		
		Iterator(Iterator&& i) : impl(i.impl){
			i.impl = nullptr;
		}

		Iterator& operator++(){
			impl->next();
			return *this;
		}

		Iterator operator++(int){
			Iterator tmp(impl->copy());
			impl->next();
			return tmp;
		}	

		bool operator==(const Iterator& i){
			return impl->equals(*i.impl);
		}
		
		bool operator!=(const Iterator& i){
			return !impl->equals(*i.impl);
		}

		kv operator*(){
			return kv(impl);
		}

		Iterator& operator=(const Iterator& i){
			if(&i != this){
				impl->assign(*i.impl);
			}
			return *this;
		}

		const K& key(){
			return impl->key();
		}

		const V& value(){
			return impl->value();
		}

		~Iterator(){
			delete impl;
		}
	
	private:
		friend class SymbolTable;
		Iterator(IteratorImpl *impl) : impl(impl){
		}
		IteratorImpl *impl;
	};
	/**
	 * Put a key-value pair into the symbol table.
	 * @param k The key.
	 * @param v The value.
	 * @note If {@var k} is already associated with a value,
	 * 	{@var v} will replace the old value.
	 */
	virtual void put(const K& k, const V& v) = 0;
	/**
	 * get the value assciated with a key.
	 * @param k The key.
	 * @return The pointer to the value associated with {@var k}.
	 * @note If {@var k} is not contained in the symbol table,
	 * 	{@code nullptr} should be returned. 
	 */
	virtual const V* get(const K& k) const  = 0;
	/**
	 * Delete a key and its associated value from the symbol table, if the key
	 * 	is contained in the symbol table.
	 * @param k The key to be deleted. 
	 */
	virtual void remove(const K& k) = 0;
	/**
	 * Test whether a key is contained in the symbol table.
	 * @param k The key to be tested.
	 * @return {@code true} is the {@var k} is contained in the symbol table,
	 * 	otherwise {@code false}. 
	 */
	virtual bool contains(const K& key) const {
		return get(key) != nullptr;
	}
	/**
	 * @return {@code true} if the symbol table is empty, {@code false} otherwise.
	 */
	virtual bool isEmpty() const {
		return size() == 0;
	}
	/**
	 * @return Number of key-value pairs in this symbol table. 
	 */
	virtual unsigned int size() const  = 0;

	/**
	 * remove all the elements in the symbol table.
	 */
	virtual void clear() = 0;

	Iterator begin() const {
		return Iterator(implBegin());
	}

	Iterator end() const {
		return Iterator(implEnd());
	}

	virtual ~SymbolTable(){
	}
};

}
#endif //MY_SYMBOL_TABLE_H
