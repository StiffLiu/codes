#ifndef MY_ORDERED_SYMBOL_TABLE_H
#define MY_ORDERED_SYMBOL_TABLE_H
#include "my_symbol_table.h"
namespace my_lib{
/**
 * An ordered symbol table is a symbol table 
 * where there is a total order relation between any two keys.
 */
template<class K, class V>
class OrderedSymbolTable : public SymbolTable<K, V>{
	typedef SymbolTable<K, V> Super;
public:

	bool remove(const K& k) override = 0;

	/**
	 * @return Pointer to the minimum key, {@code nullptr} doesn't exist.
	 */
	virtual const K* min() const = 0;
	
	/**
	 * @return Pointer to the maximum key, {@code nullptr) doesn't exist.
	 */
	virtual const K* max() const = 0;

	/**
	 * @return Largest key less than or equal to {@var k}.
	 */
	virtual const K* floor(const K& k) const = 0;
	
	/**
	 * @return Smallest key greater than or equal to {@var k}.
	 */
	virtual const K* ceil(const K& k) const = 0;
	
	/**
	 * @return Number of keys less than {@var k}.
	 */
	virtual unsigned int rank(const K& k) const = 0;

	/**
	 * @return the {@var index}th key in the ordered symbol table.
	 */
	virtual const K* select(unsigned int index) const = 0;

	/**
	 * Remove the minimum key.
	 * @return true if the min element is removed
	 */
	virtual bool removeMin(){
		const K* k = min();
		return nullptr == k ? false : remove(*k);
	}

	/**
	 * Remove the maximum key.
	 * @return true if the max element is removed
	 */
	virtual bool removeMax(){
		const K* k = max();
		return nullptr == k ? false : remove(*k);
	}

	unsigned int size() const override = 0;

	/**
	 * Number of keys in the range [{@var l}, {@var h}].
	 */
	virtual unsigned int size(const K& l, const K& h) const = 0;

	virtual bool isValid() const{
		return true;
	}
};
}
#endif //MY_ORDERED_SYMBOL_TABLE_H
