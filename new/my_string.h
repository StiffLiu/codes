#ifndef MY_LIB_STRING_H
#define MY_LIB_STRING_H
#include <cstring>
#include <iostream>
#include <cassert>
#include <unordered_map>
#include <type_traits>
#include <memory>
#include <vector>
#include <set>
#include <map>
#include <exception>
#include <queue>
#include <algorithm>

namespace my_lib{

/*
 * Most significant digit string sort algorithm.
 * */
template<class Item, class ItemTraits>
class MSD{
public:
	static void sort(Item items[], unsigned int index, unsigned int n, Item auxiliary[]){
		if(n <= ItemTraits::cutoff()){
			for (unsigned int i = 0; i < n; i++)
				for (unsigned int j = i; j > 0 && ItemTraits::compare(items[j], items[j-1], index) < 0; j--)
					ItemTraits::swap(items[j], items[j - 1]);
			return;
		}

		unsigned int itemsCount[ItemTraits::radix() + 2];
		std::memset(itemsCount, 0, sizeof(itemsCount));
		for(unsigned int i = 0;i < n;++ i){
			auto ch = ItemTraits::index(items[i], index);
			assert(ch < ItemTraits::radix() || ch == static_cast<unsigned int>(-1));
			++ itemsCount[ch + 2];
		}
		for(unsigned int i = 1;i < ItemTraits::radix() + 2;++ i)
			itemsCount[i] += itemsCount[i - 1];
		for(unsigned int i = 0;i < n;++ i){
			auto ch = ItemTraits::index(items[i], index);
			assert(ch < ItemTraits::radix() || ch == static_cast<unsigned int>(-1));
			auxiliary[itemsCount[ch + 1]++] = items[i];
		}
		for(unsigned int i = 0;i < n;++ i)
			items[i] = auxiliary[i];
		for(unsigned int i = 0;i < ItemTraits::radix();++ i){
		       	sort(items + itemsCount[i], index + 1, itemsCount[i + 1] - itemsCount[i], auxiliary);
		}
	}
public:
	static void sort(Item items[], unsigned int n){
		auto auxiliary = new Item[n];
		sort(items, 0, n, auxiliary);
		delete[] auxiliary;
	}
};

struct CStrTraits{

	static unsigned int radix(){
		return 256;
	}

	static unsigned int cutoff(){
		return 15;
	}

	template<class Str>
	static int compare(Str str1, Str str2, unsigned int index){
		return std::strcmp(reinterpret_cast<const char*>(str1 + index),
			           reinterpret_cast<const char*>(str2 + index));
	}

	template<class Str>
	static void swap(Str& str1, Str& str2){
		std::swap(str1, str2);
	}

	template<class Str>
	static unsigned int index(Str& str, unsigned int i){
		assert(str[i] >= 0);
		return static_cast<unsigned int>(str[i]) <= 0 ? -1 : str[i];
	}

	//assuming that i won't go out of scope of the c-string
	template<class Str>
	static auto addr(Str& str, unsigned int i) -> decltype(&str[i]){
		return str[i] == 0 ? nullptr : &str[i];
	}

	template<class Str>
	static void swap(Str array[], unsigned int i, unsigned int j){
		std::swap(array[i], array[j]);
	}

	static const char *toKey(const std::string& str){
		return str.c_str();
	}

	static char to(unsigned int index){
		return static_cast<char>(index);
	}
};

/*
 * Just an interesting one, with all non-digits considered as equal and less than digit.
 * */
struct NumCStrTraits : public CStrTraits{

	static unsigned int radix(){
		return 11;
	}

	static unsigned int cutoff(){
		return 15;
	}

	template<class Str>
	static int compare(Str str1, Str str2, unsigned int index){
		str1 += index; str2 += index;
		while(*str1 && *str2){
			if(*str1 != *str2){
				if(*str1 >= '0' && *str1 <= '9'){
					return *str2 >= '0' && *str2 <= '9'&& *str2 > *str1 ? -1 : 1;
				}else if(*str2 >= '0' && *str2 <= '9')
					return -1;
			}
			++ str1;
			++ str2;
		}
		if(*str2 >= '0' && *str2 <= '9')
			return -1;
		if(*str1 >= '0' && *str1 <= '9')
			return 1;
		if(*str1 < *str2)
			return -1;
		return *str1 > *str2 ? 1 : 0;
	}

	template<class Str>
	static unsigned int index(Str& str, unsigned int i){
		char ch = str[i];
		if(ch >= '0' && ch <= '9')
			return ch - '0' + 1;
		return ch == 0 ? -1 : 0;
	}
};

using CStrMSD = MSD<const char*, CStrTraits>;

template<class T>
bool compareForTest(T t1, T t2){
	return t1 < t2;
}

template<>
bool compareForTest(const char *s1, const char *s2){
	return std::strcmp(s1, s2) < 0;
}


template<class Array, class SizeType>
bool isSorted(Array array, SizeType low, SizeType high){
	for(SizeType i = low + 1;i < high;++ i){
		if(compareForTest(array[i], array[i - 1])){
			for(SizeType j = low;j < high;++ j)
				std::cout << j << "th : " << array[j] << std::endl;
		 	return false;
		}
	}
	return true;
}

template<class ArrayTraits>
class Quick3WaySort{
template<class Array, class SizeType>
static void sort(Array array, SizeType low, SizeType high, SizeType index){
begin_of_sort:
	if(high == low + 2){
		auto gt = high - 1;
		auto v1 = ArrayTraits::addr(array[low], index);
		auto v2 = ArrayTraits::addr(array[gt], index);
		while(decltype(v1)() != v1 && decltype(v2)() != v2 && *v1 == *v2){
			++index;
			v1 = ArrayTraits::addr(array[low], index);
			v2 = ArrayTraits::addr(array[gt], index);
		}
		if(decltype(v2)() != v2 && (decltype(v1)() == v1 || *v2 < *v1)){
			ArrayTraits::swap(array, low, gt);
		}
		assert(isSorted(array, low, high));
		return;
	}

	/*std::cout << "sorting array.. : " << std::endl;
	for(SizeType j = low;j < high;++ j)
		std::cout << j << "th : " << array[j] << std::endl;*/

	auto lt = low;
	auto gt = high - 1;
	auto i = low;++ i;
	auto v = ArrayTraits::addr(array[low], index);
	if(decltype(v)() != v){
		while (i <= gt){
			auto u = ArrayTraits::addr(array[i], index);
			if(decltype(u)() == u || *u < *v) ArrayTraits::swap(array, lt++, i++);
			else if(*u > *v) ArrayTraits::swap(array, i, gt--);
			else ++i;
		}
		/*for(SizeType i = low;i < lt;++ i) 
			assert(*ArrayTraits::addr(array[i], index) < *ArrayTraits::addr(array[lt], index));
		for(SizeType i = lt + 1;i <= gt;++ i) 
			assert(*ArrayTraits::addr(array[lt], index) == *ArrayTraits::addr(array[i], index));
		for(SizeType i = gt + 1;i < high;++ i) 
			assert(*ArrayTraits::addr(array[lt], index) < *ArrayTraits::addr(array[i], index));*/
		//lt - low > gt + 1 - lt
		if(2 * lt > gt + 1 + low){

			//gt + 1 - lt > 1
			if(gt > lt) sort(array, lt, gt + 1, index + 1);
			//assert(isSorted(array, lt, gt + 1));

			//high - (gt + 1) > lt - low
			if(high + low > lt + gt + 1){
				//lt - low > 1
				if(lt > low + 1) sort(array, low, lt, index);
				// sort(array, gt + 1, high, index);
				//assert(isSorted(array, low, lt));
				
				//high - (gt + 1) < 2
				if(high < gt + 3) return;
				//avoid recursive call : sort(array, gt + 1, high, index)
				low = gt + 1;
				goto begin_of_sort;
			}

			//high - (gt + 1) > 1
			if(high > gt + 2) sort(array, gt + 1, high, index);
			//assert(isSorted(array, gt + 1, high));
			
			//lt - low < 2
			if(lt < low + 2) return;
			//avoid recursiv call : sort(array, low, lt, index);
			high = lt;
			goto begin_of_sort;
		}else{
			if(lt > low + 1) sort(array, low, lt, index);
			//assert(isSorted(array, low, lt));

			// high - (gt + 1) > (gt + 1) - lt
			if(high + lt > 2 * gt + 2){
				//gt + 1 - lt > 1
				if(gt > lt) sort(array, lt, gt + 1, index + 1);
				//assert(isSorted(array, lt, gt + 1));

				//high - (gt + 1) < 2
				if(high < gt + 3) return;
				//avoid recursive call : sort(array, gt + 1, high, index)
				low = gt + 1;
				goto begin_of_sort;
			}

			//high - (gt + 1) > 1
			if(high > gt + 2) sort(array, gt + 1, high, index);
			//assert(isSorted(array, gt + 1, high));

			//gt + 1 - lt < 2
			if(gt < lt + 1) return;
			//avoid recursive call : sort(array, lt, gt + 1, index + 1);
			low = lt; high = gt + 1; ++index;
			goto begin_of_sort;
		}
	}else{
		while (i <= gt)
			if(decltype(v)() != ArrayTraits::addr(array[i], index)) ArrayTraits::swap(array, i, gt--);
			else i++;
		//avoid recursive call : sort(array, gt, high, index);
		//assert(isSorted(array, low, gt + 1));
		low = gt + 1;
		if(high > low + 1) goto begin_of_sort;
		//assert(isSorted(array, gt, high));
	}

}
public:
	template<class Array, class SizeType>
	static void sort(Array array, SizeType n){
		if (n <= 1) return;
		sort(array, SizeType(), n, SizeType());
	}
};

using CStrQ3WS = Quick3WaySort<CStrTraits>;

template<class Key, class Value, class KeyTraits>
class Trie{
	struct Node{
		Value value = Value();
		Node *node = nullptr;
		unsigned int size = 0;
		bool hasValue = false;
		unsigned int validateSize(){
			unsigned int sum = 0;
			if(node != nullptr){
				for(unsigned int i = 0;i < KeyTraits::radix();++ i){
					auto one = node[i].validateSize();
					if(one == static_cast<unsigned int>(-1))
						return static_cast<unsigned int>(-1);
					sum += one;
				}

				//It can be proved that sum must be greater than zero,
				//if we reach here. But do a sanity check.
				if(sum == 0)
					return -1;
			}
			return (hasValue ? size == sum + 1 : size == sum) ? size : -1;
		}
		bool isValid(){
			return validateSize() > 0;
		}

		void deleet(){
			if(node != nullptr)
				for(unsigned int i = 0;i < KeyTraits::radix();++ i)
					node[i].deleet();
			delete[] node;
			node = nullptr;
		}
		void removeChild(){
			deleet();
			size = (hasValue ? 1 : 0);
		}

		template<class Seq, class Prefix>
		void getAllKeys(Seq& seq, const Prefix& prefix){
			if(hasValue)
				seq.push_back(prefix);
			if(node != nullptr)
				for(unsigned int i = 0;i < KeyTraits::radix();++ i)
					node[i].getAllKeys(seq, prefix + KeyTraits::to(i));
		}

		unsigned int numNodes(){
			unsigned int count = 1;
			if(node != nullptr){
				for(unsigned int i = 0;i < KeyTraits::radix();++ i)
					count += node[i].numNodes();
			}
			return count;
		}

	};


	Node* getNode(const Key& key, bool exactMatch = true) const{
		unsigned int index = 0;
		auto subtree = root;
		//Loop invariant:
		//	subtree is the node where the index-length prefix string of key resides.
		auto i = KeyTraits::index(key, index);
		while(i != static_cast<unsigned int>(-1) && subtree->node != nullptr){
			assert(i < KeyTraits::radix());
			subtree = &subtree->node[i];
			i = KeyTraits::index(key, ++ index);
		}
		return i == static_cast<unsigned int>(-1) || !exactMatch ? subtree : nullptr;
	}

	Node *root = nullptr;
public:
	Trie(){
	}
	~Trie(){
		if(root != nullptr)
			root->deleet();
		delete root;
	}

	Trie(const Trie&) = delete;
	Trie& operator=(const Trie&) = delete;

	bool put(const Key& key, const Value& value){
		if(root == nullptr)
			root = new Node();
		//std::cout << "insert key : " << key << "\nvalue : " << value << std::endl;
		unsigned int index = 0;
		auto subtree = root;
		//Loop invariant:
		//	subtree is the node where the index-length prefix string of key resides.
		unsigned int i = KeyTraits::index(key, index);
		while(i != static_cast<unsigned int>(-1)){
			assert(i < KeyTraits::radix());
			if(subtree->node == nullptr){
				subtree->node = new Node[KeyTraits::radix()];
			}
			subtree = &subtree->node[i];
			++index;
			i = KeyTraits::index(key, index);
		}
		bool isReplace = subtree->hasValue;
		subtree->value = value;
		subtree->hasValue = true;

		//If not replace, increment the size of each node along the way till the node which the key resides.
		if(!isReplace){
			//std::cout << "key is new" << std::endl;
			++ subtree->size;
			index = 0;
			i = KeyTraits::index(key, index);
			subtree = root;
			while(i != static_cast<unsigned int>(-1)){
				assert(i < KeyTraits::radix());
				assert(subtree != nullptr);
				++subtree->size;
				subtree = &subtree->node[i];
				++index;
				i = KeyTraits::index(key, index);
			}
		}
		assert(root->isValid());
		return !isReplace;
	}

	const Value* get(const Key& key) const{
		auto node = getNode(key);
		return node != nullptr && node->hasValue ? &node->value : nullptr;
	}

	double usage(){
		return root->size / (double)root->numNodes();
	}

	bool remove(const Key& key){
		auto subtree = getNode(key);
		if(subtree == nullptr || !subtree->hasValue)
			return false;

		unsigned int index = 0;
		subtree->hasValue = false;
		//clear the value.
		subtree->value = Value();
		assert(subtree->size > 0);

		auto i = KeyTraits::index(key, index);
		assert(i != static_cast<unsigned int>(-1) || subtree == root);
		subtree = root;
		while(true){//i != static_cast<unsigned int>(-1)){
			assert(subtree != nullptr);
			assert(subtree->size > 0);
			//when the children node contains no value
			//remove them.
			if(subtree->size == 1 || (subtree->size == 2 && subtree->hasValue)){
				assert(subtree->size == 1 || subtree->hasValue);
				subtree->removeChild();
				break;
			}
			--subtree->size;
			if(i == static_cast<unsigned int>(-1))
				break;
			assert(i < KeyTraits::radix());
			subtree = &subtree->node[i];
			++index;
			i = KeyTraits::index(key, index);
		}
		if(root->size == 0){
			delete root;
			root = nullptr;
		}
		assert(root == nullptr || root->isValid());
		return true;
	}

	bool contains(const Key& key){
		return get(key) != nullptr;
	}

	bool isEmpty(){
		return root == nullptr;
	}

	unsigned int size(){
		return root == nullptr ? 0 : root->size;
	}

	template<class Seq, class Prefix>
	void keysWithPrefix(Seq& seq, const Prefix& prefix){
		auto node = getNode(KeyTraits::toKey(prefix), false);
		if(node != nullptr)
			node->getAllKeys(seq, prefix);
	}

	template<class Prefix>
	bool longestPrefixOf(const Key& s, Prefix& key){
		if(root == nullptr) return false;
		auto node = root;
		unsigned int index = 0;
		unsigned int i = 0;
		unsigned int longestPos = -1;
		if(root->hasValue) longestPos = 0;
		while((index = KeyTraits::index(s, i)) != static_cast<unsigned int>(-1) && nullptr != node->node){
			++i;
			if(node->node[index].hasValue)	longestPos = i;
			node = &node->node[index];
		}
		if(longestPos == static_cast<unsigned int>(-1)) return false;
		for(i = 0;i < longestPos;++ i) key.push_back(KeyTraits::to(KeyTraits::index(s, i)));
		return true;
	}

};

template<class Char, class Value, class Comparator = std::less<Char> >
class TernarySearchTrie{
	struct Node{
		Node *l = nullptr, *r = nullptr, *m = nullptr;
		unsigned int size = 0;
		Value value = Value();
		Char ch = Char();
		bool hasValue = false;
		bool isValid(){
			if(size <= 0){
				std::cerr << "size is zero" << std::endl;
			       	return false;
			}
			unsigned int count = 0;
			if(l != nullptr){
				if(!l->isValid()) return false;
				count += l->size;
			}
			if(r != nullptr){
				if(!r->isValid()) return false;
				count += r->size;
			}
			if(m != nullptr){
				if(!m->isValid()) return false;
				count += m->size;
			}
			if(hasValue) ++ count;
			if(count != size){
				std::cerr << "size mismatch, real " << size << ", counted " << count << std::endl;
			       	return false;
			}
			return true;
		}
		template<class Seq, class Prefix>
		void getAllKeys(Seq& seq, Prefix& prefix) const{
			prefix.push_back(ch);
			if(hasValue) seq.push_back(prefix);
			if(nullptr != m) m->getAllKeys(seq, prefix);
			prefix.pop_back();
			if(nullptr != l) l->getAllKeys(seq, prefix);
			if(nullptr != r) r->getAllKeys(seq, prefix);
		}

	};
	static void deleteNode(Node* node){
		if(node->l != nullptr) deleteNode(node->l);
		if(node->r != nullptr) deleteNode(node->r);
		if(node->m != nullptr) deleteNode(node->m);
		delete node;
	}
	template<class Key, class TST>
	static auto getNode(const Key& key, TST *tst) ->decltype(tst->root) {
		unsigned int i = 0;
		auto node = key.size() > 0 ? (tst->root != nullptr ? tst->root->m : nullptr) : tst->root;
		while(node != nullptr){
			while(node != nullptr){
				if(tst->comparator(key[i], node->ch)){
					node = node->l;
				}else if(tst->comparator(node->ch, key[i])){
					node = node->r;
				}else{
					++ i;
					if(i >= key.size())
						return node;
					node = node->m;
				}
			}
		}

		return node;
	}
	Node *root = nullptr;
	Comparator comparator;
public:
	TernarySearchTrie(Comparator comparator = Comparator()) : comparator(comparator){
	}
	TernarySearchTrie(const TernarySearchTrie&) = delete;
	TernarySearchTrie& operator=(const TernarySearchTrie&) = delete;
	TernarySearchTrie(TernarySearchTrie&& tst){
		root = tst.root;
		tst.root = nullptr;
	}
	~TernarySearchTrie(){
		clear();
	}
	template<class Key>
	const Value* get(const Key& key) const{
		auto node = getNode(key, this);
		return node == nullptr ? nullptr : (node->hasValue ? &(node->value) : nullptr);
	}
	bool isValid(){
		return root == nullptr || root->isValid();
	}
	bool isEmpty(){
		return root == nullptr;
	}
	unsigned int size(){
		return root == nullptr ? 0 : root->size;
	}
	void clear(){
		if(root != nullptr){
		 	deleteNode(root);
			root = nullptr;
		}
	}
	template<class Key>
	bool contains(const Key& key) const{
		return get(key) != nullptr;
	}
	template<class Key>
	bool put(const Key& key, const Value& value){
		if(root == nullptr)
			root = new Node();
		if(key.size() == 0){
			root->value = value;
			if(root->hasValue) return false;
			root->hasValue = true;
			++root->size;
			return true;
		}
		auto node = root->m;
		if(node == nullptr){
			root->m = new Node;
			node = root->m;
			node->ch = key[0];
		}
		unsigned int i = 0;
		while(i < key.size()){
			while(true){
				if(comparator(key[i], node->ch)){
					if(node->l == nullptr){
						node->l = new Node;
						node = node->l;
						node->ch = key[i];
						break;
					}
					node = node->l;
				}else if(comparator(node->ch, key[i])){
					if(node->r == nullptr){
						node->r = new Node;
						node = node->r;
						node->ch = key[i];
						break;
					}
					node = node->r;
				}else
					break;
			}
			if(node->m == nullptr){
				++ i;
				while(i < key.size()){
					node->m = new Node;
					node = node->m;
					node->ch = key[i];
					++ i;
				}
				break;
			}
			++ i;
			node = node->m;
		}
		node->value = value;

		bool isReplace = node->hasValue;
		if(!isReplace){
			node->hasValue = true;
			++ root->size;
			node = root->m;
			i = 0;
			while(i < key.size()){
				++ node->size;
				if(comparator(key[i], node->ch)){
					node = node->l;
				}else if(comparator(node->ch, key[i])){
					node = node->r;
				}else{
					node = node->m;
					++ i;
				}
			}

		}

		assert(isValid());
		return !isReplace;
	}
	template<class Key>
	bool remove(const Key& key){
		auto node = getNode(key, this);
		if(node == nullptr) return false;
		if(root->size == 1){
			deleteNode(root);
			root = nullptr;
			return true;
		}
		--root->size;
		if(node == root){
			root->hasValue = false;
			return true;
		}

		unsigned int i = 0;
		Node **link = &root->m;
		node = root->m;
		i = 0;
		while(node != nullptr){
			if(node->size == 1){
				deleteNode(node);
				*link = nullptr;
				break;
			}
			--node->size;
			if(comparator(key[i], node->ch)){
				link = &node->l;
				node = node->l;
			}else if(comparator(node->ch, key[i])){
				link = &node->r;
				node = node->r;
			}else{
				link = &node->m;
				node = node->m;
				++ i;
			}
		}
		assert(isValid());
		return true;
	}

	template<class Seq, class Prefix>
	void keysWithPrefix(Seq& seq, const Prefix& prefix) const{
		auto node = getNode(prefix, this);
		if(nullptr == node) return;

		Prefix mutablePrefix = prefix;
		if(node->hasValue) seq.push_back(mutablePrefix);
		node = node->m;
		if(nullptr != node) node->getAllKeys(seq, mutablePrefix);
	}

	template<class Key>
	const Node* longestPrefixOf(const Key& s, unsigned int& index) const {
		if(nullptr == root) return nullptr;
		const Node* value = nullptr;
		if(root->hasValue){
			value = root;
		 	index = 0;
		}

		if(s.size() > 0){
			unsigned int i = 0;
			auto node = root->m;

			while(nullptr != node){
				if(comparator(s[i], node->ch)) node = node->l;
				else if(comparator(node->ch, s[i]))	node = node->r;
				else{
					++ i;
					if(node->hasValue){
						value = node;
						index = i;
					}
					if(i >= s.size()) break;
					node = node->m;
				}
			}
		}
		return value;
	}

	template<class Key>
	void longestPrefixOf(const Key& s, Value& value, unsigned int& index) const {
		auto node = longestPrefixOf(s, index);
		if(nullptr != node) value = node->value;
	}

	template<class Key, class Prefix>
	bool longestPrefixOf(const Key& s, Prefix& key) const {
		unsigned int index = -1;
		longestPrefixOf(s, index);
		if(index == static_cast<unsigned int>(-1)) return false;
		key = Prefix();
		for(unsigned int i = 0;i < index;++ i) key.push_back(s[i]);
		return true;
	}
};

template<class Str, class DFA>
unsigned int buildDFA(const Str& str, unsigned int n, DFA& dfa){
	if(n == 0) return static_cast<unsigned int>(-1);

	unsigned int restartState = 0;
	dfa[0][str[0]] = 1;
	for(unsigned int i = 1;i < n;++ i){
		dfa[i] = dfa[restartState];
		dfa[i][str[i]] = i + 1;
		restartState = dfa[restartState][str[i]];
	}
	return restartState;
}

template<class Str>
unsigned int bruteForceSearch(const Str& src, unsigned int n,
  const Str& pattern, unsigned int m){
	for (unsigned int i = 0;i + m <= n;++ i){
		unsigned int j = 0;
		for (;j < m;++ j)
			if (pattern[j] != src[i + j]) break;
		if (j == m) return i;
	}
	return n;
}

template<class Str, class Algo, class All>
void searchAllBase(const Str& src, unsigned int n, const Str& pattern, unsigned int m,
  All& all, Algo algo){
	unsigned int index = 0;
	unsigned int count = n;
	while((index = algo(src + (n - count), count, pattern, m)) != count){
		all.insert(n + index - count);
		count = count - index;
		if (count <= 0) break;
		--count;
	}
}

template<class Str>
unsigned int bruteForceSearch(const Str& src, const Str& pattern){
	return bruteForceSearch(src, src.size(), pattern, pattern.size());
}

template<class Str, class All>
void bruteForceSearch(const Str& src, unsigned int n,
  const Str& pattern, unsigned int m, All& all){
	unsigned int (*func)(const Str&, unsigned int, const Str&, unsigned int) = bruteForceSearch;
	searchAllBase(src, n, pattern, m, all, func);
}

template<class Str, class All>
void bruteForceSearch(const Str& src, const Str& pattern, All& all){
	if (src.size() == 0 || pattern.size() == 0) return;
	return bruteForceSearch(&src[0], src.size(), &pattern[0], pattern.size(), all);
}

template<class Str>
unsigned int kmpSearch(const Str& src, unsigned int n,
	const Str& pattern, unsigned int m){
	typedef typename std::remove_cv<typename std::remove_reference<decltype(src[0])>::type >::type Char;
	std::unordered_map<unsigned int,
		std::unordered_map<Char, unsigned int> > dfa;
	buildDFA(pattern, m, dfa);

	// For empty string, always match.
	if (m == 0) return 0;
	unsigned int j = 0;
	for (unsigned int i = 0;i < n && n + j >= m + i;++ i){
		auto pos = dfa.find(j);
		if (pos == dfa.end()) j = 0;
		else{
			auto pos1 = pos->second.find(src[i]);
			j = ((pos1 == pos->second.end()) ? 0 : pos1->second);
		}
		if (j == m) return i + 1- m;
	}
	return n;
}

template<class Str, class All>
void kmpSearch(const Str& src, unsigned int n,
	const Str& pattern, unsigned int m, All& all){
	typedef typename std::remove_cv<typename std::remove_reference<decltype(src[0])>::type >::type Char;
	std::unordered_map<unsigned int,
		std::unordered_map<Char, unsigned int> > dfa;
	unsigned int restartState = buildDFA(pattern, m, dfa);

	// For empty string, always match.
	if (m == 0) return;
	unsigned int j = 0;
	for (unsigned int i = 0;i < n && n + j >= m + i;++ i){
		auto pos = dfa.find(j);
		if (pos == dfa.end()) j = 0;
		else{
			auto pos1 = pos->second.find(src[i]);
			j = ((pos1 == pos->second.end()) ? 0 : pos1->second);
		}
		if (j == m){
			all.insert(i + 1- m);
			j = restartState;
		}
	}
}

template<class Str>
unsigned int kmpSearch(const Str& src, const Str& pattern){
	return kmpSearch(src, src.size(), pattern, pattern.size());
}

template<class Str, class All>
void kmpSearch(const Str& src, const Str& pattern, All& all){
	if (src.size() == 0 || pattern.size() == 0) return;
	return kmpSearch(&src[0], src.size(), &pattern[0], pattern.size(), all);
}

template<class Str>
unsigned int boyerMooreSearch(const Str& src, unsigned int n,
	const Str& pattern, unsigned int m){
	if (m == 0) return 0;
	if (n < m) return n;

	typedef typename std::remove_cv<typename std::remove_reference<decltype(src[0])>::type >::type Char;
	std::unordered_map<Char, unsigned int> largestIndices;
	for(unsigned int i = 0;i < m;++ i)
		largestIndices[pattern[i]] = i;
	for(unsigned int i = 0;i + m <= n;){
		unsigned int j = m - 1;
		for(;j != static_cast<unsigned int>(-1) && src[i + j] == pattern[j];-- j);
		if (j == static_cast<unsigned int>(-1)) return i;
		auto pos = largestIndices.find(src[i + j]);
		if (pos == largestIndices.end()) i += j + 1;
		else if (pos->second < j) i += j - pos->second;
		else ++i;
	}
	return n;
}

template<class Str>
unsigned int boyerMooreSearch(const Str& src, const Str& pattern){
	return boyerMooreSearch(src, src.size(), pattern, pattern.size());
}

template<class Str, class All>
void boyerMooreSearch(const Str& src, unsigned int n,
  const Str& pattern, unsigned int m, All& all){
	unsigned int (*func)(const Str&, unsigned int, const Str&, unsigned int) = boyerMooreSearch;
	searchAllBase(src, n, pattern, m, all, func);
}

template<class Str, class All>
void boyerMooreSearch(const Str& src, const Str& pattern, All& all){
	if (src.size() == 0 || pattern.size() == 0) return;
	return boyerMooreSearch(&src[0], src.size(), &pattern[0], pattern.size(), all);
}

template<class Str>
unsigned int rabinKarpSearch(const Str& src, unsigned int n,
	const Str& pattern, unsigned int m,
	unsigned int radius = 256, unsigned long long prime = 429496729uLL/*Not neccessarily be a prime number*/){
	if (m == 0) return 0;
	if (n < m) return n;

	unsigned long long rmHash = 1;
	unsigned long long patHash = 0;
	unsigned long long srcHash = 0;
	for(unsigned int i = 0;i < m - 1;++ i){
		rmHash = radius * rmHash % prime;
		patHash = (pattern[i] + patHash * radius) % prime;
		srcHash = (src[i] + srcHash * radius) % prime;
	}
	patHash = (pattern[m - 1] + patHash * radius) % prime;
	srcHash = (src[m - 1] + srcHash * radius) % prime;


	unsigned int index = m;
	while(index < n){
		if(srcHash == patHash){
			// probably do a real match here.
			return index - m;
		}
		srcHash = (prime - (src[index - m] * rmHash % prime) + srcHash) % prime;
		srcHash = (srcHash * radius + src[index]) % prime;
		++ index;
	}

	if(srcHash == patHash){
			// probably do a real match here.
		return index - m;
	}
	return n;
}

template<class Str>
unsigned int rabinKarpSearch(const Str& src, const Str& pattern){
	return rabinKarpSearch(src, src.size(), pattern, pattern.size());
}

template<class Str, class All>
void rabinKarpSearch(const Str& src, unsigned int n,
  const Str& pattern, unsigned int m, All& all){
	searchAllBase(src, n, pattern, m, all,
	  [](const Str& str, unsigned int n, const Str& pat, unsigned int m){
	  	return rabinKarpSearch(str, n, pat, m);
	  }
	);
}

template<class Str, class All>
void rabinKarpSearch(const Str& src, const Str& pattern, All& all){
	if (src.size() == 0 || pattern.size() == 0) return;
	return rabinKarpSearch(&src[0], src.size(), &pattern[0], pattern.size(), all);
}

template<class TDArraySrc, class TDArrayPat>
std::pair<unsigned int, unsigned int> bruteForceSearch(const TDArraySrc& src, unsigned int m, unsigned int n,
	const TDArrayPat& pat, unsigned int h, unsigned int v){
	if (h <= m && v <= n && m > 0 && n > 0){
		for(unsigned int i = 0;i <= m - h;++ i)
			for(unsigned int j = 0;j <= n - v;++ j){
				for(unsigned int k = 0;k < h;++ k)
					for(unsigned int s = 0;s < v;++ s)
						if (pat[k][s] != src[i + k][j + s])
							goto NotFound;
				return {i, j};
			NotFound:;
			}
	}
	return {static_cast<unsigned int>(-1), static_cast<unsigned int>(-1)};
}
template<class TDArraySrc, class TDArrayPat>
std::pair<unsigned int, unsigned int> rabinKarpSearch(const TDArraySrc& src, unsigned int m, unsigned int n,
	const TDArrayPat& pat, unsigned int h, unsigned int v, unsigned int radius = 256, unsigned long long prime = 429496729uLL){
	if (h <= m && v <= n && m > 0 && n > 0){
		unsigned long long patHash = 0, srcHash = 0;
		unsigned int row = m - h + 1;
		std::unique_ptr<unsigned long long> srcHashColsPtr(new unsigned long long[n]);//, [](unsigned long long *p){delete[] p;});
		unsigned long long *srcHashCols = srcHashColsPtr.get();
		const unsigned int radiusCol = radius;
		const unsigned long long primeCol = prime;
		unsigned long long rvHash = 1;
		unsigned long long rhHash = 1;
		for (unsigned int i = 0;i < v - 1;++ i){
			rvHash = rvHash * radiusCol % primeCol;
		}
		for (unsigned int i = 0;i < h - 1;++ i){
			rhHash = rhHash * radius % prime;
		}
		for (unsigned int i = 0;i < v;++ i){
			unsigned long long patHashCol = 0;
			srcHashCols[i] = 0;
			for(unsigned int j = 0;j < h;++ j){
				patHashCol = (patHashCol * radius + pat[j][i]) % prime;
				srcHashCols[i] = (srcHashCols[i] * radius + src[j][i]) % prime;
			}
			patHash = (patHash * radiusCol + patHashCol) % primeCol;
			srcHash = (srcHash * radiusCol + srcHashCols[i]) % primeCol;
		}
		unsigned int i = 0;
		for (;i < row;++ i){
			if(i != 0){
				srcHash = 0;
				for(unsigned int j = 0;j < v;++ j){
					srcHashCols[j] = (prime - (src[i-1][j] * rhHash % prime) + srcHashCols[j]) % prime;
					srcHashCols[j] = (srcHashCols[j] * radius + src[i + h - 1][j]) % prime;
					srcHash = (srcHash * radiusCol + srcHashCols[j]) % primeCol;
				}
			}
			if(patHash == srcHash){
				// probably do a real match here
				return {i, 0};
			}
			for(unsigned int j = v;j < n;++ j){
				if(i == 0){
					srcHashCols[j] = 0;
					for(unsigned int k = 0;k < h;++ k){
						srcHashCols[j] = (srcHashCols[j] * radius + src[k][j]) % prime;
					}
				}else{
					srcHashCols[j] = (prime - (src[i-1][j] * rhHash % prime) + srcHashCols[j]) % prime;
					srcHashCols[j] = (srcHashCols[j] * radius + src[i + h - 1][j]) % prime;
				}
				srcHash = (primeCol - (srcHashCols[j - v] * rvHash % primeCol) + srcHash) % primeCol;
				srcHash = (srcHash * radiusCol + srcHashCols[j]) % primeCol;
				if(patHash == srcHash){
					// probably do a real match here
					return {i, j - v + 1};
				}

			}
		}
	}
	return {static_cast<unsigned int>(-1), static_cast<unsigned int>(-1)};
}

class RegExpr{
	static const char Char = '\0';
	struct State{
		char type = 0;
		char ch = 0;
		unsigned int start = 0, end = 0;
		State(char ch, char type = 0, unsigned int start = 0, unsigned int end = 0)
			: type(type), ch(ch), start(start), end(end){
		}
	};
	std::vector<State> pattern;
	std::vector<unsigned int> adjVertices;
	void calcStates(std::set<unsigned int>& states){
		std::vector<unsigned int> allStates(states.begin(), states.end());
		for (size_t i = 0;i < allStates.size();++ i){
			if(allStates[i] < pattern.size()){
				const auto& state = pattern[allStates[i]];
				for(unsigned int j = state.start;j < state.end;++ j)
					if (states.find(adjVertices[j]) == states.end()){
						allStates.push_back(adjVertices[j]);
						states.insert(adjVertices[j]);
					}
			}
		}
	}
	void printGraph(){
		for(size_t i = 0;i < pattern.size();++ i){
			std::cout << i << " : ";
			for(unsigned int j = pattern[i].start;j < pattern[i].end;++ j){
				std::cout << adjVertices[j] << ' ';
			}
			std::cout << std::endl;
		}
	}
public:
	RegExpr(const char* expression){
		std::string expr("(");
		expr += expression;
		expr.push_back(')');

		std::vector<std::pair<size_t, size_t> > ops;
		//std::set<std::pair<size_t, size_t> > edges;
		std::multimap<size_t, size_t> edges;

		for (size_t i = 0;i < expr.size();++ i){
			if (expr[i] == '\\'){
				++ i;
				if(i >= expr.size()) break;
				pattern.push_back(State(expr[i]));
			}else{
				if (expr[i] == '('){
					ops.push_back({i, pattern.size()});
					edges.insert({pattern.size(), pattern.size() + 1});
					pattern.push_back(State('\0', '('));
					if (i + 1 < expr.size() && expr[i + 1] == '*'){
						//throw;
						++ i;
						std::cerr << "'*' at position " << i << " will be ignored." << std::endl;
					}
				}else if (expr[i] == '.'){
					if (i + 1 < expr.size() && expr[i + 1] == '*'){
						// * will be push_back to pattern in next iteration
						edges.insert({pattern.size(), pattern.size() + 1});
						edges.insert({pattern.size() + 1, pattern.size()});
					}
					pattern.push_back(State('\0', '.'));
				}else if (expr[i] == '*'){
					edges.insert({pattern.size(), pattern.size() + 1});
					pattern.push_back(State('\0', '*'));
					if (i + 1 < expr.size() && expr[i + 1] == '*'){
						//throw;
						++ i;
						std::cerr << "'*' at position " << i << " will be ignored." << std::endl;
					}
				}else if (expr[i] == '|'){
					ops.push_back({i, pattern.size()});
					pattern.push_back(State('\0', '|'));
					if (i + 1 < expr.size() && expr[i + 1] == '*'){
						//throw;
						++ i;
						std::cerr << "'*' at position " << i << " will be ignored." << std::endl;
					}
				}else if (expr[i] == ')'){
					edges.insert({pattern.size(), pattern.size() + 1});
					size_t j = ops.size() - 1;
					for(;j != static_cast<size_t>(-1) && expr[ops[j].first] != '(';-- j);
					if (j == static_cast<size_t>(-1)) throw;
					for(decltype(j) k = j + 1;k < ops.size();++ k){
						assert(expr[ops[k].first] == '|');
						edges.insert({ops[j].second, ops[k].second + 1});
						edges.insert({ops[k].second, pattern.size()});
					}
					ops.erase(ops.begin() + j, ops.end());

					if (i + 1 < expr.size() && expr[i + 1] == '*'){
						// * will be push_back to pattern in next iteration
						edges.insert({ops[j].second, pattern.size() + 1});
						edges.insert({pattern.size() + 1, ops[j].second});
					}
					pattern.push_back(State('\0', ')'));
				}else{
					if (i + 1 < expr.size() && expr[i + 1] == '*'){
						// * will be push_back to pattern in next iteration
						edges.insert({pattern.size(), pattern.size() + 1});
						edges.insert({pattern.size() + 1, pattern.size()});
					}
					pattern.push_back(State(expr[i]));
				}
			}
		}

		if (!ops.empty()){
		 	//throw;
			std::cerr << "There're possiblely unmatched braces." << std::endl;
		}
		for(auto& kvp : edges){
			++ pattern[kvp.first].end;
		}
		unsigned int count = 0;
		for(auto& state : pattern){
			state.start = count;
			count += state.end;
			state.end = state.start;
		}
		adjVertices.resize(count);
		for(auto& kvp : edges){
			auto& state = pattern[kvp.first];
			adjVertices[state.end] = kvp.second;
			++state.end;
		}
		printGraph();
	}
	bool recognize(const std::string& src){
		if(pattern.empty()){
			return true;
		}
		std::set<unsigned int> indices;
		indices.insert(0);
		calcStates(indices);
		for (size_t i = 0;i < src.size();++ i){
			decltype(indices) newIndices;
			for (auto index : indices){
				const auto& state = pattern[index];
				if (state.type == '.' || (state.type == Char && src[i] == state.ch)){
					newIndices.insert(index + 1);
				}
			}
			if (newIndices.empty()) return false;
			indices.swap(newIndices);
			calcStates(indices);
		}
		return (indices.find(pattern.size()) != indices.end());
	}
};

template<class ForwardIterator>
void copy(ForwardIterator begin, ForwardIterator end, std::vector<bool>& bitArray){
	typedef typename std::remove_reference<typename std::remove_cv<decltype(*begin)>::type>::type Element;
	while(begin != end){
		Element e = *begin;
		for(unsigned int i = 0;i < sizeof(Element) * 8;++ i){
			bitArray.push_back(e % 2 == 1);
			e /= 2;
		}
		++begin;
	}
}

template<class Seq, class Element>
void copy(const std::vector<bool>& bitArray, Seq& seq){
	for(size_t i = 0;i < bitArray.size();++ i){
		Element e = Element();
		unsigned int j = std::min(bitArray.size(), i + sizeof(Element));
		for(;j > i;--j){
			e *= 2;
			if (bitArray[j - 1]) e += 1;
		}
		i += sizeof(Element);
	}
}

template<class BitArray, class CounterArray>
class RunLengthT{
public:
	static void compress(const BitArray& input, CounterArray& output, unsigned int maxLength = 255){
		unsigned int count = 0;
		bool lastValue = false;
		for(const auto& oneBit: input){
			if(oneBit != lastValue){
				output.push_back(count);
				count = 1;
				lastValue = !lastValue;
			}else if(count < maxLength){
				++count;
			}else{
				output.push_back(count);
				output.push_back(0);
				count = 1;
			}
		}
		output.push_back(count);
	}
	static void decompress(const CounterArray& input, BitArray& output){
		bool currentValue = false;
		for(auto count : input){
			while(count > 0){
				output.push_back(currentValue);
				--count;
			}
			currentValue = !currentValue;
		}
	}
};

using RunLength = RunLengthT<std::vector<bool>, std::vector<unsigned char> >;

struct BASplitor{
	BASplitor(const char separator = ' ', unsigned int block = 8,
		 unsigned int start = 0, unsigned int count = static_cast<unsigned int>(-1))
		: separator(separator), block(block), start(start), count(count){
	}

	friend BASplitor operator<<(std::ostream& os, BASplitor splitor){
		splitor.os = &os;
		return splitor;
	}

	template<class BitArray>
	friend std::ostream& operator<<(BASplitor splitor, const BitArray& bitArray){
		unsigned int size = bitArray.size();
		if(static_cast<unsigned int>(-1) != splitor.count) size = splitor.start + splitor.count;
		if(size > bitArray.size()) size = bitArray.size();
		for(unsigned int i = splitor.start;i < size;++ i){
			*splitor.os << bitArray[i];
			if((i + 1) % splitor.block == 0) *splitor.os << splitor.separator;
		}
		return *splitor.os;
	}
private:
	const char separator = ' ';
	unsigned int block = 8;
	unsigned int start = 0;
	unsigned int count = static_cast<unsigned int>(-1);
	std::ostream *os = nullptr;
};

class Huffman{
public:

	/*
	 * frequency : value to its frequency map, value is the index in an input array.
	 * index : assume 0 <= i < n, i and index[i] should equal each other.
	 * */
	template<class Frequency, class Index>
	static void computeFrequency(unsigned int n, Frequency& frequency, Index& index){
		for(unsigned int i = 0;i < n;++ i){
			auto pos = frequency.find(i);
			if(pos == frequency.end()){
				index[i] = i;
				frequency.insert({i, 1});
			}else{
				index[i] = pos->first;
				++pos->second;
			}
		}
	}

	template<class Frequency, class Tree, class NodeComparator>
	static void buildTree(const Frequency& frequency, Tree& treeCreator, NodeComparator comparator){
		typedef decltype(treeCreator.create(*frequency.begin())) Node;
		std::priority_queue<Node, std::vector<Node>, NodeComparator> pq(comparator);
		for(const auto& one : frequency) pq.push(treeCreator.create(one));
		while(!pq.empty()){
			auto left = pq.top();
			pq.pop();
			if(pq.empty()){
				treeCreator.root(left);
				return;
			}
			auto right = pq.top();
			pq.pop();
			pq.push(treeCreator.create(left, right));
		}
	}

	/*
	 * Build prefix free code
	 * */
	template<class Tree, class PFCode>
	static void buildPFCode(const Tree& tree, PFCode& pfCode){
		if(tree.empty()) return;
		typedef typename std::remove_cv<decltype(pfCode[tree.value(tree.root())])>::type NoCVType;
		typedef typename std::remove_reference<NoCVType>::type Code;
		std::queue<std::pair<decltype(tree.root()), Code> > codes;
		codes.push({tree.root(), Code()});
		while(!codes.empty()){
			bool isLeave = true;
			auto& node = codes.front();
			//std::cout << node.first << std::endl;
			if(tree.hasLeft(node.first)){
				isLeave = false;
				codes.push({tree.left(node.first), node.second});
				codes.back().second.push_back(false);
			}
			if(tree.hasRight(node.first)){
				isLeave = false;
				codes.push({tree.right(node.first), node.second});
				codes.back().second.push_back(true);
			}
			if(isLeave){
				pfCode[tree.value(node.first)] = node.second;
			}
			codes.pop();
		}
	}

	template<class Tree, class BitArray, class Node>
	static void saveSubTree(const Tree& tree, Node subTree, BitArray& bitArray){
		if(!tree.hasLeft(subTree) && !tree.hasRight(subTree)){
			bitArray.push_back(true);
			tree.saveValue(subTree, bitArray);
			return;
		}
		bitArray.push_back(false);
		if(tree.hasLeft(subTree)){
			saveSubTree(tree, tree.left(subTree), bitArray);
		}
		if(tree.hasRight(subTree)){
			saveSubTree(tree, tree.right(subTree), bitArray);
		}
	}

	template<class Tree, class BitArray>
	static void saveTree(const Tree& tree, BitArray& bitArray){
		if(tree.empty()) return;
		saveSubTree(tree, tree.root(), bitArray);
	}

	template<class Tree, class BitArray, class Node>
	static unsigned int readSubTree(const BitArray& bitArray,
		unsigned int start, Tree& tree, Node& n){
		if(start >= bitArray.size()) return static_cast<unsigned int>(-1);
		if(bitArray[start]) return tree.readValue(bitArray, start + 1, n);
		Node left = Node();
		auto nextStart = readSubTree(bitArray, start + 1, tree, left);
		if(static_cast<unsigned int>(-1) == nextStart)
			return static_cast<unsigned int>(-1);
		Node right = Node();
		nextStart = readSubTree(bitArray, nextStart, tree, right);
		if(static_cast<unsigned int>(-1) == nextStart)
			return static_cast<unsigned int>(-1);
		n = tree.create(left, right);
		return nextStart;
	}

	template<class BitArray, class Tree>
	static unsigned int readTree(const BitArray& bitArray, Tree& tree){
		typedef decltype(tree.root()) Node;
		Node root = Node();
		unsigned int current = readSubTree(bitArray, 0, tree, root);
		if(static_cast<unsigned int>(-1) != current){
		 	tree.root(root);
		}
		return current;
	}

	struct Node{
		unsigned int right, left;
		Node(unsigned int right, unsigned int left = -1) : right(right), left(left){
		}
		template<class Stream>
		friend Stream& print(Stream& os, Node n){
			if(static_cast<unsigned int>(-1) == n.left){
				os << "v=" << n.right;
			}else{
				os << "l=" << n.left << ", r= " << n.right;
			}
			return os;
		}
		template<class Stream>
		friend Stream& operator<<(Stream& os, Node n){
			return (os << "(").print(n) << ")";
		}
	};
	struct NodeWithFrequency : public Node{
		unsigned int frequency;
		NodeWithFrequency(unsigned int right, unsigned int frequency,
			unsigned int left = static_cast<unsigned int>(-1))
			: Node(right, left), frequency(frequency){
		}
		template<class Stream>
		friend Stream& operator<<(Stream& os, NodeWithFrequency n){
			return print(os << "(", n) << ", f=" << n.frequency << ")";
		}
	};
	template<class Tree>
	struct NodeComparator{
		const Tree& tree;
		NodeComparator(const Tree& tree) : tree(tree){
		}
		bool operator()(unsigned int n1, unsigned int n2){
			return tree.nodes[n2].frequency < tree.nodes[n1].frequency;
		}
	};

	template<class Node>
	struct Tree{

		std::vector<Node> nodes;

		bool empty() const {
			return nodes.empty();
		}

		unsigned int create(std::pair<unsigned int, unsigned int> data){
			auto index = nodes.size();
			nodes.push_back(Node(data.first, data.second));
			return index;
		}

		unsigned int create(unsigned int left, unsigned int right){
			assert(left < nodes.size() && right < nodes.size());
			auto index = nodes.size();
			nodes.push_back(Node(right, left));
			return index;
		}

		void root(unsigned int root){
			assert(nodes.size() - 1 == root);
		}

		unsigned int root() const {
			return nodes.size() - 1;
		}

		bool hasLeft(unsigned int n) const {
			return static_cast<unsigned int>(-1) != nodes[n].left;
		}

		bool hasRight(unsigned int n) const {
			return hasLeft(n);
		}

		unsigned int left(unsigned int n) const {
			return nodes[n].left;
		}

		unsigned int right(unsigned int n) const {
			return nodes[n].right;
		}

		unsigned int value(unsigned int n) const {
			return nodes[n].right;
		}
	};

	struct ToBit{
		template<class Value, class BitArray>
		void convert(Value value, BitArray& bitArray) const {
			static_assert(Value(-1) > Value(0), "Expecting an unsigned type");
			for(unsigned int i = 0;i < sizeof(Value) * 8;++ i){
				bitArray.push_back(value % 2 == 0 ? false : true);
				value = value / 2;
			}
		}
		template<class BitArray>
		void convert(bool value, BitArray& bitArray) const {
			bitArray.push_back(value);
		}
	};

	template<class BitArrayInput, class ValueToBitConverter>
	struct TreeWithFrequencyBase : public Tree<NodeWithFrequency>{
		typedef Tree<NodeWithFrequency> Super;
		using Super::create;
		const BitArrayInput& input;
		ValueToBitConverter converter;
		TreeWithFrequencyBase(const BitArrayInput& input, const ValueToBitConverter& converter)
			: input(input), converter(converter){
		}
		template<class BitArray>
		void saveValue(unsigned int n, BitArray& bitArray) const {
			converter.convert(input[Super::nodes[n].right], bitArray);
		}
	};

	template<class BitArrayInput, class ValueToBitConverter>
	struct TreeWithFrequency : public TreeWithFrequencyBase<BitArrayInput, ValueToBitConverter> {
		typedef TreeWithFrequencyBase<BitArrayInput, ValueToBitConverter> Super;
		using Super::create;
		unsigned int blockSize;
		TreeWithFrequency(const BitArrayInput& input, unsigned int blockSize,
			const ValueToBitConverter& converter = ValueToBitConverter()) : Super(input, converter), blockSize(blockSize){
		}
		unsigned int create(unsigned int left, unsigned int right){
			assert(left < Super::nodes.size() && right < Super::nodes.size());
			auto index = Super::nodes.size();
			Super::nodes.push_back(NodeWithFrequency(right,
					Super::nodes[left].frequency + Super::nodes[right].frequency, left));
			return index;
		}

		template<class BitArray>
		void saveValue(unsigned int n, BitArray& bitArray) const {
			assert(!Super::hasLeft(n) && !Super::hasRight(n));
			auto start = Super::nodes[n].right * blockSize;
			auto end = start + blockSize;
			//std::cout << Super::nodes[n].right << '\t';
			for(;start < end;++ start){
				Super::converter.convert(Super::input[start], bitArray);
				//std::cout << input[start];
			}
			//std::cout << std::endl;
		}
	};

	struct NodeWithCode : public Node{
		std::vector<bool> code;
		NodeWithCode(unsigned int left = static_cast<unsigned int>(-1),
			 unsigned int right = static_cast<unsigned int>(-1)) : Node(left, right){
		}

		template<class Stream>
		friend Stream& operator<<(Stream& os, const NodeWithCode& n){
			if(static_cast<unsigned int>(-1) == n.left){
				os << "(v=" << BASplitor() << n.code;
			}else{
				os << "(l=" << n.left << ", r= " << n.right;
			}
			return os << ")";
		}
	};

	struct TreeReader : public Tree<NodeWithCode>{
		typedef Tree<NodeWithCode> Super;
		unsigned int blockSize;
		TreeReader(unsigned int blockSize) : blockSize(blockSize){
		}

		template<class BitArray>
		unsigned int readValue(const BitArray& bitArray, unsigned int start, unsigned int& n){
			if(start + blockSize > bitArray.size()) return static_cast<unsigned int>(-1);
			auto end = start + blockSize;
			n = Super::nodes.size();
			Super::nodes.push_back(NodeWithCode());

			auto& code = Super::nodes.back().code;
			for(;start < end;++ start) code.push_back(bool(bitArray[start]));
			return end;
		}

		const std::vector<bool>& value(unsigned int n) const {
			return Super::nodes[n].code;
		}
	};

	struct TreeReaderShort : public Tree<Node>{
		typedef Tree<Node> Super;
		unsigned int blockSize;
		TreeReaderShort(unsigned int blockSize) : blockSize(blockSize){
			assert(blockSize <= 8 * sizeof(unsigned int));
		}

		template<class BitArray>
		unsigned int readValue(const BitArray& bitArray, unsigned int start, unsigned int& n){
			if(start + blockSize > bitArray.size()) return static_cast<unsigned int>(-1);
			auto end = start + blockSize;
			n = Super::nodes.size();
			Super::nodes.push_back(Node(0));

			auto& code = Super::nodes.back().right;
			unsigned int p2 = 1;
			for(;start < end;++ start){
				if(bool(bitArray[start]))
					code |= p2;
				p2 <<= 1;
			}
			return end;
		}

		struct Bits{
			struct BitsIterator{
				unsigned int value, current;
				BitsIterator(unsigned int value, unsigned int current) : value(value), current(current){
				}
				bool operator *(){
					return ((1 << current) & value) != 0;
				}
				bool operator ==(const BitsIterator bi) const {
					return current == bi.current;
				}
				bool operator !=(const BitsIterator bi) const {
					return current != bi.current;
				}
				BitsIterator& operator ++() {
					++current;
					return *this;
				}
				BitsIterator operator ++(int) {
					BitsIterator bi = *this;
					++current;
					return bi;
				}
			};
			unsigned int value, blockSize;
			BitsIterator begin(){ return {value, 0}; }
			BitsIterator end(){ return {value, blockSize}; }
			Bits(unsigned int value, unsigned int blockSize) : value(value), blockSize(blockSize){
			}
		};

		Bits value(unsigned int n) const {
			return Bits(Super::nodes[n].right, blockSize);
		}
	};

	/*
	 * Number of bits in input should be divisible by blockSize.
	 * */
	template<class BitArrayInput, class BitArrayOutput>
	static bool compress(const BitArrayInput& input, BitArrayOutput& output, unsigned int blockSize = 8){
		assert(input.size() % blockSize == 0);
		if(input.size() == 0 || blockSize == 0 || input.size() % blockSize != 0) return false;
		auto func =
			[&input, &blockSize](unsigned int i, unsigned int j){
				i *= blockSize;
				j *= blockSize;
				for(unsigned int k = 0;k < blockSize;++ k){
					auto vi = input[i + k], vj = input[j + k];
					if(vi < vj) return true;
					if(vj < vi) return false;
				}
				return  false;
			};

		std::map<unsigned int, unsigned int, decltype(func) > frequency(func);
		auto count = input.size() / blockSize;
		std::vector<unsigned int> index(count);
		std::map<unsigned int, std::vector<bool> > pfCode;
		{
			typedef TreeWithFrequency<BitArrayInput, ToBit> TWF;
			TWF tree(input, blockSize);
			computeFrequency(count, frequency, index);
			//std::cout << "frequency : " << frequency.size() << std::endl;
			buildTree(frequency, tree, NodeComparator<TWF>(tree));
			saveTree(tree, output);
			if(!tree.hasLeft(tree.root()) && !tree.hasRight(tree.root())){
				for(unsigned int i = 0;i < count;++ i){
					output.push_back(false);
				}
				return true;
			}
			buildPFCode(tree, pfCode);
		}
		for(unsigned int i = 0;i < count;++ i){
			assert(pfCode[index[i]].begin() != pfCode[index[i]].end());
			for(auto bit : pfCode[index[i]]) output.push_back(bool(bit));
		}
		return true;
	}

	template<class Tree, class Node, class BitArrayInput, class BitArrayOutput>
	static unsigned int decode(const Tree& tree, Node node, const BitArrayInput& input, unsigned int start, BitArrayOutput& output){
		if (!tree.hasLeft(node) && !tree.hasRight(node)){
			for(auto bit : tree.value(node)) output.push_back(bit);
			return start;
		}
		if (input[start]){
			assert(tree.hasRight(node));
			return decode(tree, tree.right(node), input, start + 1, output);
		}
		assert(tree.hasLeft(node));
		return decode(tree, tree.left(node), input, start + 1, output);
	}

	template<class BitArrayInput, class BitArrayOutput>
	static bool decompress(const BitArrayInput& input, BitArrayOutput& output, unsigned int blockSize=8){
			TreeReaderShort tree(blockSize);
			unsigned int current = readTree(input, tree);
			if(tree.empty()) return false;

			if(!tree.hasLeft(tree.root()) && !tree.hasRight(tree.root())){
				if(static_cast<unsigned int>(-1) == current) return false;
				while(current < input.size()){
					for(auto bit : tree.value(tree.root())) output.push_back(bit);
					++current;
				}
				return true;
			}
			while(static_cast<unsigned int>(-1) != current && current != input.size())
				current = decode(tree, tree.root(), input, current, output);
			return static_cast<unsigned int>(-1) != current;

	}
};

class LZW{
public:

	template<class BitArray>
	struct OffsetSegmentBitArray{
		unsigned int offset, segment, count;
		OffsetSegmentBitArray(unsigned int offset, unsigned int segment, unsigned int count)
		 	: offset(offset), segment(segment), count(count){
		}
		unsigned int operator[](unsigned int index) const {
			return offset + index * segment;
		}
		unsigned int size() const {return count;}
	};

	struct Codes{
		std::vector<bool> currentCode;
		std::vector<bool> codes;
		unsigned int codeSize;
		unsigned int current = 0;

		Codes(unsigned int codeSize) : codeSize(codeSize){
			currentCode.resize(codeSize);
			currentCode[0] = true;
		}

		bool hasNext(){
			while(current * codeSize >= codes.size()){
				unsigned int i = codeSize - 1;
				for(;i > 0 && currentCode[i];--i);
				if(i != 0){
					currentCode[i] = true;
					for(++i;i < codeSize;++ i) currentCode[i] = false;
					for(i = 0;i < codeSize;++ i) codes.push_back(currentCode[i]);
					continue;
				}
				return false;
			}
			return true;
		}

		unsigned int next(){
			if(!hasNext()) current = 0;
			return current++;
		}

		void reset(){
			current = 0;
		}

		template<class BitArray>
		void copy(unsigned int code, BitArray& bitArray){
			unsigned int start = code * codeSize;
			unsigned int end = start + codeSize;
			for(;start < end;++ start) bitArray.push_back(codes[start]);
		}

		template<class BitArray>
		void end(BitArray& bitArray){
			bitArray.push_back(true);
			for(unsigned int i = 1;i < codeSize;++ i) bitArray.push_back(false);
		}

		template<class BitArray>
		bool isCode(const BitArray& bitArray, unsigned int index){
			return bitArray[index];
		}

		template<class BitArray>
		bool isEnd(const BitArray& bitArray, unsigned int index){
			if(bitArray[index]){
				for(unsigned int i = 1;i < codeSize;++ i) if(bitArray[index + i]) return false;
				return true;
			}
			return false;
		}

	};

	template<class BitArrayInput, class BitArrayOutput>
	static bool compress(const BitArrayInput& input, 
		BitArrayOutput& output, unsigned int blockSize = 8, unsigned int codeSize = 9){
		if(input.size() % blockSize != 0 || codeSize <= blockSize) return false;

		typedef OffsetSegmentBitArray<BitArrayInput> OSBA;
		OSBA bitArray(0, blockSize, input.size() / blockSize);
		auto func =
			[&input, &blockSize](unsigned int i, unsigned int j){
				for(unsigned int k = 0;k < blockSize;++ k){
					auto vi = input[i + k], vj = input[j + k];
					if(vi < vj) return true;
					if(vj < vi) return false;
				}
				return  false;
			};
		TernarySearchTrie<unsigned int, unsigned int, decltype(func) > tst(func);
		Codes codes(codeSize);

		while(bitArray.size() > 0){
			unsigned int index = static_cast<unsigned int>(-1), value = 0;
			tst.longestPrefixOf(bitArray, value, index);
			if(static_cast<unsigned int>(-1) == index){
				for(unsigned int i = 0;i + blockSize < codeSize;++ i){
					output.push_back(false);
				}
				for(unsigned int i = 0;i < blockSize;++ i){
					output.push_back(input[bitArray.offset + i]);
				}
				index = 1;
			}else{
				assert(index > 1);
				codes.copy(value, output);
			}
			if (bitArray.size() >= index + 1){
				if(!codes.hasNext()){
					//std::cout << "restarting..." << std::endl;
					codes.end(output);
					codes.reset();
					tst.clear();
				}else{
					auto next = codes.next();
					//std::cout << "code for " << BASplitor(' ', blockSize, bitArray.offset, blockSize * (index + 1))
					//	<< input << " is " << BASplitor(' ', codeSize, next * codeSize, codeSize) << codes.codes << std::endl;
					tst.put(OSBA(bitArray.offset, blockSize, (index + 1)), next);
				}
			}
			bitArray.offset += index * blockSize;
			bitArray.count -= index;
		}
		//std::cout << "codes : " << BASplitor(' ', codeSize) << codes.codes << std::endl;
		return true;
	}

	template<class BitArrayInput, class BitArrayOutput>
	static bool decompress(const BitArrayInput& input, 
		BitArrayOutput& output, unsigned int blockSize = 8 , unsigned int codeSize = 9){
		if(input.size() % codeSize != 0 || codeSize <= blockSize) return false;
		typedef OffsetSegmentBitArray<BitArrayInput> OSBA;
		OSBA bitArray(0, blockSize, input.size() / blockSize);
		auto func =
			[&input, &codeSize](unsigned int i, unsigned int j){
				for(unsigned int k = 0;k < codeSize;++ k){
					auto vi = input[i + k], vj = input[j + k];
					if(vi < vj) return true;
					if(vj < vi) return false;
				}
				return  false;
			};

		auto func1 =
			[&output, &blockSize](unsigned int i, unsigned int j){
				for(unsigned int k = 0;k < blockSize;++ k){
					auto vi = output[i + k], vj = output[j + k];
					if(vi < vj) return true;
					if(vj < vi) return false;
				}
				return  false;
			};

		std::map<unsigned int, std::pair<unsigned int, unsigned int>, decltype(func) > code2Str(func);
		TernarySearchTrie<unsigned int, unsigned int, decltype(func1) > tst(func1);
		Codes codes(codeSize);

		unsigned int lastStart = 0;
		for(unsigned int i = 0;i < input.size();i += codeSize){
			unsigned int lastEnd = output.size();
			//std::cout << "current :" << BASplitor(' ', codeSize, i, codeSize) << input << std::endl;
			if(codes.isCode(input, i)){
				if(codes.isEnd(input, i)){
					tst.clear();
					code2Str.clear();
					codes.reset();
					continue;
				}
				auto pos = code2Str.find(i);
				if(pos == code2Str.end()){
					if(lastStart == lastEnd){
						// Wrong stream to be decompressed.
						return false;
					}
					for(unsigned int j = lastStart;j < lastEnd;++ j) output.push_back(output[j]);
					for(unsigned int j = 0;j < blockSize;++ j) output.push_back(output[lastStart + j]);
					//std::cout << "here... " << std::endl;
				}else{
					for(unsigned int j = pos->second.first;j < pos->second.second;++ j) output.push_back(output[j]);
				}
			}else{
				for(unsigned int j = i + codeSize - blockSize, k = 0;k < blockSize;++ k)
					output.push_back(input[j + k]);
			}
			if(codes.hasNext()){
				lastStart = lastEnd;
				//std::cout << "output :" << BASplitor(' ', blockSize) << output << std::endl;
				assert(lastEnd + blockSize <= output.size());
				tst.put(OSBA(lastStart, blockSize, (lastEnd - lastStart) / blockSize + 1), 
					codes.next());
				code2Str[i] = {lastStart, lastEnd + blockSize};
			}
		}
		return true;
	}
};

template<class String, class Less, class SizeType = size_t >
class SuffixArray{
	SizeType* sufficies = nullptr;
	String str = String();
	SizeType length = SizeType();
	Less less;
	friend class Comparator;
	struct Comparator{
		const SuffixArray *sa = nullptr;
		Comparator(const SuffixArray* sa) : sa(sa){}
		bool operator()(SizeType i, SizeType j){
			return sa->less(sa->str, i, sa->str, j);
		}
		bool operator()(SizeType i, const String& k){
			return sa->less(sa->str, i, k, SizeType());
		}
		bool operator()(const String& k, SizeType i){
			return sa->less(k, SizeType(), sa->str, i);
		}
	};
public:
	SuffixArray(const typename std::remove_cv<String>::type& str, SizeType length, 
		const Less& less) : str(str), length(length), less(less){
		if(SizeType() < length){
			sufficies = new SizeType[length];
			for(SizeType i = SizeType();i < length;++ i) sufficies[i] = i;
		}	
		std::sort(sufficies, sufficies + length, Comparator(this));
	}

	SuffixArray(const SuffixArray&) = delete;
	SuffixArray& operator=(const SuffixArray&) = delete;
	~SuffixArray(){delete[] sufficies;}

	//length of the string
	SizeType size() const {return length;}

	//index of the i-th smallest suffix string.
	SizeType index(SizeType i) const {
		return sufficies[i];
	}

	//length of longest common prefix of the i-th smallest and (i+1)-th smallest suffix string.
	SizeType lcp(SizeType index) const {
		SizeType len = SizeType();
		if(index < length){
			SizeType i = sufficies[index];
			++ index;
			if(index < length){
				SizeType j = sufficies[index];
				while(i < length && j < length && str[i] == str[j]){
					++len;++i;++j;
				}
			}
		}
		return len;
	}

	//number of suffixes less than key
	SizeType rank(const String& key) const {
		if(nullptr == sufficies) return SizeType();
		return std::lower_bound(sufficies, sufficies + length, key, Comparator(this)) - sufficies;
	}
};

}
#endif
