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

template<class Item, class ItemTraits>
class Quick3WaySort{
static void sort(Item array[], unsigned int low, unsigned int high, unsigned int index){
	if (high <= low)
	       	return;

	auto lt = low, gt = high - 1;
	auto i = low + 1;
	auto v = ItemTraits::addr(array[low], index);
	if(v != nullptr){
		while (i <= gt){
			auto u = ItemTraits::addr(array[i], index);
			if(u == nullptr || *u < *v) ItemTraits::swap(array, lt++, i++);
			else if(*u > *v) ItemTraits::swap(array, i, gt--);
			else ++i;
		}
		sort(array, low, lt, index);
		sort(array, lt, gt + 1, index+1);
	}else{
		while (i <= gt)
			if(ItemTraits::addr(array[i], index) != nullptr) ItemTraits::swap(array, i, gt--);
			else i++;
	}
	sort(array, gt+1, high, index);
}
public: 
	static void sort(Item array[], unsigned int n){
		sort(array, 0, n, 0);
	}
};

using CStrQ3WS = Quick3WaySort<const char*, CStrTraits>;

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

template<class Char, class Value>
class TernarySearchTrie{
	struct Node{
		Char ch = Char();
		Node *l = nullptr, *r = nullptr, *m = nullptr;
		Value value = Value();
		bool hasValue = false;
		unsigned int size = 0;
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
				if(key[i] < node->ch){
					node = node->l;
				}else if(node->ch < key[i]){
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
public:
	TernarySearchTrie(){
	}
	TernarySearchTrie(const TernarySearchTrie&) = delete;
	TernarySearchTrie& operator=(const TernarySearchTrie&) = delete;
	TernarySearchTrie(TernarySearchTrie&& tst){
		root = tst.root;
		tst.root = nullptr;
	}
	~TernarySearchTrie(){
		if(root != nullptr) deleteNode(root);
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
				if(key[i] < node->ch){
					if(node->l == nullptr){
						node->l = new Node;
						node = node->l;
						node->ch = key[i];
						break;
					}
					node = node->l;
				}else if(node->ch < key[i]){
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
				if(key[i] < node->ch){
					node = node->l;
				}else if(node->ch < key[i]){
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
			if(key[i] < node->ch){
				link = &node->l;
				node = node->l;
			}else if(node->ch < key[i]){
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

	template<class Key, class Prefix>
	bool longestPrefixOf(const Key& s, Prefix& key){
		if(nullptr == root) return false;
		unsigned int index = -1;
		if(root->hasValue) index = 0;
		
		if(s.size() > 0){
			unsigned int i = 0;
			auto node = root->m;
			
			while(nullptr != node){
				if(s[i] < node->ch) node = node->l;
				else if(node->ch < s[i])	node = node->r;
				else{
					++ i;
					if(node->hasValue) index = i;
					if(i >= s.size()) break;
					node = node->m;
				}
			}
		}
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
		unsigned int row = m - h + 1, col = n - v + 1;
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
		char ch = 0;
		char type = 0;
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

}
#endif
