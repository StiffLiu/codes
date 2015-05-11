#include "my_rand_string_generator.h"
#include "my_list_symbol_table.h"
#include "my_bin_search_symbol_table.h"
#include "my_access_counted_object.h"
#include "my_binary_search_tree.h"
#include "my_hash_table.h"
#include "my_test_base.h"
#include <set>
#include <ctime>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <iterator>
#include <vector>
#include <fstream>
#include <sstream>
#include <queue>
#include <cmath>
#include <mutex>
#include <cassert>
#include <map>

namespace{
using namespace std;
using namespace my_lib;
#define ARSIZE(ar) (sizeof (ar) / sizeof *(ar))
class TestPlot : public my_lib::StatPlot<TestPlot&>{
	typedef my_lib::StatPlot<TestPlot&> Super;
	unsigned int start = 100;
	double totalTime1 = 0;
	double totalOp1 = 0;
	double totalTime2 = 0;
	double totalOp2 = 0;
public:
	TestPlot() : Super(4, *this){
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
			0.0, 0.0, 1.0,
			1.0, 1.0, 0.0,
		};
		colors.assign(c, c + ARSIZE(c));
	}	
	bool operator()(double *values){
		vector<string> strs;
		auto func = [&strs](size_t i1, size_t i2){
			return strs[i1] < strs[i2];
		};
		set<unsigned int, decltype(func)> strSet(func);
		RandStringGenerator gen(2, 50);
		
		while(strs.size() < start){
			strs.push_back(string());	
			gen.randStringAlphaNum(strs.back());
			if(!strSet.insert(strs.size() - 1).second)
				strs.pop_back();
			
		}

		clock_t begin = 0;
		clock_t end = 0;
		ListSymbolTable<string, bool> listST;
		BinSearchSymbolTable<string, bool> binSearchST;
		unsigned int half = start / 2; 
		double ops = 0;
		double duration = 0;

		begin = clock();
		for(unsigned int i = 0;i < half;++ i){
			listST.put(strs[i], false);	
		}
		for(auto i = 0;i < 10;++ i){
			for(unsigned int j = 0;j < half;++ j){
				listST.get(strs[j]);
				listST.get(strs[j + half]);
			}
		}
		end = clock();
		
		duration = end - begin;
		ops = 21 * half; 
		totalTime1 += duration;
		totalOp1 += ops;
		values[0] = start;
		values[1] = duration / ops;
		values[2] = start;
		values[3] = totalTime1 / totalOp1;

		begin = clock();
		for(unsigned int i = 0;i < half;++ i){
			binSearchST.put(strs[i], false);	
		}
		for(auto i = 0;i < 10;++ i){
			for(unsigned int j = 0;j < half;++ j){
				binSearchST.get(strs[j]);
				binSearchST.get(strs[j + half]);
			}
		}
		end = clock();
		
		duration = end - begin;
		ops = 21 * half; 
		totalTime2 += duration;
		totalOp2 += ops;
		values[4] = start;
		values[5] = duration / ops;
		values[6] = start;
		values[7] = totalTime2 / totalOp2;

		start += 2;
		return true;
	}

};

template<class Key, class Comparator>
struct KeyComparator{
	Comparator comparator;
	KeyComparator(const Comparator comparator = Comparator()) : comparator(comparator){
	}
	int operator()(const Key& k1, const Key& k2) const {
		if(comparator(k1, k2))
			return -1;
		if(comparator(k2, k1))
			return 1;
		return 0;	
	}
};

/**
 * symbol table implementation comes from the book "Algorithms 4 edition"
 * For experimental purpose
 */
template<class Key, class Value, class Comparator = KeyComparator<Key, less<Key> > >
class BST : public OrderedSymbolTable<Key, Value>{
   typedef OrderedSymbolTable<Key, Value> Super;
public:

    BST(const Comparator& comparator = Comparator()) : comparator(comparator){
    }

    // return number of key-value pairs in BST
    unsigned int size() const override {
        return size(root);
    }

    void put(const Key& key, const Value& val) override {
        root = put(root, key, val);
        assert(check());
    }

    const Value* get(const Key& key) const override {
        return get(root, key);
    }
    
    bool remove(const Key& key) override {
	auto ret = get(key);
        root = deleet(root, key);
        assert(check());
	return ret != nullptr;
    }

    bool isEmpty() const override {
	    return root == nullptr;
    }

    void clear() override {
	    clear(root);
    }

    const Key* min() const override {
        if (isEmpty()) return nullptr;
        return &min(root)->key;
    } 

    const Key* max() const override {
        if (isEmpty()) return nullptr;
        return &max(root)->key;
    } 

    const Key* floor(const Key& key) const override {
        const Node* x = floor(root, key);
        if (x == nullptr) return nullptr;
        else return &x->key;
    } 

    const Key* ceil(const Key& key) const override {
        const Node* x = ceiling(root, key);
        if (x == nullptr) return nullptr;
        else return &x->key;
    }
    
    unsigned int rank(const Key& key) const override {
        return rank(key, root);
    } 

    const Key* select(unsigned int k) const override{
        if (k < 0 || k >= size())  return nullptr;
        const Node* x = select(root, k);
        return &x->key;
    }

    bool removeMin() override {
        if (isEmpty()) return false; 
        root = deleteMin(root);
        assert(check());
	return true;
    }

    bool removeMax() override {
        if (isEmpty()) return false; 
        root = deleteMax(root);
        assert(check());
	return true;
    }

    unsigned int size(const Key& lo, const Key& hi) const override {
        if (comparator(hi, lo) > 0) return 0;
        if (Super::contains(hi)) return rank(hi) - rank(lo) + 1;
        else              return rank(hi) - rank(lo);
    }

    // height of this BST (one-node tree has height 0)
    int height() const { return height(root); }

    const Comparator& getComparator(){
	return comparator;
    }

    ~BST(){
	    clear();
    }

private:
    struct Node {
        Key key;           // sorted by key
        Value val;         // associated data
        Node *left = nullptr, *right = nullptr;  // left and right subtrees
        unsigned int N;             // number of nodes in subtree

        Node(Key key, Value val, unsigned int N) {
            this->key = key;
            this->val = val;
            this->N = N;
        }
    };

    Comparator comparator;
    Node *root = nullptr;             // root of BST

    void clear(Node *x){
	    if(x != nullptr){
		    clear(x->left);
		    clear(x->right);
	    }
	    delete x;
    }

    // return number of key-value pairs in BST rooted at x
    unsigned int size(Node* x) const {
        if (x == nullptr) return 0;
        else return x->N;
    }

    const Value* get(const Node* x, const Key& key) const {
        if (x == nullptr) return nullptr;
        int cmp = comparator(key, x->key);
        if      (cmp < 0) return get(x->left, key);
        else if (cmp > 0) return get(x->right, key);
        else              return &x->val;
    }


    Node* put(Node* x, const Key& key, const Value& val) {
        if (x == nullptr) return new Node(key, val, 1);
        int cmp = comparator(key, x->key);
        if      (cmp < 0) x->left  = put(x->left,  key, val);
        else if (cmp > 0) x->right = put(x->right, key, val);
        else              x->val   = val;
        x->N = 1 + size(x->left) + size(x->right);
        return x;
    }

    Node* deleteMin(Node* x) {
        if (x->left == nullptr){
		Node *node= x->right;
		delete x;
	      	return node; 
	}
        x->left = deleteMin(x->left);
        x->N = size(x->left) + size(x->right) + 1;
        return x;
    }

    Node* deleteMax(Node* x) {
        if (x->right == nullptr){
		Node *node = x->left;
		delete x;
	       	return node; 
	}
        x->right = deleteMax(x->right);
        x->N = size(x->left) + size(x->right) + 1;
        return x;
    }

    Node* deleet(Node* x, const Key& key) {
        if (x == nullptr) return nullptr;
        int cmp = comparator(key, x->key);
        if      (cmp < 0) x->left  = deleet(x->left,  key);
        else if (cmp > 0) x->right = deleet(x->right, key);
        else { 
            if (x->right == nullptr) return x->left;
            if (x->left  == nullptr) return x->right;
            Node *t = x;
            x = min(t->right);
            x->right = deleteMin(t->right);
            x->left = t->left;
        } 
        x->N = size(x->left) + size(x->right) + 1;
        return x;
    } 

    static Node* min(const Node* x) { 
        if (x->left == nullptr) return const_cast<Node*>(x); 
        else                return min(x->left); 
    } 

    static Node* max(const Node* x) { 
        if (x->right == nullptr) return const_cast<Node*>(x); 
        else                 return max(x->right); 
    } 

    const Node* floor(const Node* x, const Key& key) const {
        if (x == nullptr) return nullptr;
        int cmp = comparator(key, x->key);
        if (cmp == 0) return x;
        if (cmp <  0) return floor(x->left, key);
        const Node* t = floor(x->right, key); 
        if (t != nullptr) return t;
        else return x; 
    } 

    const Node* ceiling(const Node* x, const Key& key) const {
        if (x == nullptr) return nullptr;
        int cmp = comparator(key, x->key);
        if (cmp == 0) return x;
        if (cmp < 0) { 
            const Node* t = ceiling(x->left, key); 
            if (t != nullptr) return t;
            else return x; 
        } 
        return ceiling(x->right, key); 
    } 

    // Return key of rank k. 
    const Node* select(const Node* x, unsigned int k) const {
        if (x == nullptr) return nullptr; 
        unsigned int t = size(x->left); 
        if      (t > k) return select(x->left,  k); 
        else if (t < k) return select(x->right, k-t-1); 
        else            return x; 
    } 

    // Number of keys in the subtree less than key.
    unsigned int rank(const Key& key, const Node* x) const {
        if (x == nullptr) return 0; 
        int cmp = comparator(key, x->key); 
        if      (cmp < 0) return rank(key, x->left); 
        else if (cmp > 0) return 1 + size(x->left) + rank(key, x->right); 
        else              return size(x->left); 
    } 

    int height(const Node* x) const {
        if (x == nullptr) return -1;
	int lH = height(x->left) + 1;
	int rH = height(x->right) + 1;
	return (lH > rH ? lH : rH);
    }


    bool check() {
	return true;
        if (!isBST()){
    		cout << "Not in symmetric order" << endl;
		return false;
	}
        if (!isSizeConsistent()){
		cout << "Subtree counts not consistent" << endl;
		return false;
	}
        if (!isRankConsistent()){
		cout << "Ranks not consistent" << endl;
		return false;
	}
        return true; 
    }

    // does this binary tree satisfy symmetric order?
    // Note: this test also ensures that data structure is a binary tree since order is strict
    bool isBST() const {
        return isBST(root, nullptr, nullptr);
    }

    // is the tree rooted at x a BST with all keys strictly between min and max
    // (if min or max is null, treat as empty constraint)
    // Credit: Bob Dondero's elegant solution
    bool isBST(const Node* x, const Key* min, const Key* max) const {
        if (x == nullptr) return true;
        if (min != nullptr && comparator(x->key, *min) <= 0) return false;
        if (max != nullptr && comparator(x->key, *max) >= 0) return false;
        return isBST(x->left, min, &x->key) && isBST(x->right, &x->key, max);
    } 

    // are the size fields correct?
    bool isSizeConsistent() const { return isSizeConsistent(root); }
    bool isSizeConsistent(const Node* x) const {
        if (x == nullptr) return true;
        if (x->N != size(x->left) + size(x->right) + 1) return false;
        return isSizeConsistent(x->left) && isSizeConsistent(x->right);
    } 

    // check that ranks are consistent
    bool isRankConsistent() const {
        for (unsigned int i = 0; i < size(); i++){
	    const Key* key = select(i);
            if (key == nullptr || i != rank(*key)) return false;
	}
        //for (Key key : keys())
        //    if (key.compareTo(select(rank(key))) != 0) return false;
        return true;
    }

protected:
    typename Super::IteratorImpl* implBegin() const override {
	    return nullptr;
    }

    typename Super::IteratorImpl* implEnd() const override{
	    return nullptr;
    }
};

template<template<class U, class V> class T>
class TestOne  : public my_lib::StatPlot<TestOne<T>&> {
	typedef my_lib::StatPlot<TestOne<T>&> Super;
	unsigned int start = 100;
	unsigned int totalIteration = 10;
	unsigned int currentIteration = 0;
	double totalTime = 0;
	double totalOp = 0;
public:
	TestOne() : Super(2, *this){
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
		};
		Super::colors.assign(c, c + ARSIZE(c));
		Super::maxNumPoints = 0;
		Super::interval = 1;
	}
	bool operator()(double *values){
		vector<string> strs;
		auto func = [&strs](size_t i1, size_t i2){
			return strs[i1] < strs[i2];
		};
		set<unsigned int, decltype(func)> strSet(func);
		RandStringGenerator gen(2, 50);
		
		while(strs.size() < start){
			strs.push_back(string());	
			gen.randStringAlphaNum(strs.back());
			if(!strSet.insert(strs.size() - 1).second)
				strs.pop_back();
			
		}

		clock_t begin = 0;
		clock_t end = 0;
		T<string, bool> st;
		unsigned int half = start / 2; 
		double ops = 0;
		double duration = 0;

		begin = clock();
		for(unsigned int i = 0;i < half;++ i){
			st.put(strs[i], false);	
		}
		for(unsigned int j = 0;j < half;++ j){
			st.get(strs[j]);
			st.get(strs[j + half]);
		}
		end = clock();

		duration = end - begin;
		++currentIteration;
		if(currentIteration >= totalIteration){
			currentIteration = 0;
			totalTime = 0;
			totalOp = 0;
			++ start;
		}else{
			if(!Super::graphs[1].empty()){
				Super::graphs[1].pop_back();
				Super::graphs[1].pop_back();
			}
		}

		ops = 3 * half;
		values[0] = start;
		values[1] = duration / ops;
		totalTime += duration;
		totalOp += ops;
		values[2] = start;
		values[3] = totalTime / totalIteration / totalOp;
		return true;
	}
};
template<template<class U, class V> class T>
class AverageCost : public my_lib::StatPlot<AverageCost<T>&> {
	typedef my_lib::StatPlot<AverageCost<T>&> Super;
	unsigned int start = 100;
	unsigned int totalIteration = 10;
	unsigned int currentIteration = 0;
	double totalTime = 0;
	double totalOp = 0;
	size_t current = 0;
	std::vector<clock_t> times;
public:
	AverageCost() : Super(2, *this){
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
		};
		Super::colors.assign(c, c + ARSIZE(c));

		if(false){
			std::vector<std::string> vec;
			if(false){
				std::ifstream in("/development/documents/books/algos4/algs4.cs.princeton.edu/31elementary/tale.txt");
				std::copy(std::istream_iterator<std::string>(in), 
					std::istream_iterator<std::string>(),
					std::back_inserter(vec));
			}else{
				auto count = (1 << 16);
				RandStringGenerator gen(5, 10);
				cout << "generating " << count << " strings with length between 5 and 10 " << endl;
				for(auto i = 0;i < count;++ i){
					vec.push_back(string());
					gen.randStringPrintable(vec.back());
					//cout << vec.back() << endl;
				}
				sort(vec.begin(), vec.end());
				reverse(vec.begin(), vec.end());
			}

			T<std::string, bool> st;
			times.resize(vec.size());
			//clock_t start = clock();
			for(size_t i = 0;i < vec.size();++ i){
				st.put(vec[i], false);
				times[i] = st.getComparator().getCounter();//clock() - start;
			}
			std::cout << "total operations : " << vec.size() << std::endl;
			std::cout << "total keys : " << st.size() << std::endl;
			std::cout << "total time consumed : " << times.back() << std::endl;
			Super::maxNumPoints = vec.size();
		}else{
			std::vector<unsigned int> vec;
			auto count = (1 << 16);
			T<unsigned int, bool> st;
			for(auto i = 0;i < count;++ i)
				vec.push_back(i);
			//random_shuffle(vec.begin(), vec.end());
			reverse(vec.begin(), vec.end());
			times.resize(vec.size());
			//clock_t start = clock();
			for(size_t i = 0;i < vec.size();++ i){
				st.put(vec[i], false);
				times[i] = st.getComparator().getCounter();//clock() - start;
			}
			std::cout << "total operations : " << vec.size() << std::endl;
			std::cout << "total keys : " << st.size() << std::endl;
			std::cout << "total time consumed : " << times.back() << std::endl;
			Super::maxNumPoints = vec.size();

		}
		Super::interval = 1;
	}
	bool operator()(double *values){
		auto i = current;
		if(i >= times.size())
			return false;
		values[0] = values[2] = i;
		values[1] = (i == 0 ? times[i] : times[i] - times[i - 1]);
		values[3] = times[i] / (double)(i + 1);
		++ current;
		return true;
	}
};
	
template<class K, class V>
class DefaultListST : public ListSymbolTable<K, V, 
	std::list<std::pair<K, V> >, AccessCountedFunctor<unsigned long, std::less<K> > >{
};
template<class K, class V>
class DefaultBinST : public BinSearchSymbolTable<K, V>{
};
template<class K, class V>
class DefaultBST : public BinarySearchTree<K, V>{
};

template<class K, class V, template<class K, class V, class C> class BSTTree = BinarySearchTree, class C = less<K> >
class ForTestBST: public BSTTree<K, V, C>{
	typedef BSTTree<K, V, C> Super;
	typedef typename Super::Node Node;
public:
	typedef typename Super::BSTNodeTraits BSTNodeTraits;
	typedef typename Super::NodePtr NodePtr;
	NodePtr getRoot(){
		return Super::root;
	}
};

template<class T>
struct NodeTraitsAdaptor{
	typedef typename T::NodePtr NodePtr;
	bool left(NodePtr n, NodePtr& c) const{
		c = T::left(n);
		return c != NodePtr();
	}
	bool right(NodePtr n, NodePtr& c) const{
		c = T::right(n);
		return c != NodePtr();
	}
};
struct Node{
	Node *l = nullptr, *r = nullptr;
};

struct NodeAdaptor{
	bool left(Node *n, Node*& l){
		if(n->l == nullptr)
			return false;
		l = n->l;
		return true;
	}
	bool right(Node *n, Node*& r){
		if(n->r == nullptr)
			return false;
		r = n->r;
		return true;
	}
};

Node *randTree(unsigned int n){
	if(n == 0)
		return nullptr;
	if(n == 1)
		return new Node();
	n -= 1;
	unsigned int t = rand() % n;
	Node *node = new Node();
	node->l = randTree(t);
	node->r = randTree(n - t);
	return node;
}

template<class T>
class BSTTreePlot : public TreePlot{
	typedef ForTestBST<unsigned int, unsigned int, AVLTree> BST;
	typedef typename BST::NodePtr NodePtr;
	BST bst;
	mutable string str;
	vector<unsigned int> values;
public:
	BSTTreePlot(const T& t){
		for(unsigned int i = 0;i < t.count;++ i)
			bst.put(t.values[i], i);
		vector<NodePtr> nodes;
		if(bst.getRoot() != NodePtr()){
			nodes.push_back(bst.getRoot());
			for(unsigned int i = 0;i < bst.size();++ i){
				NodePtr c = BST::BSTNodeTraits::left(nodes[i]);
				values.push_back(*BST::BSTNodeTraits::key(nodes[i]));
				if(c != NodePtr())
					nodes.push_back(c);
				c = BST::BSTNodeTraits::right(nodes[i]);
				if(c != NodePtr())
					nodes.push_back(c);
			}
		}

		calculate(bst.getRoot(),
			NodeTraitsAdaptor<BST::BSTNodeTraits>());
		cout << "height : " << bst.height() << endl;
	}
protected:
	const char* getString(unsigned int index) const override{
		str.clear();

		ostringstream os;
		if(index < values.size()){
			os << values[index];
			str = os.str();
		}
		
		if(str.empty())
			return nullptr;
		return str.c_str();
	}
};

class AveTreeHeightPlot : public my_lib::StatPlot<AveTreeHeightPlot&>{
	typedef my_lib::StatPlot<AveTreeHeightPlot&> Super;
	typedef BST<double, double, AccessCountedFunctor<long, KeyComparator<double, less<double> > > > BSTtree;
	unsigned int start = 100;
public:
	AveTreeHeightPlot() : Super(3, *this){
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
			0.0, 0.0, 1.0,
		};
		colors.assign(c, c + ARSIZE(c));
	}	
	bool operator()(double *values){
		vector<unsigned int> nums;
		nums.resize(start);
		
		for(unsigned int i = 0;i < start;++ i)
			nums[i] = i;

		const unsigned int iteration = 100;
		double accumulatedHeight = 0;
		double accumulatedCounts = 0;
		for(unsigned int i = 0;i < iteration;++ i){
			BSTtree bst;
			random_shuffle(nums.begin(), nums.end());
			for(unsigned int j = 0;j < start;++ j)
				bst.put(nums[j], j);
			accumulatedHeight += bst.height();

			long countsInBuild = bst.getComparator().getCounter();
			for(unsigned int j = 0;j < 10;++ j)
				bst.contains(nums[rand() % start]);
			accumulatedCounts += bst.getComparator().getCounter() - countsInBuild;
		}
		values[4] = values[2] = values[0] = start;
		values[1] = accumulatedHeight / iteration;
		values[3] = 2 * log(start) - 2;
		values[5] = accumulatedCounts / (iteration * 10);

		start += 100;
		return true;
	}

};


struct TreeVec{
	unsigned int num = 0;
	vector<int> counts;
	friend ostream& operator<<(ostream& os, const TreeVec& tree){
		os << tree.num << ": ";
		copy(tree.counts.begin(), tree.counts.end(), ostream_iterator<int>(os, ", "));
		return os;
	}
	TreeVec(unsigned int num, vector<int>&& counts) : num(num), counts(counts){
	}
	TreeVec(unsigned int num, vector<int>& counts) : num(num), counts(counts){
	}
	TreeVec(){
	}
};

class AllTrees{
	vector<vector<TreeVec> > allTrees;	
	static unsigned int combination(unsigned int n, unsigned int m){
		static unsigned int total[31][31] = {{0}};
		if(n > 30 || m > n)
			return 0;
		if(total[0][0] == 0){
			for(unsigned int i = 0;i < 31;++ i)
				total[i][0] = 1;
			for(unsigned int i = 1;i < 31;++ i)
				for(unsigned int j = i + 1;j < 31;++ j)
					total[i][j] = 0;

			for(unsigned int i = 1;i < 31;++ i)
				for(unsigned int j = 1;j <= i;++ j)
					total[i][j] = total[i - 1][j] + total[i - 1][j - 1];
		}
		return total[n][m];
	}
public:
	AllTrees(unsigned int n){
		if(n >= 1){
			allTrees.resize(n + 1);
			allTrees[1].push_back(TreeVec(1, {1}));
			for(unsigned int i = 2;i <= n;++ i){

				vector<TreeVec>& trees = allTrees[i];
				for(unsigned int j = 1;j < i - 1;++ j){
					unsigned int k = i - j - 1;
					vector<TreeVec>& leftTrees = allTrees[j];
					vector<TreeVec>& rightTrees = allTrees[k];
					size_t leftNum = leftTrees.size();
					size_t rightNum = rightTrees.size();
					unsigned int count = combination(j + k, j); 

					for(size_t l = 0;l < leftNum;++ l)
						for(size_t r = 0;r < rightNum;++ r){
							trees.push_back(TreeVec());

							TreeVec& left = leftTrees[l];
							TreeVec& right = rightTrees[r];
							TreeVec& tree = trees.back();
							tree.counts.push_back(i);
							tree.counts.insert(tree.counts.end(), left.counts.begin(), left.counts.end());
							tree.counts.insert(tree.counts.end(), right.counts.begin(), right.counts.end());
							tree.num = left.num * right.num * count;	
						}
					
				}

				vector<TreeVec>& subTrees = allTrees[i - 1];
				size_t subNum = 2 * subTrees.size();
				size_t start = trees.size();

				trees.resize(start + subNum);
				for(size_t s = 0;s < subNum;s += 2){
					TreeVec& subTree = subTrees[s / 2];
					TreeVec& tree1 = trees[start + s];
					TreeVec& tree2 = trees[start + s + 1];

					tree1.counts.push_back(i);
					tree1.num = subTree.num;
					tree1.counts.insert(tree1.counts.end(), subTree.counts.begin(), subTree.counts.end());
					
					tree2.counts.push_back(i);
					tree2.num = subTree.num;
					tree2.counts.insert(tree2.counts.end(), subTree.counts.begin(), subTree.counts.end());
					tree2.counts[1] = -tree2.counts[1];
				}
			}

		}
	}
	const vector<vector<TreeVec> >& get() const{
		return allTrees;
	}	
	void print()const {	
		for(const auto& i : allTrees) for(const auto& j : i) cout << j << endl;
	}
};

class PlotAllTrees : public TreePlot{
	AllTrees allTrees{6};
	unsigned int current = 0;
	string str;
	struct GetNodeChildren{
		bool left(const int *r, const int *&l){
			int c = (*r < 0 ? -*r : *r);
			if(c <= 1)
				return false;
			if(*(r + 1) < 0)
				return false;
			l = r + 1;
			return true;
		}
		bool right(const int *r, const int *&right){
			int c = (*r < 0 ? -*r : *r);
			if(c <= 1)
				return false;
			int rc = *(r + 1);
			if(rc >= 0){
				if(rc == c - 1)
					return false;
				right = r + 1 + rc;
			}else{
				right = r + 1;	
			}
			return true;
		}
	};
	void populate(){
		if(current < allTrees.get().back().size()){
			auto& tree = allTrees.get().back()[current];
			auto& counts = tree.counts;	
			copy(counts.begin(), counts.end(), ostream_iterator<int>(cout, ", "));
			cout << endl;
			ostringstream os;
			os << tree;
			str = os.str();
			const int *root = &counts[0];
			calculate(root, GetNodeChildren());
		}

	}
protected:
	virtual const char *getString(unsigned int index) const{
		if(index == 0)
			return str.c_str();
		return nullptr; 
	}
public:
	PlotAllTrees(){
		populate();
	}
	void keyboard(unsigned char key, int x, int y) override {
		next();
		redisplay();
	}
	void next(){
		++current;
		populate();
	}
};

struct UIntValues{
	unsigned int count;
	unsigned int *values;
	enum Type{RANDOM, ORDERED, OPTIMAL};
	UIntValues(unsigned int count, unsigned int *values, Type type)
		:count(count), values(values){
		if(type != OPTIMAL){
			for(unsigned int i = 0;i < count;++ i)
				values[i] = i;
			if(type == RANDOM)
				std::random_shuffle(values, values + count);
		}else{
			std::queue<unsigned int> intervals;
			intervals.push(0);
			intervals.push(count);
			for(size_t i = 0;!intervals.empty();++ i){
				unsigned int b = intervals.front();
				intervals.pop();
				unsigned int e = intervals.front();
				intervals.pop();
				unsigned int m = (b + e) / 2;
				values[i] = m;
				if(m + 1 < e){
					intervals.push(m + 1);
					intervals.push(e);
				}
				if(b < m){
					intervals.push(b);
					intervals.push(m);
				}
			}
		}
	}

};

void optimalBST(double *p, double *q, unsigned int n, unsigned int *r){
	double *e = new double [n * n];
	double *w = new double [n * n];
	for(unsigned int i = 0;i < n;++ i){
		unsigned int index = i * n + i;
		double cost = q[i] + p[i] + q[i + 1];
		e[index] = cost;
		w[index] = cost;
		r[index] = i;
	}

	for(unsigned int l = 1;l < n;++ l){
		unsigned int i = 0;
		unsigned int j = l;
		while(j < n){
			double weight = q[i] + p[i] + w[(i + 1) * n + j];
			double cost = q[i] + e[(i + 1) * n + j] + weight; 
			double minCost = cost;
			unsigned int root = i;

			cost = e[i * n + j - 1] + q[j + 1] + weight;
			if(cost < minCost){
				minCost = cost;
				root = j;
			}

			for(unsigned int k = i + 1;k < j;++ k){
				cost = e[i * n + k - 1] + e[(k + 1) + j] + weight;
				if(cost < minCost){
					minCost = cost;
					root = k;
				}
			}

			r[i * n + j] = root;
			e[i * n + j] = minCost;
			w[i * n + j] = weight;
			++ i;
			j = i + l;
		}
	}
	delete[] e;
	delete[] w;
}
void printOptimalBST(ostream &os, unsigned int *r, unsigned int n){
	struct Print{
		static void func(ostream& os, unsigned int i, unsigned int j, unsigned int n, unsigned int *r){
			unsigned int root = r[i * n + j];
			if(root != i){
				os << "k" << r[i * n + (root - 1)] << " is the left child of k" << root << endl;
				func(os, i, root - 1, n, r);
			}else{
				os << "d" << root << " is the left child of k" << i << endl;
			}
			if(root != j){
				os << "k" << r[(root + 1) * n + j] << " is the right child of k" << root << endl;
				func(os, root + 1, j, n, r);
			}else{
				os << "d" << (root + 1) << " is the right child of k" << root << endl;
			}
		}
	};
	Print::func(os, 0, n - 1, n, r);
}
class HashChainSizePlot : public my_lib::StatPlot<HashChainSizePlot&> {
	typedef my_lib::StatPlot<HashChainSizePlot&> Super;
	my_lib::SeparateChainingHashTable<unsigned int, unsigned int> hashTable;
	std::unordered_map<unsigned int, unsigned int> unorderedMap;
	std::vector<unsigned int> keys;
	double ave1 = 0, var1 = 0, ave2 = 0, var2 = 0;
	unsigned int numKeysToDelete = 0;
	bool paused = false, nonDelete = false;
public:
	HashChainSizePlot() : Super(2, *this), hashTable(1.0){
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
		};
		Super::colors.assign(c, c + ARSIZE(c));
	}

	bool getTitle(RenderInfo& renderInfo) const override{
		Super::getTitle(renderInfo);
		std::ostringstream os;
		os << ", ave1 : " << ave1 << ", var1 : " << var1 << ", chain count : " << hashTable.getChainCount();
		os << ", ave2 : " << ave2 << ", var2 : " << var2 << ", bucket count : " << unorderedMap.bucket_count();

		renderInfo.title += os.str();
		renderInfo.charSize = 0.0001;
		return true;
	}

	void keyboard(unsigned char key, int x, int y) override{
		if(key == 'p')
			paused = !paused;
		else if(key == 'n')
			nonDelete = !nonDelete;
	}

	bool operator()(double *values){
		if(paused)
			return false;
		if(nonDelete || numKeysToDelete == 0 || keys.empty()){
			auto key = rand(), value = rand();
			if(hashTable.get(key) == nullptr)
				keys.push_back(key);
			hashTable.put(key, value);
			unorderedMap[key] = value;
			if(rand() % 500 == 1){
				numKeysToDelete = rand() % hashTable.size();
			}
		}else{
			size_t index = rand() % keys.size();
			swap(keys[index], keys[keys.size() - 1]);
			hashTable.remove(keys.back());
			unorderedMap.erase(keys.back());
			keys.pop_back();

			--numKeysToDelete;
		}

		unsigned int num1 = hashTable.getChainCount();
		unsigned int num2 = unorderedMap.bucket_count();

		if(num1 <= 0 || num2 <= 0)
			return false;
		
		maxNumPoints = std::max(num1, num2);
		
		std::vector<double>& graph0 = graphs[0];
		std::vector<double>& graph1 = graphs[1];
		graph0.clear();
		graph1.clear();

		ave1 = hashTable.size() / (double)num1;
		var1 = 0;
		for(unsigned int i = 0;i + 1 < num1;++ i){
			graph0.push_back(i);
			graph0.push_back(hashTable.getChainSize(i));

			double tmp = graph0.back() - ave1;
			var1 += tmp * tmp;
		}
		values[0] = num1 - 1;
		values[1] = hashTable.getChainSize(num1 - 1);
		var1 += (values[1] - ave1) * (values[1] - ave1);
		var1 /= num1;

		ave2 = unorderedMap.size() / (double)num2;
		var2 = 0;
		for(unsigned int i = 0;i + 1 < num2;++ i){
			graph1.push_back(i);
			graph1.push_back(unorderedMap.bucket_size(i));

			double tmp = graph1.back() - ave2;
			var2 += tmp * tmp;
		}
		values[2] = num2 - 1;
		values[3] = unorderedMap.bucket_size(num2 - 1);
		var2 += (values[3] - ave2) * (values[3] - ave2);
		var2 /= num2;
		return true;
	}

};
class EmptyChainNumPlot: public my_lib::StatPlot<EmptyChainNumPlot&> {
	typedef my_lib::StatPlot<EmptyChainNumPlot&> Super;
	my_lib::SeparateChainingHashTable<unsigned int, unsigned int> hashTable;
	std::unordered_map<unsigned int, unsigned int> unorderedMap;
	std::vector<unsigned int> keys;
	double ave1 = 0, var1 = 0, ave2 = 0, var2 = 0;
	unsigned int numKeysToDelete = 0;
	bool nonDelete = true;
public:
	EmptyChainNumPlot() : Super(4, *this), hashTable(1.0){
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
			1.0, 1.0, 0.0,
			0.0, 1.0, 1.0,
		};
		maxNumPoints  = 1000000;
		Super::colors.assign(c, c + ARSIZE(c));
	}

	bool operator()(double *values){
		if(nonDelete || numKeysToDelete == 0 || keys.empty()){
			auto key = rand(), value = rand();
			if(hashTable.get(key) == nullptr)
				keys.push_back(key);
			hashTable.put(key, value);
			unorderedMap[key] = value;
			if(rand() % 500 == 1){
				numKeysToDelete = rand() % hashTable.size();
			}
		}else{
			size_t index = rand() % keys.size();
			swap(keys[index], keys[keys.size() - 1]);
			hashTable.remove(keys.back());
			unorderedMap.erase(keys.back());
			keys.pop_back();

			--numKeysToDelete;
		}
		
		unsigned int count = hashTable.getChainCount();
		unsigned int emptyCount = 0;
		for(unsigned int i = 0;i < count;++ i)
			if(hashTable.getChainSize(i) == 0)
				++emptyCount;
		values[0] = hashTable.size() / (double)count;
		values[1] = emptyCount / (double)count;
		values[2] = values[0];
		values[3] = exp(-values[0]);

		count = unorderedMap.bucket_count();
		emptyCount = 0;
		for(unsigned int i = 0;i < count;++ i)
			if(unorderedMap.bucket_size(i) == 0)
				++emptyCount;
		values[4] = unorderedMap.size() / (double)count;
		values[5] = emptyCount / (double)count;
		values[6] = values[4];
		values[7] = exp(-values[4]);
		return true;
	}
};
class LinearProbingPlot : public my_lib::StatPlot<LinearProbingPlot&> {
	typedef my_lib::StatPlot<LinearProbingPlot&> Super;
	my_lib::LinearProbingHashTable<unsigned int, unsigned int> hashTable;
	std::vector<unsigned int> keys;
	mutable std::mutex m;
	unsigned int numKeysToDelete = 0;
	bool paused = false, nonDelete = false;
	unsigned int count = 0, maxLen = 0, aveHit = 0, aveMiss = 0;
public:
	LinearProbingPlot() : Super(3, *this), hashTable(0.5){
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
			1.0, 1.0, 0.0,
		};
		Super::colors.assign(c, c + ARSIZE(c));
	}	

	bool getTitle(RenderInfo& renderInfo) const override{
		std::unique_lock<std::mutex> lk(m);
		std::ostringstream os;

		auto& usage = hashTable.getSlotUsage();

		os << "size : " << hashTable.size() << ", slot size : " << usage.size() << ", count : "
		   << count << ", max len : " << maxLen << ", load factor : "
		   << hashTable.size() / (double)usage.size() << ", ave len : "
		   << hashTable.size() / (double)(count > 0 ? count : 1)
		   << ", ave hit : " << aveHit << ", ave Miss : " << aveMiss;

		renderInfo.charSize = 0.0001;
		renderInfo.title = os.str();
		
		return true;
	}

	void keyboard(unsigned char key, int x, int y) override{
		if(key == 'p')
			paused = !paused;
		else if(key == 'n')
			nonDelete = !nonDelete;
	}

	bool operator()(double *values){
		if(paused)
			return false;
		std::unique_lock<std::mutex> lk(m);
		count = 0;
		maxLen = 0;
		if(nonDelete || numKeysToDelete == 0 || keys.empty()){
			auto key = rand(), value = rand();
			if(hashTable.get(key) == nullptr)
				keys.push_back(key);
			hashTable.put(key, value);
			if(rand() % 500 == 1){
				numKeysToDelete = rand() % hashTable.size();
			}
		}else{
			size_t index = rand() % keys.size();
			swap(keys[index], keys[keys.size() - 1]);
			hashTable.remove(keys.back());
			keys.pop_back();

			--numKeysToDelete;
		}

		auto& usage = hashTable.getSlotUsage();
		unsigned int side = (unsigned int)sqrt(usage.size()) + 1;
		unsigned int i = 0;
		while(usage[i] && i < usage.size())
			++ i;
		if(i >= usage.size())
			return false;

		unsigned int k = i;
		std::vector<double>& graph1 = graphs[1];
		graph1.clear();

		for(unsigned int j = 0;j < side;++ j){
			graph1.push_back(j);
			graph1.push_back(0.0);
		}

		aveHit = aveMiss = 0;
		unsigned int maxLenIndex = 0;
		for(unsigned int j = 0;j < usage.size();++ j, i = (i + 1) % usage.size()){
			while(!usage[i] && j < usage.size()){
				++ j;
				i = (i + 1) % usage.size();
			}
			if(j < usage.size()){
				++count;
				unsigned int len = 0;
				unsigned int len1= 0;
				unsigned int oldIndex = i;
				while(usage[i] && j < usage.size()){
					//for our case the hash value is the unsigned integer
					//itself.
					if(*hashTable.getKeyAtSlot(i) % usage.size() == i){
						aveHit += len1 * (len1 + 1) / 2;
						len1 = 0;
					}
					++ j;
					i = (i + 1) % usage.size();
					++ len;
					++ len1;
				}
				if(len >= side)
					graph1[2 * side - 1] = graph1[2 * side - 1] + 1;
				else
					graph1[2 * len - 1] = graph1[2 * len - 1] + 1;
				aveMiss += len * (len + 1) / 2;
				aveHit += len1 * (len1 + 1) / 2;
				if(len > maxLen){
					maxLen = len;
					maxLenIndex = oldIndex;
				}
			}
		}

		double mxCount = 0;
		for(unsigned int j = 0;j < side;++ j){
			mxCount = std::max(mxCount, graph1[2 * j + 1]);
		}
		for(unsigned int j = 0;j < side;++ j){
			graph1[2 * j + 1] = graph1[2 * j + 1] / mxCount * side;
		}


		std::vector<double>& graph2 = graphs[0];
		std::vector<double>& graph3 = graphs[2];
		graph2.clear();
		graph3.clear();
		i = k;

		unsigned int maxLenEndIndex = (maxLenIndex + maxLen) % usage.size();
		for(unsigned int j = 0;j < usage.size();++ j, i = (i + 1) % usage.size())
			if(usage[i]){
				unsigned int x = j / side;
				unsigned int y = j % side;
				if(x % 2 == 1)
					y = side - 1 - y;
				if(maxLenEndIndex > maxLenIndex ? (i >= maxLenIndex && i < maxLenEndIndex) 
					: (i < maxLenEndIndex || i >= maxLenIndex)){
					graph3.push_back(x);
					graph3.push_back(y);
				}else{
					graph2.push_back(x);
					graph2.push_back(y);
				}
			}

		assert(!graph1.empty());
		assert(!graph3.empty());
		if(graph2.empty()){
			graph2.push_back(0);
			graph2.push_back(0);
			values[0] = values[1] = 0;
		}else{
			values[0] = graph2[graph2.size() - 2];
			values[1] = graph2.back();
			graph2.pop_back();
			graph2.pop_back();
		}
		values[2] = graph1[graph1.size() - 2];
		values[3]= graph1.back();
		graph1.pop_back();
		graph1.pop_back();
		values[4] = graph3[graph3.size() - 2];
		values[5] = graph3.back();
		graph3.pop_back();
		graph3.pop_back();

		return true;
	}
};
class HashFunctionCompare: public my_lib::StatPlot<HashFunctionCompare&>{
	typedef my_lib::StatPlot<HashFunctionCompare&> Super;
	static int hashCode1(const std::string& str){
		int hash = 0;
		size_t count = str.size();
		for(size_t i = 0;i < count;++ i)
			hash = hash * 31 + str[i];
		return hash;
	}
	static int hashCode2(const std::string& str){
		int hash = 0;
		size_t count = str.size();
		size_t skip = std::max((size_t)1, count / 8);
		for(size_t i = 0;i < count;i += skip)
			hash = hash * 37 + str[i];
		return hash;
	}
	my_lib::RandStringGenerator gen;
	std::map<int, int> counts1;
	std::map<int, int> counts2;
	unsigned int count = 0;
public:
	HashFunctionCompare() : Super(2, *this), gen(8, 50){
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
		};
		maxNumPoints  = 10000000;
		Super::colors.assign(c, c + ARSIZE(c));
	}
	void cal1(double* values){
		std::string str;
		gen.randStringPrintable(str);
		++count;
		++ counts1[hashCode1(str)];
		++ counts2[hashCode2(str)];

		auto begin = counts1.begin(), end = counts1.end();
		values[1] = 0;
		while(begin != end){
			if(begin->second > 1)
				values[1] += begin->second;
			++begin;
		}

		begin = counts2.begin();
		end = counts2.end();
		values[3] = 0;
		while(begin != end){
			if(begin->second > 1)
				values[3] += begin->second;
			++begin;
		}
		values[0] = values[2] = count;
	}
	void cal2(double* values){
		std::string str;
		gen.randStringPrintable(str);
		unsigned int h1 = (unsigned int)hashCode1(str);
		unsigned int h2 = (unsigned int)hashCode2(str);
		double test = 65536;
		unsigned int mask = 65535;
		values[0] = (h1 & mask) / test;
		values[1] = (h1 >> 16) / test;
		values[2] = (h2 & mask) / test;
		values[3] = (h2 >> 16) / test;
	}
	bool operator()(double *values){
		cal1(values);
		return true;
	}
};

class CompareWithStdHashMap: public my_lib::StatPlot<CompareWithStdHashMap&>{
	typedef my_lib::StatPlot<CompareWithStdHashMap&> Super;
	my_lib::SeparateChainingHashTable<unsigned int, unsigned int> hashTable{2.0};
	std::unordered_map<unsigned int, unsigned int> stdHashMap;
	std::random_device rd;
	std::uniform_int_distribution<unsigned int> intGenerator;
	double times[2] = {0.0, 0.0};
	unsigned int current = 0;
public:
	bool getTitle(RenderInfo& renderInfo) const override{
		Super::getTitle(renderInfo);
		std::ostringstream os;
		os << "my : (red, " << hashTable.getChainCount() << ", " << (double)hashTable.size() / hashTable.getChainCount() << "), std : (green, " 
		   << stdHashMap.bucket_count() << ", " << stdHashMap.load_factor() << ")";

		renderInfo.title = os.str();
		renderInfo.charSize = 0.0001;
		return true;
	}
	CompareWithStdHashMap() : Super(2, *this){
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
		};
		stdHashMap.max_load_factor(2.0);
		Super::colors.assign(c, c + ARSIZE(c));
	}
	bool operator()(double* values){
		unsigned int key = intGenerator(rd);
		clock_t start = 0;
		//avoid the influence of hardware cache.
		if(current % 2 == 0){
			start = clock();
			hashTable.put(key, key);
			times[0] += (clock() - start);
			start = clock();
			stdHashMap[key] = key;
			times[1] += (clock() - start);
		}else{
			start = clock();
			stdHashMap[key] = key;
			times[1] += (clock() - start);
			start = clock();
			hashTable.put(key, key);
			times[0] += (clock() - start);
		}
		hashTable.put(key, key);

		values[0] = values[2] = current;
		values[1] = times[0], values[3] = times[1];
		++ current;
		return true;
	}
};

class ChainLengthStatistic: public my_lib::StatPlot<ChainLengthStatistic&>{
	typedef my_lib::StatPlot<ChainLengthStatistic&> Super;
	vector<unsigned int> keys;
	std::random_device rd;
	std::uniform_int_distribution<unsigned int> intGenerator;
public:
	ChainLengthStatistic() : Super(1, *this){
		double c[]={
			1.0, 0.0, 0.0,
		};
		Super::colors.assign(c, c + ARSIZE(c));
	}
	bool operator()(double* values){
		keys.push_back(intGenerator(rd));

		vector<unsigned int> lengths;

		lengths.resize(keys.size() / 100);
		if(lengths.size() == 0)
			lengths.resize(1);
		for(size_t i = 0;i < keys.size();++ i)
			++ lengths[keys[i] % lengths.size()];

		unsigned int mn = keys.size();
		unsigned int mx = 0;
		for(size_t i = 0;i < lengths.size();++ i){
			mn = min(mn, lengths[i]);
			mx = max(mx, lengths[i]);
		}
		values[0] = values[2] = keys.size();
		values[1] = mn;
		values[3] = mx;

		return true;
	}
};


template<class K, class V, class C>
static unsigned int count(const SeparateChainingHashTable<K, V, C>& val){
	return val.getChainCount();
}
template<class K, class V, class C>
static unsigned int count(const LinearProbingHashTable<K, V, C>& val){
	return val.getSlotUsage().size();
}

typedef AccessCountedFunctor<Incrementor<unsigned long long>, std::equal_to<unsigned int> > EqualComparisonCounter;
typedef AccessCountedFunctor<Incrementor<unsigned long long>, std::less<unsigned int> > LessComparisonCounter;
template<class SymbolTable, int N>
class AccumulatedTimeComparison: public my_lib::StatPlot<AccumulatedTimeComparison<SymbolTable, N>&> {
	static_assert(N >= 1, "number of symbol tables must be at least 1");
protected:
	typedef my_lib::StatPlot<AccumulatedTimeComparison<SymbolTable, N>&> Super;
	SymbolTable* symbolTables[N];
	double times[N];
	unsigned long long counter[N];
	unsigned int current = 0;
	bool paused = false, isRand = false;
	int cmd = 'i';
	unsigned int maxValue = RAND_MAX;
	mutable std::mutex m;
public:
	AccumulatedTimeComparison() : Super(N, *this){
		keyboard('c', 0, 0);
	}
	void keyboard(unsigned char key, int x, int y) override{
		switch(key){
			case 'p': paused = !paused; break;
			case 's': isRand = !isRand; break;
			case 'c':
				for(auto& t : times) t = 0;
				for(auto& c : counter) c = 0;
				break;
			case 'g': case 'i': case 'd': 
				isRand = false;
				cmd = key;
				break;
			case 'm': cin >> maxValue; break;
			case 'r':{
				std::unique_lock<std::mutex> lk(m);
				current = 0;
				for(auto& ht : symbolTables) ht->clear();
				for(auto& g : Super::graphs) g.clear();
				for(auto& t : times) t = 0.0;
				for(auto& c : counter) c = 0;
				cmd = 'i';
				break;
			}
			 
		}
	}
	bool getTitle(typename Super::RenderInfo& renderInfo) const override{
		std::unique_lock<std::mutex> lk(m);
		renderInfo.charSize = 0.00008;
		renderInfo.titleY = 0.98;
		std::ostringstream os;
		for(unsigned int i = 0;i < N;++i){
			os << "(" << i << ", "
			   << symbolTables[i]->size() << ", " << counter[i] << ", "
			   << (current == 0 ? 0 : times[i] / current) << ") ";
		}
		renderInfo.title = os.str();
		return true;
	}
	bool operator()(double* values){
		if(paused)
			return false;
		std::unique_lock<std::mutex> lk(m);
		clock_t start = 0;
		vector<unsigned int> keys;
		int patchSize = 100;
		keys.resize(patchSize);
		for(auto& key : keys) key = rand() % maxValue;
		if(isRand){
			cmd = "gid"[rand()%3];
		}
		//The symbol table used later in this loop may use less time
		//then the hash tables used ealiear due to hardware memory cache.
		//This may make the comparision inaccurate.
		//So we access these tables forward and backward interchangably.
		for(int i = (current % 2 == 0 ? 0 : N - 1);i < N && i >= 0;i += (current % 2 == 0 ? 1 : -1)){
			switch(cmd){
			case 'i':
				start = clock();
				for(const auto& key : keys) symbolTables[i]->put(key, key);
				times[i] += (clock() - start);
				break;
			case 'g':
				start = clock();
				for(const auto& key : keys) symbolTables[i]->get(key);
				times[i] += (clock() - start);
				break;
			case 'd':
				start = clock();
				for(const auto& key : keys) symbolTables[i]->remove(key);
				times[i] += (clock() - start);
				break;
			}
			values[2 * i + 1] = times[i];
			values[2 * i] = current;
		}
		++current;
		return true;
	}
	

};

template<class HashTable>
class HashTableLoadFactorCompare : public AccumulatedTimeComparison<HashTable, 5>{
	typedef AccumulatedTimeComparison<HashTable, 5> Super;
	HashTable hashTables[5];
public:
	HashTableLoadFactorCompare(const double *lfs) 
		: hashTables{lfs[0], lfs[1], lfs[2], lfs[3], lfs[4]}{
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
			1.0, 1.0, 0.0,
			1.0, 0.0, 1.0,
			0.0, 1.0, 1.0,
		};
		for(unsigned int i = 0;i < 5;++ i){
			Super::symbolTables[i] = hashTables + i;
			hashTables[i].getComparator().setCounter(Super::counter + i);
		}
		Super::Super::colors.assign(c, c + ARSIZE(c));
	}
	bool getTitle(typename Super::Super::RenderInfo& renderInfo) const override{
		renderInfo.charSize = 0.00008;
		renderInfo.titleY = 0.98;
		std::ostringstream os;
		static const char* colorNames[] = {"red", "green", "yellow", "meganta", "cyan"};
		for(unsigned int i = 0;i < 5;++i){
			os << "(" << colorNames[i] << ", " << count(hashTables[i]) << ", " 
			   << hashTables[i].getLoadFactor() << ", " << hashTables[i].size() 
			   << ", " << Super::counter[i]<< ") ";
		}
		renderInfo.title = os.str();
		return true;
	}
};

class SymbolTableCompare: public AccumulatedTimeComparison<SymbolTable<unsigned int, unsigned int>, 5>{
	typedef AccumulatedTimeComparison<SymbolTable<unsigned int, unsigned int>, 5> Super;
	my_lib::LinearProbingHashTable<unsigned int, unsigned int, EqualComparisonCounter> lpHashTable{0.5};
	my_lib::SeparateChainingHashTable<unsigned int, unsigned int, EqualComparisonCounter> scHashTable{2.0};
	my_lib::RBTree<unsigned int, unsigned int, LessComparisonCounter> rbTree;
	my_lib::AVLTree<unsigned int, unsigned int, LessComparisonCounter> avlTree;
	my_lib::BinarySearchTree<unsigned int, unsigned int, LessComparisonCounter> bsTree;
public:
	SymbolTableCompare(){
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
			1.0, 1.0, 0.0,
			1.0, 0.0, 1.0,
			0.0, 1.0, 1.0,
		};
		Super::symbolTables[0] = &lpHashTable;
	        Super::symbolTables[1] = &scHashTable;
	        Super::symbolTables[2] = &rbTree;
	        Super::symbolTables[3] = &avlTree;
	        Super::symbolTables[4] = &bsTree;
		lpHashTable.getComparator().setCounter(Super::counter);
		scHashTable.getComparator().setCounter(Super::counter + 1);
		rbTree.getComparator().setCounter(Super::counter + 2);
		avlTree.getComparator().setCounter(Super::counter + 3);
		bsTree.getComparator().setCounter(Super::counter + 4);
		Super::Super::colors.assign(c, c + ARSIZE(c));
	}
	bool getTitle(typename Super::Super::RenderInfo& renderInfo) const override{
		renderInfo.charSize = 0.00008;
		renderInfo.titleY = 0.98;
		std::ostringstream os;
		static const char* colorNames[] = {"red", "green", "yellow", "meganta", "cyan"};
		static const char* stNames[] = {"lp", "sc", "rb", "avl", "bs"};
		for(unsigned int i = 0;i < 5;++i){
			os << "(" << stNames[i] << ", " << colorNames[i] << ", " 
			   << Super::symbolTables[i]->size() << ", " << Super::counter[i]<< ") ";
		}
		renderInfo.title = os.str();
		return true;
	}
};

struct AnotherHashFunction{
	unsigned int operator()(unsigned int key) const{
		return (unsigned int)(key * 31 + 2047);
	}
};

class DiffChainTypeCompare : public AccumulatedTimeComparison<SymbolTable<unsigned int, unsigned int>, 3>{
	typedef AccumulatedTimeComparison<SymbolTable<unsigned int, unsigned int>, 3> Super;
	my_lib::SeparateChainingHashTable<unsigned int, unsigned int, EqualComparisonCounter> vecChain{2.0};
	my_lib::SeparateChainingHashTable<unsigned int, unsigned int, LessComparisonCounter, my_lib::TypeConverter<unsigned int>, 
		my_lib::RBTree<unsigned int, unsigned int, LessComparisonCounter> > rbTreeChain{2.0};
	my_lib::DoubleHashSeparateChaining<unsigned int, unsigned int, AnotherHashFunction, EqualComparisonCounter> dhChain{2.0};
public:
	DiffChainTypeCompare(){
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
			1.0, 1.0, 0.0,
		};
	        Super::symbolTables[0] = &vecChain;
	        Super::symbolTables[1] = &rbTreeChain;
	        Super::symbolTables[2] = &dhChain;
		vecChain.getComparator().setCounter(Super::counter);
		rbTreeChain.getComparator().setCounter(Super::counter + 1);
		dhChain.getComparator().setCounter(Super::counter + 2);
		Super::Super::colors.assign(c, c + ARSIZE(c));
	}
	
};
class LinearProbingCompareEstimate: public my_lib::StatPlot<LinearProbingCompareEstimate&> {
	typedef my_lib::StatPlot<LinearProbingCompareEstimate&> Super;
	my_lib::LinearProbingHashTable<unsigned int, unsigned int, EqualComparisonCounter> lpHashTable{2.0};
	unsigned long long counter = 0;
	unsigned int current = 0;
public:
	LinearProbingCompareEstimate() : Super(2, *this){
		double c[]={
			1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
		};
		lpHashTable.getComparator().setCounter(&counter);
		Super::colors.assign(c, c + ARSIZE(c));
	}
	bool operator()(double *values){
		auto key = rand();
		static const double constant = sqrt(asin(1.0));
		lpHashTable.put(key, key);
		values[0] = values[2] = current;
		values[1] = counter;;
		values[3] = constant * pow(lpHashTable.size(), 1.5);
		++current;
		return true;
	}
};

int test(int argc, char *argv[]){
	//TestOne<DefaultBST> testPlot;
	//AverageCost<DefaultBST> testPlot;	
	/*TreePlot testPlot;
	Node *tree = randTree(100);
	testPlot.calculate(tree, NodeAdaptor());*/
	unsigned int count = 100;
	unsigned int nums[count];
	UIntValues vals(count, nums, UIntValues::RANDOM);
	//plot.print();
	BSTTreePlot<UIntValues> testPlot(vals);
	
	//HashChainSizePlot  testPlot;
	//LinearProbingPlot testPlot;
	//EmptyChainNumPlot testPlot;
	//AveTreeHeightPlot testPlot;
	//HashFunctionCompare testPlot;
	//SeperateChainingStatistic testPlot;
	//double lfs[5] = {0.1, 0.5, 0.9, 2.0, 4.0};
	//HashTableLoadFactorCompare<my_lib::SeparateChainingHashTable<unsigned int, unsigned int, EqualComparisonCounter> >  testPlot(lfs);
	//double lfs[5] = {0.1, 0.2, 0.5, 0.7, 0.9};
	//HashTableLoadFactorCompare<my_lib::LinearProbingHashTable<unsigned int, unsigned int, EqualComparisonCounter> >  testPlot(lfs);
	//SymbolTableCompare testPlot;
	//CompareWithStdHashMap testPlot;
	//DiffChainTypeCompare testPlot;
	//LinearProbingCompareEstimate testPlot;
	testPlot.run(argc, argv);
	return 0;
}

int test1(int argc, char *argv[]){
	const unsigned int n = 5;
	double p[n] = {0.15, 0.10, 0.05, 0.10, 0.20};
	double q[n + 1] = {0.05, 0.10, 0.05, 0.05, 0.05, 0.10};
	unsigned int r[n * n];
	optimalBST(p, q, n, r);
	printOptimalBST(cout, r, n);
	return 0;
}

}
int main(int argc, char *argv[]){
	return test(argc, argv);
}
