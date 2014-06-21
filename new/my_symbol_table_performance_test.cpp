#include "my_rand_string_generator.h"
#include "my_list_symbol_table.h"
#include "my_bin_search_symbol_table.h"
#include "my_access_counted_object.h"
#include "my_binary_search_tree.h"
#include "my_test_base.h"
#include <set>
#include <ctime>
#include <vector>
#include <iostream>
#include <iterator>
#include <vector>
#include <fstream>
#include <sstream>
#include <queue>
#include <cmath>
#include <cassert>

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
	KeyComparator(const Comparator = Comparator()) : comparator(comparator){
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
    
    void remove(const Key& key) override {
        root = deleet(root, key);
        assert(check());
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

    void removeMin() override {
        if (isEmpty()) return; 
        root = deleteMin(root);
        assert(check());
    }

    void removeMax() override {
        if (isEmpty()) return; 
        root = deleteMax(root);
        assert(check());
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
        for (int i = 0; i < size(); i++){
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
	unsigned int height(){
		vector<NodePtr> parents;
		unsigned int height = 0;
		if(Super::root == NodePtr())
			return 0;

		parents.push_back(Super::root);

		while(!parents.empty()){
			vector<NodePtr> childs;
			for(size_t i = 0;i < parents.size();++ i){
				NodePtr p = parents[i];
				NodePtr r = BSTNodeTraits::right(p), l = BSTNodeTraits::left(p);
				if(r != NodePtr())
					childs.push_back(r);
				if(l != NodePtr())
					childs.push_back(l);

			}
			parents.swap(childs);
			++ height;
		}
		return height - 1;
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
	typedef ForTestBST<unsigned int, unsigned int, RBTree> BST;
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
	UIntValues(unsigned int count, unsigned int *values, bool isRandom = true)
		:count(count), values(values){
		if(isRandom){
			for(unsigned int i = 0;i < count;++ i)
				values[i] = i;
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
}
int main(int argc, char *argv[]){
	//TestOne<DefaultBST> testPlot;
	//AverageCost<DefaultBST> testPlot;	
	/*TreePlot testPlot;
	Node *tree = randTree(100);
	testPlot.calculate(tree, NodeAdaptor());*/
	unsigned int count = 100;
	unsigned int nums[count];
	UIntValues vals(count, nums, false);
	//plot.print();
	BSTTreePlot<UIntValues> testPlot(vals);
	
	//AveTreeHeightPlot testPlot;
	testPlot.run(argc, argv);
	return 0;
}
