#ifndef MY_BINARY_SEARCH_TREE_H
#define MY_BINARY_SEARCH_TREE_H
#include "my_ordered_symbol_table.h"
#include <utility>
#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <memory>
#include <exception>

namespace my_lib{
template<class K, class V, class BSTNodeTraitsParam, class C = std::less<K> > 
class BinarySearchTreeBase : public OrderedSymbolTable<K, V>{
public:
	typedef BSTNodeTraitsParam BSTNodeTraits;
	typedef typename BSTNodeTraits::Node Node;
	typedef typename BSTNodeTraits::NodePtr NodePtr;
	BinarySearchTreeBase(C comparator = C()) : comparator(comparator){
	}
	BinarySearchTreeBase(BinarySearchTreeBase&& bst){
		root = bst.root;
		comparator = bst.comparator;
		bst.root = NodePtr();
	}
	BinarySearchTreeBase(const BinarySearchTreeBase& bst){
		root = (bst.root == NodePtr() ? NodePtr() : BSTNodeTraits::copyNode(bst.root));
		comparator = bst.comparator;
	}
	BinarySearchTreeBase& operator=(const BinarySearchTreeBase& bst){
		if(this == &bst)
			return *this;
		clear();
		root = (bst.root == NodePtr() ? NodePtr() : BSTNodeTraits::copyNode(bst.root));
		comparator = bst.comparator;
		return *this;
	}

	void put(const K& k, const V& v) override {
		putInner(k, v);
	}

	bool remove(const K&k) override {
		NodePtr node = removeInner(k);
		BSTNodeTraits::deleteNode(node);
		return node != NodePtr();
	}

	const V* get(const K& k) const override{
		NodePtr node = getNode(k);
		if(node == NodePtr())
			return nullptr;
		return BSTNodeTraits::value(node);
	}

	const K* min() const override{
		if(root == NodePtr())
			return nullptr; 
		NodePtr leftMost = BSTNodeTraits::lmc(root);
		return BSTNodeTraits::key(leftMost);
	}

	const K* max() const override{
		if(root == NodePtr())
			return nullptr;
		NodePtr rightMost = BSTNodeTraits::rmc(root);
		return BSTNodeTraits::key(rightMost);
	}

	const K* floor(const K& k) const override{
		NodePtr subtree = root;
		while(subtree != NodePtr()){
			if(comparator(*BSTNodeTraits::key(subtree), k)){
				if(BSTNodeTraits::right(subtree) == NodePtr())
					break;
				subtree = BSTNodeTraits::right(subtree);
			}else if(comparator(k, *BSTNodeTraits::key(subtree))){
				if(BSTNodeTraits::left(subtree) == NodePtr()){
					subtree = BSTNodeTraits::parent(BSTNodeTraits::rmp(subtree));
					break;
				}
				subtree = BSTNodeTraits::left(subtree);
			}else{
				break;
			}
		}
		if(subtree == NodePtr())
			return nullptr;
		return BSTNodeTraits::key(subtree);
	}

	const K* ceil(const K& k) const override{
		NodePtr subtree = root;
		while(subtree != NodePtr()){
			if(comparator(*BSTNodeTraits::key(subtree), k)){
				if(BSTNodeTraits::right(subtree) == NodePtr()){
					subtree = BSTNodeTraits::parent(BSTNodeTraits::lmp(subtree));
					break;
				}
				subtree = BSTNodeTraits::right(subtree);
			}else if(comparator(k, *BSTNodeTraits::key(subtree))){
				if(BSTNodeTraits::left(subtree) == NodePtr())
					break;
				subtree = BSTNodeTraits::left(subtree);
			}else{
				break;
			}
		}
		if(subtree == NodePtr())
			return nullptr; 
		return BSTNodeTraits::key(subtree);
	}

	unsigned int rank(const K& k) const override{
		NodePtr subtree = root;
		unsigned int count = 0;
		while(subtree != NodePtr()){
			if(comparator(*BSTNodeTraits::key(subtree), k)){
				count += BSTNodeTraits::safeCount(BSTNodeTraits::left(subtree)) + 1;
				if(BSTNodeTraits::right(subtree) == NodePtr()){
					break;
				}
				subtree = BSTNodeTraits::right(subtree);
			}else if(comparator(k, *BSTNodeTraits::key(subtree))){
				if(BSTNodeTraits::left(subtree) == NodePtr()){
					break;
				}
				subtree = BSTNodeTraits::left(subtree);
			}else{
				count += BSTNodeTraits::safeCount(BSTNodeTraits::left(subtree));
				break;
			}
		}
		return count;
	}

	const K* select(unsigned int index) const override{
		NodePtr n = selectInner(index);
		if(n == NodePtr())
			return nullptr;
		return BSTNodeTraits::key(n);			
	}

	bool leftRotate(NodePtr node){
		if(node == root){
			if(BSTNodeTraits::leftRotate(node)){
				root = BSTNodeTraits::parent(root);
				return true;
			}
			return false;
		}
		return BSTNodeTraits::leftRotate(node);
	}

	bool rightRotate(NodePtr node){
		if(node == root){
			if(BSTNodeTraits::rightRotate(node)){
				root = BSTNodeTraits::parent(root);
				return true;
			}
			return false;
		}
		return BSTNodeTraits::rightRotate(node);
	}
	bool removeMin() override {
		NodePtr node = removeMinInner();
		BSTNodeTraits::deleteNode(node);
		return node != NodePtr();
	}

	bool removeMax() override {
		NodePtr node = removeMaxInner();
		BSTNodeTraits::deleteNode(node);
		return node != NodePtr();
	}

	bool isEmpty() const override{
		return root == NodePtr();
	}

	unsigned int size() const override{
		return BSTNodeTraits::safeCount(root);
	}

	virtual int height() const{
		if(root == NodePtr())
			return -1;

		std::vector<NodePtr> parents;
		int height = 0;
		parents.push_back(root);

		while(!parents.empty()){
			std::vector<NodePtr> childs;
			for(size_t i = 0;i < parents.size();++ i){
				auto p = parents[i];
				auto r = BSTNodeTraits::right(p), l = BSTNodeTraits::left(p);
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

	unsigned int size(const K& l, const K& h) const override{
		unsigned int i = rank(l), j = rank(h);
		if(i > j)
			return 0;
		if(getNode(h) != NodePtr())
			return j - i + 1;
		return j - i;
	}

	void clear() override{
		if(root != NodePtr()){
			unsigned int count = BSTNodeTraits::count(root);
			std::vector<NodePtr> ps(count);
			ps[0] = root;
			for(unsigned int i = 0, j = 1;i < count;++ i){
				NodePtr n = ps[i];
				NodePtr l = BSTNodeTraits::left(n), r = BSTNodeTraits::right(n);
				if(l != NodePtr()){
					ps[j] = l;
					++ j;
				}	
				if(r != NodePtr()){
					ps[j] = r;
					++ j;
				}
				BSTNodeTraits::deleteNode(n);
			}	
			root = NodePtr();
		}
	}
	/**
	 * @return true if the data structure is valid.
	 */
	bool isValid() const override {
		assert(root == NodePtr() || BSTNodeTraits::parent(root) == NodePtr());
		return root == NodePtr() || 
			(BSTNodeTraits::parent(root) == NodePtr() && BSTNodeTraits::isValid(root));
	}

	const C& getComparator(){
		return comparator;
	}

	~BinarySearchTreeBase(){
		clear();
	}
private:
	typedef OrderedSymbolTable<K, V> Super;
	typedef typename Super::IteratorImpl IteratorImpl;
protected:
	NodePtr putInner(const K& k, const V& v){
		if(root == NodePtr()){
			root = BSTNodeTraits::createNode(k, v, NodePtr(), NodePtr(), NodePtr(), 1);
			return root;
		}
		NodePtr subtree = root;
		NodePtr newNode = NodePtr();
		//Loop invariant :
		//	subtree is the root of the sub tree, in which k lies.
		while(subtree != NodePtr()){
			if(comparator(*BSTNodeTraits::key(subtree), k)){
				if(BSTNodeTraits::right(subtree) == NodePtr()){
					newNode = BSTNodeTraits::createNode(k, v, subtree, NodePtr(), NodePtr(), 1);
					BSTNodeTraits::right(subtree, newNode);
					break;					
				}
				subtree = BSTNodeTraits::right(subtree);
			}else if(comparator(k, *BSTNodeTraits::key(subtree))){
				if(BSTNodeTraits::left(subtree) == NodePtr()){
					newNode = BSTNodeTraits::createNode(k, v, subtree, NodePtr(), NodePtr(), 1);
					BSTNodeTraits::left(subtree, newNode);
					break;
				}
				subtree = BSTNodeTraits::left(subtree);
			}else{
				BSTNodeTraits::value(subtree, v);
				return NodePtr();
			}
		}

		while(subtree != NodePtr()){
			BSTNodeTraits::count(subtree, BSTNodeTraits::count(subtree) + 1);
			subtree = BSTNodeTraits::parent(subtree);
		}
		return newNode;
	}

	NodePtr removeInner(const K& k){
		NodePtr node = getNode(k);
		if(node != NodePtr()){
			if(BSTNodeTraits::left(node) != NodePtr() && BSTNodeTraits::right(node) != NodePtr()){
				NodePtr successor = BSTNodeTraits::right(node);
				while(BSTNodeTraits::left(successor) != NodePtr())
					successor = BSTNodeTraits::left(successor);
				if(node == root)
					root = successor;
				if(BSTNodeTraits::parent(successor) != NodePtr()){
					if(successor == BSTNodeTraits::left(BSTNodeTraits::parent(successor)))
						BSTNodeTraits::left(BSTNodeTraits::parent(successor), node);
					else
						BSTNodeTraits::right(BSTNodeTraits::parent(successor), node);
				}
				if(BSTNodeTraits::parent(node) != NodePtr()){
					if(node == BSTNodeTraits::left(BSTNodeTraits::parent(node)))
						BSTNodeTraits::left(BSTNodeTraits::parent(node), successor);
					else
						BSTNodeTraits::right(BSTNodeTraits::parent(node), successor);
				}

				BSTNodeTraits::swap(node, successor);
				
				if(BSTNodeTraits::left(node) != NodePtr())
					BSTNodeTraits::parent(BSTNodeTraits::left(node), node);
				if(BSTNodeTraits::right(node) != NodePtr())
					BSTNodeTraits::parent(BSTNodeTraits::right(node), node);
				if(BSTNodeTraits::left(successor) != NodePtr())
					BSTNodeTraits::parent(BSTNodeTraits::left(successor), successor);
				if(BSTNodeTraits::right(successor) != NodePtr())
					BSTNodeTraits::parent(BSTNodeTraits::right(successor), successor);
			}
			//At least one of node's children is null
			NodePtr c = BSTNodeTraits::left(node);
			if(BSTNodeTraits::right(node) != NodePtr())
				c = BSTNodeTraits::right(node);
			if(node != root){
				NodePtr pNode = BSTNodeTraits::parent(node);
				(node == BSTNodeTraits::left(pNode)) ? BSTNodeTraits::left(pNode, c) : BSTNodeTraits::right(pNode, c);
				if(c != NodePtr())
					BSTNodeTraits::parent(c, pNode);
				
				while(pNode != NodePtr()){
					BSTNodeTraits::count(pNode, BSTNodeTraits::count(pNode) - 1);
					pNode = BSTNodeTraits::parent(pNode);
				}
			}else{
				root = c;
				if(root != NodePtr())
					BSTNodeTraits::parent(root, NodePtr());
			}
		}
		return node;
	}

	NodePtr selectInner(unsigned int index) const{
		if(root == NodePtr())
			return NodePtr();
		if(index >= BSTNodeTraits::count(root))
			return NodePtr();

		NodePtr n = root;
		unsigned int tmp = 0;
		if(BSTNodeTraits::left(n) != NodePtr())
			tmp = BSTNodeTraits::count(BSTNodeTraits::left(n));
		while(index != tmp){
			if(index > tmp){
				n = BSTNodeTraits::right(n);
				index -= (tmp + 1);
			}else{
				n = BSTNodeTraits::left(n);
			}
			tmp = BSTNodeTraits::safeCount(BSTNodeTraits::left(n));
		}
		return n;
	}

	NodePtr removeMinInner(){
		if(root != NodePtr()){
			NodePtr leftMost = root;
			while(BSTNodeTraits::left(leftMost) != NodePtr()){
				BSTNodeTraits::count(leftMost, BSTNodeTraits::count(leftMost) - 1);
				leftMost = BSTNodeTraits::left(leftMost);
			}
			if(leftMost != root){
				BSTNodeTraits::left(BSTNodeTraits::parent(leftMost), BSTNodeTraits::right(leftMost));
				if(BSTNodeTraits::right(leftMost) != NodePtr())
					BSTNodeTraits::parent(BSTNodeTraits::right(leftMost), BSTNodeTraits::parent(leftMost));
			}else{
				root = BSTNodeTraits::right(root);
				BSTNodeTraits::parent(root, NodePtr());
			}
			return leftMost;
		}
		return NodePtr();
	}

	NodePtr removeMaxInner(){
		if(root != NodePtr()){
			NodePtr rightMost = root;
			while(BSTNodeTraits::right(rightMost) != NodePtr()){
				BSTNodeTraits::count(rightMost, BSTNodeTraits::count(rightMost) - 1);
				rightMost = BSTNodeTraits::right(rightMost);
			}
			if(rightMost != root){
				BSTNodeTraits::right(BSTNodeTraits::parent(rightMost), BSTNodeTraits::left(rightMost));
				if(BSTNodeTraits::left(rightMost) != NodePtr())
					BSTNodeTraits::parent(BSTNodeTraits::left(rightMost), BSTNodeTraits::parent(rightMost));
			}else{
				root = BSTNodeTraits::left(root);
				BSTNodeTraits::parent(root, NodePtr());
			}
			return rightMost;
		}
		return NodePtr();
	}

	NodePtr getNode(const K& k) const{
		NodePtr subtree = root;
		//Loop invariant :
		//	subtree is the root of the sub tree, in which k may lies.
		while(subtree != NodePtr()){
			if(comparator(*BSTNodeTraits::key(subtree), k)){
				if(BSTNodeTraits::right(subtree) == NodePtr())
					return NodePtr();
				subtree = BSTNodeTraits::right(subtree);
			}else if(comparator(k, *BSTNodeTraits::key(subtree))){
				if(BSTNodeTraits::left(subtree) == NodePtr())
					return NodePtr();
				subtree = BSTNodeTraits::left(subtree);
			}else{
				return subtree;
			}
		}
		return NodePtr();
	}
	struct BinarySearchTreeIterator : public IteratorImpl{
		void  next() override{
			if(node != NodePtr())
				node = BSTNodeTraits::suc(node);
		}
		
		bool equals(const IteratorImpl& i) const override{
			const BinarySearchTreeIterator *itor = dynamic_cast<const BinarySearchTreeIterator*>(&i);
			if(itor != nullptr)
				return node == itor->node;
			return false;
		}

		bool operator==(const BinarySearchTreeIterator& itor){
			return node == itor.node;
		}

		bool operator!=(const BinarySearchTreeIterator& itor){
			return node != itor.node;
		}

		void assign(const IteratorImpl& i) override{
			const BinarySearchTreeIterator *itor = dynamic_cast<const BinarySearchTreeIterator*>(&i);
			if(itor != nullptr)
				this->node = itor->node;
		}

		const K& key() const override{
			return *BSTNodeTraits::key(node);
		}

		const V& value() const override{
			return *BSTNodeTraits::value(node);
		}

		BinarySearchTreeIterator* copy() const override{
			return new BinarySearchTreeIterator(node);
		}

		BinarySearchTreeIterator(NodePtr node) : node(node){
		}
	private:
		NodePtr node;
		
	};

	BinarySearchTreeIterator *implBegin() const override{
		return new BinarySearchTreeIterator(root == NodePtr() ? NodePtr() : BSTNodeTraits::lmc(root));
	}

	BinarySearchTreeIterator *implEnd() const override{
		return new BinarySearchTreeIterator(NodePtr());
	}	
protected:
	NodePtr root = NodePtr();
	C comparator;
};

template<class N, class NodeBase>
struct NodeTraits{
	typedef N Node;
	typedef Node* NodePtr;
	typedef const Node* ConstNodePtr;

	/**
	 * @return Left most child of {@var node}
	 */
	template<class Ptr>
	static Ptr lmc(Ptr leftMost){
		while(left(leftMost) != NodePtr())
			leftMost = left(leftMost);
		return leftMost;
	}
	/**
	 * @return Right most child of {@var node}
	 */
	template<class Ptr>
	static Ptr rmc(Ptr rightMost){
		while(right(rightMost) != NodePtr())
			rightMost = right(rightMost);
		return rightMost;
	}
	/**
	 * @return Left most parent of {@var node}
	 */
	template<class Ptr>
	static Ptr lmp(Ptr leftMost){
		Ptr p = parent(leftMost);
		while(p != NodePtr() && right(p) == leftMost){
			leftMost = p;
			p = parent(leftMost);
		}
		return leftMost;
	}
	/**
	 * @return Right most parent of this node
	 */
	template<class Ptr>
	static Ptr rmp(Ptr rightMost){
		Ptr p = parent(rightMost);
		while(p != NodePtr() && left(p) == rightMost){
			rightMost = p;
			p = parent(rightMost);
		}
		return rightMost;
	}
	/**
	 * @return Successor of {@var node}
	 */
	template<class Ptr>
	static Ptr suc(Ptr node){
		if(right(node) != NodePtr())
			return lmc(right(node));
		node = lmp(node);
		if(parent(node) != NodePtr() && node == left(parent(node)))
			return parent(node);
		return NodePtr();
	}
	/**
	 * @return Predecessor of this node
	 */
	template<class Ptr>
	static Ptr pre(Ptr node){
		if(left(node) != NodePtr())
			return left(node)->rmc();
		node = rmp(node);
		if(parent(node) != NodePtr() && node == right(parent(node)))
			return parent(node);
		return NodePtr();
	}
	/**
	 * @return left child of {@var node}
	 */
	template<class Ptr>
	static Ptr left(Ptr node){
		return NodeBase::left(node);
	}
	template<class Ptr>
	static void left(Ptr node, Ptr l){
		NodeBase::left(node, l);
	}
	/**
	 * @return right child of {@var node}
	 */
	template<class Ptr>
	static Ptr right(Ptr node){
		return NodeBase::right(node);
	}
	template<class Ptr>
	static void right(Ptr node, Ptr r){
		NodeBase::right(node, r);
	}
	/**
	 * @reurn parent of {@var node}
	 */
	template<class Ptr>
	static Ptr parent(Ptr node){
		return NodeBase::parent(node);
	}
	template<class Ptr>
	static void parent(Ptr node, Ptr p){
		NodeBase::parent(node, p);
	}	
	/**
	 *
	 * @return Number of nodes rooted at {@var node}
	 */
	template<class Ptr>
	static unsigned int count(Ptr node){
		return NodeBase::count(node);
	}
	template<class Ptr>
	static void count(Ptr node, unsigned int c){
		NodeBase::count(node, c);
	}
	template<class Ptr>
	static unsigned int safeCount(Ptr node){
		return node == NodePtr() ? 0 : count(node);
	}

	static void swap(NodePtr node1, NodePtr node2){
		NodeBase::swap(node1, node2);
	}

	static NodePtr copyNode(NodePtr node){
		auto parent = new Node(*node);
		auto left = NodeBase::left(parent);
		if(left != NodePtr()){
			left = copyNode(left);
			NodeBase::parent(left, parent);
			NodeBase::left(parent, left);
		}
		auto right = NodeBase::right(parent);
		if(right != NodePtr()){
			right = copyNode(right);
			NodeBase::parent(right, parent);
			NodeBase::right(parent, right);
		}
		return parent;
	}

	/**
	 * @return {@code true} if the data structure is valid.
	 */
	template<class Ptr>
	static bool isValid(const Ptr node){
		if(parent(node) == node){
			assert(false);
			return false;
		}
		unsigned int tmp = 0;
		if(left(node) != NodePtr()){
			if(parent(left(node)) != node || !isValid(left(node)))
				assert(false);//return false;
			tmp = count(left(node));
		}
		if(right(node) != NodePtr()){
			if(parent(right(node)) != node || !isValid(right(node)))
				assert(false);//return false;
			tmp += count(right(node));
		}
		assert(tmp + 1 == count(node));
		return tmp + 1 == count(node);
	}

	template<class Ptr>
	static bool leftRotate(Ptr node){
		NodePtr r = right(node);
		if(r != NodePtr()){
			NodePtr p = parent(node);
			NodePtr b = left(r);
			if(p != NodePtr()){
				if(node == left(p))
					left(p, r);
				else
					right(p, r);
			}
			parent(r, p);
			parent(node, r);
			left(r, node);
			right(node, b);
			if(b != NodePtr())
				parent(b, node);
			count(node, safeCount(b) + safeCount(left(node)) + 1);
			count(r, count(node)+ safeCount(right(r)) + 1);
			return true;
		}
		return false;
	}

	template<class Ptr>
	static bool rightRotate(Ptr node){
		NodePtr l = left(node);
		if(l != NodePtr()){
			NodePtr p = parent(node);
			NodePtr b = right(l);
			if(p != NodePtr()){
				if(node == left(p))
					left(p, l);
				else
					right(p, l);
			}
			parent(l, p);
			parent(node, l);
			right(l, node);
			left(node, b); 
			if(b != NodePtr())
				parent(b, node);
			count(node,  safeCount(b) + safeCount(right(node)) + 1);
			count(l, count(node) + safeCount(left(l)) + 1);
			return true;
		}
		return false;
	}
	
};
template<class K, class V, class BSTNode, class NodeBase>
struct BSTNodeTraits : public NodeTraits<BSTNode, NodeBase>{
private:
public:
	typedef typename NodeTraits<BSTNode, NodeBase>::Node Node;
	typedef typename NodeTraits<BSTNode, NodeBase>::NodePtr NodePtr;

	static NodePtr createNode(const K& k, const V& v){
		return new Node(k, v);		
	}
	static NodePtr createNode(const K& k, const V& v, NodePtr p, NodePtr l, NodePtr r){
		return new Node(k, v, p, l, r);
	}
	static NodePtr createNode(const K& k, const V& v, NodePtr p, NodePtr l, NodePtr r, unsigned int count){
		return new Node(k, v, p, l, r, count);
	}
	template<class Ptr>
	static const K* key(Ptr p){
		return NodeBase::template key<Ptr, K>(p);
	}
	template<class Ptr>
	static const V* value(Ptr node){
		return NodeBase::template value<Ptr, V>(node);
	}
	template<class Ptr>
	static void value(Ptr node, const V& v){
		NodeBase::value(node, v);
	}
	
	template<class Ptr>
	static void deleteNode(Ptr node){
		delete node;
	}
};
template<class Node>
struct BSTNodeBase{
	typedef Node* NodePtr;
	template<class Ptr>
	static Ptr left(Ptr node){
		return node->l;
	}
	static void left(NodePtr node, NodePtr l){
		node->l = l;
	}
	template<class Ptr>
	static Ptr right(Ptr node){
		return node->r;
	}
	static void right(NodePtr node, NodePtr r){
		node->r = r;
	}
	template<class Ptr>
	static Ptr parent(Ptr node){
		return node->p;
	}
	static void parent(NodePtr node, NodePtr p){
		node->p = p;
	}
	template<class Ptr>
	static unsigned int count(Ptr node){
		return node->c;
	}
	template<class Ptr, class K>
	static const K* key(Ptr p){
		return &p->value.first;
	}
	template<class Ptr, class V>
	static const V* value(Ptr node){
		return &node->value.second;
	}
	template<class Ptr, class V>
	static void value(Ptr node, const V& v){
		node->value.second = v;
	}
	static void count(NodePtr node, unsigned int count){
		node->c= count;
	}
	static void swap(NodePtr node1, NodePtr node2){
		std::swap(node1->p, node2->p);
		std::swap(node1->l, node2->l);
		std::swap(node1->r, node2->r);
		std::swap(node1->c, node2->c);
	}
};

template<class K, class V>
struct BSTNode{
	typedef BSTNode* NodePtr;
	typedef std::pair<K, V> Pair;
	NodePtr p, l, r;
	unsigned int c;
	Pair value;
	BSTNode(const K& k, const V& v) : BSTNode(k, v, NodePtr(), NodePtr(), NodePtr(), 1){
	}
	BSTNode(const K& k, const V& v, NodePtr p, NodePtr l, NodePtr r, unsigned int count)
		:p(p), l(l), r(r), c(count), value(Pair(k, v)){
	}
};
template<class K, class V, class C = std::less<K> >
using BinarySearchTree = BinarySearchTreeBase < K, V, BSTNodeTraits<K, V, BSTNode<K, V>, BSTNodeBase<BSTNode<K, V> > >, C > ;

template<class RBNode>
struct RBNodeBase : public BSTNodeBase<RBNode>{
	typedef typename BSTNodeBase<RBNode>::NodePtr NodePtr;
	template<class Ptr>
	static typename RBNode::Color color(Ptr node){
		return node == Ptr() ? RBNode::BLACK : node->color;
	}
	static void color(NodePtr node, typename RBNode::Color c){
		assert(node != NodePtr());
		node->color = c;
	}

	static void swap(NodePtr node1, NodePtr node2){
		BSTNodeBase<RBNode>::swap(node1, node2);
		std::swap(node1->color, node2->color);
	}
};

template<class K, class V, class RBNode, class RBNodeBase>
struct RBNodeTraits : public BSTNodeTraits<K, V, RBNode, RBNodeBase>{
	template<class Ptr>
	static typename RBNode::Color color(Ptr node){
		return RBNodeBase::color(node);
	}
	template<class Ptr>
	static void color(Ptr node, typename RBNode::Color c){
		RBNodeBase::color(node, c); 
	}

	static typename RBNode::Color black(){
		return RBNode::BLACK;
	}

	static typename RBNode::Color red(){
		return RBNode::RED;
	}
};

template<class K, class V>
struct RBNode{
	typedef RBNode* NodePtr;
	typedef std::pair<K, V> Pair;
	typedef bool Color;

	static const Color RED = true;
	static const Color BLACK = false;
	RBNode(const K& k, const V& v) : RBNode(k, v, NodePtr(), NodePtr(), NodePtr(), 1){
	}
	RBNode(const K& k, const V& v, NodePtr p, NodePtr l, NodePtr r, unsigned int count)
		:value(k, v), p(p), l(l), r(r), c(count), color(RED){
	}
private:
	friend class RBNodeBase<RBNode<K, V> >;
	friend class BSTNodeBase<RBNode<K, V> >;
	Pair value;
	NodePtr p, l, r;
	unsigned int c;
	bool color = RED;
};

template<class K, class V, class RBNodeTraitsParam, class C = std::less<K> >
class RBTreeBase : public BinarySearchTreeBase<K, V, RBNodeTraitsParam, C>{
	typedef BinarySearchTreeBase<K, V, RBNodeTraitsParam, C> Super;
public:
	typedef typename Super::BSTNodeTraits RBNodeTraits;
	typedef typename Super::NodePtr NodePtr;
	RBTreeBase(const C& comparator = C()) : Super(comparator){
	}
	void put(const K& k, const V& v) override {
		NodePtr node = Super::putInner(k, v);
		if(node != NodePtr())
			fixInsert(node);
	}

	bool remove(const K&k) override {
		NodePtr node = Super::removeInner(k);
		if(node != NodePtr()){
			assert(RBNodeTraits::left(node) == NodePtr() || RBNodeTraits::right(node) == NodePtr()); 
			fixRemove(node);
			RBNodeTraits::deleteNode(node);
			return true;
		}
		return false;
	}

	bool removeMin() override {
		NodePtr node = Super::removeMinInner();
		if(node != NodePtr()){
			assert(RBNodeTraits::left(node) == NodePtr() || RBNodeTraits::right(node) == NodePtr()); 
			fixRemove(node);
			RBNodeTraits::deleteNode(node);
			return true;
		}
		return false;
	}

	bool removeMax() override {
		NodePtr node = Super::removeMaxInner();
		if(node != NodePtr()){
			assert(RBNodeTraits::left(node) == NodePtr() || RBNodeTraits::right(node) == NodePtr()); 
			fixRemove(node);
			RBNodeTraits::deleteNode(node);
			return true;
		}
		return false;
	}

	bool isValid() const override {
		if(!Super::isValid())
			return false;
		if(Super::root != NodePtr()){
			if(RBNodeTraits::color(Super::root) != RBNodeTraits::black())
				//assert(false);
				return false;

			std::vector<NodePtr> nodes;
			std::vector<unsigned int> bhs;
			unsigned int bh = 0;
			bool isBhSet = false;
			nodes.push_back(Super::root);
			bhs.push_back(1);

			for(size_t i = 0;i < nodes.size();++ i){
				NodePtr node = nodes[i];
				NodePtr l = RBNodeTraits::left(node);
				NodePtr r = RBNodeTraits::right(node);
				if(RBNodeTraits::color(node) == RBNodeTraits::red() && (RBNodeTraits::color(r) != RBNodeTraits::black() || 
					RBNodeTraits::color(l) != RBNodeTraits::black()))
					//assert(false);
					return false;
				if(l != NodePtr()){
					bhs.push_back(RBNodeTraits::color(l) == RBNodeTraits::black() ? (bhs[i] + 1) : bhs[i]);
					nodes.push_back(l);
				}else if(!isBhSet){
					bh = bhs[i];
					isBhSet = true;
				}else if(bhs[i] != bh)
					//assert(false);
					return false;
				if(r != NodePtr()){
					bhs.push_back(RBNodeTraits::color(r) == RBNodeTraits::black() ? (bhs[i] + 1) : bhs[i]);
					nodes.push_back(r);
				}else if(!isBhSet){
					bh = bhs[i];
					isBhSet = true;
				}else if(bhs[i] != bh)
					//assert(false);
					return false;
				
			}
		}
		return true;
	}
protected:
	void fixInsert(NodePtr node){
		assert(RBNodeTraits::color(node) == RBNodeTraits::red());
		assert(RBNodeTraits::left(node) == NodePtr());
		assert(RBNodeTraits::right(node) == NodePtr());
		NodePtr p = RBNodeTraits::parent(node);
		while(RBNodeTraits::color(p) == RBNodeTraits::red()){
			NodePtr pp = RBNodeTraits::parent(p);
			NodePtr uncle = NodePtr();
			if(p == RBNodeTraits::left(pp)){
				uncle = RBNodeTraits::right(pp);
				if(RBNodeTraits::color(uncle) == RBNodeTraits::red()){
					RBNodeTraits::color(p, RBNodeTraits::black());
					RBNodeTraits::color(uncle, RBNodeTraits::black());
					RBNodeTraits::color(pp, RBNodeTraits::red());
					node = pp;
				}else{
					if(node == RBNodeTraits::right(p)){
						Super::leftRotate(p);
					}else{
						node = p;
					}
					RBNodeTraits::color(node, RBNodeTraits::black());
					RBNodeTraits::color(pp, RBNodeTraits::red());
					Super::rightRotate(pp);
					break;	
				}
			}else if(p == RBNodeTraits::right(pp)){
				assert(p == RBNodeTraits::right(pp));
				uncle = RBNodeTraits::left(pp);
				if(RBNodeTraits::color(uncle) == RBNodeTraits::red()){
					RBNodeTraits::color(p, RBNodeTraits::black());
					RBNodeTraits::color(uncle, RBNodeTraits::black());
					RBNodeTraits::color(pp, RBNodeTraits::red());
					node = pp;
				}else{
					if(node == RBNodeTraits::left(p)){
						Super::rightRotate(p);
					}else{
						node = p;
					}
					RBNodeTraits::color(node, RBNodeTraits::black());
					RBNodeTraits::color(pp, RBNodeTraits::red());
					Super::leftRotate(pp);
					break;
				}
			}else{
				assert("corrupted data structure");
			}
			p = RBNodeTraits::parent(node);

		}
		RBNodeTraits::color(Super::root, RBNodeTraits::black());
		//assert(isValid());
	}
	void fixRemove(NodePtr node){
		if(node != NodePtr() && RBNodeTraits::color(node) != RBNodeTraits::red()){
			NodePtr p = RBNodeTraits::parent(node);
			node = (RBNodeTraits::left(node) == NodePtr() ? RBNodeTraits::right(node): RBNodeTraits::left(node));
			while(node != Super::root && RBNodeTraits::color(node) == RBNodeTraits::black()){
				if(node == RBNodeTraits::left(p)){
					NodePtr sib = RBNodeTraits::right(p);
					if(RBNodeTraits::color(sib) == RBNodeTraits::red()){
						RBNodeTraits::color(sib, RBNodeTraits::black());
						RBNodeTraits::color(p, RBNodeTraits::red());
						sib = RBNodeTraits::left(sib);
						Super::leftRotate(p);
					}
					if(sib == NodePtr()){
						assert(RBNodeTraits::color(p) == RBNodeTraits::red());
						node = p;
						break;
					}
					if(RBNodeTraits::color(RBNodeTraits::left(sib)) == RBNodeTraits::red()){
						RBNodeTraits::color(RBNodeTraits::left(sib), RBNodeTraits::color(p));
						RBNodeTraits::color(p, RBNodeTraits::black());
						Super::rightRotate(sib);
						Super::leftRotate(p);
						break;
					}else if(RBNodeTraits::color(RBNodeTraits::right(sib)) == RBNodeTraits::red()){
						RBNodeTraits::color(RBNodeTraits::right(sib), RBNodeTraits::color(p));
						Super::leftRotate(p);
						break;
					}else{
						RBNodeTraits::color(sib, RBNodeTraits::red());
						node = p;
					}

				}else{
					NodePtr sib = RBNodeTraits::left(p);
					if(RBNodeTraits::color(sib) == RBNodeTraits::red()){
						RBNodeTraits::color(sib, RBNodeTraits::black());
						RBNodeTraits::color(p, RBNodeTraits::red());
						sib = RBNodeTraits::right(sib);
						Super::rightRotate(p);
					}
					if(sib == NodePtr()){
						assert(RBNodeTraits::color(p) == RBNodeTraits::red());
						node = p;
						break;
					}
					if(RBNodeTraits::color(RBNodeTraits::right(sib)) == RBNodeTraits::red()){
						RBNodeTraits::color(RBNodeTraits::right(sib), RBNodeTraits::color(p));
						RBNodeTraits::color(p, RBNodeTraits::black());
						Super::leftRotate(sib);
						Super::rightRotate(p);
						break;
					}else if(RBNodeTraits::color(RBNodeTraits::left(sib)) == RBNodeTraits::red()){
						RBNodeTraits::color(RBNodeTraits::left(sib), RBNodeTraits::color(p));
						Super::rightRotate(p);
						break;
					}else{
						RBNodeTraits::color(sib, RBNodeTraits::red());
						node = p;
					}
				}				
				p = RBNodeTraits::parent(node);
			}
			if(node != NodePtr()) 
				RBNodeTraits::color(node, RBNodeTraits::black());
			//assert(isValid());
		}
	}
};

template<class K, class V, class C = std::less<K> > 
using RBTree = RBTreeBase < K, V, RBNodeTraits<K, V, RBNode<K, V>, RBNodeBase<RBNode<K, V> > >, C > ;

template<class AVLNode>
struct AVLNodeBase : public BSTNodeBase<AVLNode>{
	typedef typename BSTNodeBase<AVLNode>::NodePtr NodePtr;
	template<class Ptr>
	static int height(Ptr node){
		if(node == Ptr())
			return -1;
		return node->h;
	}
	static void height(NodePtr node, int h){
		assert(node != NodePtr());
		node->h= h;
	}

	static void swap(NodePtr node1, NodePtr node2){
		BSTNodeBase<AVLNode>::swap(node1, node2);
		std::swap(node1->h, node2->h);
	}
};

template<class K, class V, class AVLNode, class AVLNodeBase>
struct AVLNodeTraits : public BSTNodeTraits<K, V, AVLNode, AVLNodeBase>{
	template<class Ptr>
	static int height(Ptr node){
		return AVLNodeBase::height(node);
	}
	template<class Ptr>
	static void height(Ptr node, int h){
		AVLNodeBase::height(node, h); 
	}

};

template<class K, class V>
struct AVLNode{
	typedef AVLNode* NodePtr;
	typedef std::pair<K, V> Pair;

	AVLNode(const K& k, const V& v) : AVLNode(k, v, NodePtr(), NodePtr(), NodePtr(), 1){
	}
	AVLNode(const K& k, const V& v, NodePtr p, NodePtr l, NodePtr r, unsigned int count)
		:value(k, v), p(p), l(l), r(r), c(count), h(0){
	}
private:
	friend class AVLNodeBase<AVLNode<K, V> >;
	friend class BSTNodeBase<AVLNode<K, V> >;
	Pair value;
	NodePtr p, l, r;
	unsigned int c;
	int h;
};

template<class K, class V, class AVLNodeTraitsParam, class C = std::less<K> >
class AVLTreeBase : public BinarySearchTreeBase<K, V, AVLNodeTraitsParam, C>{
	typedef BinarySearchTreeBase<K, V, AVLNodeTraitsParam, C> Super;
public:
	typedef typename Super::BSTNodeTraits AVLNodeTraits;
	typedef typename Super::NodePtr NodePtr;
	AVLTreeBase(const C& comparator = C()) : Super(comparator){
	}
	void put(const K& k, const V& v) override {
		NodePtr node = Super::putInner(k, v);
		if(node != NodePtr())
			fixInsert(node);
	}

	bool remove(const K&k) override {
		NodePtr node = Super::removeInner(k);
		if(node != NodePtr()){
			assert(AVLNodeTraits::left(node) == NodePtr() || AVLNodeTraits::right(node) == NodePtr()); 
			fixRemove(node);
			AVLNodeTraits::deleteNode(node);
			return true;
		}
		return false;
	}

	bool removeMin() override {
		NodePtr node = Super::removeMinInner();
		if(node != NodePtr()){
			assert(AVLNodeTraits::left(node) == NodePtr() || AVLNodeTraits::right(node) == NodePtr()); 
			fixRemove(node);
			AVLNodeTraits::deleteNode(node);
			return true;
		}
		return false;
	}

	bool removeMax() override {
		NodePtr node = Super::removeMaxInner();
		if(node != NodePtr()){
			assert(AVLNodeTraits::left(node) == NodePtr() || AVLNodeTraits::right(node) == NodePtr()); 
			fixRemove(node);
			AVLNodeTraits::deleteNode(node);
			return true;
		}
		return false;
	}

	int height() const override {
		if(Super::root == NodePtr())
			return -1;
		return AVLNodeTraits::height(Super::root);
	}

	bool isValid() const override {
		if(!Super::isValid())
			return false;
		if(Super::root != NodePtr()){
			std::vector<NodePtr> nodes;
			nodes.push_back(Super::root);

			for(size_t i = 0;i < nodes.size();++ i){
				NodePtr node = nodes[i];
				NodePtr l = AVLNodeTraits::left(node);
				NodePtr r = AVLNodeTraits::right(node);
				if(l != NodePtr()){
					nodes.push_back(l);
				}
				if(r != NodePtr()){
					nodes.push_back(r);
				}
				auto lH = AVLNodeTraits::height(l);
				auto rH = AVLNodeTraits::height(r);
				if(AVLNodeTraits::height(node) != std::max(lH, rH) + 1){
					//std::cout << "key : " << *AVLNodeTraits::key(node) << std::endl;
					assert(false);
					return false;
				}
				if(lH > rH + 1 || rH > lH + 1){
					assert(false);
					return false;
				}
			}
		}
		return true;
	}
protected:
	void fixInsert(NodePtr node){
		assert(AVLNodeTraits::left(node) == NodePtr());
		assert(AVLNodeTraits::right(node) == NodePtr());
		assert(AVLNodeTraits::height(node) == 0);
		auto p = AVLNodeTraits::parent(node);
		auto h = AVLNodeTraits::height(node) + 1;
		if(p != NodePtr() && h > AVLNodeTraits::height(p)){
			AVLNodeTraits::height(p, h);
			node = p;
			p = AVLNodeTraits::parent(node);
			while(p != NodePtr()){
				if(node == AVLNodeTraits::left(p)){
					auto sib = AVLNodeTraits::right(p);
					auto deltaH = AVLNodeTraits::height(node) - AVLNodeTraits::height(sib);
					if(deltaH < -1 || deltaH > 1){
						assert(deltaH == 2 || deltaH == -2);
						if(AVLNodeTraits::height(AVLNodeTraits::left(node)) + 1 == AVLNodeTraits::height(node)){
							assert(AVLNodeTraits::height(AVLNodeTraits::right(node)) + 2 == AVLNodeTraits::height(node) ||
								AVLNodeTraits::height(AVLNodeTraits::right(node)) + 1 == AVLNodeTraits::height(node));
							Super::rightRotate(p);
							h = AVLNodeTraits::height(AVLNodeTraits::left(p)) + 1;
							if(h < AVLNodeTraits::height(p)){
								AVLNodeTraits::height(p, h);
								break;
							}
							assert(h == AVLNodeTraits::height(p));
							AVLNodeTraits::height(node, h + 1);
						}else{
							auto r = AVLNodeTraits::right(node);
							h = AVLNodeTraits::height(r);
							assert(AVLNodeTraits::height(AVLNodeTraits::left(node)) + 2 == AVLNodeTraits::height(node));
							assert(h + 1 == AVLNodeTraits::height(node));
							Super::leftRotate(node);
							Super::rightRotate(p);
							AVLNodeTraits::height(r, h + 1);
							AVLNodeTraits::height(node, h);
							AVLNodeTraits::height(p, h);
							break;
						
						}
					}else{
						h = AVLNodeTraits::height(node) + 1;
						if(h <= AVLNodeTraits::height(p)){
							assert(h == AVLNodeTraits::height(p));
							break;
						}
						AVLNodeTraits::height(p, h);
						node = p;
					}
					
				}else{
					assert(node == AVLNodeTraits::right(p) && "corrupted data structure" != nullptr);
					auto sib = AVLNodeTraits::left(p);
					auto deltaH = AVLNodeTraits::height(node) - AVLNodeTraits::height(sib);
					if(deltaH < -1 || deltaH > 1){
						assert(deltaH == 2 || deltaH == -2);
						if(AVLNodeTraits::height(AVLNodeTraits::right(node)) + 1 == AVLNodeTraits::height(node)){
							assert(AVLNodeTraits::height(AVLNodeTraits::left(node)) + 2 == AVLNodeTraits::height(node) ||
								AVLNodeTraits::height(AVLNodeTraits::left(node)) + 1 == AVLNodeTraits::height(node));
							Super::leftRotate(p);
							h = AVLNodeTraits::height(AVLNodeTraits::right(p)) + 1;
							if(h < AVLNodeTraits::height(p)){
								AVLNodeTraits::height(p, h);
								break;
							}
							AVLNodeTraits::height(node, h + 1);
						}else{
							auto l = AVLNodeTraits::left(node);
							h = AVLNodeTraits::height(l);
							assert(AVLNodeTraits::height(AVLNodeTraits::right(node)) + 2 == AVLNodeTraits::height(node));
							assert(h + 1 == AVLNodeTraits::height(node));
							Super::rightRotate(node);
							Super::leftRotate(p);
							AVLNodeTraits::height(l, h + 1);
							AVLNodeTraits::height(node, h);
							AVLNodeTraits::height(p, h);
							break;
						}
					}else{
						h = AVLNodeTraits::height(node) + 1;
						if(h <= AVLNodeTraits::height(p)){
							assert(h == AVLNodeTraits::height(p));
							break;
						}
						AVLNodeTraits::height(p, h);
						node = p;
					}
				}
				p = AVLNodeTraits::parent(node);
			}
		}
		//assert(isValid());
	}

	void fixRemove(NodePtr node){
		auto p = AVLNodeTraits::parent(node);
		if(p != NodePtr()){
			{
				NodePtr l = AVLNodeTraits::left(p);
				NodePtr r = AVLNodeTraits::right(p);
				if(AVLNodeTraits::height(l) > AVLNodeTraits::height(r) + 1){
					node = r;
				}else if(AVLNodeTraits::height(r) > AVLNodeTraits::height(l) + 1){
					node = l;
				}else{
					auto h = std::max(AVLNodeTraits::height(r), AVLNodeTraits::height(l)) + 1;
					if(h >= AVLNodeTraits::height(p)){
						assert(h == AVLNodeTraits::height(p));
						assert(isValid());
						return;
					}
					AVLNodeTraits::height(p, h);
					node = p;
					p = AVLNodeTraits::parent(node);
				}
			}
			while(p != NodePtr()){
				auto h = AVLNodeTraits::height(node);
				/*if(node != NodePtr())
					std::cout << "key to adjust : " << *AVLNodeTraits::key(node) << std::endl;
				std::cout << "parent key : " << *AVLNodeTraits::key(p) << std::endl;
				if(AVLNodeTraits::parent(p) != NodePtr()){
					std::cout << "parent of parent key : " << *AVLNodeTraits::key(AVLNodeTraits::parent(p)) << std::endl;
				}*/
				if(node == AVLNodeTraits::left(p)){
					auto sib = AVLNodeTraits::right(p);
					auto sibH = AVLNodeTraits::height(sib);
					auto deltaH = h - sibH;
					/*if(sib != NodePtr())
						std::cout << "right sib key : " << *AVLNodeTraits::key(sib) << std::endl;*/
					if(deltaH < -1 || deltaH > 1){
						assert(deltaH == 2 || deltaH == -2);
						/*if(AVLNodeTraits::right(sib) != NodePtr())
							std::cout << "right child of right sib key : " << *AVLNodeTraits::key(AVLNodeTraits::right(sib)) << std::endl;
						if(AVLNodeTraits::left(sib) != NodePtr())
							std::cout << "left child of right sib key : " << *AVLNodeTraits::key(AVLNodeTraits::left(sib)) << std::endl;*/
						if(AVLNodeTraits::height(AVLNodeTraits::right(sib)) + 1 == sibH){
							Super::leftRotate(p);

							if(h + 1 == AVLNodeTraits::height(AVLNodeTraits::right(p))){
								AVLNodeTraits::height(p, h + 2);
								AVLNodeTraits::height(sib, h + 3);
								break;
							}
							AVLNodeTraits::height(p, h + 1);
							AVLNodeTraits::height(sib, h + 2);
							node = sib;
						}else{
							Super::rightRotate(sib);
							Super::leftRotate(p);
							AVLNodeTraits::height(sib, h + 1);
							AVLNodeTraits::height(p, h + 1);
							node = AVLNodeTraits::parent(p);
							AVLNodeTraits::height(node, h + 2);
						}
					}else{
						h = std::max(h, sibH) + 1;
						if(h >= AVLNodeTraits::height(p)){
							assert(h == AVLNodeTraits::height(p));
							break;
						}
						AVLNodeTraits::height(p, h);
						node = p;
					}
					
				}else{
					assert(node == AVLNodeTraits::right(p));
					auto sib = AVLNodeTraits::left(p);
					auto sibH = AVLNodeTraits::height(sib);
					auto deltaH = h - sibH;
					/*if(sib != NodePtr())
						std::cout << "left sib key : " << *AVLNodeTraits::key(sib) << std::endl;*/
					if(deltaH < -1 || deltaH > 1){
						assert(deltaH == 2 || deltaH == -2);
						/*if(AVLNodeTraits::left(sib) != NodePtr())
							std::cout << "left child of left sib key : " << *AVLNodeTraits::key(AVLNodeTraits::left(sib)) << std::endl;
						if(AVLNodeTraits::right(sib) != NodePtr())
							std::cout << "right child of left sib key : " << *AVLNodeTraits::key(AVLNodeTraits::right(sib)) << std::endl;*/
						if(AVLNodeTraits::height(AVLNodeTraits::left(sib)) + 1 == sibH){
							Super::rightRotate(p);

							if(h + 1 == AVLNodeTraits::height(AVLNodeTraits::left(p))){
								AVLNodeTraits::height(p, h + 2);
								AVLNodeTraits::height(sib, h + 3);
								break;
							}
							AVLNodeTraits::height(p, h + 1);
							AVLNodeTraits::height(sib, h + 2);
							node = sib;
						}else{
							Super::leftRotate(sib);
							Super::rightRotate(p);
							AVLNodeTraits::height(sib, h + 1);
							AVLNodeTraits::height(p, h + 1);
							node = AVLNodeTraits::parent(p);
							AVLNodeTraits::height(node, h + 2);
						}
					}else{
						h = std::max(h, sibH) + 1;
						if(h >= AVLNodeTraits::height(p)){
							assert(h == AVLNodeTraits::height(p));
							break;
						}
						AVLNodeTraits::height(p, h);
						node = p;
					}
				}
				p = AVLNodeTraits::parent(node);
			}
		}
		assert(isValid());
	}

};

template<class K, class V, class C = std::less<K> > 
using AVLTree = AVLTreeBase < K, V, AVLNodeTraits<K, V, AVLNode<K, V>, AVLNodeBase<AVLNode<K, V> > >, C > ;

template<class NodePtr>
void output(NodePtr node, unsigned int depth, const char *ch = " "){
	for(unsigned int i = 0;i < depth;++ i) std::cout << ch;
	auto count = node->childCount();
	std::cout << node->kvpCount() << ':';
	for(unsigned int i = 0;i < count;++ i)
		std::cout << '(' << node->key(i) << ',' << node->value(i) << ')' << ch;
	std::cout << std::endl;
	if(!node->isLeaf())
		for(unsigned int i = 0;i <= count;++ i) output(node->child(i), depth + 1, ch);
}

template<class K, class V, class KVPArray = std::vector<std::pair<K, V> > >
class BTreeNode{
public:
	BTreeNode(){
	}
	BTreeNode(const BTreeNode&) = delete;
	BTreeNode& operator=(const BTreeNode&) = delete;
	
	size_t childCount() const {
		return items.size();
	}

	size_t itemsCount() const {
		return items.size();
	}

	size_t kvpCount() const {
		return count;
	}

	bool isLeaf() const {
		return nullptr == children;
	}

	template<class NodeTraits>
	const V* get(const K& k, NodeTraits traits) const {
		auto node = getNode(k, traits);
		return nullptr == node.first || static_cast<size_t>(-1) == node.second ? 
			nullptr : &node.first->items[node.second].second;
	}

	const BTreeNode* child(size_t index) const {
		return nullptr == children || index > items.size() ? nullptr : children[index];
	}

	const K& key(size_t index) const {
		return items[index].first;
	}

	const V& value(size_t index) const{
		return items[index]. second;
	}

	const K* min() const {
		if(items.size() == 0) return nullptr;
		if(nullptr == children) return &key(0);
		return children[0]->min();
	}

	template<class NodeTraits>
	void removeMin(NodeTraits traits) {
		if(count == 0) return;
		if(nullptr == children) items.erase(items.begin());
		else ensureNum(0, traits)->removeMin(traits);
		--count;
	}

	template<class NodeTraits>
	void removeMax(NodeTraits traits) {
		if(count == 0) return;
		if(nullptr == children) items.pop_back();
		else ensureNum(items.size(), traits)->removeMax(traits);
		--count;
	}

	const BTreeNode* minNode() const {
		if(nullptr == children) return this;
		return children[0]->minNode();
	}

	const K* max() const{
		if(items.size() == 0) return nullptr;
		if(nullptr == children) return &key(items.size() - 1);
		return children[items.size()]->max();
	}

	const BTreeNode* maxNode() const {
		if(nullptr == children) return this;
		return children[items.size()]->max();
	}

	template<class NodeTraits>
	static const K* floor(const BTreeNode* node, const K& k, NodeTraits traits) {
		auto pos = std::lower_bound(node->items.begin(), node->items.end(), k, traits);
		size_t index = pos - node->items.begin();
		if(pos != node->items.end() && !traits(k, *pos)) return &pos->first;
		if(nullptr != node->children){
		 	auto ret = floor(node->children[index], k, traits);
			if(nullptr != ret) return ret;
		}
		if(0 == index) return nullptr;
	 	return &node->key(index - 1);
	}
	
	template<class NodeTraits>
	static const K* ceil(const BTreeNode* node, const K& k, NodeTraits traits) {
		auto pos = std::lower_bound(node->items.begin(), node->items.end(), k, traits);
		size_t index = pos - node->items.begin();
		if(pos != node->items.end() && !traits(k, *pos)) return &pos->first;
		if(nullptr != node->children){
		 	auto ret = ceil(node->children[index], k, traits);
			if(nullptr != ret) return ret;
		}

		if(pos == node->items.end()) return nullptr;
		return &pos->first;
	}
	
	template<class NodeTraits>
	static unsigned int rank(const BTreeNode* node, const K& k, NodeTraits traits) {
		auto pos = std::lower_bound(node->items.begin(), node->items.end(), k, traits);
		size_t index = pos - node->items.begin();
		unsigned int sum = index;
		if(nullptr != node->children){
			for(size_t i = 0;i < index;++ i) sum += node->children[i]->count;
			sum += rank(node->children[index], k, traits);
		}
		return sum;
	}

	template<class NodeTraits>
	static const K* select(const BTreeNode* node, size_t index, NodeTraits traits) {
		if(index >= node->count) return nullptr;
		if(nullptr == node->children) return &node->key(index);
		for(size_t i = 0;i <= node->childCount();++ i){
			size_t sz = node->child(i)->count;
			//std::cout << "index is " << index << " sz is " << sz << std::endl;
			if(index == sz) return &node->key(i);
			if(index < sz) return select(node->child(i), index, traits);
			index -= (sz + 1);
		}
		return nullptr;
	}	

	~BTreeNode(){
		if(nullptr != children && count > 0){
			for(size_t i = 0;i <= items.size();++ i)
				delete children[i];
		}
		delete[] children;
	}

	template<class Stream>
	Stream& printID(Stream& stream) const {
		if(items.size() > 0) stream << items[0].first;
		return stream << "....";
	}

	template<class NodeTraits>
	bool isSizeValid(NodeTraits traits) const {
			size_t sz = childCount();
			if(sz < traits.num() - 1){
				printID(std::cout) << " : size less than " << traits.num() - 1 << std::endl;
				return false;
			}
			if(sz >= 2 * traits.num()){
				printID(std::cout) << " : size not less than " << 2 * traits.num() << std::endl;
			 	return false;
			}
			return true;
	}

	template<class NodeTraits>
	bool isValid(NodeTraits traits) const {
		for(size_t i = 1;i < items.size();++ i)
			if (!traits(items[i - 1], items[i])){
				printID(std::cout) << " : " << i << "th key not in order" << std::endl;
			 	return false;
			}
		size_t sum = items.size();
		if(nullptr != children){
			for(size_t i = 0;i <= items.size();++ i){
				if(!children[i]->isSizeValid(traits) || !children[i]->isValid(traits)) return false;
				sum += children[i]->count;
			}
			for(size_t i = 0;i < items.size();++ i){
				if(!traits(*children[i]->max(), items[i])){
					printID(std::cout) << " : " << i << "th key less than left child's max" << std::endl;
				 	return false;
				}
				if(!traits(items[i], *children[i + 1]->min())){
					printID(std::cout) << " : " << i << "th key not less than right child's min" << std::endl;
				 	return false;
				}
			}
		}

		if(sum != count){
			printID(std::cout) << "recorded k-v size is : " << count <<  "but actual size is : " << sum << std::endl;
			return false;
		}
		return true;
	}

	template<class NodeTraits>
	static BTreeNode* put(BTreeNode* root, const K& k, const V& v, NodeTraits traits){
		// Enforce that the number of children for a BTreeNode must be at least 3
		if(traits.num() <= 2) return nullptr;
		BTreeNode* newNode = root->put(k, v, traits);
		if(nullptr != newNode && newNode != root){
			auto ret = new BTreeNode;
			ret->children = new BTreeNode*[2 * traits.num()];
			ret->items.push_back(root->items.back());
			root->items.pop_back();
			ret->children[0] = root;
			ret->children[1] = newNode;
			ret->count = root->count + newNode->count + 1;
			root = ret;
		}
		return root;
	}

	template<class NodeTraits>
	static BTreeNode* copy(BTreeNode* node, NodeTraits traits){
		if(nullptr == node) return nullptr;
		return new BTreeNode(*node, traits);
	}

	template<class NodeTraits>
	static bool remove(BTreeNode*& node, const K& k, NodeTraits traits){
		if(removeInner(node, k, traits)){
			while(node->childCount() == 0){
				if(nullptr == node->children){
					node = nullptr;
					break;
				}

				BTreeNode* toDelete = node;
				node = node->children[0];

				//avoid child link being deleted.
				toDelete->count = 0;
				delete toDelete;
			}
			return true;
		}
		return false;
	}
private:
	/*
	 * Merge the {@param index}-th child and {@param index + 1}-th child 
	 * 	and the key at {@param index} into one node.
	 * Assuming that both the left and right child have {@value num - 1} keys.
	 * */
	template<class NodeTraits>
	BTreeNode* mergeAt(size_t index, NodeTraits traits){
		BTreeNode* lc = children[index];
		BTreeNode* rc = children[index + 1];
		assert(lc->items.size() == traits.num() - 1);
		assert(rc->items.size() == traits.num() - 1);

		// adjust children's link first
		for(size_t i = index + 1;i < items.size();++ i) children[i] = children[i + 1];
		if(nullptr != lc->children){
			assert(nullptr != rc->children);

			std::copy(rc->children, rc->children + rc->items.size() + 1, 
				lc->children + lc->items.size() + 1);
		}

		//move items after link is adjusted.
		lc->items.push_back(items[index]);
		items.erase(items.begin() + index);
		lc->items.insert(lc->items.end(), rc->items.begin(), rc->items.end());
		
		//adjust counts
		lc->count += (rc->count + 1);
		assert(lc->items.size() == 2 * traits.num() - 1);

		//avoiding child link being destructed
		rc->items.clear();
		rc->count = 0;
		delete rc;
		return lc;
	}

	/*
	 * Insert the {@var index}-th key to the first key of the {@param index + 1}-th child
	 * And move the last key of {@var index}-th child to the {@param index)-th key.
	 * */
	template<class NodeTraits>
	BTreeNode* l2r(size_t index, NodeTraits traits){
		BTreeNode* lc = children[index];
		BTreeNode* rc = children[index + 1];
		if(nullptr != rc->children){
			assert(nullptr != lc->children);
			for(size_t i = rc->items.size() + 1;i > 0;-- i) rc->children[i] = rc->children[i - 1];
			rc->children[0] = lc->children[lc->items.size()];
		}

		rc->items.insert(rc->items.begin(), items[index]);
		items[index] = lc->items.back();
		lc->items.pop_back();

		if(nullptr != rc->children){
			assert(nullptr != lc->children);
			lc->count -= (rc->children[0]->count + 1);
			rc->count += (rc->children[0]->count + 1);
		}else{
			assert(nullptr == lc->children);
			-- lc->count;
			++ rc->count;
		}
		return lc;
	}

	/*
	 * Insert the {@var index}-th key to the first key of the {@param index + 1}-th child
	 * And move the last key of {@var index}-th child to the {@param index)-th key.
	 * */
	template<class NodeTraits>
	BTreeNode* r2l(size_t index, NodeTraits traits){
		BTreeNode* lc = children[index];
		BTreeNode* rc = children[index + 1];
		if(nullptr != lc->children){
			assert(nullptr != rc->children);
			lc->children[lc->items.size() + 1] = rc->children[0];
			for(size_t i = 0;i < rc->items.size();++ i) rc->children[i] = rc->children[i + 1];
		}

		lc->items.push_back(items[index]);
		items[index] = *rc->items.begin();
		rc->items.erase(rc->items.begin());

		if(nullptr != lc->children){
			assert(nullptr != rc->children);
			rc->count -= (lc->children[lc->items.size()]->count + 1);
			lc->count += (lc->children[lc->items.size()]->count + 1);
		}else{
			assert(nullptr == rc->children);
			-- rc->count;
			++ lc->count;
		}
		return lc;
	}

	template<class NodeTraits>
	BTreeNode* ensureNum(size_t index, NodeTraits traits){
		BTreeNode* lc = children[index];
		if(lc->items.size() >= traits.num()) return lc;
		if(index == items.size()){
				--index;
				if(children[index]->items.size() < traits.num()) return mergeAt(index, traits);
				l2r(index, traits);
		}else{
				if(children[index + 1]->items.size() < traits.num()) return mergeAt(index, traits);
				r2l(index, traits);
		}
		assert(lc->items.size() >= traits.num());
		return lc;
	}

	template<class NodeTraits>
	static bool removeInner(BTreeNode* node, const K& k, NodeTraits traits){
		auto pos = std::lower_bound(node->items.begin(), node->items.end(), k, traits);
		size_t index = pos - node->items.begin();
		if(pos != node->items.end() && !traits(k, *pos)){
			return node->removeAt(index, traits);
		}
		//output(node, 0);
		if(nullptr != node->children && removeInner(node->ensureNum(index, traits), k, traits)){
			//pintID(std::cout);
			--node->count;
			return true;
		}
		return false;
	}

	template<class NodeTraits>
	BTreeNode(const BTreeNode& node, NodeTraits traits){
		items = node.items;
		count = node.count;
		children = new BTreeNode*[traits.num() * 2];
		for(size_t i = 0;i <= items.size();++ i) children[i] = new BTreeNode(*node.children[i], traits);
	}

	template<class NodeTraits>
	bool removeAt(size_t index, NodeTraits traits){
		if(index >= items.size()) return false;
		--count;
		if(nullptr == children){
			items.erase(items.begin() + index);
		}else{
			auto lc = children[index], rc = children[index + 1];
			size_t ls = lc->items.size(), rs = rc->items.size();
			if(ls < traits.num() && rs < traits.num()){
				return mergeAt(index, traits)->removeAt(traits.num() - 1, traits);
			}else if(ls >= rs){
				items[index] = lc->maxNode()->items.back();
				lc->removeMax(traits);
			}else{
				items[index] = rc->minNode()->items[0];
				rc->removeMin(traits);
			}
		}
		return true;
	}

	BTreeNode* maxNode(){
		if(nullptr == children) return this;
		return children[items.size()]->maxNode();
	}

	BTreeNode* minNode(){
		if(nullptr == children) return this;
		return children[0]->minNode();
	}

	template<class NodeTraits>
	std::pair<const BTreeNode*, size_t> getNode(const K& k, NodeTraits traits) const {
		auto pos = std::lower_bound(items.begin(), items.end(), k, traits);
		size_t index = pos - items.begin();
		if(pos != items.end() && !traits(k, *pos)) return {this, index};
		if(nullptr == children) return {nullptr, static_cast<size_t>(-1)};
	 	return children[index]->getNode(k, traits);
	}

	template<class NodeTraits, class Iterator>
	inline BTreeNode* insertItem(Iterator pos, const K& k, const V& v, NodeTraits traits){
			if(items.size() >= 2 * traits.num() - 1){
				BTreeNode* ret = new BTreeNode;
				auto& retItems = ret->items;
				if(static_cast<size_t>(pos - items.begin()) >= traits.num()){
					retItems.insert(retItems.end(), items.begin() + traits.num(), pos);
					retItems.push_back({k, v});
					retItems.insert(retItems.end(), pos, items.end());
					items.resize(traits.num());
				}else{
					retItems.insert(retItems.end(), items.begin() + traits.num() - 1, items.end());
					items.resize(traits.num() - 1);
					items.insert(pos, {k, v});
				}
				this->count = traits.num() - 1;
				ret->count = traits.num();
				assert(this->items.size() == traits.num());
				assert(ret->items.size() == traits.num());
				return ret;
			}else{
				items.insert(pos, {k, v});
				++count;
			}
			return this;
	}

	template<class NodeTraits>
	BTreeNode* put(const K& k, const V& v, NodeTraits traits){
		//printID(std::cout) << " put: " << k << ", " << v << std::endl;
		auto pos = std::lower_bound(items.begin(), items.end(), k, traits);
		if(pos != items.end() && !traits(k, *pos)){
			//printID(std::cout) << " previous value is : " << pos->second << std::endl;
			pos->second = v;
		 	return nullptr;
		}
		if(nullptr != children){
			//printID(std::cout) << " internal node : " << children << std::endl;
			size_t index = pos - items.begin();
			auto child = children[index];
			auto newNode = child->put(k, v, traits);
			// child->printID(std::cout << "child ") << std::endl;
			if(nullptr == newNode) return nullptr;
			if(newNode == child){
				++ count;
			 	return this;
			}
			// newNode->printID(std::cout << "new ") << std::endl;

			BTreeNode* ret = insertItem(pos, child->items.back().first, 
				child->items.back().second, traits);
			child->items.pop_back();
			if(this != ret){
				ret->children = new BTreeNode*[2 * traits.num()];
				if(index == 2 * traits.num() - 1){
					ret->children[traits.num()] = newNode;
				}else{
					ret->children[traits.num()] = children[2 * traits.num() - 1];
					for(size_t i = 2 * traits.num() - 1;i > index + 1;-- i) children[i] = children[i - 1];
					children[index + 1] = newNode;
				}
				std::copy(children + traits.num(), children + 2 * traits.num(), ret->children);
				for(size_t i = 0;i < traits.num();++ i) this->count += this->children[i]->count;
				for(size_t i = 0;i <= traits.num();++ i) ret->count += ret->children[i]->count;
			}else{
				for(size_t i = items.size();i > index + 1;-- i) children[i] = children[i - 1];
				children[index + 1] = newNode;
			}
			return ret;
		}

		// printID(std::cout << "leaf node ") << std::endl;
		return insertItem(pos, k, v, traits);
	}

	KVPArray items;
	BTreeNode **children = nullptr;
	size_t count = 0;
};

template<class K, class V, class C = std::less<K> >
class BTreeBase : public OrderedSymbolTable<K, V>{
	typedef OrderedSymbolTable<K, V> Super;
public:
	BTreeBase(unsigned int num, const C& comparator = C()) : traits(num, comparator){
		if(num < 3) throw std::invalid_argument("degree must be greater than 2 ");
	}
	BTreeBase(const BTreeBase& tree){
		traits = tree.traits;
		root = Node::copy(tree.root, traits);
	}
	BTreeBase& operator=(const BTreeBase&) = delete;

	static const int NUM = 6;
	typedef BTreeNode<K, V> Node;
	typedef Node* NodePtr;
	
	void put(const K& k, const V& v) override {
		if (NodePtr() == root) root = new Node;
		root = Node::put(root, k, v, traits);
		//assert(root->isValid(traits));
	}

	const Node* getRoot() const {
		return root;
	}

	const V* get(const K& k) const override {
		return NodePtr() == root ? nullptr : root->get(k, traits);
	}

	bool remove(const K& k) override {
		if(NodePtr() == root) return false;
		return Node::remove(root, k, traits);
	}

	unsigned int size() const override {
		return NodePtr() == root ? 0 : root->kvpCount();
	}

	void clear() override {
		delete root;
		root = NodePtr();
	}

	const K* min() const override {
		return NodePtr() == root ? nullptr : root->min();
	}
	
	const K* max() const override {
		return NodePtr() == root ? nullptr : root->max();
	}

	const K* floor(const K& k) const override {
		return NodePtr() == root ? nullptr : Node::floor(root, k, traits);
	}
	
	const K* ceil(const K& k) const override {
		return NodePtr() == root ? nullptr : Node::ceil(root, k, traits);
	}
	
	unsigned int rank(const K& k) const override {
		return NodePtr() == root ? 0 : Node::rank(root, k, traits);
	}

	const K* select(unsigned int index) const override {
		return NodePtr() == root ? nullptr : Node::select(root, index, traits);
	}	

	unsigned int size(const K& l, const K& h) const override {
		unsigned int i = rank(l), j = rank(h);
		return i > j ? 0 : (nullptr != get(h) ? j - i + 1 : j - i);
	}	

	bool isValid() const override {
		return NodePtr() == root || root->isValid(traits);
	}

	bool removeMin() override {
		if(NodePtr() == root) return false;
		root->removeMin(traits);
		return true;
	}

	bool removeMax() override {
		if(NodePtr() == root) return false;
		root->removeMax(traits);
		return true;
	}

	~BTreeBase(){ delete root;}
protected:
	struct BTreeIterator : public Super::IteratorImpl{
		unsigned int index;
		const BTreeBase* tree;

		BTreeIterator(unsigned int index, const BTreeBase* tree) : index(index), tree(tree){}

		virtual void next() override{
			++index;
		}

		virtual bool equals(const typename Super::IteratorImpl& i) const override {
			const BTreeIterator* itor = dynamic_cast<const BTreeIterator*>(&i);
			return nullptr != itor && itor->index == index && itor->tree == tree;
		}

		virtual void assign(const typename Super::IteratorImpl& i) {
			const BTreeIterator* itor = dynamic_cast<const BTreeIterator*>(&i);
			if(nullptr != itor){
				index = itor->index;
				tree = itor->tree;
			}
		}

		virtual const K& key() const {
			return *tree->select(index);
		}

		virtual const V& value() const {
			return *tree->get(key());
		}

		virtual typename Super::IteratorImpl* copy() const {
			return new BTreeIterator(index, tree);
		}
	};
	typename Super::IteratorImpl *implBegin() const override{
		return new BTreeIterator(0, this);
	}

	typename Super::IteratorImpl *implEnd() const override{
		return new BTreeIterator(size(), this);
	}	
private:
	NodePtr root = NodePtr();
	struct NodeTraits{
		size_t n = 0;
		C comparator;
		NodeTraits (size_t n = NUM, const C& comparator = C()) : n(n), comparator(comparator){}
		size_t num(){return n;}
		bool operator()(const std::pair<K, V>& kvp1, const std::pair<K, V>& kvp2){
			return kvp1.first < kvp2.first;
		}
		bool operator()(const K& k, const std::pair<K, V>& kvp){
			return k < kvp.first;
		}
		bool operator()(const std::pair<K, V>& kvp, const K& k){
			return kvp.first < k;
		}
	}traits;
};
}
#endif //MY_BINARY_SEARCH_TREE_H
