#ifndef MY_BINARY_SEARCH_TREE_H
#define MY_BINARY_SEARCH_TREE_H
#include "my_ordered_symbol_table.h"
#include <utility>
#include <vector>
#include <cassert>
#include <iostream>

namespace my_lib{
template<class K, class V, class BSTNodeTraitsParam, class C = std::less<K> > 
class BinarySearchTreeBase : public OrderedSymbolTable<K, V>{
public:
	typedef BSTNodeTraitsParam BSTNodeTraits;
	typedef typename BSTNodeTraits::Node Node;
	typedef typename BSTNodeTraits::NodePtr NodePtr;
	BinarySearchTreeBase(C comparator = C()) : comparator(comparator){
	}
	BinarySearchTreeBase(const BinarySearchTreeBase&) = delete;
	BinarySearchTreeBase& operator=(const BinarySearchTreeBase&) = delete;

	void put(const K& k, const V& v) override {
		putInner(k, v);
	}

	void remove(const K&k) override {
		NodePtr node = removeInner(k);
		BSTNodeTraits::deleteNode(node);
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
	void removeMin() override {
		BSTNodeTraits::deleteNode(removeMinInner());
	}

	void removeMax() override {
		BSTNodeTraits::deleteNode(removeMaxInner());
	}

	bool isEmpty() const override{
		return root == NodePtr();
	}

	unsigned int size() const override{
		return BSTNodeTraits::safeCount(root);
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
	virtual bool isValid() const{
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
				BSTNodeTraits::parent(root, NodePtr());
				root = BSTNodeTraits::right(root);
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
				BSTNodeTraits::parent(root, NodePtr());
				root = BSTNodeTraits::right(root);
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
		return new BinarySearchTreeIterator(BSTNodeTraits::lmc(root));
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
class BinarySearchTree : public BinarySearchTreeBase<K, V, BSTNodeTraits<K, V, BSTNode<K, V>, BSTNodeBase<BSTNode<K, V> > >, C>{
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

#include <iostream>
template<class K, class V, class RBNodeTraitsParam, class C = std::less<K> >
class RBTreeBase : public BinarySearchTreeBase<K, V, RBNodeTraitsParam, C>{
	typedef BinarySearchTreeBase<K, V, RBNodeTraitsParam, C> Super;
public:
	typedef typename Super::BSTNodeTraits RBNodeTraits;
	typedef typename Super::NodePtr NodePtr;
	void put(const K& k, const V& v) override {
		NodePtr node = Super::putInner(k, v);
		if(node != NodePtr())
			fixInsert(node);
	}

	void remove(const K&k) override {
		NodePtr node = Super::removeInner(k);
		if(node != NodePtr()){
			assert(RBNodeTraits::left(node) == NodePtr() || RBNodeTraits::right(node) == NodePtr()); 
			fixRemove(node);
			RBNodeTraits::deleteNode(node);
		}
	}

	void removeMin() override {
		NodePtr node = Super::removeMinInner();
		if(node != NodePtr()){
			assert(RBNodeTraits::left(node) == NodePtr() || RBNodeTraits::right(node) == NodePtr()); 
			fixRemove(node);
			RBNodeTraits::deleteNode(node);
		}
	}

	void removeMax() override {
		NodePtr node = Super::removeMaxInner();
		if(node != NodePtr()){
			assert(RBNodeTraits::left(node) == NodePtr() || RBNodeTraits::right(node) == NodePtr()); 
			fixRemove(node);
			RBNodeTraits::deleteNode(node);
		}
	}

	bool isValid() const override {
		if(!Super::isValid())
			return false;
		if(Super::root != NodePtr()){
			if(RBNodeTraits::color(Super::root) != RBNodeTraits::black())
				assert(false);//return false;

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
					assert(false);//return false;
				if(l != NodePtr()){
					bhs.push_back(RBNodeTraits::color(l) == RBNodeTraits::black() ? (bhs[i] + 1) : bhs[i]);
					nodes.push_back(l);
				}else if(!isBhSet){
					bh = bhs[i];
					isBhSet = true;
				}else if(bhs[i] != bh)
					assert(false);//return false;
				if(r != NodePtr()){
					bhs.push_back(RBNodeTraits::color(r) == RBNodeTraits::black() ? (bhs[i] + 1) : bhs[i]);
					nodes.push_back(r);
				}else if(!isBhSet){
					bh = bhs[i];
					isBhSet = true;
				}else if(bhs[i] != bh)
					assert(false);//return false;
				
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
			}else{
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
class RBTree: public RBTreeBase<K, V, RBNodeTraits<K, V, RBNode<K, V>, RBNodeBase<RBNode<K, V> > >, C>{
};
}
#endif //MY_BINARY_SEARCH_TREE_H
