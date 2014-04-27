#ifndef MY_BINARY_SEARCH_TREE_H
#define MY_BINARY_SEARCH_TREE_H
#include "my_ordered_symbol_table.h"
#include <utility>
#include <vector>

namespace my_lib{
template<class K, class V, class C = std::less<K> > 
class BinarySearchTree : public OrderedSymbolTable<K, V>{
public:
	BinarySearchTree(C comparator = C()) : comparator(comparator){
	}
	void put(const K& k, const V& v) override{
		if(root == nullptr){
			root = new Node(nullptr, nullptr, nullptr, k, v, 1);
			return;
		}
		Node *subtree = root;
		//Loop invariant :
		//	subtree is the root of the sub tree, in which k lies.
		while(subtree != nullptr){
			if(comparator(subtree->value.first, k)){
				if(subtree->r == nullptr){
					subtree->r = new Node(subtree, nullptr, nullptr, k, v, 1);
					break;				       	
				}
				subtree = subtree->r;
			}else if(comparator(k, subtree->value.first)){
				if(subtree->l == nullptr){
					subtree->l = new Node(subtree, nullptr, nullptr, k, v, 1);
					break;
				}
				subtree = subtree->l;
			}else{
				subtree->value.second = v;
				return;
			}
		}
		while(subtree != nullptr){
			++ subtree->count;
			subtree = subtree->p;
		}
	}

	const V* get(const K& k) const override{
		Node *node = getNode(k);
		if(node == nullptr)
			return nullptr;
		return &node->value.second;
	}

	void remove(const K& k) override{
		Node *node = getNode(k);
		if(node != nullptr){
			if(node->l != nullptr && node->r != nullptr){
				Node *successor = node->r;
				while(successor->l != nullptr)
					successor = successor->l;
				if(node == root)
					root = successor;
				if(successor->p != nullptr){
					if(successor == successor->p->l)
						successor->p->l = node;
					else
						successor->p->r = node;
				}
				if(node->p != nullptr){
					if(node == node->p->l)
						node->p->l = successor;
					else
						node->p->r = successor;
				}
				std::swap(node->p, successor->p);
				std::swap(node->l, successor->l);
				std::swap(node->r, successor->r);
				std::swap(node->count, successor->count);
				if(node->l != nullptr)
					node->l->p = node;
				if(node->r != nullptr)
					node->r->p = node;
				if(successor->l != nullptr)
					successor->l->p = successor;
				if(successor->r != nullptr)
					successor->r->p = successor;
			}
			//At least one of node's children is null
			Node *c = node->l;
			if(node->r != nullptr)
				c = node->r;
			if(node != root){
				Node *pNode = node->p;
				(node == pNode->l) ? (pNode->l = c) : (pNode->r = c);
				if(c != nullptr)
					c->p = pNode;
				
				while(pNode != nullptr){
					 -- (*pNode).count; 
					pNode = pNode->p;
				}
			}else{
				root = c;
				if(root != nullptr)
					root->p = nullptr;
			}
			delete node;
		}
	}

	const K* min() const override{
		if(root == nullptr)
			return nullptr;
		Node *leftMost = root->lmc();
		return &leftMost->value.first;
	}

	const K* max() const override{
		if(root == nullptr)
			return nullptr;
		Node *rightMost = root->rmc();
		return &rightMost->value.first;
	}

	const K* floor(const K& k) const override{
		Node *subtree = root;
		while(subtree != nullptr){
			if(comparator(subtree->value.first, k)){
				if(subtree->r == nullptr)
					break;
				subtree = subtree->r;
			}else if(comparator(k, subtree->value.first)){
				if(subtree->l == nullptr){
					subtree = subtree->rmp()->p;
					break;
				}
				subtree = subtree->l;
			}else{
				break;
			}
		}
		if(subtree == nullptr)
			return nullptr;
		return &subtree->value.first;
	}

	const K* ceil(const K& k) const override{
		Node *subtree = root;
		while(subtree != nullptr){
			if(comparator(subtree->value.first, k)){
				if(subtree->r == nullptr){
					subtree = subtree->lmp()->p;
					break;
				}
				subtree = subtree->r;
			}else if(comparator(k, subtree->value.first)){
				if(subtree->l == nullptr)
					break;
				subtree = subtree->l;
			}else{
				break;
			}
		}
		if(subtree == nullptr)
			return nullptr;
		return &subtree->value.first;
	}

	unsigned int rank(const K& k) const override{
		Node *subtree = root;
		unsigned int count = 0;
		while(subtree != nullptr){
			if(comparator(subtree->value.first, k)){
				count += (subtree->l == nullptr ? 0 : subtree->l->count) + 1;
				if(subtree->r == nullptr){
					break;
				}
				subtree = subtree->r;
			}else if(comparator(k, subtree->value.first)){
				if(subtree->l == nullptr){
					break;
				}
				subtree = subtree->l;
			}else{
				count += (subtree->l == nullptr ? 0 : subtree->l->count);
				break;
			}
		}
		return count;
	}

	const K* select(unsigned int index) const override{
		if(root == nullptr)
			return nullptr;
		if(index >= root->count)
			return nullptr;

		Node *n = root;
		unsigned int tmp = 0;
		if(n->l != nullptr)
			tmp = n->l->count;
		while(index != tmp){
			if(index > tmp){
				n = n->r;
				index -= (tmp + 1);
			}else{
				n = n->l;
			}
			tmp = (n->l == nullptr ? 0 : n->l->count);
		}
		return &n->value.first;			
	}

	void removeMin() override{
		if(root != nullptr){
			Node *leftMost = root;
			while(leftMost->l != nullptr){
				--leftMost->count;
				leftMost = leftMost->l;
			}
			if(leftMost != root){
				leftMost->p->l = leftMost->r;
				if(leftMost->r != nullptr)
					leftMost->r->p = leftMost->p;
			}else{
				root->p = nullptr;
				root = root->r;
			}
			delete leftMost;
		}
	}

	void removeMax() override{
		if(root != nullptr){
			Node *rightMost = root;
			while(rightMost->r != nullptr){
				--rightMost->count;
				rightMost = rightMost->r;
			}
			if(rightMost != root){
				rightMost->p->r = rightMost->l;
				if(rightMost->l != nullptr)
					rightMost->l->p = rightMost->p;
			}else{
				root->p = nullptr;
				root = root->l;
			}
			delete rightMost;
		}
	}

	bool isEmpty() const override{
		return root == nullptr;
	}

	unsigned int size() const override{
		if(root == nullptr)
			return 0;
		return root->count; 
	}

	unsigned int size(const K& l, const K& h) const override{
		unsigned int i = rank(l), j = rank(h);
		if(i > j)
			return 0;
		if(getNode(h) != nullptr)
			return j - i + 1;
		return j - i;
	}

	void clear() override{
		if(root != nullptr){
			unsigned int count = root->count;
			std::vector<Node*> ps(count);
			ps[0] = root;
			for(unsigned int i = 0, j = 1;i < count;++ i){
				Node *n = ps[i];
				if(n->l != nullptr){
					ps[j] = n->l;
					++ j;
				}	
				if(n->r != nullptr){
					ps[j] = n->r;
					++ j;
				}
				delete n;
			}	
			root = nullptr;
		}
	}
	/**
	 * @return true if the data structure is valid.
	 */
	bool isValid() const{
		return root == nullptr || (root->p == nullptr && root->isValid());
	}

	~BinarySearchTree(){
		clear();
	}
private:
	typedef OrderedSymbolTable<K, V> Super;
	typedef typename Super::IteratorImpl IteratorImpl;
protected:
	typedef std::pair<K, V> Pair;
	struct Node{
		Node *p, *l, *r;
		unsigned int count;
		Pair value;
		Node(Node *p, Node *l, Node *r, const K& k, const V& v, unsigned int count)
			:p(p), l(l), r(r), count(count), value(Pair(k, v)){
		}
		/**
		 * @return Left most child of this node
		 */
		Node *lmc(){
			Node *leftMost = this;
			while(leftMost->l != nullptr)
				leftMost = leftMost->l;
			return leftMost;
		}
		/**
		 * @return Right most child of this node
		 */
		Node *rmc(){
			Node *rightMost = this;
			while(rightMost->r != nullptr)
				rightMost = rightMost->r;
			return rightMost;
		}
		/**
		 * @return Left most parent of this node
		 */
		Node *lmp() {
			Node *leftMost = this;
			Node *p = leftMost->p;
			while(p != nullptr && p->r == leftMost){
				leftMost = p;
				p = leftMost->p;
			}
			return leftMost;
		}
		/**
		 * @return Right most parent of this node
		 */
		Node *rmp(){
			Node *rightMost = this;
			Node *p = rightMost->p;
			while(p != nullptr && p->l == rightMost){
				rightMost = p;
				p = rightMost->p;
			}
			return rightMost;
		}
		/**
		 * @return Successor of this node
		 */
		Node *suc(){
			if(r != nullptr)
				return r->lmc();
			Node *node = lmp();
			if(node->p != nullptr && node == node->p->l)
				return node->p;
			return nullptr;
		}
		/**
		 * @return Predecessor of this node
		 */
		Node *pre(){
			if(l != nullptr)
				return l->rmc();
			Node *node = rmp();
			if(node->p != nullptr && node == node->p->r)
				return node->p;
			return nullptr;
		}

		/**
		 * @return {@code true} if the data structure is valid.
		 */
		bool isValid() const{
			if(p == this){
				return false;
			}
			unsigned int tmp = 0;
			if(l != nullptr){
				if(l->p != this || !l->isValid())
					return false;
				tmp = l->count;
			}
			if(r != nullptr){
				if(r->p != this || !r->isValid())
					return false;
				tmp += r->count;
			}
			return tmp + 1 == count;
		}
	};
	Node *getNode(const K& k) const{
		Node *subtree = root;
		//Loop invariant :
		//	subtree is the root of the sub tree, in which k may lies.
		while(subtree != nullptr){
			if(comparator(subtree->value.first, k)){
				if(subtree->r == nullptr)
					return nullptr;
				subtree = subtree->r;
			}else if(comparator(k, subtree->value.first)){
				if(subtree->l == nullptr)
					return nullptr;
				subtree = subtree->l;
			}else{
				return subtree;
			}
		}
		return nullptr;
	}
	struct BinarySearchTreeIterator : public IteratorImpl{
		void  next() override{
			if(node != nullptr)
				node = node->suc();
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
			return node->value.first;
		}

		const V& value() const override{
			return node->value.second;
		}

		BinarySearchTreeIterator* copy() const override{
			return new BinarySearchTreeIterator(node);
		}

		BinarySearchTreeIterator(Node *node) : node(node){
		}
	private:
		Node *node;
		
	};

	BinarySearchTreeIterator *implBegin() const override{
		return new BinarySearchTreeIterator(root->lmc());
	}

	BinarySearchTreeIterator *implEnd() const override{
		return new BinarySearchTreeIterator(nullptr);
	}	
protected:
	Node *root = nullptr;
	C comparator;
};
}
#endif //MY_BINARY_SEARCH_TREE_H

