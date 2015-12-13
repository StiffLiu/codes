#include <cassert>
#include <iostream>
namespace my_lib{

class Node234{
	bool _isLeaf{};
protected:
	Node234(bool isLeaf_ = true) : _isLeaf(isLeaf_){}
public:
	bool isLeaf() const {return _isLeaf;}
};

template<class Key> unsigned int height(const Node234* node_);
template<class Key> bool isValid(const Node234* node_);
template<class Key> const Key& deleet(const Node234* node_);
template<class Key> Node234* unite(Node234* n1, Node234* n2_);
template<class Stream, class Key>
void print(Stream& os_, const Node234* node_, const char *prefix_ = " ",
	const char* keySep_ = " ", const char *childSep_ = "\n", unsigned int depth = 0);

template<class Key>
struct Leaf234 : public Node234{
	Key _key{};
public:
	Leaf234(const Key& key_ = Key{}) : _key(key_){}
	const Key& key() const {return _key;}
	bool isValid(){return _isLeaf;}
};

template<class Key> const Leaf234<Key>* minLeaf(const Node234* node_);
template<class Key>
const Key& min(const Node234* node_){return minLeaf<Key>(node_)->key();}

template<class Key>
struct Internal234 : public Node234{
	typedef Leaf234<Key> LeafNode;
	const LeafNode* _keys[3]{nullptr, nullptr, nullptr};
	Node234* _children[4]{nullptr, nullptr, nullptr, nullptr};
	Internal234* unite(unsigned int h, Internal234* node_){
		if(0 == h){
			assert(!_children[0]->isLeaf());
			if(nullptr == _keys[1] && nullptr == node_->_keys[1]){
				// node_ and this can be merged into one node.
				_children[2] = node_->_children[0];
				_children[3] = node_->_children[1];
				_keys[2] = node_->_keys[0];
				recalcMin(1);
				delete node_;
				return this;
			}
			if(nullptr != _keys[2]) return node_; //this node is full, return node_ directly.

			//this node is not full, moving one key from node_ from this node.
			auto index1 = (nullptr != _keys[1] ? 3 : 2);
			auto index2 = node_->keyCount();
			_children[index1] = node_->_children[index2];
			recalcMin(index1 - 1);

			node_->removeKey(index2 - 1);
			//if node_ has more than one key, after moving one from node_ to this
			//it still has keys, so return it.
			//std::cout << "trying to merge this and node_ " << index2 << std::endl;
			if(index2 > 1) return node_;

			//there's no keys in node_ any more.
			// delete it.
			delete node_;
			return this;
		}

		auto index = select();
		//std::cout << "at line : " << __LINE__ << " selected is : " << index << std::endl;
		auto child = ((Internal234*)_children[index])->unite(h - 1, node_);

		//std::cout << "at line : " << __LINE__ << " recal min " << index << std::endl;
		if(index > 0) recalcMin(index - 1);
		if(index < 3 && nullptr != _children[index + 1]) recalcMin(index);
		//std::cout << "at line : " << __LINE__ << " recal min done" << index << std::endl;

		if(child != _children[index]){
			//A new node is returned.
			if(nullptr != _keys[2]){
				//if this node is full, split it.
				//std::cout << "at line : " << __LINE__ << " new node created "<< index << std::endl;
				if(index < 2){
					
					auto root = new Internal234;
					root->_children[1] = _children[3];
					root->_children[0] = _children[2];
					root->_keys[0] = _keys[2];

					removeKey(2);
					if(index == 1){
					 	_children[2] = child; 
						recalcMin(1);
					}else{
						_children[2] = _children[1];
						_children[1] = child;
						_keys[1] = _keys[0];
						recalcMin(0);
					}
					return root;
				}

				Internal234* root = (3 == index ? new Internal234{_children[3], child} : 
						(2 == index ? new Internal234{child, _children[3]} : nullptr));
				assert(nullptr != root);
				removeKey(2);
				return root;
			}

			assert(index <=2);

			//this node is not full, the new node can be inserted to the end.
			//std::cout << "at line : " << __LINE__ << " append " << index << std::endl;
			if(2 == index){
				_children[3] = child;
				recalcMin(2);
			}else if(1 == index){
				_children[3] = _children[2];
				_children[2] = child;
				recalcMin(1);
				if(nullptr != _children[3]) recalcMin(2);
			}else if(0 == index){
				_children[3] = _children[2];
				_children[2] = _children[1];
				_children[1] = _children[0];
				_children[0] = child;
				_keys[2] = _keys[1];
				_keys[1] = _keys[0];
				recalcMin(0);
			}else assert(false);
			return this;
		}

		/*std::cout << "at line : " << __LINE__ << " recal min " << index << std::endl;
		if(index > 0) recalcMin(index - 1);
		if(index < 3 && nullptr != _children[index + 1]) recalcMin(index);
		std::cout << "at line : " << __LINE__ << " recal min done" << index << std::endl;*/
		return this;
	}
public:
	Internal234(const Internal234&) = delete;
	Internal234& operator=(const Internal234&) = delete;
	Internal234(Node234* n1_, Node234* n2_) : Node234(false){
		_children[0] = n1_;
		_children[1] = n2_;
		recalcMin(0);
	}
	Internal234() : Node234(false){}
	bool isValid() const {
		if(nullptr == _keys[0]) /*throw __LINE__;*/return false;
		if(nullptr == _children[0] || nullptr == _children[1]) /*throw __LINE__;*/return false;
		if(nullptr != _keys[1] && nullptr == _children[2]) /*throw __LINE__;*/return false;
		if(nullptr != _keys[2] && nullptr == _children[3]) /*throw __LINE__;*/return false;
		
		const auto& mmin = my_lib::minLeaf<Key>;
		const auto& mheight = my_lib::height<Key>;
		const auto& misValid = my_lib::isValid<Key>;
		auto h = mheight(_children[0]);
		//printKeys(std::cout);
		//std::cout << "\nExpected height is : " << h << " child 1 height is : " << mheight(_children[1]) << std::endl;
		if(h != mheight(_children[1]) || !misValid(_children[1])) /*throw __LINE__;*/return false;
		if(nullptr != _children[2] && (h != mheight(_children[2]) || !misValid(_children[2]))) /*throw __LINE__;*/return false;
		if(nullptr != _children[3] && (h != mheight(_children[3]) || !misValid(_children[3]))) /*throw __LINE__;*/return false;
		if(mmin(_children[0]) != _keys[0] && mmin(_children[1]) != _keys[0]) /*throw __LINE__;*/return false;
		if(nullptr != _keys[1] && mmin(_children[1]) != _keys[1] && mmin(_children[2]) != _keys[1]) /*throw __LINE__;*/return false;
		if(nullptr != _keys[2] && mmin(_children[2]) != _keys[2] && mmin(_children[3]) != _keys[2]) /*throw __LINE__;*/return false;
		return true;
	}

	const Key& min() const {
		return minLeaf()->key();
	}

	const LeafNode* minLeaf() const {
		auto min = _keys[0];
		if(nullptr != _keys[1]){
			if(_keys[1]->key() < min->key()) min = _keys[1];
			if(nullptr != _keys[2] && _keys[2]->key() < min->key()) min = _keys[2];
		}
		return min;
	}

	unsigned int height() const {
		return 1 + my_lib::height<Key>(_children[0]);
	}

	void removeKey(unsigned int index){
		_keys[index] = nullptr;
		_children[index + 1] = nullptr;
	}

	Internal234* insertInner(LeafNode *leaf_){
		assert((Node234*)leaf_ != (Node234*)this);
		if(_children[0]->isLeaf()){
			//std::cout << "children is leaf" << std::endl;
			if(nullptr != _keys[2]){
				auto root = new Internal234(_children[3], leaf_);
				removeKey(2);
				return root;
			}
			auto index = (nullptr != _keys[1] ? 2 : 1);
			_children[index + 1] = leaf_;
			recalcMin(index);
			return this;
		}

		auto index = select();
		//std::cout << "selected is : " << index << std::endl;
		Internal234* child = ((Internal234*)_children[index])->insertInner(leaf_);
		if(child != _children[index]){
			//std::cout << "a new node is created" << std::endl;
			assert(index == keyCount());
			// child is splited.
			// make room for a new node.

			if(nullptr == _keys[2]){
				assert(index != 3);
				_children[index + 1] = child;
				if(index > 0) recalcMin(index - 1);
				recalcMin(index);
				return this;
			}

			//std::cout << "this node will be splitted" << std::endl;
			assert(index == 3);
			auto root = new Internal234{_children[index], child};
			removeKey(2);
			return root;
		}

		//std::cout << "recalculating min value" << std::endl;
		//assert(index != 3);
		if(index > 0) recalcMin(index - 1);
		if(index < 3 && nullptr != _keys[index]) recalcMin(index);
		return this;
	}

	Internal234* insert(LeafNode *leaf_){
		auto newNode = insertInner(leaf_);
		return newNode != this ? new Internal234{this, newNode} : this;
	}

	Internal234* recalcMin(unsigned int index_){
		const auto& mmin = my_lib::minLeaf<Key>;
		auto m1 = mmin(_children[index_]), m2 = mmin(_children[index_ + 1]);
		_keys[index_] = (m1->key() < m2->key() ? m1 : m2);
		return this;
	}

	unsigned int select(){
		auto index = 3;
		if(nullptr == ((Internal234*)_children[0])->_keys[2]) index = 0;
		else if(nullptr == _children[2] || nullptr == ((Internal234*)_children[1])->_keys[2]) index = 1;
		else if(nullptr == _children[3] || nullptr == ((Internal234*)_children[2])->_keys[2]) index = 2;
		return index;
	}

	unsigned int keyCount(){
		return (nullptr != _keys[2]) ? 3 : ((nullptr != _keys[1]) ? 2 : 1);
	}

	// remove the minimum, not implemented.
	// to make the time bound O(log(n)), a pointer to the parent is needed in each node.
	// void extractMin(){
	// }

	// delete a leaf from the heap, not implemented.
	// to make the time bound O(log(n)), a pointer to the parent is needed in each node.
	// void remove(LeafNode* leaf_){
	// }

	// set value for a given leaf node, not implemented.
	// to make the time bound O(log(n)), a pointer to the parent is needed in each node.
	// void set(LeafNode* leaf_, const Key& key){
	// }

	Internal234* unite(Internal234* another_){
		assert(another_ != this);
		auto thisHeight = height();
		auto anotherHeight = another_->height();
		if(thisHeight == anotherHeight) return new Internal234(this, another_);

		if(thisHeight > anotherHeight){
		 	auto newNode = unite(thisHeight - anotherHeight, another_);
			return newNode != this ? new Internal234(this, newNode) : this;
		}
		auto newNode = another_->unite(anotherHeight - thisHeight, this);
		return newNode != another_ ? new Internal234(newNode, another_) : another_;
	}

	Internal234* unite(Node234* node_){
		return node_->isLeaf() ? insert((LeafNode*)node_) : unite((Internal234*)node_);
	}

	void deleet(){
		if(_children[0]->isLeaf()){
			delete _children[0];
			delete _children[1];
			delete _children[2];
			delete _children[3];
		}else{
		  ((Internal234*)_children[0])->deleet();
		  ((Internal234*)_children[1])->deleet();
		  if(nullptr != _children[2]) ((Internal234*)_children[2])->deleet();
		  if(nullptr != _children[3]) ((Internal234*)_children[3])->deleet();
		}
		delete this;
	}

	template<class Stream>
	void printKeys(Stream& os_, const char* keySep_ = " ") const {
		os_ << _keys[0] << " : " << _keys[0]->key();
		if(nullptr != _keys[1]) os_ << keySep_ << _keys[1] << " : " << _keys[1]->key();
		if(nullptr != _keys[2]) os_ << keySep_ << _keys[2] << " : " << _keys[2]->key();
	}
	template<class Stream>
	void printChildren(Stream& os_, const char *prefix_ = " ", const char* keySep_ = " ",
		const char *childSep_ = "\n", unsigned int depth_ = 0) const{
		os_ << childSep_;
		const auto& mprint = my_lib::print<Stream,Key>;
		mprint(os_, _children[0], prefix_, keySep_, childSep_, depth_ + 1);
		os_ << childSep_;
		mprint(os_, _children[1], prefix_, keySep_, childSep_, depth_ + 1);
		if(nullptr != _children[2]){
			os_ << childSep_;
			mprint(os_, _children[2], prefix_, keySep_, childSep_, depth_ + 1);
		}
		if(nullptr != _children[3]){
			os_ << childSep_;
			mprint(os_, _children[3], prefix_, keySep_, childSep_, depth_ + 1);
		}
	}
	template<class Stream>
	void print(Stream& os_, const char *prefix_ = " ", const char* keySep_ = " ",
		const char *childSep_ = "\n", unsigned int depth_ = 0) const{
		for(unsigned int i = 0;i < depth_;++ i) os_ << prefix_;
		printKeys(os_, keySep_);
		printChildren(os_, prefix_, keySep_, childSep_, depth_);
	}

	template<class Stream>
	friend Stream& operator<<(Stream& os_, const Internal234& node_){
		node_.print(os_);
		return os_;
	}
};

template<class Key>
unsigned int height(const Node234* node_){
	return (node_->isLeaf()) ? 0 : ((const Internal234<Key>*)node_)->height();
};

template<class Key>
bool isValid(const Node234* node_){
	return (node_->isLeaf()) || ((const Internal234<Key>*)node_)->isValid();
}

template<class Key>
const Leaf234<Key>* minLeaf(const Node234* node_){
	return (node_->isLeaf()) ? (const Leaf234<Key>*)node_ : ((const Internal234<Key>*)node_)->minLeaf();
}

template<class Key>
void deleet(Node234* node_){
	if(nullptr == node_ || node_->isLeaf()) delete (Leaf234<Key>*)node_;
	else ((Internal234<Key>*)node_)->deleet();
}

template<class Key>
Node234* unite(Node234* n1_, Node234* n2_){
	typedef Leaf234<Key> Leaf;
	typedef Internal234<Key> Internal;
	if(n1_->isLeaf()) 
		return n2_->isLeaf() ? new Internal{n1_, n2_} : ((Internal*)n2_)->insert((Leaf*)n1_);
	return n2_->isLeaf() ? ((Internal*)n1_)->insert((Leaf*)n2_) : ((Internal*)n1_)->unite((Internal*)n2_);
}

template<class Stream, class Key>
void print(Stream& os_, const Node234* node_, const char *prefix_ = " ",
	const char* keySep_ = " ", const char *childSep_ = "\n", unsigned int depth_ = 0){
	if(node_->isLeaf()){
		for(unsigned int i = 0;i < depth_;++ i) os_ << prefix_;
		os_ << node_ << " : " << ((Leaf234<Key>*)node_)->key();
	}else ((Internal234<Key>*)node_)->print(os_, prefix_, keySep_, childSep_, depth_);
}

template<class Key>
class Heap234{
	Node234 *_root{};
	size_t _size{};
public:
	Heap234(const Heap234&) = delete;
	Heap234& operator=(const Heap234&) = delete;
	Heap234(){}
	const Key& min(){ return my_lib::min<Key>(_root);}
	bool empty(){return nullptr != _root;}
	size_t size(){return _size;}
	void clear(){
		my_lib::deleet<Key>(_root);
		_root = nullptr;
		_size = 0;
	}

	void push(const Key& key_){
		auto leaf = new Leaf234<Key>(key_);
		if(nullptr == _root) _root = leaf;
		else _root = unite<Key>(_root, leaf);
		++ _size;
	}

	bool isValid(){return nullptr == _root || my_lib::isValid<Key>(_root);}

	Heap234& merge(Heap234& another_){
		if(this == &another_) return another_;
		if(nullptr == _root) _root = another_._root;
		else if(nullptr != another_._root) _root = unite<Key>(_root, another_._root);
		_size += another_._size;
		another_._root = nullptr;
		another_._size = 0;
		return *this;
	}

	template<class Stream>
	friend Stream& operator<<(Stream& os_, const Heap234& heap_){
		if(nullptr != heap_._root) print<Stream, Key>(os_, heap_._root);
		return os_;
	}
	~Heap234(){clear();}
};
}

