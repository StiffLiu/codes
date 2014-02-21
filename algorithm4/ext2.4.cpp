/**
 * heap with explict links.
 */
template<class T>
class LinkedHeap {
public:
	template<class U>
	struct Link {
		T val;
		U p, l, r;
		Link(const T& val, U p, U l, U r) :
				val(val), p(p), l(l), r(r) {
		}
	};
	struct PointerLink: public Link<PointerLink*> {
		PointerLink(T val, PointerLink* p, PointerLink* l, PointerLink* r) :
				Link<PointerLink*>(val, p, l, r) {
		}
	};
	void insert(const T& val) {
		if (last == 0) {
			last = new PointerLink(val, 0, 0, 0);
			return;
		}
		PointerLink *l = last;
		while (l->parent != 0 && l == l->parent->r)
			;
		l = l->parent;
		if (l->parent != 0) {
			if (l->right == 0) {
				last = new PointerLink(val, l, 0, 0);
				l->right = last;
				return;
			}
			l = l->right;
		}
		while (l->l != 0)
			l = l->l;
		last = new PointerLink(val, l, 0, 0);
		l->l = last;
	}
	LinkedHeap() :
			last(0) {
	}
	~LinkedHeap() {
		deleteNode();
	}
	int height() {
		return height(last);
	}
	LinkedHeap(const LinkedHeap&) = delete;
	LinkedHeap& operator=(const LinkedHeap&) = delete;
private:
	int height(PointerLink* link) {
		if (link == 0 || (link->l == 0 && link->r == 0))
			return 0;
		auto hl = height(link->l);
		auto hr = height(link->r);
		if (hl > hr)
			return 1 + hl;
		return 1 + hr;
	}
	void deleteNode() {
		PointerLink *root = last;
		while (root->p != 0)
			root = root->p;
		deleteNode(root);
		last = 0;
	}
	void deleteNode(PointerLink *l) {
		if (l != 0) {
			if (l->l != 0)
				deleteNode(l->l);
			if (l->r != 0)
				deleteNode(l->r);
		}
		delete l;
	}
	PointerLink *last;
};
#include <queue>
#include <map>
#include <iostream>
#include <iomanip>

template<class Key, class Value>
struct KeyValue {
	Key key;
	Value value;
	KeyValue(const Key& k, const Value& v) : key(k), value(v){
	}
	bool operator==(const KeyValue& keyValue) const {
		return key == keyValue.key;
	}
	bool operator<(const KeyValue& keyValue) const {
		return key < keyValue.key;
	}
	bool operator>(const KeyValue& keyValue) const {
		return key > keyValue.key;
	}
};
/**
 * This function find all the pair of integers (a, b) and (c, d) such that :
 * 		1. 0 <= a, b, c, d <= n
 * 		2. a^3 + b^3 = c^3 + d^3
 * @param n
 */
void cubesum(unsigned int n) {
	typedef KeyValue<unsigned int, std::pair<unsigned int , unsigned int>> One;
	std::priority_queue<One, std::vector<One>, std::greater<One> > sums;
	for(unsigned int i = 0;i <= n;++ i)
		sums.push({i * i * i, {i, 0}});

	decltype(One::key) sum = 0;
	std::cout << std::setw(5) << sum;
	while(!sums.empty()){
		auto one = sums.top();
		if(one.key != sum){
			sum = one.key;
			std::cout << std::endl << std::setw(5) << sum;
			sum = one.key;
		}
		std::cout <<std::setw(10) << '{' << one.value.first << ',' << one.value.second << '}' << ' ';
		if(one.value.second < n){
			auto i = one.value.first, j = one.value.second + 1;
			sums.push({i * i * i + j * j * j, {i, j}});
		}
		sums.pop();
	}
}

int ext24Main(int argc, char *argv[]){
	cubesum(100);
	return 0;
}
