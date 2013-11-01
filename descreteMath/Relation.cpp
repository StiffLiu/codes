#include <iostream>
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <map>
#include <algorithm>
#include <set>
#include <vector>
using namespace std;

class Relation{
	unsigned int *relation;
	unsigned int count;
	void set(unsigned int n, unsigned int *rel = NULL){
		delete []relation;
		count = n;
		if(n == 0)
			relation = NULL;
		else{
			n *= count;		
			unsigned int rows = ((n / unitBits) + (n % unitBits == 0 ? 0 : 1));
			relation = new unsigned int[rows];
			set(rel);
		}
	}
public:
	static const int unitBits = sizeof(unsigned int) * 8;
	Relation(unsigned int n = 0){
		relation = NULL;
		set(n);
	}
	Relation(const Relation& rel){
		relation = NULL;		
		set(rel.count, rel.relation);
	}
	Relation& operator =(const Relation& rel){
		if(this == &rel)
			return *this;		
		set(rel.count, rel.relation);
		return *this;
	}
	unsigned int getCount() const{
		return count;
	}
	bool isReflexive(){
		for(unsigned int i = 0;i < count;++ i)
			if(!get(i, i))
				return false;
		return true;
	}
	bool isIrreflexive(){
		for(unsigned int i = 0;i < count;++ i)
			if(get(i, i))
				return false;
		return true;
	}
	bool isTransitive(){
		Relation tmp = *this;
		for(unsigned int i = 1;i < count;++ i){
			tmp *= *this;
			if(!isSubset(tmp))
				return false;
		}
		return true;
	}
	bool isSymmetric(){
		for(unsigned int i = 1;i < count;++ i)
			for(unsigned int j = 0;j < i;++ j)
				if(get(i, j) != get(j, i))
					return false;
		return true;
	}
	bool isAntisymmetric(){
		for(unsigned int i = 1;i < count;++ i)
			for(unsigned int j = 0;j < i;++ j)
				if(get(i, j) && get(j, i))
					return false;
		return true;
	}
	bool isSubset(const Relation& rel){
		assert(count == rel.count);
		unsigned int n = count * rel.count;
		unsigned int rows = ((n / unitBits) + (n % unitBits == 0 ? 0 : 1));		
		for(unsigned int i = 0;i < rows;++ i)
			if(rel.relation[i] & ~relation[i])
				return false;
		return true;
	}
	void reflexiveClosure(){
		for(unsigned int i = 0;i < count;++ i)
			set(i, i);
	}
	void symmetricClosure(){
		for(unsigned int i = 1;i < count;++ i)
			for(unsigned int j = 0;j < i;++ j)
				if(get(i, j) != get(j, i)){
					set(i, j);
					set(j, i);
				}
	}
	void set(unsigned int *rel = NULL){
		unsigned int n = count * count;		
		unsigned int rows = ((n / unitBits) + (n % unitBits == 0 ? 0 : 1));			
		if(rel != NULL)
			for(unsigned int i = 0;i < rows;++ i)
				relation[i] = rel[i];
		else
			for(unsigned int i = 0;i < rows;++ i)
				relation[i] = 0;
	}
	void transitiveClosureUsingWarshall(){
		for(unsigned int k = 0;k < count;++ k)
			for(unsigned int i = 0;i < count;++ i)
				if(get(i, k)){
					for(unsigned int j = 0;j < count;++ j)
						if(!get(i, j) && get(k, j))
							set(i, j);
				}
	}
	void transitiveClosureUsingMatrixJoin(){
		Relation tmp = *this;
		for(unsigned int i = 1;i < count;++ i){
			tmp *= *this;
			*this |= tmp;
		}
	}
	void transitiveClosure(){
		bool warshall = true;
		if(warshall){
			transitiveClosureUsingWarshall();
		}else{
			transitiveClosureUsingMatrixJoin();
		}
	}
	void equivalenceClosure(){
		reflexiveClosure();
		symmetricClosure();
		transitiveClosure();
	}
	void set(unsigned int i, unsigned int j){
		assert(i < count && j < count);
		unsigned int nBit = i * count + j;
		unsigned int row = nBit / unitBits;
		unsigned int col = nBit % unitBits;
		relation[row] |= (1 << col);
	}
	void unset(unsigned int i, unsigned int j){
		assert(i < count && j < count);
		unsigned int nBit = i * count + j;
		unsigned int row = nBit / unitBits;
		unsigned int col = nBit % unitBits;
		relation[row] &= ~(1 << col);
	}
	bool isPartialOrder(){
		return isReflexive() && isAntisymmetric() && isTransitive();
	}
	bool isTotalOrder(){
		if(isReflexive()){
			for(unsigned int i = 1;i < count;++ i)
				for(unsigned int j = 0;j < i;++ j)
					if(get(i, j) == get(j, i))
						return false;
			return true;
		}
		return false;
	}
	unsigned int findLub(unsigned int i, unsigned j){
		if(get(i, j))
			return j;
		unsigned int *indices = new unsigned int[count];
		unsigned int current = 0;
		for(unsigned int k = 0;k < count;++ k){
			if(get(i, k) && get(j, k)){
				indices[current] = k;
				++current;
			}
		}
		for(unsigned int k = 0;k < current;++ k){
			unsigned int t = 0;
			for(;t < current;++ t)
				if(!get(indices[k], indices[t]))
					break;
			if(t == current)
				return indices[k];
		}
		return -1;
	}
	unsigned int findGlb(unsigned int i, unsigned j){
		if(get(i, j))
			return i;
		unsigned int *indices = new unsigned int[count];
		unsigned int current = 0;
		for(unsigned int k = 0;k < count;++ k){
			if(get(k, i) && get(k, j)){
				indices[current] = k;
				++current;
			}
		}
		for(unsigned int k = 0;k < current;++ k){
			unsigned int t = 0;
			for(;t < current;++ t)
				if(!get(indices[t], indices[k]))
					break;
			if(t == current)
				return indices[k];
		}
		return -1;
	}
	bool isLattice(){
		if(isPartialOrder()){
			for(unsigned int i = 0;i < count;++ i)
				for(unsigned int j = i + 1;j < count;++ j)
					if(findLub(i, j) == -1 || findGlb(i, j) == -1)
						return false;
			return true;
		}
		return false;
	}
	bool get(unsigned int i, unsigned int j) const{
		assert(i < count && j < count);
		unsigned int nBit = i * count + j;
		unsigned int row = nBit / unitBits;
		unsigned int col = nBit % unitBits;
		return relation[row] & (1 << col);
	}
	bool linearize(){
		if(isPartialOrder()){
			if(isTotalOrder())
				return true;
			bool *isUsed = new bool[count];
			unsigned int i = 0;
			for(unsigned int j = 0;j < count;++ j)
				isUsed[j] = false;
			while(i < count){
				for(i = 0;i < count;++ i){
					if(!isUsed[i]){
						unsigned int j = 0;
						for(;j < count;++ j)
							if(j != i && !isUsed[j] && get(j, i))
								break;
						if(j >= count){
							isUsed[i] = true;
							for(j = 0;j < count;++ j)
								if(!isUsed[j])
									set(i, j);
							i = 0;
							break;
						}
					}
				}
			}
			delete []isUsed;
			return true;
		}
		assert(false);
		return false;		
	}
	void power(unsigned int p){
		if(p == 0){
			unsigned int nBits = count * count;
			unsigned int row = nBits / unitBits;
			if(nBits % unitBits != 0)
				++ row;
			for(unsigned int i = 0;i < row;++ i)
				relation[i] = 0;	
			for(unsigned int i = 0;i < count;++ i)
				set(i, i);
		}else{
			Relation cp = *this;
			p >>= 1;
			while(p != 0){
				cp *= cp;
				if(p % 2 == 1) (*this) *= cp;
				p >>= 1;
			}
		}
	}
	Relation& operator |=(const Relation& rel){
		assert(count == rel.count);
		unsigned int nBits = count * count;
		unsigned int row = nBits / unitBits;
		if(nBits % unitBits != 0)
			++ row;
		for(unsigned int i = 0;i < row;++ i)
			relation[i] |= rel.relation[i];		
		return *this;
	}
	Relation& operator &=(const Relation& rel){
		assert(count == rel.count);
		unsigned int nBits = count * count;
		unsigned int row = nBits / unitBits;
		if(nBits % unitBits != 0)
			++ row;
		for(unsigned int i = 0;i < row;++ i)
			relation[i] &= rel.relation[i];		
		return *this;
	}
	Relation& operator *=(const Relation& rel){
		assert(count == rel.count);
		Relation tmp = *this;
		for(unsigned int i = 0;i < count;++ i){
			for(unsigned int j = 0;j < count;++ j){
				bool val = false;
				for(unsigned int k = 0;k < count;++ k){
					val |= (tmp.get(i, k) & rel.get(k, j));
				}
				if(val)
					set(i, j);
				else
					unset(i, j);
			}
		}
		return *this;
	}
	void moveToNext(){
		unsigned int n = count * count;
		unsigned int vSize = n / unitBits;
		unsigned int carry = 1;
		for(unsigned int i = 0;i < vSize;++ i){
			unsigned int oldCarry = carry;
			if(relation[i] == -1)
				carry = 1;
			else
				carry = 0;
			relation[i] += oldCarry;
		}
		if(n % unitBits != 0){
			relation[vSize] += carry;
		}
	}
	friend ostream& operator<<(ostream& os, const Relation& rel){
		for(unsigned int i = 0;i < rel.count;++ i){
			for(unsigned int j = 0;j < rel.count;++ j)
				os << (rel.get(i, j) ? '1' : '0') << (j == rel.count - 1 ? '\n' : ' ');
		}
		return os;
	}
	template<class T>
	friend ostream& display(ostream& os, const Relation& rel, T vals){
		os << '{';
		for(unsigned int i = 0;i < rel.count;++ i){
			for(unsigned int j = 0;j < rel.count;++ j)
				if(rel.get(i, j))
					os << '(' << vals[i] << ',' << vals[j] << ") ";
		}
		os << '}';
		return os;
	}
};
template<class T>
class Combinator{
	T count;
	T *combinations;
	Combinator(const Combinator&);
	Combinator& operator=(const Combinator&);
public:
	Combinator(T count){
		T n = count + 1;
		this->count = count;
		combinations = new T[n * n];
		for(T i = 0;i <= count;++ i)
			combinations[i * n + i] = 1;
		for(T i = 0;i <= count;++ i)
			combinations[i * n] = 1;
		for(T i = 2;i <= count;++ i){
			for(int j = 1;j < i;++ j)
				combinations[i * n + j] = combinations[(i - 1) * n + j - 1] + combinations[(i - 1) * n + j];
			for(int j = i + 1;j <= count;++ j)
				combinations[i * n + j] = 0;
		}
		
	}
	T *operator[](T i){
		return &combinations[i * (count + 1)];
	}
	~Combinator(){
		delete []combinations;
	}
};
class EquivalenceRelationIterator{
	int *tmpVals;
	int count;
	vector<int*> pts;
	EquivalenceRelationIterator(const EquivalenceRelationIterator&);
	EquivalenceRelationIterator& operator=(const EquivalenceRelationIterator&);
public:
	EquivalenceRelationIterator(unsigned int count){
		this->count = count;
		tmpVals = NULL;
		if(count != 0){
			tmpVals = new int[count];
			for(int i = 0;i < count;++ i)
				tmpVals[i] = i;
		}
	}
	void nextCombination(int *vals, int m, int n){
		assert(m >= n);
		if(m == n){
			pts.push_back(vals + n);
			if(vals + n == tmpVals + count){
				int count = pts.size();--count;
				cout << '{';
				for(int i = 0;i < count;++ i){
					int *start = pts[i], *end = pts[i + 1];
					if(start != end){
						cout << '{';
						while(start < end){
							cout << *start;					
							++start;
							if(start == end)
								break;
							cout << ' ';
						}
						cout << '}';
					}
				}
				cout << '}' << endl;
				
			}else{
				iterate(vals + n, tmpVals + count - vals - n);
			}
			pts.pop_back();
		}else if(n == 0){
			pts.push_back(vals + n);
			iterate(vals + n, tmpVals + count - vals - n);
			pts.pop_back();
		}else{
			nextCombination(vals + 1, m - 1, n - 1);
			int *last = vals + m -1;
			int tmp = *vals;
			*vals = *last;
			*last= tmp;
			nextCombination(vals, m - 1, n);
			*last = *vals;
			*vals = tmp;
		}
	}
	void print(){
		iterate(tmpVals, count);
	}
	void iterate(int *tmp, int cnt){
		for(int i = 1;i <= cnt;++ i){
			pts.push_back(tmp);
			nextCombination(tmp + 1, cnt - 1, i - 1);
			pts.pop_back();
		}
	}
	static void countNum(unsigned int n){
		Combinator<unsigned long long> combinator(n);
		unsigned long long *counts = new unsigned long long[n + 1];
		if(n >= 0)
			counts[0] = 1;
		if(n >= 1)
			counts[1] = 1;
		for(unsigned int i = 2;i <= n;++ i){
			counts[i] = 0;
			for(unsigned int j = 0;j < i;++ j)
				counts[i] += combinator[i - 1][j] * counts[i - 1 - j];
		}
		for(unsigned int i = 0;i <= n;++ i)
			cout << counts[i] << ' ';
		cout << endl;
	}
};

class NaryRelation{
	vector<void*> datas;
	unsigned int degree;
public:
	NaryRelation(unsigned int deg = 0){
		degree = deg;
	}
	unsigned int getDegree(){
		return degree;
	}
	vector<void*>& getDatas(){
		return datas;
	}
	bool validateFields(const vector<int>& fields){
		for(unsigned int i = 0;i < fields.size();++ i)
			if(fields[i] >= degree){
				assert(false);
				return false;
			}
		return true;
	}
	template<class C>
	NaryRelation projection(const vector<int>& fields, C comp){
		if(fields.empty() || !validateFields(fields))
			return NaryRelation();
		assert(datas.size() % degree == 0);
		int count = datas.size();
		set<int ,C> tmp(comp);
		vector<int> indices;
		for(int i = 0;i < count;i += degree)
			if(tmp.find(i) == tmp.end()){
				tmp.insert(i);
				indices.push_back(i);
			}
		assert(!tmp.empty());
		NaryRelation ret(fields.size());
		for(int i = 0;i < indices.size();++ i){
			int index = indices[i];
			for(int j = 0;j < fields.size();++ j)
				ret.datas.push_back(datas[index + fields[j]]);
		}
		return ret;
		
	}
	template<class C1, class C2>
	NaryRelation join(const vector<int>& fields1, const NaryRelation& rel, const vector<int>& fields2,
		C1 c1, C2 c2){
		if(fields1.empty() || !validateFields(fields1)
			|| fields2.empty() || !rel.validateFields(fields2) ||
			fields1.size() != fields2.size())
			return NaryRelation();
		NaryRelation ret;//(degree + ;
		return ret;		
	}
	template<class O>
	ostream& print(ostream& os, O o){
		int count = datas.size();		
		for(int i = 0;i < count;i += degree){
			for(int j = 0;j < degree;++ j){
				o.putOne(os, datas[i + j]);
				if(j < degree - 1)
					o.putFieldDelimiter(os);
			}
			if(i < count - degree)
				o.putRecordDelimiter(os);
		}
		return os;
	}
};
Relation& randRelation(Relation& rel, int max = RAND_MAX){
	int count = rel.getCount();
	unsigned int *val = NULL;	
	int byteCount = (count * count) / Relation::unitBits;
	if((count * count) % Relation::unitBits != 0)
		++byteCount;
	val = new unsigned int[byteCount];
	for(int i = 0;i < byteCount;++ i)
		val[i] = rand() % RAND_MAX;
	rel.set(val);
	delete []val;
	return rel;
}
int test1(int argc, char *argv[]){
	int vals[]={1,2,3,4,5};
	int count = sizeof vals / sizeof(*vals);
	int relCount = (1 << (count * count));
	Relation rel(count);
	for(int i = 1;i <= relCount;++ i){
		if(rel.isLattice())/*isPartialOrder())*/ /*rel.isReflexive() && rel.isTransitive() && rel.isSymmetric())*/{
			display(cout, rel, vals);
			cout << endl;
		}
		rel.moveToNext();
	}
	return 0;
}
int test2(int argc, char *argv[]){
	int count = 5;
	Relation rel1(count), rel2(count);
	if(true){
		randRelation(rel1);
	}else{
		unsigned int vals[][2]={
			{1, 2},{2, 3},{3, 4},
			{4, 5},{5, 2},{3, 2},
		};
		for(int i = 0;i < (sizeof vals / sizeof *vals);++ i){
			rel1.set(vals[i][0] - 1, vals[i][1] - 1);
		}
	
	}
	rel2 = rel1;
	cout << "original relation:" << endl;
	cout << rel1 << endl;
	rel1.transitiveClosureUsingWarshall();
	rel2.transitiveClosureUsingMatrixJoin();
	cout << "transitive closure using warshall algorithm:" << endl;
	cout << rel1 << endl;
	cout << "transitive closure using matrix join algorithm:" << endl;
	cout << rel2 << endl;	
	return 0;
}
int test3(int argc, char *argv[]){
	int count = 6;
	int iterCount = 30;
	Relation rel(count);
	clock_t start = clock();
	while((clock() - start) < CLOCKS_PER_SEC * iterCount)
		if(randRelation(rel).isPartialOrder())
			break;
	if(!rel.isPartialOrder()){
		cerr << "couldnot genetate a random partial order in " << iterCount << " seconds" << endl;
		return 1;
	}else{
		cout << "orignal relation : " << endl;
		cout << rel << endl;
		cout << "after linearization : " << endl;
		rel.linearize();
		assert(rel.isTotalOrder());
		cout << rel << endl;
	}
	return 0;
}
int test4(int argc, char *argv[]){
	EquivalenceRelationIterator tmp(10);
	tmp.print();
	cout << endl;
	EquivalenceRelationIterator::countNum(27);
}
class Comparer{
	vector<void*> *datas;
	vector<int>* flds;
public:
	Comparer(vector<void*> *datas = NULL,
		vector<int>* flds = NULL) : datas(datas), flds(flds){}
	int compare(void *p1, void *p2)const{
	}
	int compare(int index1, int index2) const{
		for(int i = 0;i < flds->size();++ i){
			int idx = (*flds)[i];
			int ret = (int)(*datas)[index1 + idx] - (int)(*datas)[index2 + idx];
			if(ret != 0)
				return ret;
		}
		return 0;
	}
	bool operator()(int index1, int index2) const{
		return compare(index1, index2) < 0;
	}
};
class Outputer{
public:
	void putOne(ostream& os, void *data){
		os << (int)data;
	}
	void putFieldDelimiter(ostream& os){
		os << '\t';
	}
	void putRecordDelimiter(ostream& os){
		os << endl;
	}
};
int test5(int argc, char *argv[]){
	int degree = 15;
	NaryRelation rel(degree);
	vector<int> tmpFields;
	vector<void*>& tmp = rel.getDatas();
	tmpFields.push_back(0);
	tmpFields.push_back(1);
	for(int i = 0;i < 100;++ i)
		for(int j = 0;j < degree;++ j)
			tmp.push_back((void*)(rand() % 10));
	rel.print(cout, Outputer());
	NaryRelation proj = rel.projection(tmpFields, Comparer(&tmp, &tmpFields));
	cout << "\nprojection:" << endl;
	proj.print(cout, Outputer());
	cout << endl;
	return 0;
}
int test6(int argc, char *argv[]){
	int n = 10;
	Relation rel(n);
	for(int i = 0;i < n * n / 8;++ i){
		int tmp = rand() % (n * n);
		rel.set(tmp / n, tmp % n);
	}
	cout << rel << endl;
	rel.transitiveClosureUsingWarshall();
	cout << rel << endl;
}
int main(int argc, char *argv[]){
	srand(time(0));
	return test6(argc, argv);
}
