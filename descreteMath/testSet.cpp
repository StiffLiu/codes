#include <vector>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <set>
#include <map>
using namespace std;

template<class T, class UniversalSet>
class Set{
	static const int bitsCount = (sizeof(unsigned int) * 8);
	const UniversalSet *u;
	unsigned int *vals;
	int getVsize() const{
		int size = u->getSize();
		int vSize = size / bitsCount;
		if(size % bitsCount != 0)
			++ vSize;
		return vSize;
	}
public:
	template<class ForwardIterator>
	Set(const UniversalSet *u, ForwardIterator begin, ForwardIterator end) : u(u){
		int vSize = getVsize();
		vals = new unsigned int[vSize];
		memset((void*)vals, 0, vSize * bitsCount / 8);
		while(begin != end){
			int index = u->getIndex(*begin);
			int sIndex = index / bitsCount;
			int pIndex = index % bitsCount;
			vals[sIndex] |= (1 << pIndex);
			++begin;
		}
	}
	Set(const Set& s) : u(s.u){
		int vSize = getVsize();
		vals = new unsigned int[vSize];
		memcpy((void*)vals, (void*)s.vals, vSize * bitsCount / 8);
	}
	friend ostream& operator<<(ostream& os, const Set& s){
		int vSize = s.getVsize();
		for(int i = 0;i < vSize - 1;++ i){
			int bIndex = i * bitsCount;
			unsigned int val = s.vals[i];
			while(val != 0){
				if((val & 1) != 0)
					os << s.u->getValue(bIndex) << ' '; 
				val >>= 1;
				++bIndex;
			}
		}	
		if(vSize > 0){
			int tmp = s.u->getSize() - bitsCount * (vSize - 1);
			unsigned int val = (1 << tmp);
			int bIndex = (vSize - 1) * bitsCount;
			val -= 1;
			val = (val & s.vals[vSize - 1]);
			while(val != 0){
				if((val & 1) != 0)
					os << s.u->getValue(bIndex) << ' '; 
				val >>= 1;
				++bIndex;
			}
		}
		return os;
	}
	Set& operator&=(const Set<T, UniversalSet>& s){
		if(u == s.u){
			int vSize = s.getVsize();
			for(int i = 0;i < vSize;++ i)
				vals[i] &= s.vals[i];
		}
		return *this;
	}
	Set& operator|=(const Set<T, UniversalSet>& s){
		if(u == s.u){
			int vSize = s.getVsize();
			for(int i = 0;i < vSize;++ i)
				vals[i] |= s.vals[i];
		}
		return *this;
	}
	Set& operator^=(const Set<T, UniversalSet>& s){
		if(u == s.u){
			int vSize = s.getVsize();
			for(int i = 0;i < vSize;++ i)
				vals[i] ^= s.vals[i];
		}
		return *this;
	}
	Set& operator-=(const Set<T, UniversalSet>& s){
		if(u == s.u){
			int vSize = s.getVsize();
			for(int i = 0;i < vSize;++ i)
				vals[i] &= ~s.vals[i];
		}
		return *this;
	}
	friend Set operator-(const Set& s1, const Set& s2){
		Set s3 = s1;
		s3 -= s2;
		return s3;
	}
	friend Set operator^(const Set& s1, const Set& s2){
		Set s3 = s1;
		s3 ^= s2;
		return s3;
	}
	friend Set operator|(const Set& s1, const Set& s2){
		Set s3 = s1;
		s3 |= s2;
		return s3;
	}
	friend Set operator&(const Set& s1, const Set& s2){
		Set s3 = s1;
		s3 &= s2;
		return s3;
	}
	void moveToNext(){
		int vSize = getVsize();
		int carry = 1;
		for(int i = 0;i < vSize;++ i){
			int oldCarry = carry;
			if(vals[i] == -1)
				carry = 1;
			else
				carry = 0;
			vals[i] += oldCarry;
		}
	}
	void flip(){
		int vSize = getVsize();
		for(int i = 0;i < vSize;++ i){
			vals[i] = ~vals[i];
		}
	}
	Set operator~(){
		Set tmp = *this;
		tmp.flip();
		return tmp;
	}
	Set& operator=(const Set& s){
		if(this != &s){
			if(u == s.u)
				memcpy((void*)vals, (void*)s.vals, getVsize() * sizeof(unsigned int));
		}	
		return *this;		
	}
	~Set(){
		delete []vals;
	}
};
template<class Functor>
bool isOneToOne(Functor& f, unsigned int n){
	vector<bool> results(n, false);
	for(int i = 1;i <= n;++ i){
		unsigned int ret = f(i);
		if(results[ret])
			return false;
		results[ret] = true;
	}
	return true;
}
template<class Functor>
bool isOneToOne(const Functor& f, unsigned int n){
	vector<bool> results(n, false);
	for(int i = 1;i <= n;++ i){
		unsigned int ret = f(i);
		if(results[ret])
			return false;
		results[ret] = true;
	}
	return true;
}
template<class Functor>
bool isOnTo(Functor& f, unsigned int n){
	return isOneToOne(f, n);
}
template<class Functor>
bool isOnTo(const Functor& f, unsigned int n){
	return isOneToOne(f, n);
}
class InverseFunc{
	map<unsigned int, unsigned int> relations;
public:
	template<class Functor>
	InverseFunc(Functor& f, unsigned int n){
	}
	template<class Functor>
	InverseFunc(const Functor& f, unsigned int n){
	}
	unsigned int operator()(unsigned int index){
		map<unsigned int, unsigned int>::iterator pos = relations.find(index);
		if(pos == relations.end())
			return 0;
		return pos->second;
	}
};
class IntUniversalSet{
	unsigned int size;
public:
	IntUniversalSet(unsigned int size) : size(size){}
	int getSize() const{
		return size;
	}
	unsigned int getValue(int index) const{
		return index;
	}
	int getIndex(unsigned int val) const{
		return val;
	}
};
unsigned int func1(unsigned int val){
	return val / 2;
}
struct TempFunc{
	unsigned int max;
	TempFunc(unsigned int lim) : max(lim){}
	unsigned int operator()(unsigned int val) const{
		return max - val + 1;
	}
};
template<class T>
set<pair<T,T> > CartesianProduct(const set<T>& s1, const set<T>& s2){
	typename set<T>::const_iterator begin = s1.begin(), end = s1.end();
	set<pair<T, T> > ret;
	while(begin != end){
		typename set<T>::iterator start = s1.begin(), term = s1.end();
		while(start != term){
			ret.insert(make_pair(*begin, *start));	
			++start;
		}
		++begin;
	}
	return ret;
}
int main(int argc, char *argv[]){
	IntUniversalSet tmp(50);
	set<unsigned int> sets1;
	set<unsigned int> sets2;
	srand(time(0));
	for(int i = 0;i < 50;++ i){
		if(rand() % 3 == 0)
			sets1.insert(i);
		if(rand() % 3 == 0)
			sets2.insert(i);
	}
	Set<unsigned int, IntUniversalSet> set1(&tmp, sets1.begin(), sets1.end());
	Set<unsigned int, IntUniversalSet> set2(&tmp, sets2.begin(), sets2.end());
	set<pair<unsigned int, unsigned int> > results = CartesianProduct(sets1, sets2);
	cout << "S1 : " << set1 << endl << "S2:" << set2 << endl;
	cout << "S1 | S2 : " << (set1 | set2) << endl;
	cout << "S1 & S2 : " << (set1 & set2) << endl;
	cout << "S1 - S2 : " << (set1 - set2) << endl;
	cout << "S1 ^ S2 : " << (set1 ^ set2) << endl;	
	cout << "~S1 : " << (~set1) << endl;
	cout << "func1 is one-to-one:" << isOneToOne(func1, 20) << endl;
	cout << "func1 is onto:" << isOnTo(func1, 20) << endl;
	cout << "func2 is one-to-one:" << isOneToOne(TempFunc(20), 20) << endl;
	cout << "func2 is onto:" << isOnTo(TempFunc(20), 20) << endl;
	cout << "----------------------------------------------" << endl;
	cout << "S1 * S2 : ";
	set<pair<unsigned int, unsigned int> >::iterator begin = results.begin(), end = results.end();
	int i = 0;
	while(begin != end){
		cout << '(' << begin->first << ',' << begin->second << ") ";
		++begin;
		++ i;
		if((i + 1) % 8 == 0)
			cout << endl;
	}
	cout << endl;
	cout << "----------------------------------------------" << endl;
	IntUniversalSet tmp0(6);
	Set<unsigned int, IntUniversalSet> set3(&tmp0, (unsigned int*)NULL, (unsigned int*)NULL);
	for(i = 0;i < 64;++ i){
		cout << '{' << set3 << '}' << endl;
		set3.moveToNext();
	}
	cout << endl;
	
}
