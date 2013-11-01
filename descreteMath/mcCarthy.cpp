#include <iostream>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <map>
#include <set>
#include <vector>
using namespace std;
map<int, int> dependence;
int m = 100, n = 1;
void set(int m, int n){
	assert(n > 0);
	assert(m >= 0);
	for(int i = 1;i <= 100;++ i){
		int val = i + n;
		assert(val >= 0);
		if(val <= 100){
			dependence[i] = val;
		}
		else{
			val -= m;
			if(val <= 100){
				dependence[i] = val;
			}		
		}
	}
}
void init(){
	dependence.clear();
	::set(m, n);
}


bool check(bool output = true){
	std::set<int> all;
	std::vector<int> indices;
	for(int i = 1;i <= 100;++ i){
		int current = i;
		map<int,int>::iterator pos;
		bool findCycle = false;
		all.clear();
		if(output){
			indices.clear();
		}
		while((pos = dependence.find(current)) != dependence.end()){
			if(pos->second == i || all.find(pos->second) != all.end()){
				findCycle = true;
				break;
			}
			all.insert(pos->second);
			if(output)
				indices.push_back(pos->second);
			current = pos->second;
		}
		if(output){
			std::vector<int>::iterator begin = indices.begin(), end = indices.end();
			cout << i << ':';
			while(begin != end){
				cout << *begin << ' ';
				++begin;
			}
			if(findCycle)
				cout << " cycular find";
			cout << endl;
		}
		if(findCycle)
			return false;
	}
	return true;
}
template<int n>
void out(bool ret[n][n]){
	for(int i = 0;i < n;++ i){
		for(int j = 0;j < n;++ j)
			cout << ret[i][j] << '\t';
		cout << endl;
	}
}
template<int n>
void outPair(bool ret[n][n]){
	int count = 0;
	for(int i = 0;i < n;++ i){
		for(int j = 0;j < n;++ j)
			if(ret[i][j]){
				cout << '(' << (i + 1) << ',' << (j + 1) << ')' << ' ';
				++count;
				if(count % 8 == 0){
					cout << endl;
				}
			}
	}
}
int main(int argc, char *argv[]){
	string str(60, '-');
	bool ret[10][10];
	for(int i = 1;i <= 10;++ i){
		for(int j = 1;j <= 10;++ j){
			m = i;
			n = j;
			init();
			ret[i - 1][j - 1] = check(false);
		}
	}
	out(ret);
	//outPair(ret);
	/*m = 10,n = 11;
	init();
	cout << check() << endl;*/
}

