#include <iostream>
#include <vector>
using namespace std;

void findPrime(vector<bool>& ret){
	int n = ret.size();
	int i = 0;
	if(n < 1)
		return;
	ret[0] = false;
	if(n < 2)
		return;
	ret[1] = true;
	while(i < n){
		while(i < n && !ret[i]) ++i;
		for(int j = 2;j * (i + 1) <= n;++ j)
			ret[j * (i + 1) - 1] = false;
		++i;
	}
}
void out(unsigned int n, vector<bool>& tmp){
	unsigned int total = tmp.size();
	int count = 0;
	for(unsigned int i = 0;n*i < total;++i)
		if(tmp[n*i]){
			cout << (n*i + 1) << ' ';
			++count;
			if(count % 8 == 0)
				cout << endl;
		}
}
void out(unsigned int n, unsigned int *vals, vector<bool>& tmp){
	unsigned int total = tmp.size();
	for(unsigned int i = 0;n*i < total;++i){
		bool isAllPrime = true;
		for(unsigned int j = 0;j < n;++ j){
			isAllPrime &= tmp[vals[j]*i];
		}
		if(isAllPrime){
			for(unsigned int j = 0;j < n;++ j){
				cout << (vals[j]*i + 1) << ((j < n-1)? ", " : "\n");
			}
		}
	}
}
int main(int argc, char *argv[]){
	unsigned int total = 1000;
	vector<bool> tmp(total, true);
	unsigned int vals[]={6, 12, 18};
	int count = 0;
	clock_t start = clock();
	findPrime(tmp);
	clock_t end = clock();
	//cout << "\n time exhuasted:" << (end - start) << ", count=" << count << endl;
	out(3, vals, tmp);	
}
