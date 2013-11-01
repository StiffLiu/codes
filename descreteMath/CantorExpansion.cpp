#include<iostream>
#include<string>
#include <vector>
#include<algorithm>
using namespace std;
typedef vector<int> RetType;
void cantorExpansion(int n, RetType& ret){
	int i = 2;
	if(n == 0){
		ret.push_back(0);
		return;
	}
	while(n != 0){
		ret.push_back(n%i);
		n /= i;
		++i;
	}
}
int factorial(int n){
	int ret = 1;
	for(int i = 2;i <= n;++ i)
		ret *= i;
	return ret;
}
void add(const RetType& a, const RetType& b, RetType& ret){
	int carry = 0;
	size_t i = 0;
	for(;i < a.size() && i < b.size();++ i){
		int current = a[i] + b[i] + carry;
		carry = current / (int)(i + 2);
		ret.push_back(current % (int)(i + 2));
	}
	for(;i < a.size();++ i){
		int current = a[i] - '0' + carry;
		carry = current / (int)(i + 2);
		ret.push_back(current % (int)(i + 2));
	}
	for(;i < b.size();++ i){
		int current = b[i] - '0' + carry;
		carry = current / (int)(i + 2);
		ret.push_back(current % (int)(i + 2));
	}
	if(carry != 0){
		ret.push_back(carry);
	}
}
ostream& out(ostream& os, RetType& ret){
	for(size_t i = 0;i < ret.size();++ i)
		cout << ret[i] << "*factorial("<<(i+1)<<")" << ((i+1)<ret.size()?"+":"");
	return os;
}
int test(int argc, char *argv[]){
	{
		RetType a,b,ret;
		cantorExpansion(1000000, a);
		cantorExpansion(10000, b);
		add(a, b, ret);
		out(cout, ret) << endl;
	}
	{
		RetType a, b, ret, c;
		for(int i = 1;i <= 6;++ i)
			a.push_back(i);
		b.push_back(1);
		add(a, b, ret);
		out(out(out(cout, a) << endl, b) << endl, ret) << endl;
		add(a, a, c);
		out(cout, c) << endl;
		cout << endl;
	}
	
}
int iteratePermutation(int n){
	int total = factorial(n);
	RetType ret;
	vector<int> values;
	values.resize(n);
	for(int i = 0;i < total;++ i){
		int count = 0;
		int j = 1;
		ret.clear();
		values[0]= 1;
		cantorExpansion(i, ret);
		count = ret.size();
		//cout << count << endl;
		while(j <= count){
			int cnts = ret[j - 1];
			values.insert(values.begin() + j - cnts, j + 1);
			++j;
		}
		//cout << count << endl;
		while(j <= n - 1){
			values[j] = (j + 1);
			++ j;
		}
		//cout << count << endl;
		for(int k = 0;k < n;++ k)
			cout << values[k] << ' ';
		cout << endl;
	}
}
int main(int argc, char *argv[]){
	/*if(argc > 1){
		RetType ret;
		int tmp = atoi(argv[1]);
		cantorExpansion(abs(tmp), ret);
		out(cout, ret) << endl;
	}*/
	iteratePermutation(4);
}
