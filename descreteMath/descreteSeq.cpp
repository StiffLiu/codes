#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <cstdlib>
using namespace std;

int toBinary(unsigned long long val, unsigned int nBit){
	if(nBit > 0){
		const unsigned int len = (sizeof (unsigned long long))*8 - 1;
		if(nBit > len + 1)
			nBit = len;
		for(unsigned int i = 0;i < nBit;++ i){
			unsigned long long tmp = 1;
			tmp <<= (len - i);
			cout << ((tmp & val) == 0 ? '0' : '1');
		}
	}
}

void calculate(unsigned int n){
	if(n == 0)
		return;
	if(n == 1){
		cout << 0 << endl << 1 << endl;
		return ;
	}
	unsigned long long *ps_n_2 = new unsigned long long[1 << n];
	unsigned long long *ps_n_1 = new unsigned long long[1 << n];
	unsigned long long *ps_n = new unsigned long long[1 << n];
	unsigned int c_n_2 = 0;
	unsigned int c_n_1 = 1;
	unsigned int c_n = 2;
	const unsigned int len = (sizeof (unsigned long long)) * 8 - 1;
	ps_n_1[0] = 0;
	ps_n[0] = 0;
	ps_n[1] = 1;
	ps_n[1] <<= len;
	for(unsigned int i = 2;i <= n;++ i){
		unsigned long long tmp = 1;
		unsigned long long *tmp_n = ps_n_2;
		ps_n_2 = ps_n_1;
		ps_n_1 = ps_n;
		ps_n = tmp_n;
		c_n_2 = c_n_1;
		c_n_1 = c_n;
		c_n = c_n_2 + c_n_1;
		tmp <<= (len - i + 2);
		for(unsigned int j = 0;j < c_n_2;++ j){
			ps_n[j] = (ps_n_2[j] | tmp);
		}
		tmp >>= 1;
		for(unsigned int j = 0;j < c_n_1;++ j){
			ps_n[j + c_n_2] =  (ps_n_1[j] | tmp);
		}
	}
	for(unsigned int i = 0;i < c_n;++ i){
		toBinary(ps_n[i], n);
		cout << endl;
	}
	delete []ps_n;
	delete []ps_n_1;
	delete []ps_n_2;	
}

void calc(unsigned long long tmp, unsigned int s, unsigned int t, 
	unsigned int n, vector<unsigned long long>* result){
	const unsigned int len = (sizeof (unsigned long long)) * 8 - 1;
	assert(s >= t);
	if(t == n){
		result->push_back(tmp);
		return;
	}
	if(s < n){
		unsigned long long one = 1;
		one <<= (len - s - t);
		calc(tmp | one, s + 1, t, n, result);
	}
	if(t < s){
		calc(tmp, s, t + 1, n, result);
	}
}
void listFormula(const string& variables){
	if(variables.empty())
		return;
	vector<unsigned long long> result;
	int n = variables.size() - 1;
	const unsigned int len = (sizeof (unsigned long long)) * 8 - 1;
	calc(0, 0, 0, n, &result);
	for(int i = 0;i < result.size();++ i){
		unsigned long long val = result[i];
		int current = 0;
		vector<string> stack;
		stack.push_back(string());
		stack.back().push_back(variables[current]);
		for(int j = 0;j < 2 * n;++ j){
			unsigned long long one = 1;
			one <<= (len - j);
			if(val & one){
				++ current;
				stack.push_back(string());
				stack.back().push_back(variables[current]);
			}else{
				string tmp;
				string ret;
				tmp.swap(stack.back());
				stack.pop_back();
				//if(tmp.size() > 1 || stack.back().size() > 1){
					ret.push_back('('); ret += stack.back();
					ret += tmp; ret.push_back(')'); 
					stack.back().swap(ret);
				/*}else{
					ret = stack.back();
					ret += tmp;
					stack.back().swap(ret);
				}*/
			}
		}
		assert(stack.size() == 1);
		toBinary(val, 2 * n);
		cout << " : " << stack.back() << endl;
	}
}
template<class T>
T calcNth(T c0, T c1, T a0, T a1, unsigned int n){
	switch(n){
		case 0: return a0;
		case 1: return a1;
		default:{
			T a_n_2 = 0;
			T a_n_1 = a0;
			T a_n = a1;
			for(unsigned int i = 2;i <= n;++ i){
				a_n_2 = a_n_1;
				a_n_1 = a_n;
				a_n = c0 * a_n_2 + c1 * a_n_1;
			}
			return a_n;
		}
	}
	assert(false);
}
int main(int argc, char *argv[]){
	for(unsigned int i = 1;i <= 5;++ i){
		calculate(i);
		cout << "-----" << endl;
	}
	string tmp("abcdef");
	listFormula(tmp);
	cout << "-----" << endl;
	cout << calcNth<long long>(1, 1, 1, 1, 100)
		<< ' ' << calcNth<long long>(1, 1, 1, 1, 500) 
		<< ' ' << calcNth<long long>(1, 1, 1, 1, 1000) << endl;
	if(argc > 1)
		cout << calcNth<long long>(1, 1, 1, 1, atoi(argv[1])) << endl;
}
