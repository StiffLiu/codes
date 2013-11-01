#include<iostream>
#include <vector>

using namespace std;
int count = 0;
unsigned long long archerman(unsigned long long m, unsigned long long n, 
	std::vector<unsigned long long>* relies = NULL){
	//A(m,n) = 2 * n, if m == 0
	if(m == 0){
		return 2 * n;
	}
	//A(m,n) = 0, if m >= 1 and n == 0
	if(n == 0){
		return 0;
	}
	//A(m,n) = 2, if m >= 1 and n == 1
	if(n == 1){
		return 2;
	}
	if(n == 2){
		return 4;
	}
	if(m == 1){
		return (1 << n);
	}
	//A(m,n)=A(m-1,A(m,n-1)), if m >= 1 and n >= 2
	if(relies != NULL){
		//relies->push_back(m - 1);
		if(n - 1 == 1){
			//cout << relies->size() << endl;
			//relies->push_back(m);relies->push_back(n - 1);
		}
	}
	unsigned long long ret = archerman(m, n - 1, relies);
	ret = archerman(m - 1, ret, NULL);
	return ret;
}

int main(int argc, char *argv[]){
	vector<unsigned long long> rets;
	for(int i = 2;i <=4;++ i){
		for(int j = 2;j <= 4;++ j)
			cout << "A(" << i << ',' << j << ")=" << archerman(i, j, NULL) << endl;
	}
	cout << archerman(3,4);
	//cout << ret << endl;
	/*for(int i = 0;i < rets.size();++ i)
		std::cout << rets[i] << ' ';*/	
}
