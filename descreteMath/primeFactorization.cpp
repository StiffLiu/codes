#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
using namespace std;

void findPrime(vector<bool>& ret){
	unsigned int n = ret.size();
	int i = 0;
	if(n < 1)
		return;
	ret[0] = false;
	if(n < 2)
		return;
	ret[1] = true;
	while(i < n){
		while(i < n && !ret[i]) ++i;
		for(unsigned int j = 2;(unsigned long long)j * (unsigned long long)(i + 1) <= (unsigned long long)n;++ j)
			ret[j * (i + 1) - 1] = false;
		++i;
	}
}
unsigned int allNumbers[105097565];
unsigned int findAll(){
	unsigned int total = -1;total >>= 1;
	vector<bool> tmp(total, true);
	unsigned int count = 0;
	findPrime(tmp);
	for(int i = 0;i < total;++ i)
		if(tmp[i]){
			allNumbers[count]=(i+1);
			++count;
		}
	return count;
}
template<class T>
int factorize(T n, unsigned int* primes, T *ret){
	int count = 0;
	unsigned long long prime = 2;
	for(T i = 1;prime * prime <= n;){
		if(n % prime == 0){
			ret[count] = prime;
			n /= prime;
			++count;
		}else{
			prime = primes[i];
			++i;
		}
	}
	ret[count] = n;
	return count + 1;
}
template<class T>
int factorize(T n, vector<bool>& primes, T *ret){
	int count = 0;
	for(T i = 2;i * i <= n;){
		if(primes[i - 1] && n % i == 0){
			ret[count] = i;
			n /= i;
			++count;
		}else
			++i;
	}
	ret[count] = n;
	return count + 1;
}
void test(int bit, unsigned int* primes){
	unsigned long long tmp = 1;
	unsigned long long ret[100];
	unsigned long long nums[1000];
	int total = (sizeof nums / sizeof *nums);
	tmp <<= bit;
	for(int i = 0;i < total;++ i){
		unsigned long long num = rand() * rand();
		num &= (tmp - 1);num |= (tmp);
		num |= 1;
		nums[i] = num;
	}
	cout << "bit:" << bit << endl;
	clock_t start = clock();
	for(int i = 0;i < total;++ i){
		int count = factorize(nums[i], primes, ret);
		if (i < 10){
			cout << nums[i] << ':';
			for(int i = 0; i < count;++ i)
				cout << ret[i] << ' ';
			cout << endl;
		}
	}
	cout << "calc time: " << (double)(clock() - start)/CLOCKS_PER_SEC << endl;
}
void test1(unsigned int *primes){
	unsigned long long ret[100];
	unsigned long long nums[100];
	unsigned long long tmp = 1;
	int total = (sizeof nums / sizeof *nums);
	for(int i = 0;i < total;++ i)
		nums[i] = 1099511627689uLL;
	cout << "prime:" << nums[0] << endl;
	clock_t start = clock();
	int count = 0;
	for(int i = 0;i < total;++ i){
		count = factorize(nums[i], primes, ret);
		/*cout << nums[i] << ':';
		cout << endl;*/
	}
	cout << "calc time: " << (double)(clock() - start)/CLOCKS_PER_SEC << endl;
	for(int i = 0; i < count;++ i)
		cout << ret[i] << ' ';
	cout << endl;
}
void test(int bit, vector<bool>& primes){
	unsigned long long tmp = 1;
	unsigned long long ret[100];
	unsigned long long nums[1000];
	int total = (sizeof nums / sizeof *nums);
	tmp <<= bit;
	for(int i = 0;i < total;++ i){
		unsigned long long num = rand() * rand();
		num &= (tmp - 1);num |= (tmp);
		num |= 1;
		nums[i] = num;
	}
	cout << "bit:" << bit << endl;
	clock_t start = clock();
	for(int i = 0;i < total;++ i){
		int count = factorize(nums[i], primes, ret);
		/*cout << nums[i] << ':';
		for(int i = 0; i < count;++ i)
			cout << ret[i] << ' ';
		cout << endl;*/
	}
	cout << "calc time: " << (double)(clock() - start)/CLOCKS_PER_SEC << endl;
}
void test1(vector<bool>&primes){
	unsigned long long ret[100];
	unsigned long long nums[100];
	unsigned long long tmp = 1;
	int total = (sizeof nums / sizeof *nums);
	for(int i = 0;i < total;++ i)
		nums[i] = ((tmp << 37) + 1);
	cout << "prime:" << nums[0] << endl;
	clock_t start = clock();
	int count = 0;
	for(int i = 0;i < total;++ i){
		count = factorize(nums[i], primes, ret);
		/*cout << nums[i] << ':';
		cout << endl;*/
	}
	cout << "calc time: " << (double)(clock() - start)/CLOCKS_PER_SEC << endl;
	for(int i = 0; i < count;++ i)
		cout << ret[i] << ' ';
	cout << endl;
}
void findGreatest(unsigned int *primes){
	unsigned long long num = 1;
	unsigned long long count = 0;
	unsigned long long ret[100];
	num <<= 40;num -= 1;
	clock_t start = clock();
	while(factorize(num, primes, ret) > 1){
		--num;++count;
	}
	cout << "calc time: " << (double)(clock() - start)/CLOCKS_PER_SEC 
		<< ", counts:" << count << ", result:" << num << endl;
}
int main(int argc, char *argv[]){
	if(false){
		clock_t start = clock();
		cout << findAll() << endl;
		cout << "calc time: " << (clock() - start) << endl;
		start = clock();
		if(true){
			ofstream out("primes");
			out.write((const char*)allNumbers, sizeof allNumbers);
		}
		cout << "write time: " << (clock() - start) << endl;
	}
	srand(time(0));
	{
		vector<bool> primes((1<<21), true);
		int primeCount = 0;
		unsigned int *primeArr = NULL;
		findPrime(primes);
		for(int i = 0;i < primes.size();++ i)
			if(primes[i])
				++ primeCount;
		primeArr = new unsigned int[primeCount];
		primeCount = 0;
		for(int i = 0;i < primes.size();++ i)
			if(primes[i])
				primeArr[primeCount++] = (i+1);
		for(int i = 0;i < 10;++ i)
			cout << primeArr[i] << ' ';
		cout << endl;
		test(20, primes);
		test(30, primes);
		test(40, primes);
		test1(primes);
		cout << "--------------------" << endl;
		test(20, primeArr);
		test(30, primeArr);
		test(40, primeArr);
		test1(primeArr);
		cout << "--------------------" << endl;
		findGreatest(primeArr);
		delete [] primeArr;
	}
}
