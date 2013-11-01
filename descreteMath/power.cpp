#include <cassert>
template<class T>
T power(T val, unsigned int p){
	T ret = 1;
	T cp = val;
	while(p != 0){
		if(p % 2 == 1) ret *= cp;
		cp *= cp;
		p >>= 1;
	}
	return ret;
}

#include <iostream>
int main(int argc, char *argv[]){
	assert(power(2, 5) == 32);
	assert(power(2, 10) == 1024);
	assert(power(3, 7) == 2187);
	std::cout << power(0.5, 10) << std::endl;
}
