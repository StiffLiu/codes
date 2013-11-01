#include <iostream>
#include <map>
using namespace std;

class PseudoRand{
	unsigned int multiplier;
	unsigned int modulus;
	unsigned int shift;
	unsigned int current;
public:
	PseudoRand(unsigned int multiplier,unsigned int modulus,unsigned int shift, unsigned int seed = 0){
		this->multiplier = (multiplier % modulus);
		this->shift = (shift % modulus);
		this->modulus = modulus;
		this->current = seed;
	}
	unsigned int next(){
		current = ((unsigned long long)multiplier * (unsigned long long)current + (unsigned long long)shift) % (unsigned long long)modulus;
		return current;
	}
};

int main(int argc, char *argv[]){
	PseudoRand rand((1 << 13) - 1, (1 << 23) - 1, 4);
	map<unsigned int, unsigned int> counts;
	for(int i = 0;i < 1000;++ i)
		for(int i = 0;i < 90;++ i)
			++counts[rand.next()%100];
	map<unsigned int, unsigned int>::iterator begin = counts.begin(), end = counts.end();
	int count = 0;
	while(begin != end){
		cout << begin->first << ':' << begin->second << endl;
		count += (begin->second - 900)*(begin->second - 900);
		++ begin;
	}
	cout << endl << double(count) / 90000;
}
