#include <iostream>
#include <cfloat>
#include <vector>
#include <cmath>
using namespace std;
void test(double val, vector<bool>& ret){
	if(val < FLT_MIN)
		return;
	unsigned int count = ret.size();
	unsigned int current = 0;
	unsigned int i = 1;
	if(count == 0)
		return;
	ret[i] = true;
	while(i > 0 && current < count){
		current = i * val;
		ret[current] = true;
		++i;
	}
}
int main(int argc, char *argv[]){
	vector<bool> ret(1000000, false);
	unsigned int count = ret.size();
	test(sqrt(0.9*0.9), ret);
	for(unsigned int i = 0;i < count;++ i)
		if(!ret[i])
			cout << i << ' ';
	cout << endl;
}
