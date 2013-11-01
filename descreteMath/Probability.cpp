#include <cstdlib>
#include <iostream>
#include <cassert>
#include <ctime>
#include <cmath>
#include <iomanip>
using namespace std;
const int count = 100;
long long combinations[count][count];

void calculate(){

	for(int i = 0;i < count;++ i)
		combinations[i][i] = 1;
	for(int i = 0;i < count;++ i)
		combinations[i][0] = 1;
	for(int i = 2;i < count;++ i)
		for(int j = 1;j < i;++ j)
			combinations[i][j] = combinations[i - 1][j - 1] + combinations[i - 1][j];
}

int binormialGen(unsigned int n, float p){
	assert(p > 0);
	unsigned int threshhold = (int)(RAND_MAX * p + 0.5);
	unsigned int k = 0;
	for(unsigned int i = 0;i < n;++ i)
		if(((unsigned int)rand()) <= threshhold)
			++k;
	return k;
}
void test2(int argc, char *argv[]){
	const int cnt = 2000;
	int result[cnt];
	int pCount = 50;
	const long long one = 1;
	const long long target = ((one << pCount) - 1);
	int t = (int)(pCount * (log(pCount) + 0.5772) + 0.5);
	srand(time(0));
	for(int i = 0;i < cnt;++ i){
		long long current = 0;
		int count = 0;
		while(current != target){
			current |= (one << (rand() % pCount));
			++count;
		}
		result[i] = count;
	}
	for(int i = 0;i < cnt;++ i)
		cout << setw(4) << (result[i] - t) << ((i + 1) % 30 == 0 ? "\n" : "");
	if(cnt % 30 != 0)
		cout << endl;
	cout << "expected : " << t << endl;
}
int test1(int argc, char *argv[]){
	const int n = 20;
	int vals[n+1];
	int trialCount = 20000000;
	float p = 0.6;
	srand(time(0));
	calculate();
	for(int i = 0;i < n + 1;++ i)
		vals[i] = 0;
	for(int i = 0;i < trialCount;++ i)
		++vals[binormialGen(n, p)];
	for(int i = 0;i < n + 1;++ i)
		cout << i << "\t:" << vals[i]/(double)trialCount << "\t\t" << combinations[n][i] * pow(p, i) * pow(1 - p, n - i) << endl;
}
int main(int argc, char *argv[]){
	test2(argc, argv);
}
