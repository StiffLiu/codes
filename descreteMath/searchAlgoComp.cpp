#include <iostream>
#include <cstdlib>
#include <ctime>
using namespace std;

bool ternarySearch(int *a, int n, int x){
	int i = 1;
	int j = n;
	while(i < j){
		//cout << i  << ',' << j << endl;
		int t1 = (2*i+j)/3;
		int t2 = (i+2*j)/3;
		if(a[t1-1] == x){
			i = t1;
			break;
		}
		if(a[t2-1] == x){
			i = t2;
			break;
		}
		if(x<a[t1-1])
			j = t1-1;
		else if(a[t2-1]<x)
			i = t2+1;
		else{
			i = t1+1;
			j = t2-1;
		}
	}
	//cout << i  << ',' << j << endl;
	return (a[i-1]==x);
}

bool quadralSearch(int *a, int n, int x){
	int i = 1;
	int j = n;
	while(i < j){
		//cout << i  << ',' << j << endl;
		int t1 = (3*i+j)/4;
		int t2 = (i+j)/2;
		int t3 = (i+3*j)/4;
		if(a[t1-1]==x){
			i = t1;
			break;
		}
		if(a[t2-1]==x){
			i = t2;
			break;
		}
		if(a[t3-1]==x){
			i = t3;
			break;
		}
		if(x<a[t1-1])
			j = t1-1;
		else if(a[t3-1]<x)
			i = t3+1;
		else if(a[t2-1]<x){
			i = t2+1;
			j = t3-1;
		}else{
			i = t1+1;
			j = t2-1;
		}
	}
	//cout << i  << ',' << j << endl;
	return (a[i-1]==x);
}

bool binarySearch(int *a, int n, int x){
	int i = 1;
	int j = n;
	while(i < j){
		//cout << i  << ',' << j << endl;
		int t1 = (i+j)/2;
		if(a[t1-1]==x){
			i = t1;
			break;
		}
	if(x<a[t1-1])
		j = t1-1;
	else
		i = t1+1;
	}
	//cout << i  << ',' << j << endl;
	return (a[i-1]==x);
}
int main(int argc, char *argv[]){
	int n = 1000000;
	int *vals = new int[n];
	int inc = 1;
	double time1 = 0, time2 = 0, time3 = 0;
	int cycle=10000;
	int round = 1000;
	for(int i = 0;i < n;++ i){
		vals[i] = 2 * i;
	}
	srand(time(0));
	for(int i = 0;i < round;++ i){
		int target = rand();
		clock_t start = clock();
		for(int j = 0;j < cycle;++j)
			binarySearch(vals, n, target);
		time1 += (clock() - start);
		start = clock();
		for(int j = 0;j < cycle;++ j)
			ternarySearch(vals, n, target);
		time2 += (clock() - start);
		start = clock();
		for(int j = 0;j < cycle;++ j)
			quadralSearch(vals, n, target);
		time3 += (clock() - start);
	}
	cout << "binary average: " << time1 / (round * cycle) << endl;
	cout << "ternary average: " << time2 / (round * cycle) << endl;
	cout << "quadral average: " << time3 / (round * cycle) << endl;
}
