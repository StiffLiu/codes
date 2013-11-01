#include <iostream>
using namespace std;

void hanoiMove(int start, int n, char s, char d, char m){
	if(n == 1){
		cout << (start + n) << ' ' << s << ' ' << d << endl;
		return;
	}
	hanoiMove(start, n - 1, s, m, d);
	cout << (start + n) << ' ' << s << ' ' << d << endl;
	hanoiMove(start, n - 1, m, d, s);
}

void frameStewart(int n, char s, char d, char m1, char m2){
	if(n == 1){
		cout << n << ' ' << s << ' ' << d << endl;
		return;
	}
	int k = n / 2;
	frameStewart(n - k, s, m1, m2, d);
	hanoiMove(n - k, k, s, d, m2);
	frameStewart(n - k, m1, d, s, m2);	
}

int main(int argc, char *argv[]){
	hanoiMove(0, 4, 'a', 'c', 'b');
	cout << "----------------------------" << endl;
	frameStewart(4, 'a', 'd', 'b', 'c');
}
