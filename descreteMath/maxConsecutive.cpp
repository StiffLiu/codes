#include <iostream>
using namespace std;

template<class T>
T maxConsecutive(T *vals, int n, int& startIndex, int& endIndex){
	if(n < 1)
		return T();
	T curMax = vals[0];
	T zero = T();
	T curVal = vals[0];
	startIndex = 0;
	endIndex = 0;
	for(int i = 1;i < n;++ i){
		if(curVal < 0){
			startIndex = i;
			curVal = vals[i];
		}
		else
			curVal += vals[i];
		if(curVal > curMax){
			endIndex = i;
			curMax = curVal;
		}
	}
	return curMax;
}

int main(int arg, char *argv[]){
	int vals1[] = {-2, 4, -1, 3, 5, -6, 1, 2};
	int vals2[] = {4, 1, -3, 7, -1, -5, 3, -2};
	int vals3[] = {-1, 6, -3, -4, -5, 8, -1, 7};
	int cnt1 = (sizeof vals1 / sizeof *vals1);
	int cnt2 = (sizeof vals2 / sizeof *vals2);
	int cnt3 = (sizeof vals3 / sizeof *vals3);
	int s1, e1, s2, e2, s3, e3;
	int max1 = maxConsecutive(vals1, cnt1, s1, e1),
		max2 = maxConsecutive(vals2, cnt2, s2, e2),
		max3 = maxConsecutive(vals3, cnt3, s3, e3);
	cout << max1 << " : " << s1 << "---" << e1 << endl
		<< max2 << " : " << s2 << "---" << e2 << endl
		<< max3 << " : " << s3 << "---" << e3 << endl;
	return 0;
}
