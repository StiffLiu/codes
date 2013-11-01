#include <algorithm>
#include <iostream>
#include <cassert>
#include <iterator>
using namespace std;

template<class T>
bool nextPermutation(T *vals, int n){
	int i = n - 2;
	for(;i >= 0 && vals[i] >= vals[i + 1];-- i);
	if(i < 0){
		reverse(vals + i + 1, vals + n);
		return false;
	}

	int k = n - 1;
	for(;vals[k] <= vals[i];-- k);
	swap(vals[k], vals[i]);	
	reverse(vals + i + 1, vals + n);
	return true;
}
template<class T>
bool nextCombination(T *vals, int m, int n){
	assert(m <= n);
	int i = m - 1;
	while(i >= 0 && vals[i] == (n - m + i + 1)){
		assert(vals[i] <= n);
		if(i > 0)
			assert(vals[i] > vals[i - 1]);
		--i;
	}
	if(i < 0){
		for(int i = 0;i < m;++ i)
			vals[i] = (i + 1);
		return false;
	}
	T tmp = vals[i];++tmp;
	while(i < m){
		vals[i] = tmp;
		++tmp;
		++ i;
	}
	return true;
}
template<class T>
bool nextCombinationWithRepetition(T *vals, int m, int n){
	assert(m <= n);
	int i = m - 1;
	while(i >= 0 && vals[i] == n){
		assert(vals[i] <= n);
		if(i > 0)
			assert(vals[i] >= vals[i - 1]);
		--i;
	}
	if(i < 0){
		for(int i = 0;i < m;++ i)
			vals[i] = 1;
		return false;
	}
	T tmp = vals[i];++tmp;
	while(i < m){
		vals[i] = tmp;
		++ i;
	}
	return true;
}
template<class T>
bool nextRPermutation(T *vals, int m, int n){
	if(nextPermutation(vals, m))
		return true;
	return nextCombination(vals, m, n);
}
void listAllSolutions(int varSize, int result){
	if(varSize <= 0 || result < 0)
		return;
	if(varSize == 1){
	}
	int n = varSize + result - 1;
	int *tmp = new int[varSize - 1];
	int current = 0;
	for(int i = 1;i <= varSize;++ i)
		tmp[i - 1] = i;
	do{
		int prev = 0;
		for(int i = 0;i < varSize - 1;++ i){
			cout << (tmp[i] - prev - 1) << " ";
			prev = tmp[i];
		}
		cout << (n - prev) << " ";cout << "\t";
		++ current;
		if(current % 8 == 0)
			cout << endl;
	}while(nextCombination(tmp, varSize - 1, n));
	cout << (current % 8 != 0 ? "\n" : "") << "total: " << current << endl;
	return;
}
int countWin(int m){
	if(m <= 0)
		return 0;
	struct Temp{
		static void fun(int m, int fWin, int sWin, int& total){
			assert(fWin <= m && sWin <= m);
			if(fWin == m){
				++total;
				return;
			}
			if(sWin == m){
				++total;
				return;
			}
			fun(m, fWin + 1, sWin, total);
			fun(m, fWin, sWin + 1, total);
		}
	};
	int ret = 0;
	Temp::fun(m, 0, 0, ret);
	return ret;
}

void output(unsigned int val, int cnt){
	for(int i = 0;i < cnt; ++ i)
		cout << ((val & (1 << i)) != 0 ? '1' : '0');
}
int test1(int argc, char *argv[]){
	int allVals[] = {1, 2, 3, 4};
	int count = (sizeof allVals / sizeof *allVals);
	int current = 0;
	cout << "-------------permutation without repitition : P(4,4)-----------------" << endl;
	do{
  		copy (allVals, allVals + count, ostream_iterator<int>(cout," ")); cout << '\t';
		++ current;
		if(current % 8 == 0)
			cout << endl;
	}while(nextPermutation(allVals, count));
	cout << (current % 8 != 0 ? "\n" : "") << "total: " << current << endl;
	/**************************************************************/
	int allVals1[] = {1, 1, 2, 2, 3};
	int count1 = (sizeof allVals1 / sizeof *allVals1);
	cout << "-------------permutation with repitition : P(5;2,2,1)-----------------" << endl;
	current = 0;
	do{
  		copy (allVals1, allVals1 + count1, ostream_iterator<int>(cout," ")); cout << '\t';
		++ current;
		if(current % 8 == 0)
			cout << endl;
	}while(nextPermutation(allVals1, count1));
	cout << (current % 8 != 0 ? "\n" : "") << "total: " << current << endl;
	/**************************************************************/
	int comVals1[] = {1, 2, 3};
	int count2 = (sizeof comVals1 / sizeof *comVals1);
	cout << "-------------combination without repitition : C(7,3)-----------------" << endl;
	current = 0;
	do{
  		copy (comVals1, comVals1 + count2, ostream_iterator<int>(cout," ")); cout << '\t';
		++ current;
		if(current % 8 == 0)
			cout << endl;
	}while(nextCombination(comVals1, count2, 7));
	cout << (current % 8 != 0 ? "\n" : "") << "total: " << current << endl;
	/**************************************************************/
	int comVals2[] = {1, 1, 1, 1};
	int count3 = (sizeof comVals2 / sizeof *comVals2);
	cout << "-------------combination with repitition : C(4 + 4 - 1, 4 - 1)-----------------" << endl;
	current = 0;
	do{
  		copy (comVals2, comVals2 + count3, ostream_iterator<int>(cout," ")); cout << '\t';
		++ current;
		if(current % 8 == 0)
			cout << endl;
	}while(nextCombinationWithRepetition(comVals2, count3, 4));
	cout << (current % 8 != 0 ? "\n" : "") << "total: " << current << endl;
	/**************************************************************/
	int comVals3[] = {1, 2, 3};
	int count4 = (sizeof comVals3 / sizeof *comVals3);
	cout << "-------------r permutation without repitition : P(6, 3)-----------------" << endl;
	current = 0;
	do{
  		copy (comVals3, comVals3 + count4, ostream_iterator<int>(cout," ")); cout << '\t';
		++ current;
		if(current % 8 == 0)
			cout << endl;
	}while(nextRPermutation(comVals3, count4, 6));
	cout << (current % 8 != 0 ? "\n" : "") << "total: " << current << endl;
	/**************************************************************/
	cout << "-------------all solutions to : x1+x2+x3+x4=7-----------------" << endl;
	listAllSolutions(4, 7);	
	/**************************************************************/
	cout << "-------------winning count:1-9-----------------" << endl;
	for(int i = 1;i < 9; ++i)
		cout << i << ':' << countWin(i) << ' ';
	cout << endl;
	/**************************************************************/
	int count5 = 8;
	int comVals4Temp[] = {1, 2, 3, 4, 5, 6};
	int cnts = (sizeof comVals4Temp / sizeof *comVals4Temp);
	current = 0;
	do{
  		unsigned int tmp = -1;
		for(int i = 0;i < cnts;i += 2){
			for(int j = comVals4Temp[i] - 1;j < comVals4Temp[i + 1] - 1;++ j)
				tmp &= ~(1 << j);
		}
		output(tmp, count5);cout << ' ';
		++ current;
		if(current % 8 == 0)
			cout << endl;
	}while(nextCombination(comVals4Temp, cnts, count5 + 1));
	cout << (current % 8 != 0 ? "\n" : "") << "total: " << current << endl;
	return 0;
}
#include <sstream>
void inclusionExclusion(int n){
	int *vals = new int[n];
	ostringstream formula;
	for(int i = 1;i <= n;++ i){
		for(int j = 0;j < i;++ j)
			vals[j] = (j + 1);
		do{
			if(i % 2 == 1){
				if(!formula.str().empty())
					formula << '+';
			}else{
				formula << '-';
			}
			formula << '|';
			for(int j = 0;j < i;++ j){
				if(j != 0)
					formula << 'U';
				formula << 'A' << vals[j];
			}
			formula << '|';
		}while(nextCombination(vals, i, n));
		formula << endl;
	}
	delete []vals;
	cout << formula.str() << endl;
}
int test2(int argc, char *argv[]){
	inclusionExclusion(6);
	return 0;
}

int main(int argc, char *argv[]){
	return test2(argc, argv);
}
