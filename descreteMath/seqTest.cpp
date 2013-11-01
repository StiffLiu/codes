#include <iostream>
#include <iomanip>
using namespace std;
const int count = 2000;
long long combinations[count][count];
bool combinationsOddity[count][count];
int combinationsModulee[count][count];
int modulee = 5;
void calculate(){
	for(int i = 0;i < count;++ i)
		combinations[i][i] = 1;
	for(int i = 0;i < count;++ i)
		combinations[i][0] = 1;
	for(int i = 2;i < count;++ i)
		for(int j = 1;j < i;++ j)
			combinations[i][j] = combinations[i - 1][j - 1] + combinations[i - 1][j];
}

void calculateOddity(){
	for(int i = 0;i < count;++ i)
		combinationsOddity[i][i] = true;
	for(int i = 0;i < count;++ i)
		combinationsOddity[i][0] = true;
	for(int i = 2;i < count;++ i)
		for(int j = 1;j < i;++ j)
			combinationsOddity[i][j] = combinationsOddity[i - 1][j - 1] ^ combinationsOddity[i - 1][j];
}

void calculateModulee(){
	for(int i = 0;i < count;++ i)
		combinationsModulee[i][i] = 1;
	for(int i = 0;i < count;++ i)
		combinationsModulee[i][0] = 1;
	for(int i = 2;i < count;++ i)
		for(int j = 1;j < i;++ j){
			long long tmp = combinationsModulee[i - 1][j - 1];
			tmp += combinationsModulee[i - 1][j];
			combinationsModulee[i][j] = tmp % modulee;
		}
}
void out(int lim){
	for(int i = 0;i < lim;++ i){
		for(int j = 0;j <= i;++ j)
			cout << setw(8) << combinations[i][j];
		cout << endl;
	}
	cout << combinations[72][34] << endl;
}
void outOddity(int lim){
	for(int i = 0;i < lim;++ i){
		cout << setw(3) << i << ':';
		for(int j = 0;j <= i;++ j)
			cout << setw(3) << combinationsOddity[i][j];
		cout << endl;
	}
	cout << combinations[72][34] << endl;
}
int main(int argc, char *argv[]){
	calculate();
	calculateOddity();
	calculateModulee();
	outOddity(70);
	for(int i = 0;i < 10;++ i)
		cout << combinations[3 * i][i] << ' ';
	cout << endl;
	for(int i = 1;i <= count / 2;++ i)
		if(combinationsModulee[2 * i - 1][i - 1] != 0)
			cout << i /*<< ':' << combinationsModulee[2 * i - 1][i - 1]*/ << ' ';
	cout << endl;
}
