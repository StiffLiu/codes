#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <algorithm>
#include <iterator>
using namespace std;

class Chipset{
	vector<bool> condition;
public:
	Chipset(unsigned int num, unsigned int badnum = 0) : condition(num, true){
		if(badnum == 0)
			badnum = rand() % (num / 2);
		
		cout << "badnum : " << badnum << endl;
		for(unsigned int i = 0;i < badnum;++ i){
			unsigned int index = rand() % (num - i);
			unsigned int k = 0;
			for(unsigned int j = 0;j < index;++ k){
				if(condition[k])
					++ j;
			}
			while(!condition[k])
				++ k;
			assert(k < num);
			condition[k] = false;
		}
	}
	unsigned int chipCount()const{
		return condition.size();
	}
	bool isGood(unsigned int i)const{
		return condition[i];
	}
	bool test(unsigned int i, unsigned int j)const{
		if(condition[i])
			return condition[j];
		return rand() % 2 == 0;
	}
	friend ostream& operator<<(ostream& os, const Chipset& chipSet){
		copy(chipSet.condition.begin(), chipSet.condition.end(), ostream_iterator<bool>(os, " "));
		return os;
	}
};
class ChipTester{
	const Chipset& chipSet;
	mutable unsigned int iteration;

	unsigned int goodCount(unsigned int* indices, unsigned int n) const{
		unsigned int cnt = 0;
		for(unsigned int i = 0;i < n;++ i)
			if(chipSet.isGood(indices[i]))
				++cnt;
		if(2 * cnt <= n){
			cout << cnt << ',' << n << endl;
			cout << chipSet << endl;
			for(unsigned int i = 0;i < n;++ i)
				cout << indices[i] << ' ';
			cout << endl;
		}
		return cnt;
	}
	unsigned int find(unsigned int* indices, unsigned int n)const{
		assert(goodCount(indices, n) * 2 > n);
		if(n < 1){
			++iteration;
			assert(false);
			return - 1;
		}
		if(n == 2){
			++iteration;
			unsigned int i1 = indices[0];
			unsigned int i2 = indices[1];
			if(chipSet.test(i1, i2))
				return i2;
		}
		unsigned int startIndex = n % 2;
		unsigned int count = n / 2;
		unsigned int k = startIndex;
		for(unsigned int i = 0;i < count;++ i){
			unsigned int index = startIndex + 2 * i;
			unsigned int i1 = indices[index];
			unsigned int i2 = indices[index + 1];
			++iteration;
			if(chipSet.test(i1, i2)){				
				indices[k] = i2;
				++k;
			}
			else if(chipSet.test(i2, i1)){
				indices[k] = i1;
				++k;
			}
		}
		return find(indices, k);
	}
public:
	ChipTester(const Chipset& chipSet) : chipSet(chipSet), iteration(0){
	}
	unsigned int getIteration(){
		return iteration;
	}
	bool test() const{
		unsigned int chipCount = chipSet.chipCount();
		vector<bool> testResult(chipCount, true);
		unsigned int oneGood = 0;
		unsigned int *indices = new unsigned int[chipCount];
		for(unsigned int i = 0;i < chipCount;++ i)
			indices[i] = i;
		iteration = 0;
		oneGood = find(indices, chipCount);
		bool ret = true;
		for(unsigned int i = 0;i < chipCount;++ i){
			if(i != oneGood){
				testResult[i] = chipSet.test(oneGood, i);
			}
			ret = (ret && (testResult[i]  == chipSet.isGood(i)));
		}
		iteration += chipCount;
		delete [] indices;
		if(!ret)
			copy(testResult.begin(), testResult.end(), ostream_iterator<bool>(cout, " "));
		return ret;
	}
};
int main(int argc, char *argv[]){	
	srand(time(0));
	for(int i = 1;i <= 100;++ i){
		unsigned int count = rand() % 80 + 20;
		unsigned int badnum = badnum = rand() % (count / 2);
		if(rand() % 2 == 0){
			badnum = count / 2;
			if(count % 2 == 0 && badnum > 0)
				badnum -= 1;
		}
		cout << "count " << count << endl;
		Chipset chipSet(count, badnum);
		ChipTester tester(chipSet);
		//cout << "exptected :\n\t" << chipSet << endl;
		//cout << "test result :\n\t";
		if(!tester.test())
			cout << endl << "exptected :\n\t" << chipSet << endl;
		cout << "count : " << count << ", badnum : " << badnum << ", ratio : " << (tester.getIteration() / (double)count) << '\t';
		if(i % 8 == 0)
			cout << endl;
	}
	cout << endl;
}
