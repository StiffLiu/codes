#include <iostream>
#include <vector>
#include <ctime>
#include <cstring>
#include <set>
#include <map>
using namespace std;

class Decompose{
	static const unsigned int maxLen = 600000;
	unsigned int nums[maxLen];
	unsigned int limit;
public:
	Decompose(unsigned int limit){
		if(limit > maxLen)
			limit = maxLen;
		this->limit = limit;
		//nums = new unsigned int[limit];
		memset(nums, 0, limit * sizeof(unsigned int));
		for(unsigned int i = 1;i <= limit;++ i){
			unsigned int st = i * i;
			if(st <= limit)
				nums[st - 1] = 1;
			if(st + 1 <= limit)
				nums[st] = 2;
			if(st + 4 <= limit)
				nums[st + 3] = 2;
		}
		for(unsigned int i = 3;i <= limit;++ i){
			if(nums[i - 1] == 0){
				nums[i - 1] = i;
				for(int j = 1;j < i;++ j){
					unsigned int tmp1 = nums[j - 1] + nums[i - j - 1];
					if(tmp1 == 2){
						nums[i - 1] = tmp1;
						break;
					}
					if(tmp1 < nums[i - 1])
						nums[i - 1] = tmp1;
				}
			}
		}
	}
	~Decompose(){
		//delete [] nums;
	}
	int outGreaterThan(unsigned int val){
		int count = 0;
		for(int i = 1;i <= limit;++ i)
			if(nums[i - 1] > val){
				++count;
				cout << i << ' ';
			}
		return count;
	}
	void test(){
		unsigned int i = 1;
		unsigned int cubic;
		set<unsigned int> cubicSet;
		vector<bool> results(limit, false);
		while((cubic = i * i * i) <= limit){
			cubicSet.insert(cubic);
			++ i;
		}
		cubicSet.insert(0);
		set<unsigned int>::iterator begin = cubicSet.begin(), end = cubicSet.end();
		int count = 0;
		while(begin != end){
			int start = *begin;		
			for(int i = start + 1;i <= limit;++ i){			
				if(nums[i - 1] <= 2){
					results[i - start - 1] = true;
				}
			}
			++begin;
		}
		for(int i = 1;i <= limit;++ i){
			if(!results[i - 1]){
				cout << i << ' ';
				++ count;
				if(count % 10 == 0)
					cout << endl;
				if(count % 100 == 0)
					cout << endl;
			}
		}
		cout << "\n total:" << count << endl;
	}
	static void calc(unsigned int lim){
		vector<unsigned int> all;
		unsigned int cube = 0;
		unsigned int i = 1;
		map<unsigned int, unsigned int> counts;
		while((cube = i * i * i) <= lim){
			all.push_back(cube);
			++i;
		}
		unsigned int count = all.size();
		for(unsigned int j = 0;j < count;++ j){
			for(unsigned int k = j;k < count;++ k)
				++ counts[all[j] + all[k]];
		}
		map<unsigned int, unsigned int>::iterator begin = counts.begin(), end =counts.end();
		count = 0;
		while(begin != end){
			if(begin->second > 1){
				cout << begin->first << ',' << begin->second << ' ';
				++count;
				if(count % 10 == 0)
					cout << endl;
				if(count % 100 == 0)
					cout << endl;
			}
			++begin;
		}
		cout << endl << endl;
		/*begin = counts.begin(), end =counts.end();
		count = 0;
		while(begin != end){
			if(begin->second <= 1){
				cout << begin->first << ',' << begin->second << ' ';
				++count;
			if(count % 10 == 0)
				cout << endl;
			}
			++begin;
		}*/
	}
};
void test1(){
	clock_t start = clock();
	cout << "begin calculate" << endl;
	Decompose test(100000);
	test.test();
	cout << "\ndone" << endl;
	cout << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
}
void test2(){
	clock_t start = clock();
	cout << "begin calculate" << endl;
	Decompose::calc(500000000);
	cout << "\ndone" << endl;
	cout << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
}
int main(int argc,char * argv[]){
	test2();
}
