#include <iostream>
#include <vector>
#include <cassert>
#include <iterator>
#include <algorithm>
#include <stack>
#include <queue>

/*
 * A test of git hub
 */
using namespace std;
/*
 * Illustartion of exercise 1.3.49 in book Algorithms Edition 4
 * This is not a solution of the problem, just a research.
 * The problem is not solved by me now.
 */
int test100(int argc, char *argv[]) {
	const unsigned int start = 3, end = 12;
	for (unsigned int i = start; i <= end; ++i) {
		vector<vector<unsigned int> > stacks;
		vector<unsigned int > seq;
		stacks.resize(3);
		seq.resize(3 * i);
		for (unsigned int j = 0; j < 3; ++j) {
			int index = j * i;
			for(unsigned int k = 0;k < i;++ k)
				stacks[j].push_back(index + k);
		}

		unsigned j = 0;
		unsigned cnt = 0;
		while(cnt < 3 * i){
			vector<unsigned int>& current = stacks[j];
			vector<unsigned int>& next = stacks[(j + 1) % 3];
			vector<unsigned int>& nextNext = stacks[(j + 2) % 3];
			bool needBack = false;
			assert(!current.empty());
			seq[current.back()] = cnt;
			++cnt;
			current.pop_back();
			if(!next.empty()){
				current.push_back(next.back());
				next.pop_back();
				needBack = true;
			}
			if(!next.empty()){
				current.push_back(next.back());
				next.pop_back();
			}
			if(!nextNext.empty()){
				next.push_back(nextNext.back());
				nextNext.pop_back();
			}
			if(needBack){
				next.push_back(current.back());
				current.pop_back();
			}
			j = (j + 1) % 3;
			for(unsigned int k = 0;stacks[j].empty() && k < 3;++ k, j = (j + 1) % 3);
		}
		for (unsigned int j = 0; j < 3; ++j) {
			unsigned int index = j * i;
			for(unsigned int k = 0;k < i;++ k)
				cout << seq[index + k] << '\t';
			cout << endl;
		}
		cout << endl;
	}
	return 0;
}
static bool couldOverflow(bool *ops, unsigned int n){
	int c = 0;
	unsigned int i = 0;
	for(;i < n && c >= 0;++ i)
		c += (ops[i] ? 1 : -1);
	return i != n;
}
static bool isGenerable(int *seq, unsigned int n){
	if(n <= 2)
		return true;
	stack<int> vals;
	for(unsigned int i = 0, k = 0;i < n;){
		if(vals.empty() || vals.top() < seq[i]){
			if(k >= n)
				return false;
			vals.push(k);
			++ k;
		}else if(vals.top() == seq[i]){
			vals.pop();
			++ i;
		}else
			return false;			
	}
	assert(vals.empty());
	return true;
}
static void generate(int *seq, bool *ops, unsigned int n){
	stack<int> vals;
	unsigned int j = 0;
	unsigned int k = 0;
	for(unsigned int i = 0;i < n;++ j){
		if(vals.empty() || rand() % 2 == 0){
			ops[j] = true;
			vals.push(i);
			++ i;
		}else{
			ops[j] = false;
			seq[k] = vals.top();
			vals.pop();
			++ k;
		}		
	}
	while(!vals.empty()){
		ops[j] = false;
		seq[k] = vals.top();
		vals.pop();
		++ j;
		++ k;
	}
}
static int test101(int argc, char *argv[]){
	const int n = 20;
	const int iteration = 10;
	bool ops[2 * n];
	int nums[n];

	for(int k = 0;k < iteration;++ k){
		for(int i = 0;i < n;++ i){
			nums[i] = i;
			ops[i] = true;
			ops[i + n] = false;
		}
		generate(nums, ops, n);
		//std::random_shuffle(nums, nums + n);
		//std::random_shuffle(ops, ops + 2 * n);
		copy(ops, ops + 2 * n, ostream_iterator<int>(cout, " "));
		cout << endl << "could over flow : " << couldOverflow(ops, 2 * n) << endl;
		copy(nums, nums + n, ostream_iterator<int>(cout, " "));
		cout << endl << "could be generated : " << isGenerable(nums, n) << endl;
	}
}
static int test102(int argc, char *argv[]){
	const int m = 100, n = 10;
	std::queue<int> nums;
	for(int i = 0;i < n;++ i)
		 nums.push(i);
	while(nums.size() > 1){
		for(int j = 0;j < m - 1;++ j){
			nums.push(nums.front());
			nums.pop();
		}
		cout << nums.front() << ' ';
		nums.pop();
	}
	cout << "\nsurvied : " << nums.back() << endl;
}
int main(int argc, char *argv[]){
	return test102(argc, argv);
}
