#include <iostream>
#include <vector>
using namespace std;

class Chomp{
	unsigned int rows;
	vector<bool> cookies;
	vector<vector<int> > allSteps;
public:
	Chomp(unsigned int _rows, unsigned int _cols) : rows(_rows), cookies(_rows * _cols, false){}
	void calc(){
		vector<int> steps;
		calc(steps);
	}
	void out(){
		int stepCount = 0;
		for(int i = 0;i < allSteps.size();++ i){
			vector<int>& tmp = allSteps[i];
			if(tmp.size() < allSteps[stepCount].size()){
				stepCount = i;
			}
			for(int j = 0;j < tmp.size();++ j)
				cout << tmp[j] << ' ';
			cout << endl;
		}
		cout << "----------------all possible steps:" << allSteps.size() << "----------------" << endl;
		cout << "----------------one of the min steps----------------" << endl;
		vector<int>& tmp = allSteps[stepCount];
		for(int j = 0;j < tmp.size();++ j)
			cout << tmp[j] << ' ';
		cout << endl;
	}
private:
	void calc(vector<int>& steps){
		unsigned int cols = (unsigned int)cookies.size() / rows;
		bool has = false;
		//cout << cols << ',' << rows << endl;
		for(unsigned int i = 0;i < rows - 1;++ i)
			for(unsigned int j = 0;j < cols - 1;++ j){
				int cur = i * cols + j;
				if(cur != 0 && !cookies[cur]){
					int right = cur + 1;
					int down = cur + cols;
					bool rOld = cookies[right];
					bool dOld = cookies[down];
					cookies[cur] = true;
					cookies[right] = true;
					cookies[down] = true;
					steps.push_back(cur);
					calc(steps);
					steps.pop_back();
					cookies[cur] = false;
					cookies[right] = rOld;
					cookies[down] = dOld;					 
					has = true;
				}
			}
		if(rows > 0 && cols > 0){
			int start = (rows - 1) * cols;
			for(unsigned int j = 0;j < cols - 1;++ j){
				int cur = start + j;
				if(cur != 0 && !cookies[cur]){
					int right = cur + 1;
					bool rOld = cookies[right];
					cookies[cur] = true;
					cookies[right] = true;
					steps.push_back(cur);
					calc(steps);
					steps.pop_back();
					cookies[cur] = false;
					cookies[right] = rOld;				 
					has = true;
				}
			}			
			for(unsigned int j = 0;j < rows - 1;++ j){
				int cur = j * cols + cols - 1;
				if(cur != 0 && !cookies[cur]){
					int down = cur + cols;
					bool dOld = cookies[down];
					cookies[cur] = true;
					cookies[down] = true;
					steps.push_back(cur);
					calc(steps);
					steps.pop_back();
					cookies[cur] = false;
					cookies[down] = dOld;				 
					has = true;
				}
			}
			int cur = cols * rows - 1;
			if(cur != 0 && !cookies[cur]){
				cookies[cur] = true;
				steps.push_back(cur);
				calc(steps);
				steps.pop_back();
				cookies[cur] = false;			 
				has = true;
			}
		}
		if(!has)
			allSteps.push_back(steps);
	}
};
int main(int argc, char *argv[]){
	Chomp temp(3, 4);
	temp.calc();
	temp.out();
}
