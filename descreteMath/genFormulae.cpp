#include <iostream>
#include <vector>
#include <string>
using namespace std;
vector<string> symbols;
vector<int> indices;
char *opers=NULL;
char *unary = "!";
char *arith = "+-*/";
char *logic = "&|>=";
void genNew(const string& op1, char op, const string& op2, string& ret){
	ret.push_back('(');ret += op1;ret.push_back(')');
	ret.push_back(op);
	ret.push_back('(');ret += op2;ret.push_back(')');
}
void initArith(){
	opers = arith;
	//Step First
	symbols.push_back("x");
	symbols.push_back("y");
	symbols.push_back("z");
	indices.push_back(0);
}
void initLogic(){
	opers = logic;
	//Step First
	symbols.push_back("p");
	symbols.push_back("q");
	symbols.push_back("F");
	symbols.push_back("T");
	indices.push_back(0);
}
void generate(int n){
	if(symbols.empty())
		return;
	//Step Second to Fifth
	for(int i = 1;i < n;++ i){
		int startIndex = indices[i - 1];
		int endIndex = symbols.size();
		indices.push_back(symbols.size());
		for(int j = startIndex;j < endIndex;++ j){
			if(symbols[j].empty())
				continue;
			for(int k = 0;k < startIndex;++ k){
				if(symbols[k].empty())
					continue;
				for(char * s = opers;*s;++ s){
					symbols.push_back(string());
					genNew(symbols[j], *s, symbols[k], symbols.back());
					symbols.push_back(string());
					genNew(symbols[k], *s, symbols[j], symbols.back());
				}
				symbols.push_back(string());
			}
		}
		for(int j = startIndex;j < endIndex;++ j){
			if(symbols[j].empty())
				continue;
			if(unary != NULL){
				for(char * s = unary;*s;++ s){
					symbols.push_back(string());
					symbols.back().push_back(*s);symbols.back().push_back('(');
					symbols.back() += symbols[j];symbols.back().push_back(')');
				}
			}
			for(int k = startIndex;k < endIndex;++ k){
				if(symbols[k].empty())
					continue;
				for(char * s = opers;*s;++ s){
					symbols.push_back(string());
					genNew(symbols[j], *s, symbols[k], symbols.back());
				}
				symbols.push_back(string());
			}
		}		
	}
}
int main(int argc, char *argv[]){
	int start = 0;
	int count = 0;
	initLogic();
	generate(4);
	for(int i = 0;i < indices.size();++ i){
		int start = indices[i];
		int end = symbols.size();
		if(i < indices.size() - 1){
			end = indices[i + 1];
		}
		cout << "--------------generation " << (i + 1) << ", total " << (end - start) << "--------------" << endl;
		for(int j = start;j < end;++ j){
			if(symbols[j].empty()){
				cout << "\n" << endl;
				continue;
			}
			cout << symbols[j] << (j == end - 1 ? "" : ",");
			++count;
		}
		cout << endl;
		start = end;
	}
	cout << "total " << count << endl;
}
