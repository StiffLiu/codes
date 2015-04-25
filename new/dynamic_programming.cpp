#include <iostream>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <vector>
#include <ctime>
#include <cstring>
#include <map>
#include <algorithm>
#include <sstream>
#include <string>

namespace{
using namespace std;

/**
 * Assembly-line scheduling problem.
 */
template<class T>
T fastest_way(T *a0, T *a1, T *t0, T *t1, 
	const T& e0, const T& e1, const T& x0, const T& x1, 
	unsigned int n, unsigned int *o){
	if(n == 0)
		return T();
	unsigned int *l = new unsigned int[2 * n];
	unsigned int *l0 = l;
	unsigned int *l1 = l0 + n;
	T f0 = e0 + a0[0];
	T f1 = e1 + a1[0];

	for(unsigned int i = 1;i < n;++ i){
		T b0 = f0 + a0[i];
		T c0 = f1 + t1[i - 1] + a0[i];
		T b1 = f1 + a1[i];
		T c1 = f0 + t0[i - 1] + a1[i];
		if(b0 > c0){
			f0 = c0;
			l0[i - 1] = 1;
		}else{
			f0 = b0;
			l0[i - 1] = 0;
		}
		if(b1 > c1){
			f1 = c1;
			l1[i - 1] = 0;
		}else{
			f1 = b1;
			l1[i - 1] = 1;
		}
	
	}
	
	f0 = f0 + x0;
	f1 = f1 + x1;

	T w;
	if(f0 > f1){
		o[n - 1] = 1;
		w = f1;
	}else{
		o[n - 1] = 0;
		w = f0;
	}	
	for(unsigned int i = n - 1;i > 0;-- i)
		o[i - 1] = l[n * o[i] + i - 1];
	delete[] l;
	return w;
}

template<class T, class V>
auto matrixChainOrder(const unsigned int *p, unsigned int n, 
	T& costOutput, V& indexOutput)->decltype(costOutput[0][0]){
	for(unsigned int i = 0;i < n;++ i){
		costOutput[i][i] = 0;
		indexOutput[i][i] = i;
	}

	for(unsigned int k = 1;k < n;++ k){
		unsigned int end = n - k;
		for(unsigned int i = 0;i < end;++ i){
			unsigned int j = i + k;
			unsigned int minIndex = i;
			unsigned int minCost = costOutput[minIndex + 1][j] + p[i] * p[minIndex + 1] * p[j + 1];
			for(unsigned int s = i + 1;s < j;++ s){
				unsigned int cost = costOutput[i][s] + costOutput[s + 1][j] + p[i] * p[s + 1] * p[j + 1];
				if(cost < minCost){
					minCost = cost;
					minIndex = s;	
				}
			}
			indexOutput[i][j] = minIndex;
			costOutput[i][j] = minCost;
		}
	}
	return costOutput[0][n - 1];
}
template<class T, class V>
auto matrixChainMeoized(const unsigned int *p, unsigned int n,
	T& costOutput, V& indexOutput)->decltype(costOutput[0][0]){
	for(unsigned int i = 0;i < n;++ i){
		costOutput[i][i] = 0;
		indexOutput[i][i] = i;		
	}
	for(unsigned int i = 0;i < n;++ i){
		for(unsigned int j = i + 1;j < n;++ j)
			costOutput[i][j] = -1;
	}
	typedef decltype(costOutput[0][0]) Type;
	struct Functor{
		static Type func(const unsigned int *p, unsigned int i, unsigned int j,
			T& costOutput, V& indexOutput){
			Type& val = costOutput[i][j];
			if(val != -1)
				return val;

			unsigned int minIndex = i;
			unsigned int minCost = func(p, minIndex + 1, j, costOutput, indexOutput) + p[i] * p[minIndex + 1] * p[j + 1];
			for(unsigned int s = i + 1;s < j;++ s){
				unsigned int cost = p[i] * p[s + 1] * p[j + 1];
				if(cost >= minCost)
					continue;
				cost += func(p, i, s, costOutput, indexOutput);
				
				if(cost >= minCost)
					continue;
				cost += func(p, s + 1, j, costOutput, indexOutput);
				if(cost < minCost){
					minCost = cost;
					minIndex = s;
				}
			}
			indexOutput[i][j] = minIndex;
			val = minCost;
			return val;
		}
	};
	return Functor::func(p, 0, n - 1, costOutput, indexOutput);
}

template<class T, class V>
auto matrixChainGreedy(const unsigned int *p, unsigned int n,
	T& costOutput, V& indexOutput)->decltype(costOutput[0][0]){
	for(unsigned int i = 0;i < n;++ i){
		costOutput[i][i] = 0;
		indexOutput[i][i] = i;		
	}
	typedef decltype((costOutput[0][0])) Type;
	struct Functor{
		static Type func(const unsigned int *p, unsigned int i, unsigned int j,
			T& costOutput, V& indexOutput){
			Type& val = costOutput[i][j];
			if(i == j)
				return val;
			unsigned int minIndex =  i;
			unsigned int minCost = p[i] * p[minIndex + 1] * p[j + 1];
			for(unsigned int s = i + 1;s < j;++ s){
				unsigned int cost = p[i] * p[s + 1] * p[j + 1];
				if(cost < minCost){
					minCost = cost;
					minIndex = s;
				}
			}
			indexOutput[i][j] = minIndex;
			val = func(p, i, minIndex, costOutput, indexOutput) + 
				func(p, minIndex + 1, j, costOutput, indexOutput) + minCost;

			return val;			
		}
	};
	return Functor::func(p, 0, n - 1, costOutput, indexOutput);
}

template<class T>
void outputMatrixChainOrder(ostream& os, T& indexOutput, 
	unsigned int i, unsigned int j){
	if( i > j)
		return;
	if(i == j){
		os << 'A' << i;
		return;
	}
	else if (i + 1 == j){
		os << 'A' << i << " * " << 'A' << j;
		return;
	}
	
	auto index = indexOutput[i][j];
	outputMatrixChainOrder(os, indexOutput, i, index);
	os << " * ";
	if(index + 1 != j){
		os << '(';
		outputMatrixChainOrder(os, indexOutput, index + 1, j);
		os << ')';
	}else{
		outputMatrixChainOrder(os, indexOutput, index + 1, j);
	}
}

template<class S1, class S2, class T>
void longestCommonSequence(const S1& s1, unsigned int n1, 
	const S2& s2, unsigned int n2, T& lenOutput){
	if(n1 == 0 || n2 == 0)
		return;
	lenOutput[0][0] = (s1[0] == s2[0] ? 1 : 0);
	for(unsigned int i = 1;i < n2;++ i){
		lenOutput[0][i] = (s1[0] == s2[i] ? 1 : lenOutput[0][i - 1]);
	}
	for(unsigned int i = 1;i < n1;++ i){
		lenOutput[i][0] = (s1[i] == s2[0] ? 1 : lenOutput[i - 1][0]);
	}

	for(unsigned int i = 1;i < n1;++ i){
		for(unsigned int j = 1;j < n2;++ j){
			lenOutput[i][j] = (s1[i] == s2[j] ? lenOutput[i - 1][j - 1] + 1: 
				max(lenOutput[i - 1][j], lenOutput[i][j - 1]));
		}
	}
}

template<class S1, class S2, class T>
void outputLCS(ostream& os, const S1& s1, unsigned int n1,
	const S2& s2, unsigned int n2, const T& lenOutput, char ch = '*'){
	if(n1 == 0 && n2 == 0)
		return;

	unsigned int i = n1 - 1;
	unsigned int j = n2 - 1;
	vector<bool> indices1;
	vector<bool> indices2;

	indices1.resize(n1);
	indices2.resize(n2);

	while(i > 0 && j > 0 && lenOutput[i][j] != 0){
		if(s1[i] == s2[j]){
			indices1[i] = true;
			indices2[j] = true;
			-- i;
			-- j;
		}else if(lenOutput[i][j] == lenOutput[i - 1][j]){
			-- i;
		}else if(lenOutput[i][j] == lenOutput[i][j - 1]){
			-- j;
		}else{
			assert(false);
			return;
		}
		
	}

	if(s1[i] == s2[j]){
		indices1[i] = indices2[j] = true;
	}

	for(unsigned int i = 0;i < n1;++ i)
		os << (indices1[i] ? ch : s1[i]);
	os << endl;
	
	for(unsigned int i = 0;i < n2;++ i)
		os << (indices2[i] ? ch : s2[i]);
	os << endl;

	for(unsigned int i = 0;i < n2;++ i)
		if(indices2[i]) cout << s2[i];
	cout << endl;
}

/*
 * Longest Mono Increasing Sub Sequence
 */
template<class T, class C>
pair<unsigned int, unsigned int> lmiss(const T *values, unsigned int n, const C& c, unsigned int *out){
	if(n == 0)
		return {-1, -1};
	vector<unsigned int> lengths;
	unsigned int indexMax;
	unsigned int sizeMax;
	lengths.resize(n);

	out[0] = -1;
	lengths[0] = 1;
	indexMax = 0;
	sizeMax = 1;

	for(unsigned int i = 1;i < n;++ i){
		auto& o = out[i];
		auto& l = lengths[i];
		const auto& v = values[i];

		o = -1;
		l = 1;
		for(unsigned int j = i;j >= l;-- j){
			unsigned int tmpL = 0;
			if(!c(v, values[j - 1]) && (tmpL = lengths[j - 1] + 1) > l){
				o = j - 1;
				l = tmpL;
			}
		}
		if(l > sizeMax){
			sizeMax = l;
			indexMax = i;
		}
	}
	return {indexMax, sizeMax};
}

template<class T>
pair<unsigned int, unsigned int> lmiss(const T *values, unsigned int n, unsigned int *out){
	return lmiss(values, n, std::less<T>(), out);
}

template<class T, class C>
pair<unsigned int, unsigned int> lmissFaster(const T *values, unsigned int n, const C& c, unsigned int *out){
	if(n == 0)
		return {-1, -1};

	auto func=[&values, &c](unsigned int i, unsigned int j){
		return c(values[i], values[j]);
	};
	vector<unsigned int> lengths;
	lengths.push_back(0);
	out[0] = - 1;

	for(unsigned int i = 1;i < n;++ i){
		auto pos = std::upper_bound(lengths.begin(), lengths.end(), i, func);
		if(pos == lengths.end()){
			out[i] = lengths[lengths.size() - 1];
			lengths.push_back(i);
		}else{
			*pos = i;
			out[i] = (pos == lengths.begin() ? -1 : *(pos - 1));
		}
	}
	return {lengths.back(), lengths.size()}; 
}

template<class T>
pair<unsigned int, unsigned int> lmissFaster(const T *values, unsigned int n, unsigned int *out){
	return lmissFaster(values, n, std::less<T>(), out);
}

unsigned int lmissIndices(unsigned int *out, unsigned int n, 
	unsigned int index, unsigned int *indices){
	if(n == 0)
		return 0;
	unsigned int j = 0;
	while(out[index] != static_cast<unsigned int>(-1)){
		indices[j] = index;
		index = out[index];
		++ j;
	}
	indices[j] = index;
	++ j;
	return j;
}

template<class S, class T>
void lmissPrint(S& s, const T *values, unsigned int n1, unsigned int *indices, unsigned int n2){
	vector<bool> isIn;
	isIn.resize(n1);

	for(unsigned int i = 0;i < n2;++ i)
		isIn[indices[i]] = true;  

	for(unsigned int i = 0;i < n1;++ i){
		s << values[i];
	}
	s << "\n";
	for(unsigned int i = 0;i < n1;++ i){
		if(isIn[i])
			s << values[i];
		else
			s << ' ';
	}
}
template<class T>
struct SymetricMatrix{
	T *elements = nullptr;
	unsigned int n = 0;
	struct Row{
		unsigned int no = 0;
		const SymetricMatrix& matrix;
		Row(unsigned int no, const SymetricMatrix& matrix) 
			: no(no), matrix(matrix){
		}
		T& operator[](unsigned int i){
			unsigned int row = no;
			unsigned int col = i;
			if(col > row){
				row = i;
				col = no;
			}

			return matrix.elements[row * (row + 1) / 2 + col];
		}
	};
	explicit SymetricMatrix(unsigned int n) : n(n){
		elements = new T[n * (n + 1) / 2];
	}
	~SymetricMatrix(){
		delete[] elements;
	}
	SymetricMatrix(const SymetricMatrix&) = delete;
	SymetricMatrix& operator=(const SymetricMatrix&) = delete;
	Row operator[](unsigned int i) const{
		return Row(i, *this);		
	}
	
	template<class Stream>
	friend Stream& operator<<(Stream& os, const SymetricMatrix& matrix){
		for(unsigned int i = 0;i < matrix.n;++ i){
			unsigned int start = i * (i + 1) / 2;
			for(unsigned int k = 0;k <= i;++ k)
				os << matrix.elements[start + k];
			for(unsigned int k = i + 1;k < matrix.n;++ k)
				os << matrix.elements[k * (k + 1) / 2 + i];
			os << "\n";
		}
		return os;
	}
};

template<class T>
struct TwoDArray{
	T *array;
	unsigned int row;
	TwoDArray(T * array, unsigned int row) : array(array), row(row){
	}
	T* operator[](unsigned int r) const{
		return array + r * row;
	}
};

struct MyFormatedStream{
	unsigned int len = 8;
	template<class T>
	MyFormatedStream& operator<<(const T& val){
		cout << std::setw(len) << val;
		return *this;
	}
};

int testAssemblLine(int argc, char *argv[]){
	const unsigned int n = 6;
	int a0[n] = {7, 9, 3, 4, 8, 4};
	int t0[n] = {2, 3, 1, 3, 4};
	int a1[n] = {8, 5, 6, 4, 5, 7};
	int t1[n] = {2, 1, 2, 2, 1};
	int e0 = 2, e1 = 4;
	int x0 = 3, x1 = 2;
	unsigned int o[n] = {};
	int result = fastest_way(a0, a1, t0, t1, e0, e1, x0, x1, n, o);
	for(unsigned int i = 0;i < n;++ i)
		cout << o[i] << ' ';
	cout << endl;
	cout << "fastest : " << result << endl;
	return 0;
}

int testSymetricMatrix(int argc, char *argv[]){
	SymetricMatrix<unsigned int> matrix(5);
	for(unsigned int i = 0;i < 5;++ i){
		for(unsigned int j = 0;j <= i;++ j)
			matrix[i][j] = i * (i + 1) / 2 + j;
	}
	for(unsigned int i = 0;i < 5;++ i){
		for(unsigned int j = 0;j <= i;++ j){
			assert(matrix[i][j] == i * (i + 1) / 2 + j);
			assert(matrix[i][j] == matrix[j][i]);
		}
	}
	cout << matrix;
	return 0;
}

int testMatrixChainOrder(unsigned int *p, unsigned int n){
	SymetricMatrix<int> cost(n);
	SymetricMatrix<int> index(n);
	MyFormatedStream out;
	int count = matrixChainOrder(p, n, cost, index);
	out << cost;
	cout << endl;
	out << index;
	outputMatrixChainOrder(cout, index, 0, n - 1);
	cout << "\ndynamic : " << count << endl;
	return 0;
}

int testMatrixChainGreedy(unsigned int *p, unsigned int n){
	SymetricMatrix<int> cost(n);
	SymetricMatrix<int> index(n);
	MyFormatedStream out;
	int count = matrixChainGreedy(p, n, cost, index);
	out << cost;
	cout << endl;
	out << index;
	outputMatrixChainOrder(cout, index, 0, n - 1);
	cout << "\ngreedy : " << count << endl;
	return 0;
}

int testMatrixChain(int argc, char *argv[]){
	unsigned int p1[7] = {30, 35, 15, 5, 10, 20, 25};
	unsigned int p2[7] = {5,  10, 3, 12,  5, 50,  6};
	unsigned int p3[20];
	testMatrixChainOrder(p1, 6);
	testMatrixChainGreedy(p1, 6);
	cout << endl;
	testMatrixChainOrder(p2, 6);
	testMatrixChainGreedy(p2, 6);
	cout << endl;
	srand(time(0));
	for(unsigned int j = 0;j < 10;++ j){
		for(unsigned int i = 0;i < 20;++ i)
			p3[i] = rand() % 200 + 10;
		testMatrixChainOrder(p3, 19);
		testMatrixChainGreedy(p3, 19);
		cout << endl;
	}
	return 0;
}

int testLCS(const char *str1, const char *str2){
	unsigned int n1 = strlen(str1);
	unsigned int n2 = strlen(str2);
	int array[n1 * n2];
	TwoDArray<int> cost(array, n2);
	longestCommonSequence(str1, n1, str2, n2, cost);
	outputLCS(cout, str1, n1, str2, n2, cost);
	return 0;			
}

int testLCS(unsigned int argc, char *argv[]){
	testLCS("ABCBDAB", "BDCABA");
	testLCS("10010101", "010110110");
	return 0;
}

void validateLmiss(){
	for(unsigned int i = 20;i < 1000;++ i){
		unsigned int values[i];
		unsigned int out[i];
		for(unsigned int j = 0;j < i;++ j)
			values[j] = rand() % (2 * i);

		assert(lmiss(values, i, out).second ==
			lmissFaster(values, i, out).second);
	}
}

int testLmiss(unsigned int argc, char *argv[]){
	validateLmiss();
	unsigned int count = 24;
	unsigned int values[count];
	
	for(unsigned int i = 0;i < count;++ i)
		values[i] = rand() % 100;

	unsigned int out[count];
	auto cnt = lmiss(values, count, out);
	assert(cnt.first != static_cast<unsigned int>(-1));
	assert(cnt.second!= static_cast<unsigned int>(-1));
	
	unsigned int indices[cnt.second];
	auto ret = lmissIndices(out, count, cnt.first, indices);               
	MyFormatedStream stream;
	stream.len = 4;
	for(unsigned int i = 0;i < count;++ i)
		stream << (int)out[i];
	cout << endl;
	assert(ret == cnt.second);
	lmissPrint(stream, values, count, indices, cnt.second);
	return 0;
}

template<class T>
void optimalMovement(T * values, unsigned int n, 
	T *maxValues, unsigned int* movement = nullptr){
	if(n <= 1)
		return;
	unsigned int stride = 3 * n;

	T *current = values;
	for(unsigned int i = 1;i < n;++ i, current += stride){
		T tmpVal = T();
		T *lastMaxValues = maxValues + n * (i - 1);
		T *currentMaxValues = lastMaxValues + n;
		currentMaxValues[0] = lastMaxValues[0] + current[1];
		currentMaxValues[n - 1] = lastMaxValues[n - 1] + current[3 * n - 2];
		
		if(n > 1){
			if((tmpVal = lastMaxValues[1] + current[3]) > currentMaxValues[0])
				currentMaxValues[0] = tmpVal;
			if((tmpVal = lastMaxValues[n - 2] + current[3 * n - 4]) > currentMaxValues[n - 1])
				currentMaxValues[n - 1] = tmpVal;
		}

		for(unsigned int j = 1;j < n - 1;++ j){
			unsigned int index = 3 * j;
			currentMaxValues[j] = lastMaxValues[j - 1] + current[index - 1];
			currentMaxValues[j] = max(currentMaxValues[j], lastMaxValues[j] + current[index + 1]);
			currentMaxValues[j] = max(currentMaxValues[j], lastMaxValues[j + 1] + current[index + 3]);
		}

	}

	if(movement != nullptr){
		movement[n - 1] = max_element(maxValues + (n - 1) * n, 
			maxValues + n * n) - maxValues - (n - 1) * n;

		current = values + (n - 2) * stride;
		for(unsigned int i = 0;i < n - 1;++ i, current -= stride){
			unsigned int j = n - 2 - i;
			T *currentMaxValues = maxValues + n * j;
			T *lastMaxValues = currentMaxValues + n;
			unsigned int k = movement[j + 1];
			unsigned int index = 3 * k;
			if(k > 0 && currentMaxValues[k - 1] + current[index - 1] == lastMaxValues[k])
				movement[j] = k - 1;
			else if(k < n - 1 && currentMaxValues[k + 1] + current[index + 3] == lastMaxValues[k])
				movement[j] = k + 1;
			else if(currentMaxValues[k] + current[index + 1] == lastMaxValues[k])
				movement[j] = k;
			else{
				cout << "This algorithm is buggy" << endl;
				break;
			}
		}
	}
}

struct Spliter{
	unsigned int k, n;
	char crossChar, spaceChar;
	Spliter(unsigned int k, unsigned int n, char crossChar = '+',
		char spaceChar = '-') : k(k), n(n), crossChar(crossChar), spaceChar(spaceChar){
	}
	friend ostream& operator<<(ostream& os, const Spliter& vs){
		unsigned int w = vs.k + 1;
		unsigned int total = vs.n * w + 1;
		for(unsigned int i = 0;i < total;++ i)
			os << (i % w == 0 ? vs.crossChar : vs.spaceChar);
		return os;
	}
};

struct Table{
	unsigned int s, m, k, n;
	Table(unsigned int s, unsigned int m, unsigned int k, unsigned int n)
		:s(s), m(m), k(k), n(n){
	}
	friend ostream& operator<<(ostream& os, const Table& t){
		Spliter vs(t.k, t.n);
		Spliter hs(t.k, t.n, '|', ' ');
		unsigned int h = t.s + 1;
		unsigned int total = t.m * h + 1;
		for(unsigned int i = 0;i < total;++ i)
			os << (i % h == 0 ? vs : hs) << endl;
		return os;
	}
};

int testOptimalMovement(int argc, char *argv[]){
	const unsigned int n = 6;
	int values[3 * (n - 1) * n];
	int maxValues[n * n];
	unsigned int movement[n];
	srand(time(0));
	for(unsigned int i = 0;i < 3 * (n - 1) * n;++ i)
		values[i] = rand() % 100 - 50;
	for(unsigned int i = 0;i < n;++ i)
		maxValues[i] = rand() % 100 - 50;

	const unsigned int w = 15;
	optimalMovement(values, n, maxValues, movement);
	for(unsigned int i = 0;i < n - 1;++ i){
		for(unsigned int j = 0;j < n;++ j){
			unsigned int index = 3 * (n * i + j);
			std::ostringstream ss;
			ss << (int)values[index] << ',' << values[index + 1] << ',' << values[index + 2];
			cout << std::setw(w) << ss.str();
		}
		cout << endl;
	}
	cout << endl;

	for(unsigned int i = 0;i < n;++ i){
		for(unsigned int j = 0;j < n;++ j)
			cout << std::setw(8) << maxValues[i * n + j];
		cout << endl;
	}
	cout << endl;

	for(unsigned int i = 0;i < n;++ i)
		cout << std::setw(w) << (int)movement[i];
	cout << endl;

	return 0;
}

}
int main(int argc, char *argv[]){
	testOptimalMovement(argc, argv);
	Table table(3, 6, 7, 6);
	cout << table;
	testLmiss(argc, argv);
	return 0;
}
