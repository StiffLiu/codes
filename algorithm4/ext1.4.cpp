#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <ctime>
#include <set>
#include <cassert>
#include <algorithm>
#include <iterator>
#include <fstream>

using namespace std;
/*
 * If an integer is randomly choosen from the range 1~N,
 * how many times we need to choose, before the first time 
 * a repeated integer is choosed? It is said that the expected 
 * number of chooses is about C*sqrt(N), where C is a constant.
 * This programme validate this hypothesis.
 * This programme could also answers the expected number of people we 
 * need to have so that at least a pair of them are born on the
 * same day.
 */
int birthdayProblem(int argc, char *argv[]){
	srand(time(0));
	const int start = 10, end = 1000;
	const int iteration = 100000;
	for(int i = start;i <= end;++ i){
		int count = 0;
		for(int j = 0;j < iteration;++ j){
			vector<bool> isUsed(i, false);
			int tmp = rand() % i;
			while(!isUsed[tmp]){
				isUsed[tmp] = true;
				tmp = rand() % i;
				++ count;
			}
		}

		double val1 = sqrt(i);
		double val2 = (count / (double)iteration);
		cout << val1 << "\t\t" << val2 << "\t\t" << (val2 / val1) << endl;
	}
}
/*
 * This programme validates the hypothesis that the number of integers generated
 * before all possible values are generated is ~N*lg(N).
 * Consider the following problem:
 * There are N bins. Balls are thrown into these bins with equal possibilities.
 * How many balls are expected to thrown before every bin contains at leat one ball?
 *
 * This programm answers the above problem.
 */
int couponCollectorProlem(int argc, char *argv[]){
	const int start = 10, end = 1000;
	const int iteration = 1000;
	for(int i = start;i <= end;++ i){
		int count = 0;
		for(int j = 0;j < iteration;++ j){
			vector<bool> isUsed(i, false);
			int tmp = rand() % i;
			int cnt = 0;
			while(cnt != i){
				if(!isUsed[tmp]){
					isUsed[tmp] = true;
					++ cnt;
				}
				tmp = rand() % i;
				++ count;
			}
		}

		double val1 = i * log(i);
		double val2 = (count / (double)iteration);
		cout << val1 << "\t\t" << val2 << "\t\t" << (val2 / val1) << endl;
	}
}
void generateBitonicArray(int max, int min, int n, vector<int>& result){
	int mid = rand() % n;
	int diff = max - min;
	set<int> lHalf;
	set<int, greater<int> > rHalf;
	int i = 0;
	for(;i <= mid;++ i)
		lHalf.insert(rand() % diff + min);
	lHalf.insert(min);
	lHalf.insert(max);
	for(;i < n;++ i){
		int num = rand() % diff + min;
		if(lHalf.find(num) == lHalf.end())
			rHalf.insert(num);
	}
	set<int>::iterator begin = lHalf.begin(), end = lHalf.end();
	while(begin != end){
		result.push_back(*begin);
		++begin;
	}

	set<int, std::greater<int> >::iterator start = rHalf.begin(), term = rHalf.end();
	while(start != term){
		result.push_back(*start);
		++ start;
	}
}
/*
 * a bitonic array is an array comprimised of an increasing sequence of
 * elements, followed by a decreasing sequence of elements.
 * This algorithm determines whether an integer is in a bitonic array of
 * n distinct integers. This algorithm is expected to use ~3lg(n) comparisons.
 */
int binarySearch(int *val, unsigned int n, int key, int& numOfComparison){
	++ numOfComparison;
	if(n == 0)
		return -1;
	if(n == 1){
		++ numOfComparison;
		return (val[0] == key ? 0 : -1);
	}
	unsigned int low = 0;
	unsigned int high = n - 1;
	while(low < high){
		int mid = (low + high) / 2;
		++ numOfComparison;
		if(val[mid] > key)
			high = mid;
		else 
			low = mid + 1;
	}
	++ numOfComparison;
	return (val[low] == key ? low : -1);
}
int bitonicSearch(int *val, unsigned int n, int key, int& numOfComparison){
	if(n == 0)
		return -1;
	if(n == 1){
		++ numOfComparison;
		return (val[0] == key ? 0 : -1);
	}
	unsigned int low = 0;
	unsigned int high = n - 1;
	while(low < high){
		unsigned int mid = (low + high) / 2;
		++ numOfComparison;
		if(val[mid] > val[mid + 1]){
			high = mid;
		}else{
			low = mid + 1;
		}
	}

	int ret = -1;
	if((ret = binarySearch(val, low, key, numOfComparison)) != -1)
		return ret;
	return binarySearch(val + low, n - low, key, numOfComparison);	
}
int testBitonicSearch(int argc, char *argv[]){
	const int iteration1 = 100;
	const int iteration2 = 10;
	int max = 1000000, min = 10;
	const unsigned int numOfInts = 100000;
	double expected = log(numOfInts) / log(2);
	vector<int> numbers;
	for(int i = 0;i < iteration2;++ i){
		numbers.clear();
		generateBitonicArray(max, min, numOfInts, numbers);
		int numOfComparison = 0;
		for(int j = 0;j < iteration1;++ j){
			int key = rand();
			bool ret1 = (bitonicSearch(&numbers[0], numbers.size(), key, numOfComparison) != -1);
			bool ret2 = (std::find(numbers.begin(), numbers.end(), key) != numbers.end());
			assert(ret1 == ret2);
		}
		cout << "expected : " << 3 << ", actual : " << (numOfComparison / (double)iteration1) / expected << endl;
	}
}
/*
 * A local minimum in array A of length n is an index i such that:
 * A[i] <= A[i - 1], if i > 0;
 * A[i] >= A[i + 1], if i < n;
 * It's not difficult to prove that there must exists a 
 * local minimum for an arbitrary array.
 * This algorithm finds a local minimum from an array 
 * of n DISTINCT elements with ~lg(n) comparisons.
 */
template<class T>
bool isLocalMinimum(T *val, unsigned int n, unsigned int i){
	if(i > 0 && val[i] > val[i - 1])
		return false;
	if(i < n - 1 && val[i] > val[i + 1])
		return false;
	return true;
}
template<class T>
unsigned int localMinimum(T* val, unsigned int n){
	if(n == 0)
		return -1;
	if(n == 1)
		return 0;
	if(val[0] <= val[1])
		return 0;
	if(val[n - 1] <= val[n - 2])
		return n - 1;
	unsigned int low = 1, high = n - 2;
	while(low < high){
		unsigned int mid = (low + high) / 2;
		if(val[mid] <= val[mid + 1] && val[mid] <= val[mid - 1])
			return mid;
		if(val[mid + 1] > val[mid - 1]){
			high = mid - 1;
		}else{
			low = mid + 1;
		}
	}
	return low;
}
/*
 * A local minimum in 2D m*n array A of a pair of index (i,j ), such that:
 * A[i][j] <= A[i][j-1], if j > 0;
 * A[i][j] <= A[i][j + 1], if j < n - 1;
 * A[i][j] <= A[i + 1][j], if i < m - 1;
 * A[i][j] <= A[i - 1][j], if i > 0;
 * pigeonhole principle and graph theory could be used to prove that there 
 * must exists a local minimum for an arbitrary 2D array.
 * This algorithm finds a local minimum from a m*n 2D array, the elements of 
 * which are DISTINCT, with O(m + n) comparisions.
 * The cost of the algorithm is O(m + n) follows from the recurrence relation:
 *	T(m, n) = T(m / 2, n / 2) + O(m + n)
 * Upper bound of this algorithm in worst case situation in terms of the number of comparisions
 * made is about 16 * (m + n).
 * But the average case is very good
 */
template<class T>
inline bool isLocalMinimum(const T *array2D, unsigned int m, 
	unsigned int n, unsigned int i, unsigned int j, unsigned int& comparisions){
	unsigned int index = i * n + j;
	comparisions += 1;
	if(i > 0){
		comparisions += 1;
		if(array2D[index] > array2D[index - n]){
			return false;
		}
	}

	comparisions += 1;
	if(i < m - 1){
		comparisions += 1;
		if(array2D[index] > array2D[index + n]){
			return false;
		}
	}

	comparisions += 1;
	if(j > 0){
		comparisions += 1;
		if(array2D[index] > array2D[index - 1]){
			return false;
		}
	}

	comparisions += 1;
	if(j < n - 1){
		comparisions += 1;
		if(array2D[index] > array2D[index + 1]){
			return false;
		}
	}

	return true;
}
template<class T>
inline bool isLocalMinimumFast(const T *array2D, 
	unsigned int n, unsigned int i, unsigned int j, unsigned int& comparisions){
	unsigned int index = i * n + j;
	comparisions += 1;
	if(array2D[index] > array2D[index - n])
		return false;

	comparisions += 1;
	if(array2D[index] > array2D[index + n])
		return false;

	comparisions += 1;
	if(array2D[index] > array2D[index - 1])
		return false;

	comparisions += 1;
	if(array2D[index] > array2D[index + 1])
		return false;

	return true;
}
template<class T>
void output(ostream& os, const T *vals, unsigned int n, unsigned int startRow, unsigned int endRow,
	unsigned int startCol, unsigned int endCol){
	for(unsigned int i = startRow;i <= endRow;++ i){
		for(unsigned int j = startCol;j <= endCol;++ j){
			os << vals[i * n + j] << ' ';
		}
		os << endl;
	}
}
template<class T>
void output(ostream& os, const T *vals, unsigned int m, unsigned int n){
	output(os, vals, n, 0, m - 1, 0, n - 1);
}
unsigned int calculateOdd(unsigned int n ,unsigned int total, 
	unsigned int offset, unsigned int i, unsigned int j){
	int tmp;
	switch(i % 4){
	case 0:
		if(j > 0)
			return  (total + i / 4 * 2 * (n - 1) + j - 1);
		tmp = i / 4 * 2 * (n + 1);
		break;
	case 1:
		tmp = i / 4 * 2 * (n + 1) + 1 + j;
		break;
	case 2:
		if(j < n - 1)
			return (total + (i / 4 * 2 + 1) * (n - 1) + j);
		tmp = i / 4 * 2 * (n + 1) + 1 + n;
		break;
	case 3:
		tmp = i / 4 * 2 * (n + 1) + 1 + 2 * n - j;
		break;
				
	}
	if(tmp >= offset)
		return tmp - offset;
	return total - 1 - tmp;	
}
unsigned int calculateEven(unsigned int n ,unsigned int total, 
	unsigned int offset, unsigned int i, unsigned int j){
	int tmp;
	switch(i % 4){
	case 0:
		tmp = i / 4 * 2 * (n + 1) + j;
		break;
	case 1:
		if(j < n - 1)
			return  (total + i / 4 * 2 * (n - 1) + j);
		tmp = i / 4 * 2 * (n + 1) + n;
		break;
	case 2:
		tmp = i / 4 * 2 * (n + 1) + 2 * n - j;
		break;
	case 3:
		if(j > 0)
			return (total + (i / 4 * 2 + 1) * (n - 1) + j - 1);
		tmp = i / 4 * 2 * (n + 1) + 1 + 2 * n;
		break;				
	}
	if(tmp >= offset)
		return tmp - offset;
	return total - 1 - tmp;	
}
void generateSpecialArray(int *output, unsigned int m, unsigned int n, unsigned int min, 
	unsigned int row, unsigned int col){
	if(row >= m)
		row = m - 1;
	if(col >= n)
		col = n - 1;
	unsigned int total = m / 4 * 2 * (n + 1);
	if(row % 2 != 0){
		switch(m % 4){
			case 1:
				total += 1;break;
			case 2:
				total += (n + 1);break;
			case 3:
				total += (n + 2);break;
		}
		unsigned int offset = calculateOdd(n, total, 0, row, col);
		cout << "offset " << offset << ", total " << total << ", row " << row << ", col " << col << endl;
		for(unsigned int i = 0;i < m;++ i){
			for(unsigned int j = 0;j < n;++ j){
				output[i * n + j] = min + calculateOdd(n, total, offset, i, j);
			}
		}
	}else{
		switch(m % 4){
			case 1:
				total += n;break;
			case 2:
				total += (n + 1);break;
			case 3:
				total += (2 * n + 1);break;
		}
		unsigned int offset = calculateEven(n, total, 0, row, col);
		cout << "offset " << offset << ", total " << total << ", row " << row << ", col " << col << endl;
		for(unsigned int i = 0;i < m;++ i){
			for(unsigned int j = 0;j < n;++ j){
				output[i * n + j] = min + calculateEven(n, total, offset, i, j);
			}
		}
	}
}
void generateSpecialArray(int *output, unsigned int m, unsigned int n, unsigned int min){
	generateSpecialArray(output, m, n, min, m - 1, n / 2);
}
template<class T>
pair<unsigned int, unsigned int> buggyLocalMinimum(const T * array2D, unsigned int m, unsigned int n, unsigned int& comparisions){
	typedef pair<unsigned int, unsigned int> Type;
	if(m == 0 || n == 0)
		return Type(-1, -1);
	if(m == 1)
		return Type(0, localMinimum(array2D, n));
	if(n == 1)
		return Type(localMinimum(array2D, m), 0);
	unsigned int startRow = 0, endRow = m - 1;
	unsigned int startCol = 0, endCol = n - 1;
	while(startRow != endRow && startCol != endCol){
		//cout << "--------------------"<<endl;
		//output(cout, array2D, n, startRow, endRow, startCol, endCol);
		//cout << "--------------------" << endl;
		unsigned int midRow = (startRow + endRow) / 2;
		unsigned int midCol = (startCol + endCol) / 2;
		unsigned int colMin = startCol;
		unsigned int rowMin = startRow;

		if(isLocalMinimum(array2D, m, n, midRow, colMin, comparisions))
			return Type(midRow, colMin);
		for(unsigned int i = startCol + 1;i <= endCol;++ i){
			//This test is necessary, or there would be bugs
			if(isLocalMinimum(array2D, m, n, midRow, i, comparisions))
				return Type(midRow, i);
			if(array2D[midRow * n + i] < array2D[midRow * n + colMin])
				colMin = i;
			++comparisions;
		}

		if(isLocalMinimum(array2D, m, n, rowMin, midCol, comparisions)){
			assert(isLocalMinimum(array2D, m, n, rowMin, midCol, comparisions));
			return Type(rowMin, midCol);
		}
		for(unsigned int i = startRow + 1;i <= endRow;++ i){
			//This test is necessary, or there would be bugs
			if(isLocalMinimum(array2D, m, n, i, midCol, comparisions))
				return Type(i, midCol);
			if(array2D[i * n + midCol] < array2D[rowMin * n + midCol])
				rowMin = i;
			++comparisions;
		}
		comparisions += 2;
		if(colMin == midCol && rowMin == midRow)
			return Type(midRow, midCol);
		comparisions += 3;
		if(array2D[midRow * n + colMin] < array2D[rowMin * n + midCol]){
			if(colMin > midCol){
				startCol = midCol + 1;
			}else{
				endCol = midCol;
			}
			//We could gaurantee that midRow + 1 <= endRow,
			//that is we won't access memory out of the 2d array
			if(array2D[midRow * n + colMin] > array2D[(midRow + 1) * n + colMin]){
				startRow = midRow + 1;
			}else /*This is unnecessary*//*if(midRow > 0 && array2D[midRow * n + colMin] > array2D[(midRow - 1) * n + colMin]){
				endRow = midRow;
			}else{
				return Type(midRow, colMin);
			}*/{
				endRow = midRow;
			}
		}else{
			if(rowMin > midRow){
				startRow = midRow + 1;
			}else{
				endRow = midRow;
			}
			//We could gaurantee that midCol + 1 <= endCol,
			//that is we won't access memory out of the 2d array
			if(array2D[rowMin * n + midCol] > array2D[rowMin * n + midCol + 1]){
				startCol = midCol + 1;
			}else /*This is unnecessary*//*if(midCol > 0 && array2D[rowMin * n + midCol] > array2D[rowMin * n + midCol - 1]){
				endCol = midCol;
			}else{
				return Type(rowMin, midCol);
			}*/{
				endCol = midCol;
			}
		}		
	}
	if(startRow == endRow){
		for(unsigned int j = startCol;j <= endCol;++ j)
			if(isLocalMinimum(array2D, m, n, startRow, j, comparisions))
				return Type(startRow, j);
		//ofstream out("test.txt");
		output(cout, array2D, m, n);
		cout << "startRow = " << startRow << ", startCol = " << startCol << ", endCol = " << endCol << endl;
		assert(false);
	}
	if(startCol == endCol){
		for(unsigned int j = startRow;j <= endRow;++ j)
			if(isLocalMinimum(array2D, m, n, j, startCol, comparisions))
				return Type(j, startCol);
		output(cout, array2D, m, n);
		cout << "startRow = " << startRow << ", startCol = " << startCol << ", endCol = " << endCol << endl;
		assert(false);
	}
	assert(false);
	return Type(-1, -1);
}
void testLocalMinimun(int argc, char *argv[]){
	//test for local minimum of 1D array
	{
		const unsigned int numOfInts = 40;
		int vals[numOfInts];
		for(unsigned int i = 0;i < numOfInts;++ i)
			vals[i] = i + 10;
		for(int j = 0;j < 10;++ j){
			std::random_shuffle(vals, vals + numOfInts);
			for(unsigned int i = 0;i < numOfInts;++ i)
				cout << vals[i] << ' ';
			cout << endl;
			unsigned int index = localMinimum(vals, numOfInts);
			for(unsigned int i = 0;i < index;++ i)
				cout << "   ";
			cout << "^^" << endl;
		}
	}
	//test for local minimum of 2D array
	{
		int failArray[] = {
			39, 38, 37, 36, 35, 34, 33, 32, 31,
			40, 41, 42, 43, 44, 45, 46, 47, 30,
			21, 22, 23, 24, 25, 26, 27, 28, 29,
			20, 48, 49, 50, 51, 52, 53, 54, 55,
			2, 18, 17, 16, 15, 14, 13, 12, 11,
			56, 57, 58, 59, 60, 61, 62, 63, 10,
			4, 3, 19, 1, 5, 6, 7, 8, 9,
		};
		unsigned int comparisions = 0;
		pair<unsigned int, unsigned int> index = buggyLocalMinimum(failArray, 7, 9, comparisions);
		for(unsigned int i = 0;i < 7;++ i){
			for(unsigned int j = 0;j < 9;++ j){
				if(i == index.first && j == index.second)
					cout << "*** ";
				else
					cout << failArray[i * 9 + j] << ' ';
			}
			cout << endl;
		}
		const unsigned int m = 200;
		const unsigned int n = 200;
		int *vals = new int[m * n];
		for(unsigned int i = 0;i < m * n;++ i)
			vals[i] = i + 100;
		cout << vals[m * n - 1] << endl;
		for(int j = 0;j < 1;++ j){
			//generate random 2D array with distinct numbers.
			std::random_shuffle(vals, vals + m * n);
			comparisions = 0;
			index = buggyLocalMinimum(failArray, 7, 9, comparisions);
			if(comparisions == 0)
			for(unsigned int i = 0;i < m;++ i){
				for(unsigned int j = 0;j < n;++ j){
					if(i == index.first && j == index.second)
						cout << "*** ";
					else
						cout << vals[i * n + j] << ' ';
				}
				cout << endl;
			}
			cout << "local minimum is at (" << index.first << "," << index.second << ") : "
				<< vals[index.first * n + index.second] << ", number of comparisions : "
				<< comparisions << ", ratio : " << (comparisions / (double)(m + n)) << endl;
		}
		delete []vals;
	}
}
void testSpeicalArray(unsigned int m, unsigned int n, unsigned int min){
	int *vals = new int[m * n];
	generateSpecialArray(vals, m, n, min);
	output(cout, vals, m, n);
	delete []vals;
}
void testSpeicalArray(int argc, char *argv[]){
	testSpeicalArray(23, 20, 100);
	cout << "------------------" << endl;
	testSpeicalArray(15, 10, 100);
	cout << "------------------" << endl;
	testSpeicalArray(14, 10, 100);
	cout << "------------------" << endl;
	testSpeicalArray(13, 10, 100);
	cout << "------------------" << endl;
	testSpeicalArray(12, 10, 100);
	cout << "------------------" << endl;
	testSpeicalArray(12, 5, 100);
}
int main(int argc, char *argv[]){
	testLocalMinimun(argc, argv);
}
