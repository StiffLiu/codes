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
 */
template<class T>
inline bool isLocalMinimum(T *array2D, unsigned int m, 
	unsigned int n, unsigned int i, unsigned int j, unsigned int& comparisions){
	unsigned int index = i * n + j;
	if(i > 0 && array2D[index] > array2D[index - n]){
		comparisions += 2;
		return false;
	}
	if(i < m - 1 && array2D[index] > array2D[index + n]){
		comparisions += 2;		
		return false;
	}
	if(j > 0 && array2D[index] > array2D[index - 1]){
		comparisions += 2;
		return false;
	}
	if(j < n - 1 && array2D[index] > array2D[index + 1]){
		comparisions += 2;
		return false;
	}
	return true;
}
template<class T>
void output(ostream& os, T *vals, unsigned int m, unsigned int n){
	for(unsigned int i = 0;i < m;++ i){
		for(unsigned int j = 0;j < n;++ j){
			os << vals[i * n + j] << ' ';
		}
		os << endl;
	}
}
template<class T>
pair<unsigned int, unsigned int> localMinimum(T * array2D, unsigned int m, unsigned int n, unsigned int& comparisions){
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
		unsigned int midRow = (startRow + endRow) / 2;
		unsigned int midCol = (startCol + endCol) / 2;
		unsigned int colMin = startCol;
		unsigned int rowMin = startRow;
		comparisions += (endCol - startCol);
		for(unsigned int i = startCol + 1;i <= endCol;++ i)
			if(array2D[midRow * n + i] < array2D[midRow * n + colMin])
				colMin = i;
		comparisions += (endRow - startRow);		
		for(unsigned int i = startRow + 1;i <= endRow;++ i)
			if(array2D[i * n + midCol] < array2D[rowMin * n + midCol])
				rowMin = i;
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
			}else{
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
			}else{
				endCol = midCol;
			}
		}		
	}
	if(startRow == endRow){
		for(unsigned int j = startCol;j <= endCol;++ j)
			if(isLocalMinimum(array2D, m, n, startRow, j, comparisions))
				return Type(startRow, j);
		/*ofstream out("test.txt");
		output(out, array2D, m, n);
		out << "startRow = " << startRow << ", startCol = " << startCol << ", endCol = " << endCol << endl;*/
		assert(false);
	}
	if(startCol == endCol){
		for(unsigned int j = startRow;j <= endRow;++ j)
			if(isLocalMinimum(array2D, m, n, j, startCol, comparisions))
				return Type(j, startCol);
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
		const unsigned int m = 200;
		const unsigned int n = 200;
		int vals[m * n];
		for(unsigned int i = 0;i < m * n;++ i)
			vals[i] = i + 100;
		
		for(int j = 0;j < 10;++ j){
			//generate random 2D array with distinct numbers.
			std::random_shuffle(vals, vals + m * n);
			unsigned int comparisions = 0;
			pair<unsigned int, unsigned int> index = localMinimum(vals, m, n, comparisions);
			/*for(unsigned int i = 0;i < m;++ i){
				for(unsigned int j = 0;j < n;++ j){
					if(i == index.first && j == index.second)
						cout << "*** ";
					else
						cout << vals[i * n + j] << ' ';
				}
				cout << endl;
			}*/
			//After running serveral tests we could see that, it takes less than 3*(m + n)
			//and greater than 2*(m + n) comparisions to find a local minimum of the
			//2D array, this lower bound and upper are not coincidences. It can be shown
			//by rigorous proof.
			cout << "local minimum is at (" << index.first << "," << index.second << ") : "
				<< vals[index.first * n + index.second] << ", number of comparisions : "
				<< comparisions << ", ratio : " << (comparisions / (double)(m + n)) << endl;
		}
	}
}
int main(int argc, char *argv[]){
	couponCollectorProlem(argc, argv);
}
