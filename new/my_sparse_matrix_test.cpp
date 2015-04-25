#include "my_sparse_matrix.h"
#include "my_binary_search_tree.h"
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <random>

template<class T>
void randArray(T& array){
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> generator(0, 100);

	for(int i = 0;i < array.dim();++ i)
		if(rand() % 10 < 2)
			array.set(i, generator(gen) / 100.0);
}

template<class T>
void randMatrix(T& matrix){
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> generator(0, 100);
	
	for(int i = 0;i < matrix.row();++ i)
		for(int j = 0;j < matrix.col();++ j)
			if(rand() % 10 < 2)
				matrix.set(i, j, generator(gen) / 100.0);
}

using namespace std;

int test(int argc, char *argv[]){
  using namespace my_lib;
  //typedef StdSTAdapter<std::map<unsigned int, double> > ST;
  typedef RBTree<unsigned int, double> ST;
  SparseArray<double, ST> array1(25);
  SparseArray<double, ST> array2(25, 0.4);

  randArray(array1);
  randArray(array2);

  cout << "array1          : " << array1 << endl;
  cout << "array2          : " << array2 << endl;
  cout << "array1 + 0.5    : " << array1 + 0.5 << endl;
  cout << "array2 * 0.2    : " << array2 * 0.2 << endl;
  cout << "array1 + array2 : " << array1 + array2 << endl;
  cout << "array1 - array2 : " << array1 - array2 << endl;
  cout << "array1 . array2 : " << array1.dot(array2) << endl;
  
  SparseMatrix<double, ST> matrix1(20, 20);
  SparseMatrix<double, ST> matrix2(20, 20);
  SparseMatrix<double, ST> matrix3(20, 25);

  randMatrix(matrix1);
  randMatrix(matrix2);
  randMatrix(matrix3);

  cout << "matrix1, " << matrix1.ratio() << " : \n" << matrix1 << endl;
  cout << "matrix2, " << matrix2.ratio() << " : \n" << matrix2 << endl;
  cout << "matrix3, " << matrix3.ratio() << " : \n" << matrix3 << endl;

  auto m1_p_m2 = matrix1 + matrix2;
  cout << "matrix1 + matrix2, " << m1_p_m2.ratio() << " : \n" << m1_p_m2 << endl;
  auto m1_m2 = matrix1 - matrix2;
  cout << "matrix1 - matrix2, " << m1_m2.ratio() << " : \n" << m1_m2 << endl;
  auto m2xm3 = matrix2 * matrix3;
  cout << "matrix2 * matrix3, " << m2xm3.ratio() << " : \n" << m2xm3 << endl;
  auto m1xm2xm3 = matrix1 * matrix2 * matrix3;
  cout << "matrix1 * matrix2 * matrix3, " << m1xm2xm3.ratio() << " : \n" << m1xm2xm3 << endl;
  auto m1xm1xm2xm3 = matrix1 * matrix1 * matrix2 * matrix3;
  cout << "matrix1 * matrix1 * matrix2 * matrix3, " << m1xm1xm2xm3.ratio() << " : \n" << m1xm1xm2xm3 << endl;
  auto complicated = (matrix1 + matrix2) * (matrix1 - matrix2) * matrix3;
  cout << "(matrix1 + matrix2) * (matrix1 - matrix2) * matrix3, " << complicated.ratio() << " : \n" << complicated << endl;

  cin.get();
  return 0;
}

int main(int argc, char *argv[]){
	return test(argc, argv);
}
