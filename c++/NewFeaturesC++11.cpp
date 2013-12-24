#include <iostream>
#include <vector>
#include <thread>
using namespace std;
/**
 * new features of c++11
 */
/**
 * Rvalue references and move constructors
 */

/**
 * constexpr C Generalized constant expressions
 */
template<class T>
struct SumOfArray{
	T *val;
	unsigned int n;
	constexpr SumOfArray(T *val, unsigned int n) : val(val), n(n){
	}
	constexpr T sum() const{
		return sum(0);
	}
private:
	constexpr T sum(unsigned int start) const{
		return start < n ? sum(start + 1) + val[start] : T();
	}
};
int testConstExpr(int argc, char *argv[]){
	int values[]={1, 2, 3, 4};

	//since the compiler must allocate memory
	//for an array it must know the size of the array
	//at compile time., This one shows that "constexpr" really works
	int another[SumOfArray<int>(values, 4).sum()];
	cout << (sizeof another) << endl;
	return 0;
}
class TestTLS{
	int i;
public:
	static int count;
	TestTLS(){
		i = (count ++);
		cout << i <<  " TestTLS constructor called" << endl;
	}
	~TestTLS(){
		cout << i << " TestTLS destructor called" << endl;
	}
};
int TestTLS::count = 0;
thread_local TestTLS test;
void fun(){
	cout << std::this_thread::get_id() << ":" << (void*)&test << endl;
}
thread_local std::vector<int> ints{1, 3, 4, 9, 9, 10};
int testCPP0X(int argc, char *argv[]){
	for(auto &val : ints)
		cout << val << ' ' << endl;
	std::thread t(fun);
	fun();
	t.join();
	cout << endl;
}
