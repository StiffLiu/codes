#ifndef MY_LIB_MY_MATH_H
#define MY_LIB_MY_MATH_H
namespace my_lib{
	/**
	 * This function returns the smallest integer {@var n}, such that the {@var n}th
	 * 	power of {@var e} is not less than {@var val}.
	 * {@var e} must be greater than {@code 1}, and {@var val} must not be 
	 * 	less than {@code 1}.
	 */
	template<class T>
	int ceilLog(T e, T val){
		//check the parameter. If not valid, return -1.
		if(e <= 1 || val < 1)
			return -1;
		if(val == 1)
			return 0;
		int n = 0;
		T tmp = 1;
		//loop invariant :
		//"tmp" is the n-th power of "e", and "tmp" is less than "val".
		while(tmp < val){
			++ n;
			tmp *= e;
		}
		return n;
	}
	/**
	 * This function returns the largest integer {@var n}, such that the {@var n}th
	 * 	power of {@var e} is not greater than {@var val}.
	 * {@var e} must be greater than {@code 1}, and {@var val} must not be 
	 * 	less than {@code 1}.
	 */
	template<class T>
	int floorLog(T e, T val){
		//check the parameter. If not valid, return -1.
		if(e <= 1 || val < 1)
			return -1;
		if(val == 1)
			return 0;
		int n = 1;
		T tmp = e;
		//loop invariant :
		//"tmp" is the n-th power of "e", and "tmp" is not greater than "val".
		while(tmp <= val){
			++ n;
			tmp *= e;
		}
		return n - 1;
	}
}
#endif //MY_LIB_MY_MATH_H
