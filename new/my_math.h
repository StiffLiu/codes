#ifndef MY_LIB_MY_MATH_H
#define MY_LIB_MY_MATH_H
#include <vector>
#include <cassert>
#include <random>
#include <algorithm>
#include <iostream>
#include <unordered_map>

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
	class DiscreteDistribution{
		typedef double Real;
		typedef std::vector<Real> DoubleVector;
		typedef std::random_device RandomDevice;
		typedef std::default_random_engine DefaultRandomEngine;
		typedef std::uniform_real_distribution<Real> UniformRealDistribution;
		DoubleVector accumulated;	
		RandomDevice rd;
		DefaultRandomEngine dre; 	
		UniformRealDistribution urd;	
	public:
		template<class ForwardIterator>
		DiscreteDistribution(ForwardIterator begin, ForwardIterator end){
			Real start = 0.0;
			while(begin != end){
			       start += *begin;
			       accumulated.push_back(start);
			       ++begin;
			}

			dre = DefaultRandomEngine(rd());
			if(!accumulated.empty())
				urd = UniformRealDistribution(0.0, accumulated.back());
			for(auto i = 0;i < accumulated.size();++ i)
				std::cout << accumulated[i] << std::endl;
		}
		int next(){
			if(accumulated.empty())
				return -1;
			Real tmp = urd(dre);
			int count = accumulated.size();
			int distance = std::upper_bound(accumulated.begin(), accumulated.end(), tmp) -
				accumulated.begin();
			if(distance >= count)
				return count - 1;
			return distance;
				

		}
				
	};

	/*
	 * select n numbers from the numbers: 0,1,2,.....,m-1
	 *
	 * */
	template<class Generator>
	void randomSelectDirect(unsigned int m, unsigned int n, unsigned int* result, Generator generator){
		if(n > m)
			n = m;
		//Is there some c++ 11 constructs to auto delete an array at the end of this function?
		unsigned int *numbers = new unsigned int[n];
		for(unsigned int i = 0;i < n;++ i)
			numbers[i] = i;

		for(unsigned int i = 0;i < n; ++ i){
			unsigned int index = i + generator() % (m - i);
			if(index != i){
				result[i] = numbers[index];
				result[index] = numbers[i];
			}else{
				result[i] = numbers[i];
			}
		}
		delete[] numbers;
	}

	/*
	 * select n numbers from the numbers: 0,1,2,.....,m-1
	 * more optimized solution for the case when n/m is very small
	 * but not efficient when n/m is close or greater than 1.
	 * */
	template<class Generator>
	void randomSelectOptimized(unsigned int m, unsigned int n, unsigned int* result, Generator generator){
		if(n > m)
			n = m;
		std::unordered_map<unsigned int, unsigned int> array;
		for(unsigned int i = 0;i < n; ++ i){
			unsigned int index = i + generator() % (m - i);
			auto pos = array.find(i);
			if(index != i){
				auto pos1 = array.find(index);
				if(pos1 == array.end()){
					result[i] = index;
					if(pos == array.end()){
						array[index] = i;
					}else{
						array[index] = pos->second;
						array.erase(pos);
					}
				}else if(pos == array.end()){
					//pos1 != pos->second
					result[i] = pos1->second;
					pos1->second = i;
				}else{
					result[i] = pos1->second;
					if(pos->second == index)
						array.erase(pos1);
					else
						pos1->second = pos->second;
					array.erase(pos);
				}
			}else if(pos == array.end()){
				result[i] = i;
			}else{
				result[i] = pos->second;
				array.erase(pos);
			}

		}
	}

	template<class Generator>
	void randomSelect(unsigned int m, unsigned int n, unsigned int* result, Generator generator){
		const double threshold = 0.8;
		if(n > m)
			n = m;
		if(n < m * threshold)
			randomSelectOptimized(m, n, result, generator);
		else
			randomSelectDirect(m, n, result, generator);


	}

	/*
	 * For a given integer n, this function returns the largest number a,
	 * such that a * a <= n
	 * */
	template<class T>
	T floorSqrt(T n){
		if(n == 0)
			return 0;
		T a = 1;
		while(a * a < n)
			a <<= 1;
		if(a * a == n)
			return a;
		a >>= 1;

		T b = a;
		b >>= 1;
		while(b != 0){
			T m = a + 1;
			T s = m * m;
			if(s == n)
				return m;
			if(s > n)
				return a;

			m = a + b;
			s = m * m;
			if(s == n)
				return m;
			if(s < n)
				a = m;
			b >>= 1;
		}
		return a;

	}

	/*
	 * For a given integer n, this function returns the largest number a,
	 * such that a*(a+1)/2 <=n
	 */
	template<class T>
	T floorTriSqrt0(T n){
		return (floorSqrt(8 * n + 1) - 1) / 2;
	}

	/*
	 * For a given integer n, this function returns the largest number a,
	 * such that a*(a-1)/2 <=n
	 */
	template<class T>
	T floorTriSqrt1(T n){
		return (floorSqrt(8 * n + 1) + 1) / 2;
	}

}
#endif //MY_LIB_MY_MATH_H
