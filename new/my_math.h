#ifndef MY_LIB_MY_MATH_H
#define MY_LIB_MY_MATH_H
#include <vector>
#include <cassert>
#include <random>
#include <algorithm>
#include <iostream>

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
}
#endif //MY_LIB_MY_MATH_H
