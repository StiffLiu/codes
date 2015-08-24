#include "my_string.h"
#include "my_rand_string_generator.h"
#include "my_test_base.h"
#include <cassert>
#include <cstdlib>

namespace my_lib{
class SubstrSearchCompare : public StatPlot<SubstrSearchCompare&> {
	typedef my_lib::StatPlot<SubstrSearchCompare&> Super;
	mutable unsigned int current = 10;
	class ACString : public std::string{
	public:
		mutable unsigned int count = 0;
		auto operator[](unsigned int index) const ->decltype(std::string::operator[](index)) {
			assert(index < size());
			++ count;
			return std::string::operator[](index);
		}
	};
public:
	SubstrSearchCompare() 
		: Super(5, *this){
		Super::colors.assign({
			1, 0, 1, 
			0, 1, 0, 
			0, 1, 1, 
			1, 1, 0, 
			1, 1, 1});
	}
	void generate(std::string& src, std::string& pat) const {
		bool mustMatch = rand() % 2;
		unsigned int patLength = rand() % (current / 2) + 1;
		// srcLength >= patLength;
		unsigned int srcLength = current - patLength;

		RandStringGenerator srcGen(srcLength, srcLength);

		srcGen.randStringNoWhiteSpaces(src);
		if(mustMatch){
			unsigned int startIndex = rand() % (srcLength - patLength + 1);
			pat = src.substr(startIndex, patLength);
		}else{
			RandStringGenerator patGen(patLength, patLength);
			patGen.randStringNoWhiteSpaces(pat);
		}
	}
	bool operator()(double *values) const {
		ACString src, pat;
		generate(src, pat);

		typedef unsigned int (*SearchAlgo)(const ACString&, const ACString&);
		SearchAlgo algos[] = {bruteForceSearch, kmpSearch, boyerMooreSearch, rabinKarpSearch};
		for(unsigned int i = 0;i < (sizeof(algos) / sizeof(*algos));++ i){
			ACString s(src), p(pat);
			clock_t start = clock();
			algos[i](s, p);
			values[2 * i]= current;
			// values[2 * i + 1] = (i == 3 ? s.count + p.count : 0);
			values[2 * i + 1] = (i == 3 ? (clock() - start) : 0);
		}
		values[8] = current;
		values[9] = 0;
		++current;
		return true;
	}
};
}

int main(int argc, char *argv[]){
	my_lib::SubstrSearchCompare compare;
	compare.run(argc, argv);
	return 0;
}
