#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <iterator>
#include <cstring>
#include <map>
#include <cstdlib>
#include <ctime>

/*
 * Compute the Burrows-Wheeler transform for the string: [begin, end)
 * The result is outputed to the iterator {@param output}
 * And the returned iterator is the split position in {@param output}
 *
 * Take this string for example: [begin, end)=mississippi
 * The Burrows-Wheeler transform will be: ipssm$pissii
 * 		where "$" is the imaginary end-of-file character.
 *
 * The characters will be put to {@param output} iterator in this sequence : ipssmpissii
 * And the iterator {@param output + 5} will be returned.
 * */
template<class InputIterator, class OutputIterator, class Comparator >
OutputIterator computeBWT(InputIterator begin, InputIterator end, OutputIterator output, Comparator comparator){
	if(begin == end) return output;
	//auto beginBackup = begin;

	// iterators in this vector points to the last elements of 
	// a cyclic rotation, 
	std::vector<InputIterator> suffixArray;
	auto func = [begin, end, comparator](InputIterator i1, InputIterator i2){
		if(i1 == i2) return false;
		if(i1 == end){
		 	i1 = begin;
			++ i2;
		}else if(i2 == end){
		 	i2 = begin;
			++ i1;
		}else{
			++ i1;
			++ i2;
		}
		assert(i1 != i2);
		while(i1 != end && i2 != end){
			if(comparator(*i1, *i2)) return true;
			if(comparator(*i2, *i1)) return false;
			++i1;
			++i2;
		}
		if(i1 == end){
			assert(i2 != end);
			return true;
		}
		assert(i2 == end && i1 != end);
		return false;
	};

	while(begin != end){
		suffixArray.push_back(begin);
		++begin;
	}
	suffixArray.push_back(end);
	std::sort(suffixArray.begin(), suffixArray.end(), func);
	//std::cout << "=========" << std::endl;
	OutputIterator split = output;
	for(auto& itor : suffixArray){
		/*if(itor == end){
			std::cout << beginBackup << std::endl;
		}else{
			std::cout <<(itor + 1) << std::endl;
		}*/
		if(itor != end){
			*output = *itor;
			++output;
		}else{
			split = output;
		}
	}
	return split;
}

/*
 * Compute the Burrows-Wheeler inverse transform for the string: [split, end) + [end, begin]
 * returns true  if exits, else false.
 * The result is outputed to {@param output} if exists.
 *
 * Take this input for example: [begin, end)=ipssmpissii   split=begin+5
 * That the Burrows-Wheeler tranform of the original string is ipssm$pissii
 * 		where "$" is the imaginary end-of-file character.
 * The original string should be: mississippi
 * The characters will be put to {@param output} iterator in this sequence : mississippi
 * And {@code true} will be returned.
 * */
template<class InputIterator, class OutputIterator, class Comparator>
bool computeBWI(InputIterator split, InputIterator begin, InputIterator end,
	OutputIterator output, Comparator comparator){
	typedef typename std::remove_reference<typename std::remove_cv<decltype(*begin)>::type>::type Element;
	std::map<Element, unsigned int, Comparator> counts(comparator);
	std::vector<InputIterator> iterators;

	for(auto start = begin;start != end;++start){
		if(start == split) iterators.push_back(end);
	 	++ counts[*start];
		iterators.push_back(start);
	}

	unsigned int count = 1;
	{
		for(auto& kvp : counts){
			unsigned int start = count;
			count += kvp.second;
			kvp.second = start;
		}
	}

	assert(count == iterators.size());
	std::vector<unsigned int> seq;
	{
		std::vector<unsigned int> link(count);
		for(size_t i = 0;i < iterators.size();++ i){
			if(iterators[i] != end){
				auto& start = counts[*iterators[i]];
				if(i == start){
					std::cout << i << std::endl;
				 	return false;
				}
				link[i] = start;
				++start;
			}
		}
		//std::cout << "link ";
		//for(auto index : link) std::cout << index << ' ';
		//std::cout << std::endl;
		unsigned int next = 0;
		while(iterators[next] != end){
			seq.push_back(next);
			next = link[next];
		}
		//std::cout << "seq ";
		//for(auto index : seq) std::cout << index << ' ';
		//std::cout << std::endl;
	}
	while(!seq.empty()){
		*output = *iterators[seq.back()];
		seq.pop_back();
		++output;
	}
	return true;
}

template<class InputIterator, class OutputIterator>
OutputIterator computeBWT(InputIterator begin, InputIterator end, OutputIterator output){
	typedef typename std::remove_reference<typename std::remove_cv<decltype(*begin)>::type>::type Element;
	return computeBWT(begin, end, output, std::less<Element>());
}

template<class InputIterator, class OutputIterator>
bool computeBWI(InputIterator split, InputIterator begin, InputIterator end, OutputIterator output){
	typedef typename std::remove_reference<typename std::remove_cv<decltype(*begin)>::type>::type Element;
	return computeBWI(split, begin, end, output, std::less<Element>());
}

int test_btw_bti(const char *str){
	unsigned int count = strlen(str);
	char btw[count] = {0};
	char bti[count] = {0};
	auto split = computeBWT(str, str + count, btw);
	std::cout << "original string is : " << str << std::endl;
	std::cout << "btw is : ";
	std::copy(btw, split, std::ostream_iterator<char>(std::cout, ""));
	std::cout << "$";
	std::copy(split, btw + count, std::ostream_iterator<char>(std::cout, ""));
	std::cout << std::endl;
	assert(computeBWI(split, btw, btw + count, bti));
	std::cout << "bti is : ";
	std::copy(bti, bti + count, std::ostream_iterator<char>(std::cout, ""));
	std::cout << std::endl;
	return 0;
}

int main(int argc, char *argv[]){
	//return test_btw_bti("mississippi");
	srand(time(0));
	unsigned int count = 70;
	char str[count + 1] = {0};
	for(unsigned int i = 0;i < count;++ i) str[i] = (rand() % 10) + 'a';
	return test_btw_bti(str);
}
