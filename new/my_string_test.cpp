#include "my_rand_string_generator.h"
#include "my_access_counted_object.h"
#include "my_string.h"
#include <vector>
#include <iostream>
#include <ctime>
#include <algorithm>
#include <fstream>
#include <map>
#include <set>

#define MEASURETIME(statement, desc) \
{\
	clock_t start = clock();\
	statement;\
	cout << desc << (clock() - start) << endl;\
}

void generateStrings(std::vector<std::string>& strs){
	using namespace my_lib;
	const unsigned int count = 10000;
	RandStringGenerator gen(10, 100);
	
	strs.resize(count);
	for(unsigned int i = 0;i < count;++ i) gen.randStringInt(strs[i]);
}

std::vector<std::string> generateStrings(){
	std::vector<std::string> strs;
	generateStrings(strs);
	std::ifstream is("/backup/documents/books/Algorithms 4th Edition/algs4-data/mobydick.txt");
	std::string one;
	while(is >> one)
		strs.push_back(one);
	return strs;
}

struct StrArrayTraits{

	static unsigned int radix(){
		return 256;
	}

	static unsigned int cutoff(){
		return 35;
	}

	template<class Str>
	static int compare(Str str1, Str str2, unsigned int index){
		while(str1[index] && str1[index] == str2[index])
			++index;
		return str1[index] - str2[index];
	}

	template<class Str>
	static void swap(Str& str1, Str& str2){
		std::swap(str1, str2);
	}
	
	template<class Str>
	static unsigned int index(Str& str, unsigned int i){
		return static_cast<unsigned int>(str[i]) <= 0 ? -1 : str[i];
	}

	//assuming that i won't go out of scope of the c-string
	template<class Str>
	static auto addr(Str& str, unsigned int i) -> decltype(&str[i]){
		return str[i] == 0 ? nullptr : &str[i];
	}

	template<class Str>
	static void swap(Str array[], unsigned int i, unsigned int j){
		std::swap(array[i], array[j]);
	}
};
int test_msd_sort(int argc, char *argv[]){
	using namespace my_lib;
	using namespace std;
	std::vector<std::string> strs = generateStrings();
	const unsigned int count = strs.size();
#define USE_PLAIN_ARRAY
#ifdef USE_PLAIN_ARRAY
	typedef const char * StrArray;
#else
	typedef Incrementor<unsigned long> Counter;
	typedef AccessCountedArray<Counter, const char *> StrArray;
	unsigned long msdCount = 0;
	unsigned long q3wCount = 0;
	unsigned long stdCount = 0;
	Counter msdCounter(&msdCount);
	Counter q3wCounter(&q3wCount);
	Counter stdCounter(&stdCount);
#endif
	auto cstrs = new StrArray[count];

	//cout << "---------------------before msd sort---------------------" << endl;
#ifdef USE_PLAIN_ARRAY
	for(unsigned int i = 0;i < count;++ i) /*cout << (*/cstrs[i] = strs[i].c_str()/*) << std::endl*/;
	MEASURETIME(CStrMSD::sort(cstrs, count), "msd time:")
#else
	for(unsigned int i = 0;i < count;++ i) /*cout << (*/cstrs[i] = StrArray{msdCounter, strs[i].c_str()}/*) << std::endl*/;
	MEASURETIME((MSD<StrArray, StrArrayTraits>::sort(cstrs, count)), "msd time:")
	cout << "msd count: " << msdCount << endl;
#endif
	//cout << "---------------------after  msd sort---------------------" << endl;
	//for(unsigned int i = 0;i < count;++ i) cout << cstrs[i] << std::endl;

	//cout << "---------------------before quick 3 way sort---------------------" << endl;
#ifdef USE_PLAIN_ARRAY
	for(unsigned int i = 0;i < count;++ i) /*cout << (*/cstrs[i] = strs[i].c_str()/*) << std::endl*/;
	MEASURETIME(CStrQ3WS::sort(cstrs, count), "quick 3 way:")
#else
	for(unsigned int i = 0;i < count;++ i) /*cout << (*/cstrs[i] = StrArray{q3wCounter, strs[i].c_str()}/*) << std::endl*/;
	MEASURETIME((Quick3WaySort<StrArray, StrArrayTraits>::sort(cstrs, count)), "quick 3 way:")
	cout << "quick 3 way count: " << q3wCount << endl;
#endif
	//cout << "---------------------after  quick 3 way sort---------------------" << endl;
	//for(unsigned int i = 0;i < count;++ i) cout << cstrs[i] << std::endl;

	//cout << "---------------------before std sort---------------------" << endl;
#ifdef USE_PLAIN_ARRAY
	for(unsigned int i = 0;i < count;++ i) /*cout << (*/cstrs[i] = strs[i].c_str()/*) << std::endl*/;
	MEASURETIME(std::sort(cstrs, cstrs + count, 
		[](const char *s1, const char *s2){return std::strcmp(s1, s2) < 0;}), "std sort:")
#endif
	//for(unsigned int i = 0;i < count;++ i) cout << cstrs[i] << std::endl;
	
	
	delete[] cstrs;

	return 0;
}

int test_trie(int argc, char *argv[]){
	using namespace my_lib;
	using namespace std;
	Trie<const char *, unsigned int, CStrTraits> trie;
	TernarySearchTrie<char, unsigned int> tst;
	std::map<std::string, unsigned int> strMap;
	std::vector<std::string> toBeRemoved;
	std::vector<std::string> keys;
	RandStringGenerator gen(10, 100);
	trie.put("", 0);
	tst.put(string(), 0);
	strMap[""] = 0;
	for(unsigned int i = 0;i < 100;++ i){
		std::string one;
		gen.randStringNoWhiteSpaces(one);
		strMap[one] = i;
		trie.put(one.c_str(), i);
		tst.put(one, i);
		if(i < 20)
			toBeRemoved.push_back(one);
	}
	cout << "trie usage: " << trie.usage()  << std::endl;
	assert(trie.size() == strMap.size());
	cout << "tst size : " << tst.size() << std::endl;
	assert(tst.size() == strMap.size());

	for(auto& kvp : strMap){
		auto value = trie.get(kvp.first.c_str());
		std::cout << kvp.first << std::endl;
		assert(value != nullptr);
		assert(*value == kvp.second);
		trie.put(kvp.first.c_str(), kvp.second * 2);
		value = tst.get(kvp.first);
		assert(value != nullptr);
		assert(*value == kvp.second);
		tst.put(kvp.first, kvp.second * 2);
	}
	trie.keysWithPrefix(keys, std::string(""));
	cout << "key with prefix \"\"" << std::endl;
	assert(keys.size() == strMap.size());
	for(auto& k : keys){
		assert(strMap.find(k) != strMap.end());
		std::cout << k << std::endl;
	}
	keys.clear();
	tst.keysWithPrefix(keys, std::string(""));
	assert(keys.size() == strMap.size());
	for(auto& k : keys){
		assert(strMap.find(k) != strMap.end());
	}
	keys.clear();


	trie.keysWithPrefix(keys, std::string("b"));
	cout << "\nkey with prefix \"b\"" << std::endl;
	for(auto& k : keys) std::cout << k << std::endl;
	keys.clear();
	tst.keysWithPrefix(keys, std::string("b"));
	cout << "\nkey with prefix \"b\"" << std::endl;
	for(auto& k : keys) std::cout << k << std::endl;

	assert(trie.size() == strMap.size());
	for(auto& v : toBeRemoved){
		trie.remove(v.c_str());
		assert(tst.remove(v));
		assert(trie.get(v.c_str()) == nullptr);
		assert(tst.get(v) == nullptr);
		strMap.erase(v);
	}
	for(auto& kvp : strMap){
		auto value = trie.get(kvp.first.c_str());
		assert(value != nullptr && *value == 2 * kvp.second);
		value = tst.get(kvp.first);
		assert(value != nullptr && *value == 2 * kvp.second);
	}
	assert(trie.size() == strMap.size());
	for(unsigned int i = 0;i < 20;++ i){
		std::string one;
		gen.randStringNoWhiteSpaces(one);
		//cout << "search : " << one << std::endl;
		//cout << "std map contains : " << (strMap.find(one) != strMap.end()) << std::endl;
		//cout << "trie contains : " << trie.contains(one.c_str()) << std::endl;
		if(trie.contains(one.c_str()))
			cout << "trie value is : " << *trie.get(one.c_str());
		if(tst.contains(one))
			cout << "tst value is : " << *trie.get(one.c_str());
		assert((strMap.find(one) != strMap.end()) == trie.contains(one.c_str()));
		assert((strMap.find(one) != strMap.end()) == tst.contains(one));
	}

	assert(trie.size() == strMap.size());
	assert(tst.size() == strMap.size());
	for(auto& kvp : strMap){
		trie.remove(kvp.first.c_str());
		tst.remove(kvp.first);
	}
	assert(trie.size() == 0);
	assert(trie.isEmpty());
	assert(tst.size() == 0);
	assert(tst.isEmpty());

	std::string key;
	for(const auto& str : {string("ade"), string("abc"), string("fsdfe"), string("abce"), string("abcd")}){
		trie.put(str.c_str(), 0);
		tst.put(str, 0);
		
		auto value = trie.get(str.c_str());
		assert(nullptr != value && *value == 0);
		value = tst.get(str);
		assert(nullptr != value && *value == 0);
	}
	bool exists = trie.longestPrefixOf("abcedfs", key);
	cout << "longest prefix of \"abcedfs\" is " << (!exists ? "[NON EXISTS]" : key) << endl;
	exists = tst.longestPrefixOf(std::string("abcdfsf"), key);
	cout << "longest prefix of \"abcdfsf\" is " << (!exists ? "[NON EXISTS]" : key) << endl;

	return 0;
}

int test_substr_search(int argc, char *argv[]){
	unsigned int numTests = 100;
	using namespace my_lib;
	RandStringGenerator gen(20, 100);
	for(unsigned int i = 0;i < numTests;++ i){
		std::string src;
		gen.randStringNoWhiteSpaces(src);

		bool mustMatch = rand() % 2;
		std::string pat;
		size_t matchIndex = std::string::npos;
		if(mustMatch){
			unsigned int startIndex = rand() % (src.size() - 5);
			unsigned int len = rand() % (src.size() - startIndex);
			pat = src.substr(startIndex, len);
		}else{
			gen.randStringNoWhiteSpaces(pat);
		}
		matchIndex = src.find(pat);
		std::cout << "expected match index is : " << (long long)(matchIndex == std::string::npos ? -1 : matchIndex) << std::endl;

		unsigned int result = 0;
		std::cout << "Search \n" << pat << "\n\t in \n" << src << std::endl;
		std::cout << "bruteForce : " << (result = bruteForceSearch(src, pat)) << std::endl;
		assert((matchIndex == std::string::npos && result == src.size()) || matchIndex == result);
		std::cout << "kmp : " << (result = kmpSearch(src, pat)) << std::endl;
		assert((matchIndex == std::string::npos && result == src.size()) || matchIndex == result);
		std::cout << "boyerMoore : " << (result = boyerMooreSearch(src, pat)) << std::endl;
		assert((matchIndex == std::string::npos && result == src.size()) || matchIndex == result);
		std::cout << "rabinKarp : " << (result = rabinKarpSearch(src, pat)) << std::endl;
		assert((matchIndex == std::string::npos && result == src.size())|| matchIndex == result);
		std::cout << "--------------------------------------------------" << std::endl;

	}
	{
		using namespace std;
		string src("abcabcabcabcabcab"), pat("abcab");
		set<unsigned int> all;
		bruteForceSearch(src, pat, all);
		cout << "brureForce: "; for (auto index : all) cout << index << " "; cout << endl; all.clear();
		kmpSearch(src, pat, all);
		cout << "kmp: "; for (auto index : all) cout << index << " "; cout << endl; all.clear();
		boyerMooreSearch(src, pat, all);
		cout << "boyerMoore: "; for (auto index : all) cout << index << " "; cout << endl; all.clear();
		rabinKarpSearch(src, pat, all);
		cout << "rabinKarp: "; for (auto index : all) cout << index << " "; cout << endl; all.clear();
	}

	return 0;
}

int main(int argc, char *argv[]){
	for(unsigned int i = 0;i < 1;++ i){
		test_substr_search(argc, argv);
	}
	return 0;
}
