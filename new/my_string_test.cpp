#include "my_rand_string_generator.h"
#include "my_access_counted_object.h"
#include "my_string.h"
#include <vector>
#include <iostream>
#include <ctime>
#include <algorithm>
#include <fstream>
#include <streambuf>
#include <map>
#include <set>
#include <cmath>

#define MEASURETIME(statement, desc) \
{\
	clock_t start = clock();\
	statement;\
	std::cout << desc << (clock() - start) / (double)CLOCKS_PER_SEC << std::endl;\
}

void generateStrings(std::vector<std::string>& strs){
	using namespace my_lib;
	const unsigned int count = 1000000;
	RandStringGenerator gen(10, 100);
	
	strs.resize(count);
	for(unsigned int i = 0;i < count;++ i) gen.randStringAlpha(strs[i]);
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
	assert(std::is_sorted(cstrs, cstrs + count, [](const char* s1, const char* s2){return strcmp(s1, s2) < 0;}));
	//cout << "---------------------after  msd sort---------------------" << endl;
	//for(unsigned int i = 0;i < count;++ i) cout << cstrs[i] << std::endl;

	//cout << "---------------------before quick 3 way sort---------------------" << endl;
#ifdef USE_PLAIN_ARRAY
	for(unsigned int i = 0;i < count;++ i) /*cout << (*/cstrs[i] = strs[i].c_str()/*) << std::endl*/;
	MEASURETIME(CStrQ3WS::sort(cstrs, count), "quick 3 way:")
#else
	for(unsigned int i = 0;i < count;++ i) /*cout << (*/cstrs[i] = StrArray{q3wCounter, strs[i].c_str()}/*) << std::endl*/;
	MEASURETIME((Quick3WaySort<StrArrayTraits>::sort(cstrs, count)), "quick 3 way:")
	cout << "quick 3 way count: " << q3wCount << endl;
#endif
	assert(std::is_sorted(cstrs, cstrs + count, [](const char* s1, const char* s2){return strcmp(s1, s2) < 0;}));
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
	{
		std::cout << std::endl << std::endl << std::endl;
		std::string pattern = "it is a far far better thing that i do than i have ever done";
		std::ifstream tale("F:\\documents\\books\\Algorithms 4th Edition\\algs4-data\\tale.txt");
		//std::ifstream tale("/f/documents/books/Algorithms 4th Edition/algs4-data/tale.txt");
		std::string src{std::istreambuf_iterator<char>(tale), std::istreambuf_iterator<char>()};
		unsigned int (*func[])(const std::string&, const std::string&) = 
			{bruteForceSearch, kmpSearch, boyerMooreSearch, rabinKarpSearch};
		const char *names[] = {"bruteForce ", "kmp ", "boyerMoore ", "rabinKarp "};
		for(unsigned int i = 0;i < sizeof(names) / sizeof (*names);++ i){
			unsigned int index = 0;
			MEASURETIME(index = func[i](src, pattern), names[i])
			std::cout << names[i] << "result is : " << index << std::endl;
		}
	}

	std::cout << "================================" << std::endl;
	RegExpr expr("abcd.*((f(dfdf|sdf))|cdfd|)*.*e");
	std::cout << expr.recognize("abcdcdffsdfcdfdfsdfaaae") << std::endl;
	return 0;
}

int test_2d_search(int argc, char *argv[]){
	const unsigned int m = 100, n = 50, h = 20, v = 30;
	unsigned int src[m][n], pat[h][v];
	using namespace my_lib;
	using namespace std;
	for(unsigned int i = 0;i < m;++ i)
		for(unsigned int j = 0;j < n;++ j)
			src[i][j] = rand();

	unsigned rowOffset = 10, colOffset = 6;
	for(unsigned int i = 0;i < h;++ i)
		for(unsigned int j = 0;j < v;++ j)
			pat[i][j] = src[i + rowOffset][j + colOffset];
	std::pair<unsigned int, unsigned int> result;
	result = bruteForceSearch(src, m, n, pat, h, v);
	cout << "bruteForce: (" << result.first << "," << result.second << ")" << endl;
	result = rabinKarpSearch(src, m, n, pat, h, v);
	cout << "rabinKarp: (" << result.first << "," << result.second << ")" << endl;

	return 0;
}

void compare_compress(std::vector<unsigned char>& chars, unsigned int blockSize = 8){
	using namespace my_lib;
	std::vector<bool> bitArray;
	std::vector<unsigned char> runLength;
	std::vector<bool> compressed;
	copy(chars.begin(), chars.end(), bitArray);
	RunLength::compress(bitArray, runLength);
	double runLengthRatio = (32 + 
		log(*std::max_element(runLength.begin(), runLength.end()) + 1) * runLength.size()) / bitArray.size();
	Huffman::compress(bitArray, compressed, blockSize);
	double huffmanRatio = compressed.size() / (double)bitArray.size();
	compressed.clear();
	LZW::compress(bitArray, compressed, blockSize, blockSize + 1);
	double lzwRatio = compressed.size() / (double)bitArray.size();
	std::cout << "(runLength, huffman, lzw)=(" << runLengthRatio << "," << huffmanRatio << "," << lzwRatio << ")" << std::endl;
}

int test_compress(int argc, char *argv[]){
	using namespace my_lib;
	{
		unsigned char ch = 0b01010101;
		std::vector<unsigned char> chars({ch, ch, ch, ch, ch, ch, ch, ch, ch, ch});
		std::vector<bool> bitArray;
		std::vector<unsigned char> runLength;
		std::vector<bool> decompressed;
		
		copy(chars.begin(), chars.end(), bitArray);
		BASplitor spr;
		std::cout << "bitArray : " << spr << bitArray << std::endl;
		RunLength::compress(bitArray, runLength);
		std::cout << "Run length compress : ";
		for(auto length : runLength) std::cout << (unsigned int)length << ' ';
		std::cout << std::endl;
		RunLength::decompress(runLength, decompressed);
		std::cout << "After decompress : " << spr << decompressed << std::endl;

		std::vector<bool> compressed;
		Huffman::compress(chars, compressed, 1);
		std::cout << "Huffman compress : " << spr << compressed << std::endl;
		decompressed.clear();
		Huffman::decompress(compressed, decompressed, 8);
		std::cout << "After decompress : " << spr << decompressed << std::endl;

		compressed.clear();
		decompressed.clear();
		LZW::compress(bitArray, compressed, 2, 3);
		std::cout << "LZW compress : " << spr << compressed << std::endl;
		LZW::decompress(compressed, decompressed, 2, 3);
		std::cout << "After decompress : " << spr << decompressed << std::endl;
	}
	unsigned int blockSize = 8;
	if (argc > 1) blockSize = atoi(argv[1]);
	if(blockSize < 2) blockSize = 2;
	std::cout << "block size is : " << blockSize << std::endl;
	{
		std::cout << "===========random ascii string===========" << std::endl;
		const char *alphabet = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";
		for(unsigned int i = 1000;i < 100000;i += 1000)
		{
			std::vector<unsigned char> chars;
			for(unsigned int j = 0;j < i;++j) chars.push_back(alphabet[rand() % 52]);
			compare_compress(chars, blockSize);
		}
	}
	{
		std::cout << "===========repeating string of \"ab\"===========" << std::endl;
		for(unsigned int i = 1;i < 100000;i += 1000)
		{
			std::vector<unsigned char> chars;
			for(unsigned int j = 0;j < i;++j){
			 	chars.push_back('a');
			 	chars.push_back('b');
			}
			compare_compress(chars, blockSize);
		}
	}

	return 0;
}

// longest repeated sub-string.
std::pair<size_t,size_t> lrss(const std::string& str){
	using namespace my_lib;
	auto func=[](const std::string& s1, size_t i1, const std::string& s2, size_t i2){
		return std::strcmp(s1.c_str() + i1, s2.c_str() + i2) < 0;
	};
	SuffixArray<const std::string&, decltype(func)> sa(str, str.size(), func);
	size_t len = 0, index = static_cast<size_t>(-1);
	for(size_t i = 0;i + 1 < sa.size();++ i){
		size_t idx1 = sa.index(i), idx2 = sa.index(i + 1);
		if(idx1 + len < sa.size() && idx2 + len < sa.size()){
			size_t lcp = sa.lcp(i);
			if(lcp > len){
				len = lcp;
				index = idx1;
			}
		}
	}
	return {index, len};
}

int lrss_client(int argc, char *argv[]){
	char ch;
	std::string str;
	while(std::cin.get(ch)) str.push_back(ch);
	decltype(lrss(str)) ret;
	MEASURETIME(ret = lrss(str), "longest repeated sub-string search takes(seconds) : ");
	if(static_cast<size_t>(-1) != ret.first){
		std::cout << "longest repeated sub-string is at index : " << ret.first << ", and length is :"
			<< ret.second << "\t" << str.substr(ret.first, ret.second) << std::endl;
	}
	return 0;
}
int main(int argc, char *argv[]){
	//lrss_client(argc, argv);
	test_msd_sort(argc, argv);
	return 0;
}
