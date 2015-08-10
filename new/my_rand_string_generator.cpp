#include "my_rand_string_generator.h"
#include <cctype>
#include <random>
#include <cassert>

namespace{
struct CharacterClasses{
	std::string alphas;
	std::string digits;
	std::string hexDigits;
	std::string alphaNums;
	std::string printables;
	std::string printableNoWhiteSpaces;
	CharacterClasses(){
		digits = "0123456789";
		hexDigits = "0123456789ABCDEF";
		for(int i = 0;i < 128;++ i){
			if(std::isprint(char(i))){
				printables.push_back(char(i));
				if(!std::isspace(char(i))){
					printableNoWhiteSpaces.push_back(char(i));
					if(std::isalpha(char(i)))
						alphas.push_back(char(i));

				}
			}
		}
		alphaNums = alphas + digits;
	}
	
	static void generate(unsigned int minLen, unsigned int maxLen, 
			const char *chars, unsigned int count, std::string& str){
		std::random_device rd;
		assert(count > 0);
		std::uniform_int_distribution<unsigned int> indexGenerator(0, count - 1); 
		std::uniform_int_distribution<unsigned int> lengthGenerator(minLen, maxLen);
		auto len = lengthGenerator(rd);
		decltype(len) i = str.size();
		str.resize(str.size() + len);
		len = str.size();
		for(;i < len;++ i)
			str[i] = chars[indexGenerator(rd)];
	}
}instance;
}
namespace my_lib{
void RandStringGenerator::randStringAlpha(std::string& str){
	CharacterClasses::generate(minLen, maxLen, instance.alphas.c_str(),
		instance.alphas.size(), str);
}
void RandStringGenerator::randStringInt(std::string& str){
	CharacterClasses::generate(minLen, maxLen, instance.digits.c_str(),
		instance.digits.size(), str);
}

void RandStringGenerator::randStringHex(std::string& str){
	CharacterClasses::generate(minLen, maxLen, instance.hexDigits.c_str(),
		instance.hexDigits.size(), str);
}

void RandStringGenerator::randStringBin(std::string& str){
	CharacterClasses::generate(minLen, maxLen, "01", 2, str);
}

void RandStringGenerator::randStringAlphaNum(std::string& str){
	CharacterClasses::generate(minLen, maxLen, instance.alphaNums.c_str(),
		instance.alphaNums.size(), str);
}

void RandStringGenerator::randStringPrintable(std::string& str){
	CharacterClasses::generate(minLen, maxLen, instance.printables.c_str(),
		instance.printables.size(), str);
}

void RandStringGenerator::randStringNoWhiteSpaces(std::string& str){
	CharacterClasses::generate(minLen, maxLen, instance.printableNoWhiteSpaces.c_str(),
		instance.printableNoWhiteSpaces.size(), str);
}

void RandStringGenerator::randString(std::string& str){
		std::random_device rd;
		std::uniform_int_distribution<int> charGenerator(0, 127); 
		std::uniform_int_distribution<unsigned int> lengthGenerator(minLen, maxLen);
		auto len = lengthGenerator(rd);
		for (decltype(len) i = 0; i < len; ++i)
			str.push_back(static_cast<char>(charGenerator(rd)));
}

}
