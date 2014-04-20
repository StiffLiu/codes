#ifndef MY_RAND_STRING_GENERATOR_H
#define MY_RAND_STRING_GENERATOR_H
#include <string>

namespace my_lib{
/**
 * A random string generator, that could be used to generated different kinds of strings.
 */
class RandStringGenerator{
	unsigned int minLen = 0, maxLen = 10;
public:
	/**
	 * @param minLen Minimum length of the random string.
	 * @param maxLen Maximum length of the random string.
	 */
	RandStringGenerator(unsigned int minLen = 0, unsigned int maxLen = 10){
		if(minLen > maxLen){
			this->minLen = maxLen;
			this->maxLen = minLen;
		}else{
			this->minLen = minLen;
			this->maxLen = maxLen;
		}
	}

	/**
	 * Generate a random string that contains only letters 
	 * and whose length is between {@var minLen} and {@var maxLen}.
	 * @param str The generated random string.
	 */
	void randStringAlpha(std::string& str);

	/**
	 * Generate a string that contains only decimal digits
	 * and whose length is between {@var minLen} and {@var maxLen}.
	 * @param str The generated random string.
	 */
	void randStringInt(std::string& str);

	/**
	 * Generate a string that contains only hex digits
	 * and whose length is between {@var minLen} and {@var maxLen}.
	 * @param str The generated random string.
	 */
	void randStringHex(std::string& str);

	/**
	 * Generate a random string that contains only binary digits
	 * and whose length is between {@var minLen} and {@var maxlen}.
	 * @param str The generated random string.
	 */
	void randStringBin(std::string& str);

	/**
	 * Generate a random string that contains only letters and digits.
	 * and whose length is between {@var minLen) and {@var maxLen}.
	 * @param str The generated random string.
	 */
	void randStringAlphaNum(std::string& str);

	/**
	 * Generate a random string that contains only printable characters 
	 * and whose length is between {@var minLen} and {@var maxLen}.
	 * @param str The generated random string.
	 */
	void randStringPrintable(std::string& str);

	/**
	 * Generate a random string that contains only printable characters that are
	 * not white spaces and whose length is between {@var minLen} and {@var maxLen}.
	 * @param str The generated random string.
	 */
	void randStringNoWhiteSpaces(std::string& str);

	/**
	 * Generate a random string whose length is in the range {@code [minLen, maxLen]}.
	 * @param str The generated random string.
	 */
	void randString(std::string& str);
};
}
#endif //MY_RAND_STRING_GENERATOR_H
