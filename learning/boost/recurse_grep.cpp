#include <vector>
#include <string>
#include <regex>
#include <fstream>
#include <streambuf>
#include <boost\filesystem.hpp>

namespace{
	const char *href_pat = "<a\\s+.*href\\s*=\\s*\"([^\"]+)\"[^>]*>([^<]*)";
}
template<class Func>
void grep(const std::string& text_, const std::regex& expr_, Func func_){
	std::sregex_iterator begin(text_.begin(), text_.end(), expr_),
		end = std::sregex_iterator();
	for (std::sregex_iterator i = begin; i != end; ++i) {
		func_(*i);
	}
}

std::vector<std::string> grep_href(const std::string& text_){
	std::vector<std::string> matches;
	grep(text_, std::regex(href_pat), 
		[&matches](const std::smatch& match){
			for (size_t i = 0; i < match.size(); ++i){
				matches.emplace_back(match[i].str());
			}
		}
	);
	return matches;
}

template<class Func>
void grep_file(const char *file_, const std::regex& expr_, Func func_){
	std::ifstream is(file_);
	if (is)	grep({ std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>() }, expr_, func_);
}

template<class Func>
void recurve_grep(const char *path_, const std::regex& expr_, Func func_, const char * fileNamePattern_ = nullptr){
	using namespace boost::filesystem;
	recursive_directory_iterator begin(path_), end;
	std::regex fileNamePattern(nullptr == fileNamePattern_ ? ".*" : fileNamePattern_);
	while (begin != end){
		const auto& pth = begin->path();
		if (is_regular_file(pth) && std::regex_match(pth.filename().string(), fileNamePattern)){
			grep_file(begin->path().string().c_str(), expr_, func_);
		}
		++begin;
	};
}

#include <iostream>
std::vector<std::string>  recursive_grep_href(const char *path_, const char * fileNamePattern_ = nullptr){
	std::vector<std::string> matches;
	recurve_grep(path_, std::regex(href_pat),
		[&matches](const std::smatch& match){
		for (size_t i = 0; i < match.size(); ++i){
			matches.emplace_back(match[i].str());
		}
	},
	fileNamePattern_
	);
	return matches;
}

template<class Func>
std::vector<std::string> recursive_grep_href(const char *folder_, Func func_, const char * fileNamePattern_ = nullptr){
	return recurve_grep(folder_, std::regex(href_pat), func_, fileNamePattern_);
}
std::vector<std::string>  recursive_grep_href_from_html(const char *path_){
	return recursive_grep_href(path_, R"(.*\.html)");
}