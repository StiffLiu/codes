
#include "HeaderParser.h"
#include <map>
#include <iostream>
#include <sstream>

using namespace std;
using namespace sli;
#define ARRAYSIZE(ar) ((sizeof (ar))/(sizeof (*ar)))
/*the order of these string literals should be changed with care,
so that they are consistent with the enum defitions in the class Keyword.
The first empty string represents an invalid keyword
*/
const char *Keyword::keywords[] = {
	"", "auto", "extern", "register", "static", "typedef", 
	"volatile", "const", 
	"long", "short",
	"signed", "unsigned", 
	"char", "void", "double", "int", "float",
	"break", "case", "continue", "default", "do", "else", "for", "goto",
	"if", "while", "return", "switch",
	"sizeof",
	"enum", "struct", "union",
};
struct Scope::ScopeAttribute{
		int startToken;
		int endToken;
		int declarationEndToken;
		int varListSeparator;
		bool allowFunctionDefine;
		bool allowTypeDefine;
		bool allowVarList;
		bool allowCompondType;
		bool allowInitialization;
		bool allowBitType;
		bool needDeclarationeEndToken;
		bool needNewLine;
	};

static Scope::ScopeAttribute globalScopeAttribute;
static Scope::ScopeAttribute paramListAttribute;
static Scope::ScopeAttribute compoundTypeScopeAttribute;
namespace{
	typedef bool (*StringLess)(const char *, const char *);
	static bool isLessThan(const char *str1, const char *str2){
		return strcmp(str1, str2) < 0;
	}
	class Initializer{
		typedef map<const char *, int, StringLess> Map;
		Map keywordList;
	public:

		Initializer() : keywordList(isLessThan){
			for(int i = Keyword::INVALID;i < Keyword::KEYWORD_COUNT;++ i)
				keywordList[Keyword::keywords[i]] = i;
			globalScopeAttribute.startToken = Lexer::END;
			globalScopeAttribute.endToken = Lexer::END;
			globalScopeAttribute.declarationEndToken = ';';
			globalScopeAttribute.allowFunctionDefine = true;
			globalScopeAttribute.allowTypeDefine = true;
			globalScopeAttribute.allowVarList = true;
			globalScopeAttribute.allowCompondType = true;
			globalScopeAttribute.allowInitialization = true;
			globalScopeAttribute.varListSeparator = ',';
			globalScopeAttribute.allowBitType = false;
			globalScopeAttribute.needDeclarationeEndToken = true;
			globalScopeAttribute.needNewLine = true;

			paramListAttribute.startToken =  '(';
			paramListAttribute.endToken = ')';
			paramListAttribute.declarationEndToken = ',';
			paramListAttribute.allowFunctionDefine = false;
			paramListAttribute.allowTypeDefine = false;
			paramListAttribute.allowVarList = true;
			paramListAttribute.allowCompondType = false;
			paramListAttribute.varListSeparator = ',';
			paramListAttribute.allowBitType = false;
			paramListAttribute.needDeclarationeEndToken = false;
			
			compoundTypeScopeAttribute.startToken = '{';
			compoundTypeScopeAttribute.endToken = '}';
			compoundTypeScopeAttribute.declarationEndToken = ';';
			compoundTypeScopeAttribute.allowFunctionDefine = false;
			compoundTypeScopeAttribute.allowTypeDefine = false;
			compoundTypeScopeAttribute.allowVarList = true;
			compoundTypeScopeAttribute.allowCompondType = true;
			compoundTypeScopeAttribute.varListSeparator = ',';
			compoundTypeScopeAttribute.allowBitType = true;
			compoundTypeScopeAttribute.needDeclarationeEndToken = true;
			compoundTypeScopeAttribute.needNewLine = true;
		}
		Byte getIndex(const char *name) const{
			Map::const_iterator pos = keywordList.find(name);
			if(pos == keywordList.end())
				return Keyword::INVALID;
			return pos->second;
		}
	};
	const Initializer instance;
}

Keyword::Keyword(const char *name){
	index = Keyword::INVALID;
	if(name != NULL)
		index = instance.getIndex(name);
}

const PrimType* PrimType::intType(){
	static PrimType iType(Keyword::intKeyword());
	return &iType;
}
const PrimType* PrimType::charType(){
	static PrimType cType(Keyword::charKeyword());
	return &cType;
}
const PrimType* PrimType::doubleType(){
	static PrimType dType(Keyword::doubleKeyword());
	return &dType;
}
const PrimType* PrimType::floatType(){
	static PrimType fType(Keyword::floatKeyword());
	return &fType;
}
const PrimType* PrimType::voidType(){
	static PrimType vType(Keyword::voidKeyword());
	return &vType;
}
Node::~Node(){
}
int Lexer::next(BlockIterator& itor, bool escape){
	int ch = *itor;
	if(-1 == ch)
		return EOF;
	++itor;
	if(ch == '\\' && escape && *itor == '\n'){
		++itor;
		return next(itor, escape);
	}
	return ch;
}
void Lexer::readCharSequence(int delimiter, std::string& ret){
	int ch = *itor;
	while(-1 != ch){
			if(ch == '\\'){
				++ itor;
				ch = *itor;
				switch(ch){
						case '\n':;break;
						case 'n':ret.push_back('\n');break;
						case 't':ret.push_back('\t');break;
						case 'v':ret.push_back('\v');break;
						case 'b':ret.push_back('\b');break;
						case 'r':ret.push_back('\r');break;
						case 'f':ret.push_back('\f');break;
						case 'a':ret.push_back('\a');break;
						case '0':case '1':case '2':case '3':
						case '4':case '5':case '6':case '7':{
							int val = Utility::oct2dec(ch);
							BlockIterator next = itor;
							++ next;
							if(Utility::isOctDigit(*next)){
								ch = ch * 8 + Utility::oct2dec(*next);
								itor = next;
								++next;
								if(Utility::isOctDigit(*next)){
									ch = ch * 8 + Utility::oct2dec(*next);
									itor = next;
									++next;
								}
							}
							ret.push_back((char)ch);
							break;
						}
						case 'x':{
							BlockIterator next = itor;
							++ next;
							if(Utility::isHexDigit(*next)){
								int ch = Utility::hex2dec(*next);
								itor = next;
								++next;
								if(Utility::isHexDigit(*next)){
									ch = ch * 16 + Utility::hex2dec(*next);
									itor = next;
									++next;
								}
								ret.push_back((char)ch);
								break;
							}
						}
						default:ret.push_back(ch);break;
				}					
			}else if(ch == '\n'){
				cerr << "new line in constant line no :";
				exit(__LINE__);
			}
			else if(ch != delimiter)
				ret.push_back(ch);
			else
				break;
			++ itor;
			ch = *itor;
		}
		if(ch != delimiter){
			cerr << "invalid character sequence found";
			exit(__LINE__);
		}
		++ itor;
}
void Lexer::readMultiLineComment(string& ret){
	int ch = *itor;
	while(-1 != ch){
		if(ch == '*'){
			++itor;
			if(*itor == '/'){
				++ itor;
				return;
			}
			ch = *itor;
			ret.push_back('*');
			continue;
		}
		ret.push_back(ch);
		++itor;
		ch = *itor;
	}
	cerr << "invalid comment found";
	exit(__LINE__);
}
void Lexer::skipToEndOfLine(string& ret){
	int ch;
	BlockIterator tmp = itor;
	
	while((ch = next(tmp)) != EOF){
		if(ch == '\n')
			break;
		ret.push_back(ch);
		itor = tmp;
	}
}
void Lexer::skipSpace(){
	int ch = *itor;
	while(-1 != ch && ch != '\n' && ::isspace(ch)){
		++ itor;
		ch = *itor;
	}
}

void Lexer::readIdentifier(string& ret){
	BlockIterator tmp = itor;
	int ch = next(tmp);
	while(ch != EOF && (ch == '_' || ::isalnum(ch))){
		itor = tmp;
		ret.push_back(ch);
		ch = next(tmp);
	}
}
bool Lexer::readNum(BlockIterator& itor, string& ret){
	BlockIterator tmp = itor;
	int ch = next(tmp);
	char *endParse = NULL;
	bool hasDot = false;
	while(ch != EOF){
		ret.push_back(ch);
		::strtod(ret.c_str(), &endParse);
		if(*endParse){
			ret.erase(ret.size() - 1, 1);
			break;
		}
		itor = tmp;
		if(ch == '.')
			hasDot = true;
		ch = next(tmp);
	}
	return hasDot;
}
int Lexer::peek(){
	int ch = *itor;
	if(-1 == ch){
		return EOF;
	}
	return ch;
}

int Lexer::getToken(){
	curToken.clear();
	skipSpace();
	int ch = next();
	if(ch == EOF)
		return END;
	currentBlock = itor.getBlock();
	if(ch == '_' || ::isalpha(ch)){
		string ident;
		ident.push_back(ch);
		readIdentifier(ident);
		curToken.swap(ident);
		return IDENTIFIER;
	}
	if(::isdigit(ch)){
		string num;
		char *endParse = NULL;
		BlockIterator tmp = itor;
		num.push_back(ch);
		bool hasDot = readNum(tmp, num);
		_isDoubleValue = false;
		if(hasDot || ch != '0'){
			double val = ::strtod(num.c_str(), &endParse);
			if(endParse != NULL && *endParse){
					cerr << "invalid number found at line no " << getLineNo();
					exit(__LINE__);
			}
			if(!hasDot){
				for(size_t i = 0;i < num.size();++ i)
					if(num[i] == 'e' || num[i] == 'E')
						_isDoubleValue = true;
			}else
				_isDoubleValue = true;
			doubleValue = val;
			intValue = val;
			assert(_isDoubleValue || doubleValue == intValue);
		}else{
			double val = 0;
			size_t total = num.size();
			if(total > 1){
				size_t i = 1;
				int base = 8;
				bool (*isValidDigit)(Byte) = &Utility::isOctDigit;
				int (*toDigit)(Byte) = &Utility::oct2dec;

				if(num[1] == 'x' || num[1] == 'X'){
					i = 2;
					base = 16;
					isValidDigit = &Utility::isHexDigit;
					toDigit = &Utility::hex2dec;
				}
				for(;i < total;++ i){
					if(!(isValidDigit)(num[i])){
						cerr << "invalid number found at line no " << getLineNo();
						exit(__LINE__);
					}
					val = val * base + (toDigit)(num[i]);
				}
			}
			doubleValue = val;
			intValue = val;
			assert(intValue == doubleValue);
		}
		itor = tmp;
		curToken.swap(num);
		return NUMBER;
	}
	if(ch == '"'){
		string str;
		readCharSequence('"', str);
		curToken.swap(str);
		return STRING;
	}
	if(ch == '\''){
		string str;
		readCharSequence('"', str);
		if(str.size() != 1){
			cerr << "invalid character constant found at line no " << getLineNo();
			exit(__LINE__);
		}
		curToken.swap(str);
		intValue = curToken[0];
		doubleValue = intValue;
		_isDoubleValue = false;
		return CHARACTER;
	}
	if(ch == '\\'){
		cerr << "invalid escape sequence found at line no " << getLineNo();
		exit(__LINE__);
	}
	switch(ch){
		case '{':case '}':case '[':case ']':case '(':case ')':
		case ',':case ';':case '.':case '?':case ':':
		case '~':return ch;
		case '#':case '!':case '%':case '^':case '-':case '+':case '=':
		case '<':case '>':case '|':case '/':case '&':case '*':{
			BlockIterator tmp = itor;
			int nextCh = next(tmp, ch != '/');
			if(nextCh == '='){
				// != %= ^= -= += == <= >= |= /= &= *=
				itor = tmp;
				curToken.push_back(ch);
				curToken.push_back(nextCh);
				//return OPERATOR;
				return OPERATOR;
			}
			switch(ch){
				case '/':{
					if(nextCh == '*'){
						//  /*
						itor = tmp;
						readMultiLineComment(curToken);
						return COMMENT;
					}else if(nextCh == '/'){
						// //
						itor = tmp;
						skipToEndOfLine(curToken);
						return COMMENT;
					}
				 }
				case '-':{
					if(ch == '-' && nextCh == '>'){
						// ->
						itor = tmp;
						curToken.push_back(ch);
						curToken.push_back(nextCh);
						//return OPERATOR;
						return OPERATOR;
					}
				}case '+':case '<':case '>':case '|':case '&':{
					if(nextCh == ch){
						// ++ -- << >> || &&
						itor = tmp;
						curToken.push_back(ch);
						curToken.push_back(nextCh);
						if(ch == '<' || ch == '<'){
							int thirdCh = next(tmp);
							if(thirdCh == '='){
								// <<= >>=
								itor = tmp;
								curToken.push_back(thirdCh);
							}
						}
						//return OPERATOR;
						return OPERATOR;
					}							
				}default:return ch;
			}
		}
	}
	if(ch == '\n')
		return NEWLINE;
	return BAD_TOKEN;
}


bool Modifiers::check(Modifier modifier){
	switch(modifier.category().getIndex()){
		case ModifierCategory::STORAGE:return !isStorageSet();
		case ModifierCategory::SIGN:
			if(isUnsigned())
				return modifier.getKeyword() == Keyword::unsignedKeyword();
			if(isSigned())
				return modifier.getKeyword() == Keyword::signedKeyword();
			return true;
		case ModifierCategory::LENGTH:
			if(modifier.getKeyword() == Keyword::longKeyword())
				return !isShort();
			return !isLong1() && !isLong2();
		case ModifierCategory::CV:
		case ModifierCategory::INVALID:
			return true;
	}
	return true;
}
bool Modifiers::add(Modifier modifier, bool duplicateRet){
	if(!check(modifier))
		return false;
	if(modifier.getKeyword() == Keyword::autoKeyword()){
		if(isAuto())
			return duplicateRet;
		storage = 1;
		return true;
	}
	if(modifier.getKeyword() == Keyword::staticKeyword()){
		if(isStatic())
			return duplicateRet;
		storage = 2;
		return true;
	}
	if(modifier.getKeyword() == Keyword::registerKeyword()){
		if(isRegister())
			return duplicateRet;
		storage = 3;
		return true;
	}
	if(modifier.getKeyword() == Keyword::typedefKeyword()){
		if(isTypedef())
			return duplicateRet;
		storage = 4;
		return true;
	}
	if(modifier.getKeyword() == Keyword::externKeyword()){
		if(isExtern())
			return duplicateRet;
		storage = 5;
		return true;
	}
	if(modifier.getKeyword() == Keyword::constKeyword()){
		if(isConst())
			return duplicateRet;
		constBit = 1;
		return true;
	}
	if(modifier.getKeyword() == Keyword::volatileKeyword()){
		if(isVolatile())
			return duplicateRet;
		volatileBit = 1;
		return true;
	}
	if(modifier.getKeyword() == Keyword::longKeyword()){
		if(isLong1() || isLong2())
			return duplicateRet;
		if(!isLong1())
			longBit1 = 1;
		else
			longBit2 = 1;
		return true;
	}
	
	if(modifier.getKeyword() == Keyword::shortKeyword()){
		if(isShort())
			return duplicateRet;
		shortBit = 1;
		return true;
	}
	
	if(modifier.getKeyword() == Keyword::signedKeyword()){
		if(isSigned())
			return duplicateRet;
		signedBit = 1;
		return true;
	}
	
	if(modifier.getKeyword() == Keyword::unsignedKeyword()){
		if(isUnsigned())
			return duplicateRet;
		unsignedBit = 1;
		return true;
	}
	return false;
}
bool Modifiers::check(const Type *type){
	const PrimType *prim = dynamic_cast<const PrimType*>(type);
	if(prim != NULL){
		if(prim->equals(PrimType::intType())
			|| prim->equals(PrimType::charType()))
			return true;
		return !this->isSigned() && !this->isLong1() && !this->isLong2()
			&& !this->isShort() && !this->isUnsigned();
	}
	if(dynamic_cast<const PointerType*>(type) != NULL || 
		dynamic_cast<const CompoundType*>(type) != NULL){
		return !this->isSigned() && !this->isLong1() && !this->isLong2()
			&& !this->isShort() && !this->isUnsigned();
	}
	const Typedef *defType = dynamic_cast<const Typedef*>(type);
	if(defType != NULL){
		return check(defType->getSrc());
	}
	const ModifiedType *modType = dynamic_cast<const ModifiedType*>(type);
	if(modType != NULL)
		return check(modType->getModifiers()) && check(modType->getType());
	return false;
}

bool Modifiers::check(Modifiers modifiers){
	if(modifiers.storage != this->storage && modifiers.storage * this->storage != 0)
		return false;
	if((modifiers.signedBit | this->signedBit) == 1
		&&((modifiers.unsignedBit | this->unsignedBit) == 1))
		return false;
	if((modifiers.shortBit | this->shortBit) == 1
		&& (modifiers.longBit1 | modifiers.longBit2 | this->longBit1 | this->longBit2) == 1)
		return false;
	return true;
}

Macro *GlobalScope::findMacro(const std::string& macro){
	map<string, int>::iterator pos = macros.find(macro);
	if(pos == macros.end())
		return NULL;
	return dynamic_cast<Macro*>(nodes[pos->second]);
}
std::pair<Macro*,int> GlobalScope::findMacroIndex(const std::string& macro){	
	map<string, int>::iterator pos = macros.find(macro);
	if(pos == macros.end()){
		return make_pair((Macro*)NULL, -1);
	}
	return make_pair(dynamic_cast<Macro*>(nodes[pos->second]), pos->second);
}
Type *Scope::findType(const std::string& name, unsigned int level){
	if(parent != NULL && level > 0){
		return parent->findType(name, level - 1);
	}
	for(size_t i = 0;i < nodes.size();++ i){
		Node *node = nodes[i];
		NamedCompoundType *namedType = dynamic_cast<NamedCompoundType*>(node);
		if(namedType != NULL && namedType->getName() == name)
			return namedType;
		Typedef *td = dynamic_cast<Typedef*>(node);
		if(td != NULL && td->getName() == name)
			return td;
	}
	return NULL;
}
Node *Scope::findDeclare(const std::string& name, unsigned int level){
	if(parent != NULL && level > 0){
		return parent->findDeclare(name, level - 1);
	}
	for(size_t i = 0;i < nodes.size();++ i){
		Node *node = nodes[i];
		NamedCompoundType *namedType = dynamic_cast<NamedCompoundType*>(node);
		if(namedType != NULL && namedType->getName() == name)
			return namedType;
		Declaration *dec = dynamic_cast<Declaration*>(node);
		if(dec != NULL && dec->getName() == name)
			return dec;
		Typedef *td = dynamic_cast<Typedef*>(node);
		if(td != NULL && td->getName() == name)
			return td;
	}
	return NULL;
}
Type *GlobalScope::findType(const std::string& name, unsigned int level){
	map<string, int>::iterator pos = declares.find(name);
	if(pos == declares.end())
		return NULL;
	return dynamic_cast<Type*>(nodes[pos->second]);
}
Node *GlobalScope::findDeclare(const std::string& name, unsigned int level){
	map<string, int>::iterator pos = declares.find(name);
	if(pos == declares.end())
		return NULL;
	return nodes[pos->second];
}

bool GlobalScope::removeMacro(const std::string& macro){
	map<string, int>::iterator pos = macros.find(macro);
	if(pos == macros.end())
		return false;
	delete nodes[pos->second];
	nodes[pos->second] = NULL;
	return true;
}
/****************************************************************************/
void Declaration::toString(string& str)const {
	string tmp = str;
	tmp += name;
	if(type != NULL)
		type->toString(tmp);
	str.clear();
	modifiers.toString(str);
	if(str.empty())
		str.swap(tmp);
	else{
		if(!tmp.empty() && !::isspace(tmp[0]))
			str.push_back(' ');
		str += tmp;
	}
}
void Macro::toString(std::string& str)const {
	str += "\n#define ";
	str += name;
	if(!params.empty()){
		str.push_back('(');
		for(size_t i = 0;i < params.size(); ++ i)
			str += params[i];
		str.push_back(')');
	}
	str.push_back(' ');
	str += content;
	str.push_back('\n');
}
void Comment::toString(std::string& str)const {
	str += "/*";
	str += content;
	str += "*/";
}
void Modifiers::toString(std::string& str)const{
	string& tmp = str;
	if(isAuto()){
		tmp += Keyword::autoKeyword().name();
		if(::isspace(tmp[tmp.size() - 1]))
			tmp.push_back(' ');
	}
	if(isStatic()){
		tmp += Keyword::staticKeyword().name();
		if(::isspace(tmp[tmp.size() - 1]))
			tmp.push_back(' ');
	}
	if(isRegister()){
		tmp += Keyword::registerKeyword().name();
		if(::isspace(tmp[tmp.size() - 1]))
			tmp.push_back(' ');
	}
	if(isTypedef()){
		tmp += Keyword::typedefKeyword().name();
		if(::isspace(tmp[tmp.size() - 1]))
			tmp.push_back(' ');
	}
	if(isExtern()){
		tmp += Keyword::externKeyword().name();
		if(::isspace(tmp[tmp.size() - 1]))
			tmp.push_back(' ');
	}


	if(isConst()){
		tmp += Keyword::constKeyword().name();
		if(::isspace(tmp[tmp.size() - 1]))
			tmp.push_back(' ');
	}
	if(isVolatile()){
		tmp += Keyword::volatileKeyword().name();
		if(::isspace(tmp[tmp.size() - 1]))
			tmp.push_back(' ');
	}
	if(isLong1()){
		tmp += Keyword::longKeyword().name();
		if(::isspace(tmp[tmp.size() - 1]))
			tmp.push_back(' ');
	}
	if(isLong2()){
		tmp += Keyword::longKeyword().name();
		if(::isspace(tmp[tmp.size() - 1]))
			tmp.push_back(' ');
	}
	if(isShort()){
		tmp += Keyword::shortKeyword().name();
		if(::isspace(tmp[tmp.size() - 1]))
			tmp.push_back(' ');
	}
	if(isSigned()){
		tmp += Keyword::signedKeyword().name();
		if(::isspace(tmp[tmp.size() - 1]))
			tmp.push_back(' ');
	}
	if(isUnsigned()){
		tmp += Keyword::unsignedKeyword().name();
		if(::isspace(tmp[tmp.size() - 1]))
			tmp.push_back(' ');
	}
}
void PrimType::toString(std::string& str)const {
	string tmp;
	tmp += keyword.name();
	if(!str.empty() && !::isspace(str[0]))
		tmp.push_back(' ');
	tmp += str;
	tmp.swap(str);
}
void PointerType::toString(std::string& str)const {
	string tmp;
	bool needBracket = false;
	const ModifiedType *modType = dynamic_cast<const ModifiedType*>(srcType);
	if(modType != NULL)
		needBracket = (dynamic_cast<const ArrayType*>(modType->getUnModifiedType()) != NULL);
	else if(dynamic_cast<const ArrayType*>(srcType) != NULL)
		needBracket = true;
	if(needBracket)
		tmp.push_back('(');
	tmp.push_back('*');
	tmp += str;
	if(needBracket)
		tmp.push_back(')');
	if(srcType != NULL){
		srcType->toString(tmp);
	}
	str.swap(tmp);
}
void ArrayType::toString(std::string& str)const{
	str += "[";
	str += this->arraySize;
	str += "]";
	if(srcType != NULL){
		srcType->toString(str);
	}
}
void Scope::toString(std::string& str) const{
	string tmp;
	if(att->endToken != Lexer::END)
		tmp.push_back(att->startToken);
	if(!tmp.empty() && tmp[tmp.size() - 1] != '\n' && att->needNewLine)
			tmp.push_back('\n');
	for(size_t i = 0;i < nodes.size();++ i){
		if(nodes[i] == NULL)
			continue;
		string one;
		Typedef *td = dynamic_cast<Typedef*>(nodes[i]);
		if(td != NULL)
			td->getDefinition(one);
		else{
			NamedCompoundType *namedType = dynamic_cast<NamedCompoundType*>(nodes[i]);
			if(namedType != NULL)
				namedType->getDefinition(one);
			else
				nodes[i]->toString(one);
		}
		if(dynamic_cast<Comment*>(nodes[i]) == NULL && dynamic_cast<Macro*>(nodes[i]) == NULL && 
			(i != nodes.size() - 1 || att->needDeclarationeEndToken) )
			one.push_back(att->declarationEndToken);
		tmp += one;
		if(!tmp.empty() && tmp[tmp.size() - 1] != '\n' && att->needNewLine)
			tmp.push_back('\n');
	}
	if(str.empty())
		str.swap(tmp);
	else
		str += tmp;
	if(att->endToken != Lexer::END)
		str.push_back(att->endToken);
}
void FunctionType::toString(std::string& str)const {
	string tmp;
	if(!Utility::isValidIdentifier(str.c_str())){
		tmp.push_back('(');
		tmp += str;
		tmp.push_back(')');
	}else{
		tmp += str;
	}
	if(params != NULL)
		params->toString(tmp);
	else
		tmp += "(void)";
	retType->toString(tmp);
	str.swap(tmp);
}
void Typedef::toString(std::string& str) const{
	string tmp = name;
	if(!str.empty() && !::isspace(str[0]))
		tmp.push_back(' ');
	tmp += str;
	str.swap(tmp);
}
void ModifiedType::toString(std::string& str) const{
	string tmp;
	this->modifiers.toString(tmp);
	if(!str.empty() && !::isspace(str[0]))
		tmp.push_back(' ');
	tmp += str;
	if(this->srcType != NULL)
		srcType->toString(tmp);
	str.swap(tmp);
}
void BitsType::toString(std::string& str) const{
	string tmp = "int";
	if(!str.empty() && !::isspace(str[0]))
		tmp.push_back(' ');
	tmp += str;
	tmp += bitLen;
	tmp.swap(str);
}
void AnonymousCompoundType::toString(std::string& str) const{
	string tmp("struct ");
	scope->toString(tmp);
	tmp += str;
	str.swap(tmp);
}
void NamedCompoundType::toString(std::string& str) const{
	string tmp("struct ");
	tmp += name;
	if(!str.empty() && !::isspace(str[0]))
		tmp.push_back(' ');
	tmp += str;
	str.swap(tmp);
}

void Typedef::getDefinition(std::string& str) const{
	Declaration dec(name, srcType);
	Modifiers modifiers;
	modifiers.add(Modifier(Keyword::typedefKeyword()));
	dec.setModifiers(modifiers);
	dec.toString(str);
}
void NamedCompoundType::getDefinition(std::string& str) const{
	string tmp("struct ");
	tmp += name;
	scope->toString(tmp);
	str += tmp;
}
void EnumType::getDefinition(std::string& str) const{
	string tmp("enum ");
	str += "enum ";
	tmp += name;
	tmp.push_back('{');
	str += definition;
	str.push_back('}');
}
void EnumType::toString(std::string& str)const{
	if(name.empty()){
		string tmp("enum {");
		tmp += definition;
		tmp += "}";
		tmp += str;
		str.swap(tmp);
	}else{
		string tmp("enum ");
		tmp += name;
		if(!str.empty() && !::isspace(str[0]))
			tmp.push_back(' ');
		tmp += str;
		str.swap(tmp);
	}
}
/****************************************************************************/
HeaderParser::HeaderParser(const Byte *start, const Byte *end) : 
	block(start, end), currentScope(NULL), globalScope(NULL), lexer(NULL), isParsed(false){
}
int HeaderParser::getToken(int token){
	while(token == Lexer::NEWLINE){
		token = lexer->getToken();
		if(token == '#'){
			token = preprocessInstruction();
		}
	}
	if(token == Lexer::IDENTIFIER && !lexer->getCurrentBlock()->isMacroExpansion() 
		&& expandMacro(globalScope->findMacro(lexer->getTokenContent())))
		return getToken();
	if(token == Lexer::COMMENT){
		Comment *comment = new Comment();
		comment->getContent().swap(lexer->getTokenContent());
		currentScope->add(comment);
		return getToken();
	}
	return token;
}
int HeaderParser::getToken(){
	return getToken(lexer->getToken());
}

bool CompoundTypeScope::add(Node *node){
	NamedCompoundType *namedType = dynamic_cast<NamedCompoundType*>(node);
	if(namedType != NULL){
		map<string, int>::iterator pos = declares.find(namedType->getName());
		if(pos == declares.end()){
			declares[namedType->getName()] = nodes.size();
			nodes.push_back(namedType);
			return true;
		}
	}
	Declaration *dec = dynamic_cast<Declaration*>(node);
	if(dec != NULL){
		map<string, int>::iterator pos = declares.find(dec->getName());
		if(pos == declares.end()){
			declares[dec->getName()] = nodes.size();
			nodes.push_back(dec);
			return true;
		}
	}
	AnonymousCompoundType *anonymousType = dynamic_cast<AnonymousCompoundType*>(node);
	if(anonymousType != NULL){
		nodes.push_back(node);
		return true;
	}
	Comment *comment = dynamic_cast<Comment*>(node);
	if(comment != NULL){
		nodes.push_back(node);
		return true;
	}
	assert(false);
	return false;
}
Type *CompoundTypeScope::findType(const std::string& name, unsigned int level){
	map<string, int>::iterator pos = declares.find(name);
	if(pos == declares.end())
		return level == 0 ? NULL : Scope::findType(name, level - 1);
	return dynamic_cast<Type*>(nodes[pos->second]);
}
Node *CompoundTypeScope::findDeclare(const std::string& name, unsigned int level){
	map<string, int>::iterator pos = declares.find(name);
	if(pos == declares.end())
		return level == 0 ? NULL : Scope::findDeclare(name, level - 1);
	return nodes[pos->second];
}
bool GlobalScope::add(Node *node){
	if(nodes.size() == 6227)
		cout << "this is a test" << endl;
	Macro *macro = dynamic_cast<Macro*>(node);
	if(macro != NULL)
	{
		map<string, int>::iterator pos = macros.find(macro->name);
		if(pos == macros.end()){
			macros[macro->name] = nodes.size();
			nodes.push_back(macro);
			return true;
		}
	}
	NamedCompoundType *namedType = dynamic_cast<NamedCompoundType*>(node);
	if(namedType != NULL){
		map<string, int>::iterator pos = declares.find(namedType->getName());
		if(pos == declares.end()){
			declares[namedType->getName()] = nodes.size();
			nodes.push_back(namedType);
			return true;
		}
	}
	AnonymousCompoundType *anonymousType = dynamic_cast<AnonymousCompoundType*>(node);
	if(anonymousType != NULL){
		nodes.push_back(anonymousType);
		return true;
	}
	Declaration *dec = dynamic_cast<Declaration*>(node);
	if(dec != NULL){
		map<string, int>::iterator pos = declares.find(dec->getName());
		if(pos == declares.end()){
			declares[dec->getName()] = nodes.size();
			nodes.push_back(dec);
			return true;
		}
	}
	Typedef *td = dynamic_cast<Typedef*>(node);
	if(td != NULL){
		map<string, int>::iterator pos = declares.find(td->getName());
		if(pos == declares.end()){
			declares[td->getName()] = nodes.size();
			nodes.push_back(td);
			return true;
		}
	}
	Comment *comment = dynamic_cast<Comment*>(node);
	if(comment != NULL){
		nodes.push_back(comment);
		return true;
	}
  EnumType *enumType = dynamic_cast<EnumType*>(node);
  if (enumType != NULL){
    nodes.push_back(enumType);
    return true;
  }
	assert(false);
	return false;
}

void HeaderParser::reportError(const char *msg, int exitCode){
	Position pos = lexer->getPosition();
	cerr << msg << " at line " << pos.row << ", column " << pos.col << endl;
	cin.get();
	if(exitCode != 0)
		exit(exitCode);
}

bool HeaderParser::expandMacro(Macro *macro){
	if(macro != NULL){
		if(macro->content.empty())
			return true;
		MacroExpansion macroExpansion(macro, this->globalScope, this->lexer);
		string result;
		if(!macroExpansion.expandMacro(result)){
			this->reportError("invalid macro invokation");
			return false;
		}
		Utility::trim(result);
		if(!macroExpansion.changed() && !macro->content.empty()){
			blocks.push_back(Block((const Byte*)macro->content.c_str(), (const Byte*)macro->content.c_str() + macro->content.size(), 
				(int)Block::NOT_UPDATE_POS | (int)Block::IS_MACRO_EXPANSION));
		}
		else if(!result.empty()){
			expandedMacros.push_back(std::string());
			string& str = expandedMacros.back();
			str.swap(result);
			blocks.push_back(Block((const Byte*)str.c_str(), (const Byte*)str.c_str() + str.size(), 
				(int)Block::NOT_UPDATE_POS | (int)Block::IS_MACRO_EXPANSION));
		}
		blocks.back().setStartPosition(lexer->getPosition());
		lexer->iterator().split(&blocks.back());
			return true;
	}
	return false;
}
int HeaderParser::parseMacro(){
		int token = lexer->getToken();
		if(token != Lexer::IDENTIFIER){
				reportError("macro definition expects an identifier", __LINE__);
		}
		if(globalScope->findMacro(lexer->getTokenContent()) != NULL){
				reportError("macro redefined", __LINE__);
		}
		Macro *p = new Macro;
		int ch = lexer->peek();
		p->name.swap(lexer->getTokenContent());
		if(::isspace(ch)){
				lexer->skipToEndOfLine(p->content);
				globalScope->add(p);
				return lexer->getToken();
		}
		token = lexer->getToken();
		if(token == '('){
				token = lexer->getToken();
				while(token != ')' && token != Lexer::NEWLINE){
					if(token != Lexer::IDENTIFIER){
						reportError("identifier expected, invalid macro definition", __LINE__);
					}
					p->params.push_back(string());
					p->params.back().swap(lexer->getTokenContent());
					if(p->params.end() != ::find(p->params.begin(), p->params.end(), lexer->getTokenContent())){
						reportError("duplicate parameter name  in macro", __LINE__);
					}
					token = lexer->getToken();
					if(token == ','){
						token = lexer->getToken();
						if(token != Lexer::IDENTIFIER){
							reportError("invalid macro definition", __LINE__);
						}
					}
				}
				if(token == Lexer::NEWLINE){
					reportError("')' expected, invalid macro definition", __LINE__);
				}
				lexer->skipToEndOfLine(p->content);
		}else{
				reportError("unexpected macro definition", __LINE__);
		}
		globalScope->add(p);
		return getToken();
}
void HeaderParser::includeFile(const char *fileName){
	static string baseDir = "D:\\";
	map<string, string>::iterator pos = fileContents.find(fileName);
	const char *start = NULL;
	const char *end = NULL;
	if(pos != fileContents.end()){
		const string& content = pos->second;
		start = content.c_str();
		end = start + content.size();
	}else{
		string path = baseDir;
		string& content = fileContents[fileName];
		while(*fileName && ::isspace(*fileName))
			++fileName;
		++fileName;
		path += fileName;
		path.erase(path.size() - 1, 1);
		if(!Utility::fileToString(path.c_str(), content)){
			cerr << "fail to open file " << path << endl;
			return;
		}
		start = content.c_str();
		end = start + content.size();
	}
	if(start != NULL && end != NULL){
		blocks.push_back(Block((const Byte*)start, (const Byte*)end));
		lexer->iterator().split(&blocks.back());
	}
}
int HeaderParser::skipBraces(char ch, std::string* s){
	string str;
	str.push_back(ch);
	int token = getToken();
	while(token != Lexer::END){
		switch(token){
			case '(':case '[':case '{':
				str.push_back(token);break;
			case ')':case ']':case '}':{
				char expected = Utility::getMatchedBrace(str[str.size() - 1]);
				if(expected != token){
					reportError("unmatched right brace found");
					exit(__LINE__);
				}
				str.erase(str.size() - 1, 1);
				break;
			}
		}
		if(str.empty())
			break;
		if(s != NULL)
			lexer->toString(token, *s);
		token = getToken();
	}
	if(token == Lexer::END){
		reportError("right brace expected", __LINE__);
	}
	return getToken();
}

class HeaderParser::DeclarationParser{
	ModifiedType *type;
	Modifiers storage;
	HeaderParser *headerParser;
	const Scope::ScopeAttribute *att;
	int currentToken;
	bool hasDefinition;
	void getModifiers(Modifiers& modifiers, const Type *type){
		while(currentToken == Lexer::IDENTIFIER){
			Keyword keyword(headerParser->lexer->getTokenContent().c_str());
			if(!keyword.modiferCatetory().isValid())
				break;
			Modifier modifier(keyword);
			//check modifier compatibility
			if(!modifiers.add(modifier)){
				headerParser->reportError("incompatible modifier", __LINE__);
			}
			currentToken = headerParser->getToken();
		}
		
		//check modifier compatibility with type
		if(type != NULL && !modifiers.check(type)){
			headerParser->reportError("modifiers incompatible with type", __LINE__);
		}

		//check modifier compatibility with scope attribute
	}
	ParameterList *parseParamList(){
		//parse function paramter list
		ParameterList *paramList = new ParameterList(headerParser->currentScope, &paramListAttribute);								
		headerParser->currentScope = paramList;
		this->att = headerParser->currentScope->getAttribute();		
		currentToken = headerParser->getToken();
		while(currentToken == Lexer::IDENTIFIER){
			currentToken = headerParser->parseDeclaration();
			if(currentToken == att->endToken)
				break;
			if(currentToken != att->declarationEndToken)
				headerParser->reportError("',' expected", __LINE__);
			currentToken = headerParser->getToken();
		}
		if(currentToken != att->endToken){
			headerParser->reportError("')' expected", __LINE__);
		}
		headerParser->currentScope = paramList->getParent();
		this->att = headerParser->currentScope->getAttribute();
		currentToken = headerParser->getToken();
		return paramList;
	}
	//type : modified type, function type, pointer type, bit type, array type.
	bool setNodeType(Type *type, Type *node){
		ModifiedType *modType = dynamic_cast<ModifiedType*>(type);
		if(modType != NULL)
			return setNodeType(const_cast<Type*>(modType->getType()), node);
		modType = dynamic_cast<ModifiedType*>(node);
		if(modType != NULL)
			return setNodeType(type, const_cast<Type*>(modType->getType()));
		FunctionType *funcType = dynamic_cast<FunctionType*>(type);
		if(funcType != NULL){
			assert(funcType->getRetType() == NULL);
			funcType->setRetType(node);
			return true;
		}
		PointerType *pointerType = dynamic_cast<PointerType*>(type);
		if(pointerType != NULL){
			pointerType->setSrcType(node);
			return true;
		}
		BitsType *bitsType = dynamic_cast<BitsType*>(type);
		if(bitsType != NULL){
			const PrimType *primType = dynamic_cast<const PrimType*>(node);	
			Keyword keyword;
			if(primType != NULL)
				keyword = primType->getKeyword();
			if(keyword != Keyword::intKeyword() &&
				keyword != Keyword::charKeyword()){
				headerParser->reportError("integeral type expected!");
				return false;
			}
			return true;
		}
		return false;
	}
	//innerMost, possibilities : (variable or function) declaration, function type, pointer type, bit type, array type.
	//ouuterMost, possibilities : (variable or function) declaration, function type, pointer type, bit type, array type.
	void parseInBracket(Node *&innerMost, Node *&outterMost){
		/****************start************************/		
		//pointer type
		Declaration *declaration = NULL;
		ModifiedType *declaredType = NULL;
		currentToken = headerParser->getToken();
		while(currentToken == '*'){
			Modifiers modifiers;
			PointerType *pType = new PointerType(declaredType);
			currentToken = headerParser->getToken();
			getModifiers(modifiers, pType);
			declaredType = new ModifiedType(pType);
			declaredType->setModifiers(modifiers);
			if(outterMost == NULL)
				outterMost = declaredType;
		}
		if(!isEnd()){
			Node *outterNode = NULL, *innerNode = NULL;
			if(currentToken == '('){
				parseInBracket(innerNode, outterNode);
				if(outterNode == NULL || innerNode == NULL)
					headerParser->reportError("syntax error", __LINE__);
			}else if(currentToken == Lexer::IDENTIFIER){
				const string& name = headerParser->lexer->getTokenContent();
				Keyword keyword(name.c_str());
				if(keyword.isValid()){
					headerParser->reportError("unexpected keyword", __LINE__);
				}
				declaration = new VariableDeclaration(name, NULL);
				currentToken = headerParser->getToken();
			}
			if(currentToken == '['){
				//declaration of arrays
				ArrayType *arType = new ArrayType(declaredType);
				declaredType = new ModifiedType(arType);
				//Bypass the expression that represent the size of the array.
				while((currentToken = headerParser->skipBraces('[', &arType->arraySize)) == '['){
					ArrayType *tmp = new ArrayType(	arType->getSrcType());
					arType->setSrcType(tmp);
					arType = tmp;
				}
				if(outterMost == NULL)
					outterMost = arType;
			}else if( currentToken == '('){
				//function parameter list
				if(declaration == NULL && innerNode == NULL && outterNode == NULL){
					headerParser->reportError("function declaration expects a function name", __LINE__);							
				}
				if(hasDefinition)
					headerParser->reportError("type definition is not allowed", __LINE__);
				ParameterList *paramList = parseParamList();
				FunctionType *funcType = new FunctionType(declaredType, paramList);
				declaredType = new ModifiedType(funcType);
				if(outterMost == NULL){
					assert(funcType->getRetType() == NULL);
					outterMost = declaredType;
				}
			}else if(currentToken == ':'){
				//bit type
				if(!att->allowBitType){
					headerParser->reportError("bit type not allowed in this scope!", __LINE__);
				}
				if(declaredType != NULL)
					headerParser->reportError("integeral type expected!", __LINE__);
				assert(outterMost == NULL);
				BitsType *bitsType = new BitsType();
				declaredType = new ModifiedType(bitsType);
				outterMost = declaredType;
				string str;
				bypassExpr(str);
				bitsType->bitLen.swap(str);
			}

			if(declaration != NULL){
				assert(innerNode == NULL && outterNode == NULL);
				assert(declaration->getType() == NULL && declaredType != NULL);

				const FunctionType *funcType = dynamic_cast<const FunctionType*>(declaredType->getUnModifiedType());
				if(funcType != NULL)
					declaration = new FunctionDeclaration(declaration->getName(), funcType);
				else
					declaration->setType(declaredType);
			}else{
				assert(innerNode != NULL && outterNode != NULL);
				declaration = dynamic_cast<Declaration*>(outterNode);
				if(declaration != NULL){
					assert(innerNode == outterNode);
					assert(declaration->getType() == NULL && declaredType != NULL);

					const FunctionType *funcType = dynamic_cast<const FunctionType*>(declaredType->getUnModifiedType());
					if(funcType != NULL)
						declaration = new FunctionDeclaration(declaration->getName(), funcType);
					else
						declaration->setType(declaredType);
				}else{
					Type *tmpType = dynamic_cast<Type*>(outterNode);
					if(tmpType == NULL || !setNodeType(tmpType, declaredType))
						headerParser->reportError("syntax error", __LINE__);
					
					declaration = dynamic_cast<Declaration*>(innerNode);
					if(declaration == NULL){
						tmpType = dynamic_cast<Type*>(innerNode);
						if(tmpType == NULL)
							headerParser->reportError("syntax error", __LINE__);
						declaredType = dynamic_cast<ModifiedType*>(tmpType);
						if(declaredType == NULL)
							declaredType = new ModifiedType(tmpType);
					}else{
						Type *tmpType = const_cast<Type*>(declaration->getType());
						if(tmpType == NULL)
							headerParser->reportError("syntax error", __LINE__);
						declaredType = dynamic_cast<ModifiedType*>(tmpType);
						if(declaredType == NULL)
							declaredType = new ModifiedType(tmpType);
					}
				}
			}
			if(currentToken == '=')
				headerParser->reportError("initialization not allowned", __LINE__);			
		}
		if(declaration != NULL)
			innerMost = declaration;
		else{
			assert(declaredType != NULL);
			innerMost = declaredType;
		}
		if(outterMost == NULL)
			outterMost = innerMost;
		if(currentToken != ')' || (declaredType == NULL &&  declaration == NULL))
			headerParser->reportError("syntax error", __LINE__);
		currentToken = headerParser->getToken();
		/****************end************************/
	}
	bool isEnd(){
		return currentToken == ';' || currentToken == att->declarationEndToken
			|| currentToken == att->endToken || currentToken == Lexer::END;
	}
	bool getType(){
		CompoundType *forwardDeclaration = NULL;
		Modifiers modifiers;
		{
			getModifiers(modifiers, NULL);
			if(currentToken == Lexer::IDENTIFIER){
				Keyword keyword(headerParser->lexer->getTokenContent().c_str());
				switch(keyword.category().getIndex()){
						case KeywordCategory::OPTER:
						case KeywordCategory::CONTROL:
							headerParser->reportError("unexpected keyword", __LINE__); break;
						case KeywordCategory::METATYPE:{
							//compound type
							CompoundType *compoundType = NULL;							
							currentToken = headerParser->getToken();
							if(currentToken == Lexer::IDENTIFIER){
								Keyword name(headerParser->lexer->getTokenContent().c_str());
								Type *tmpType = NULL;
								if(name.isValid()){
									headerParser->reportError("unexpected keyword", __LINE__);break;
								}								
								//find the compound type by name;
								tmpType = headerParser->currentScope->findType(headerParser->lexer->getTokenContent(), -1);
								compoundType = dynamic_cast<CompoundType*>(tmpType);
								if(tmpType != NULL && compoundType == NULL){
									//if this name is defined but is not a compound type
									headerParser->reportError("conflict definition of identifier", __LINE__);break;
								}
								if(compoundType == NULL){
									//if this compound type is not defined.
									if(!(keyword == Keyword::enumKeyword()))
										compoundType = new NamedCompoundType(headerParser->lexer->getTokenContent());
									else
										compoundType = new EnumType(headerParser->lexer->getTokenContent());
									forwardDeclaration = compoundType;
								}
								currentToken = headerParser->getToken();
							}
							if(currentToken == '{'){								
								if(compoundType == NULL){
									//anonymouse compound type
									if(!(keyword == Keyword::enumKeyword()))
										compoundType = new AnonymousCompoundType();
									else
										compoundType = new EnumType();
								}
								if(compoundType->isDefined()){
									headerParser->reportError("redefinition", __LINE__);break;
								}
								if(forwardDeclaration != NULL)
									headerParser->currentScope->add(compoundType);
								if(!(keyword == Keyword::enumKeyword())){
									//parse definition of compound type
									CompoundTypeScope *compoundTypeScope = new CompoundTypeScope(headerParser->currentScope, &compoundTypeScopeAttribute);								
									headerParser->currentScope = compoundTypeScope;
									this->att = headerParser->currentScope->getAttribute();
									compoundType->setScope(compoundTypeScope);
									while((currentToken = headerParser->getToken()) != att->endToken){
										if(currentToken == att->declarationEndToken)
											continue;
										currentToken = headerParser->parseDeclaration();
										if(currentToken != att->declarationEndToken){
											headerParser->reportError("';' expected", __LINE__);break;
										}
									}
									headerParser->currentScope = compoundTypeScope->getParent();
									this->att = headerParser->currentScope->getAttribute();
									currentToken = headerParser->getToken();
								}else{
									//headerParser->reportError("could not handle enum definitions", __LINE__);break;
									EnumType *enumType = (EnumType*)compoundType;
									currentToken = headerParser->skipBraces('{', &enumType->definition);
								}
								hasDefinition = true;
								forwardDeclaration = NULL;
							}
							if(compoundType == NULL){
								headerParser->reportError("missing name", __LINE__);break;
							}
							this->type = new ModifiedType(compoundType);
							break;
						}
						case KeywordCategory::INVALID:{
							//compound type or typedef
							//find type by name
							Type *tmpType = headerParser->currentScope->findType(headerParser->lexer->getTokenContent(), -1);
							if(tmpType != NULL){
								//identifier is a type
								this->type = new ModifiedType(tmpType);
								currentToken = headerParser->getToken();
							}
							else if(modifiers.isLong1() || modifiers.isLong2() || modifiers.isShort()
								|| modifiers.isSigned() || modifiers.isUnsigned()){
								//identifier is a variable name
								this->type = new ModifiedType(PrimType::intType());
							}
							break;
						}
						case KeywordCategory::PRIM:
							//primary type
							this->type = new ModifiedType(PrimType::get(keyword));
							assert(this->type->getType() != NULL);
							currentToken = headerParser->getToken();
							break;
				}
			}else{
				//empty type specification, check modifiers,
				//if there are "long", "short", "signed", or "unsigned", then take it as integer type.
				if(modifiers.isLong1() || modifiers.isLong2() || modifiers.isShort()
					|| modifiers.isSigned() || modifiers.isUnsigned())
					this->type = new ModifiedType(PrimType::intType());
			}
			if(type == NULL){
				headerParser->reportError("type expected", __LINE__);
			}
			getModifiers(modifiers, type->getType());
		}
		storage = modifiers;
		modifiers.clearStorage();
		storage.clearCV();
		storage.clearSign();
		storage.clearLength();
		if(type != NULL){
			if(isEnd()){
				//check if there is a forward declaration, a forward declaration could only appear in global scope;
				if(forwardDeclaration != NULL){
					//forward declaration.
					headerParser->currentScope->add(forwardDeclaration);
				}
				return false;				
			}else{				
				if(forwardDeclaration != NULL){
					//not a forward declaration
					if(!forwardDeclaration->isDefined())
						headerParser->reportError("using undefined type", __LINE__);
					else
						headerParser->currentScope->add(forwardDeclaration);
				}
				this->type->setModifiers(modifiers);
			}
		}
		return true;
	}
public:
	DeclarationParser(HeaderParser *headerParser, int currentToken)
		: headerParser(headerParser), currentToken(currentToken), att(NULL), hasDefinition(false){		
	}
	void bypassExpr(string& str){
		string s;
		Lexer *lexer = headerParser->lexer;
		s.push_back(0);
		while(currentToken != att->varListSeparator && !isEnd()){
			switch(currentToken){
				case '(':case '[':case '{':
					s.push_back(currentToken);break;
				case ')':case ']':case '}':{
					char expected = Utility::getMatchedBrace(s[s.size() - 1]);
					if(expected != currentToken){
						//headerParser->reportError("unmatched right brace found");
						//exit(__LINE__);
						return;
					}
					str.erase(str.size() - 1, 1);
					break;
				}
			}
			lexer->toString(currentToken, str);
			currentToken = headerParser->getToken();
		}
		if(currentToken == Lexer::END){
			headerParser->reportError("syntax error", __LINE__);
		}
	}
	int parse(){
		type = NULL;
		att = headerParser->currentScope->getAttribute();
		if(getType()){
			 parseDeclarator();
		}else{
			//if(!hasDefinition)
			//	headerParser->currentScope->add(type);
		}
		return currentToken;
	}
	
	void parseDeclarator(){
		while(true){
			//pointer type
			Declaration *declaration = NULL;
			ModifiedType *declaredType = type;
			while(currentToken == '*'){
				Modifiers modifiers;
				PointerType *pType = new PointerType(declaredType);
				currentToken = headerParser->getToken();
				getModifiers(modifiers, pType);
				declaredType = new ModifiedType(pType);
				declaredType->setModifiers(modifiers);
			}
			if(!isEnd()){
				Node *outterNode = NULL, *innerNode = NULL;
				if(currentToken == '('){
					this->parseInBracket(innerNode, outterNode);
					if(outterNode == NULL || innerNode == NULL)
						headerParser->reportError("syntax error", __LINE__);						
				}else if(currentToken == Lexer::IDENTIFIER){
					const string& name = headerParser->lexer->getTokenContent();
					Keyword keyword(name.c_str());
					if(keyword.isValid()){
						headerParser->reportError("unexpected keyword", __LINE__);
					}
					declaration = new VariableDeclaration(name, NULL);
					currentToken = headerParser->getToken();
				}
				
				
				if(currentToken == '['){
					//declaration of arrays					
					ArrayType *arType = new ArrayType(declaredType);
					declaredType = new ModifiedType(arType);
					//Bypass the expression that represent the size of the array.
					while((currentToken = headerParser->skipBraces('[', &arType->arraySize)) == '['){
						ArrayType *tmp = new ArrayType(arType->getSrcType());
						arType = new ArrayType(tmp);
						arType = tmp;
					}
				}else if( currentToken == '('){
					//function parameter list
					if(declaration == NULL && declaration == NULL && innerNode == NULL && outterNode == NULL){
						//function declration must have a function name
						//That is the following code is invalid
						// void *(int, int);
						//should be:
						//void *func(int, int);
						headerParser->reportError("function declaration expects a function name", __LINE__);							
					}
					if(hasDefinition)
						headerParser->reportError("type definition is not allowed", __LINE__);

					//parse function paramter list
					ParameterList *paramList = parseParamList();

					declaredType = new ModifiedType(new FunctionType(declaredType, paramList));
				}else if(currentToken == ':'){
					//bit type
					if(!att->allowBitType){
						headerParser->reportError("bit type not allowed in this scope!", __LINE__);
					}
					const PrimType *primType = dynamic_cast<const PrimType*>(declaredType->getType());	
					Keyword keyword;
					if(primType != NULL)
						keyword = primType->getKeyword();
					if(keyword != Keyword::intKeyword() &&
						keyword != Keyword::charKeyword())
						headerParser->reportError("integeral type expected!", __LINE__);
					BitsType *bitsType = new BitsType();
					string str;
					declaredType = new ModifiedType(bitsType);
					bypassExpr(str);
					bitsType->bitLen.swap(str);
				}
				if(declaration != NULL){
					assert(innerNode == NULL && outterNode == NULL);
					assert(declaration->getType() == NULL && declaredType != NULL);

					const FunctionType *funcType = dynamic_cast<const FunctionType*>(declaredType->getUnModifiedType());
					if(funcType != NULL)
						declaration = new FunctionDeclaration(declaration->getName(), funcType);
					else
						declaration->setType(declaredType);
				}else if(innerNode != NULL && outterNode != NULL){
					if(innerNode == NULL || outterNode == NULL)
						headerParser->reportError("syntax error", __LINE__);

					//assert(innerNode != NULL && outterNode != NULL);
					declaration = dynamic_cast<Declaration*>(outterNode);
					if(declaration != NULL){
						assert(innerNode == outterNode);
						assert(declaration->getType() == NULL && declaredType != NULL);

						const FunctionType *funcType = dynamic_cast<const FunctionType*>(declaredType->getUnModifiedType());
						if(funcType != NULL)
							declaration = new FunctionDeclaration(declaration->getName(), funcType);
						else
							declaration->setType(declaredType);
					}else{
						Type *tmpType = dynamic_cast<Type*>(outterNode);
						if(tmpType == NULL || !setNodeType(tmpType, declaredType))
							headerParser->reportError("syntax error", __LINE__);
						
						declaration = dynamic_cast<Declaration*>(innerNode);
						if(declaration == NULL){
							tmpType = dynamic_cast<Type*>(innerNode);
							if(tmpType == NULL)
								headerParser->reportError("syntax error", __LINE__);
							declaredType = dynamic_cast<ModifiedType*>(tmpType);
							if(declaredType == NULL)
								declaredType = new ModifiedType(tmpType);
						}else{
							Type *tmpType = const_cast<Type*>(declaration->getType());
							if(tmpType == NULL)
								headerParser->reportError("syntax error", __LINE__);
							declaredType = dynamic_cast<ModifiedType*>(tmpType);
							if(declaredType == NULL)
								declaredType = new ModifiedType(tmpType);
						}
					}
				}
				if(currentToken == '='){
					//assignment
					if(storage.isTypedef()){
						//initialization is not allowed in "typedef" qualified declration.
						//That is the following code is invalid:
						//typedef int a = 324;
						headerParser->reportError("initialization is not allowed in 'typedef' qualified declaration", __LINE__);
					}
					if(dynamic_cast<FunctionDeclaration*>(declaration)){
						headerParser->reportError("initialization is not allowed for function declaration", __LINE__);
					}
					if(!att->allowInitialization){
						//initialization not allowed in some scope
						//That is the following code is invalid under C
						//strut Test{
						// int a = 1343;
						//};
						//void func(int, int b = 3);
						headerParser->reportError("initialization is not allowed in this scope", __LINE__);
					}
					//Bypass the initialization expression.
					currentToken = headerParser->getToken();
					string str;
					bypassExpr(str);
				}

				if(declaration != NULL){
					if(headerParser->currentScope->findDeclare(declaration->getName(), 0) != NULL){
						//This is differenct from C, in C, redeclaration of an identifier is allowed
						//But redefinition is not allowed.
						//Notice that definition is different from declaration in C.
						//Definition will leads to momery allocation for variables or functions
						//Declaration is only specifications of functions or variables
						//But in my trival code, I treat them as the same.
						headerParser->reportError("redeclaration of an identifier", __LINE__);
					}
					if(storage.isTypedef()){
						Typedef *td = new Typedef();
						td->getName().swap(declaration->getName());
						td->setSrc(declaration->getType());
						headerParser->currentScope->add(td);
					}else{
						declaration->setModifiers(storage);
						headerParser->currentScope->add(declaration);
					}
				}
			}
			if(!att->allowFunctionDefine && dynamic_cast<const FunctionType*>(declaredType->getType()) != NULL){
					headerParser->reportError("function declaration not allowed in this scope", __LINE__);
			}
			if(currentToken == att->declarationEndToken || (!att->needDeclarationeEndToken && currentToken == att->endToken)){
				if(declaration == NULL){
					if(storage.isStorageSet()){
						headerParser->reportError("variable name expected", __LINE__);
					}
					headerParser->currentScope->add(declaredType);
				}
				break;
			}
			if(currentToken == att->varListSeparator){
				if(!att->allowVarList){
					// in function parameter list scope multiple definition of variables in one statement is not allowed.
					//That is the following code is invalid:
					//void func(int a, b,c);
					//But the following code is valid
					//int a, b, c;
					//struct Test{
					// int a, c;
					//};
					headerParser->reportError("multiple declarations is not allowed", __LINE__);
				}
				currentToken = headerParser->getToken();
				if(currentToken == att->declarationEndToken){
					headerParser->reportError("syntax error", __LINE__);
				}
				if(declaration == NULL){
					headerParser->reportError("syntax error", __LINE__);
				}
				continue;
			}
			headerParser->reportError("syntax error", __LINE__);
		}
	}
};
int HeaderParser::parseDeclaration(){
	DeclarationParser declarationParser(this, Lexer::IDENTIFIER);
	return declarationParser.parse();
}

struct PreprocessingExpr{
	GlobalScope *globalScope;
	string content;
	int value;
	struct Evaluator : public ExprVisitor{
		bool isValid;
		int value;
		Evaluator() : isValid(false), value(0){
		}
		virtual bool visit(NumberExpr *expr){
			isValid = !expr->isDouble;
			if(isValid)
				value = expr->number.intVal;
			return true;
		}
		virtual bool visit(FunctionExpr *expr){
			isValid = !Keyword(expr->name.c_str()).isValid();
			return true;
		}
		virtual bool visit(IdentifierExpr *expr){
			isValid = !Keyword(expr->name.c_str()).isValid();
			return true;
		}
		virtual bool visit(BinaryExpr *expr){
			Evaluator left;
			left.ExprVisitor::visit(expr->left);
			isValid = left.isValid;
			if(!isValid)
				return true;
			if(expr->optr == '?'){
				BinaryExpr *tmp = dynamic_cast<BinaryExpr *>(expr->right);
				Evaluator val;
				if(tmp == NULL || tmp->optr != ':'){
					isValid = false;
					return true;
				}
				if(left.value != 0)
					val.ExprVisitor::visit(tmp->left);
				else
					val.ExprVisitor::visit(tmp->right);
				isValid = val.isValid;
				value = val.value;
				return true;
			}
			Evaluator right;
			right.ExprVisitor::visit(expr->right);
			isValid = isValid && right.isValid;
			if(isValid)
				switch(expr->optr){
					case Expr::PLUS:value = left.value + right.value;break;
					case Expr::MINUS:value = left.value - right.value;break;
					case Expr::MULTI:value = left.value * right.value;break;
					case Expr::DIV:
						if(right.value == 0){
							cerr << "divided by zero" << endl;
							isValid = false;
							return true;
						}
						value = left.value / right.value;break;
					case Expr::MODULA:
						if(right.value == 0){
							cerr << "modular by zero" << endl;
							isValid = false;
							return true;
						}
						value = left.value % right.value;break;
					case Expr::LSHIFT:value = (left.value << right.value);break;
					case Expr::RSHIFT:value = (left.value >> right.value);break;

					case Expr::LESS:value = (left.value < right.value);break;
					case Expr::GREATER:value = (left.value > right.value);break;
					case Expr::LESSEQ:value = (left.value <= right.value);break;
					case Expr::GREATEREQ:value = (left.value >= right.value);break;
					case Expr::EQUAL:value = (left.value == right.value);break;
					case Expr::NEQUAL:value = (left.value != right.value);break;

					case Expr::AND:value = (left.value && right.value);break;
					case Expr::OR:value = (left.value || right.value);break;
					case Expr::BITAND:value = (left.value & right.value);break;
					case Expr::BITOR:value = (left.value | right.value);break;
					case Expr::BITXOR:value = (left.value ^ right.value);break;
					default:isValid = false;
				}
			return true;
		}
		virtual bool visit(UnaryExpr *expr){
			isValid = true;
			switch(expr->optr){
			case '+':ExprVisitor::visit(expr->expr);break;
			case '-':ExprVisitor::visit(expr->expr);value = -value;break;
			case '~':ExprVisitor::visit(expr->expr);value = ~value;break;
			case '!':ExprVisitor::visit(expr->expr);value = !value;break;
			default:isValid = false;
			}
			return true;
		}
	};
	PreprocessingExpr(GlobalScope *globalScope) : globalScope(globalScope), value(0){
	}
	bool evaluate(){
		string expaned;
		Block block((Byte*)&content[0], (Byte*)&content[0] + content.size());
		Lexer lexer(&block);
		int token = Lexer::END;
		while((token = lexer.getToken()) != Lexer::END){
			if(token == Lexer::IDENTIFIER){
				if(lexer.getTokenContent() == "defined"){
					bool hasBraces = false;
					bool value = false;
					token = lexer.getToken();
					if(token == '('){
						hasBraces = true;
						token = lexer.getToken();
					}
					value = globalScope->findMacro(lexer.getTokenContent()) != NULL;
					if(hasBraces && (token = lexer.getToken()) != ')'){
						cerr << "')' expected" << endl;
						return false;
					}
					if(value)
						expaned += " 1 ";
					else
						expaned += " 0 ";
				}else{
					Macro *macro = globalScope->findMacro(lexer.getTokenContent());
					if(macro != NULL){
						MacroExpansion macroExpansion(macro, this->globalScope, &lexer);
						string result;
						if(!macroExpansion.expandMacro(result)){
							cerr << "invalid macro invokation" << endl;
							return false;
						}
						expaned += result;
					}else
						lexer.toString(token, expaned);
				}
			}else
				lexer.toString(token, expaned);
		}
		if(expaned.empty()){
			return false;
		}
		ExprParser exprParser((Byte*)&expaned[0], (Byte*)&expaned[0] + expaned.size());
		Expr *expr = exprParser.parse();
		Evaluator eval;
		eval.ExprVisitor::visit(expr);
		this->value = eval.value;
		return eval.isValid;
	}
};

int HeaderParser::preprocessInstruction(){
	int token = lexer->getToken();
	if(token == Lexer::NEWLINE)
		return token;
	if(token == Lexer::IDENTIFIER){
		const string& name = lexer->getTokenContent();
		if(name == "define")
			return parseMacro();
		if(name == "undef"){
			int token = lexer->getToken();
			if(token != Lexer::IDENTIFIER){
				reportError("undef expects an identifier", __LINE__);
			}
			const string& macroName = lexer->getTokenContent();
			globalScope->removeMacro(macroName);
		}else if(name == "include"){
			string path;
			lexer->skipToEndOfLine(path);
			includeFile(path.c_str());
		}else if(name == "error"){
			string errMsg;
			lexer->skipToEndOfLine(errMsg);
			reportError(errMsg.c_str());
			exit(__LINE__);
		}else if(name == "ifdef" || name == "if" || name == "ifndef"){
			int nameSize = name.size();
			bool isBlockActive = false;
			if(nameSize == 2){
				PreprocessingExpr expr(globalScope);
				lexer->skipToEndOfLine(expr.content);
				if(!expr.evaluate())
				{					
					reportError("invalid expression found", __LINE__);
				}
				isBlockActive = (expr.value != 0);
			}else{
				int token = lexer->getToken();
				string tmp;				
				lexer->skipSpace();
				lexer->skipToEndOfLine(tmp);
				if(token != Lexer::IDENTIFIER){
					reportError("invalid pre processing instruction found", __LINE__);
				}
				if(!tmp.empty()){
					Position pos = lexer->getPosition();
					cerr << "string literal " << tmp << " is ignored at line " << pos.row << ',' << pos.col << endl;
				}
				bool macroDefined = (globalScope->findMacro(lexer->getTokenContent()) != NULL);
				isBlockActive = (macroDefined && nameSize == 5 || (!macroDefined && nameSize == 6));
			}
			if(isBlockActive){
				globalScope->blockActiveStack.push_back(true);
			}else{
				int token = lexer->getToken();
				int currentLevel = 0;
				while(token != Lexer::END){
					if(token == Lexer::NEWLINE){
						token = lexer->getToken();
						if(token == '#'){
							token = lexer->getToken();
							if(token == Lexer::IDENTIFIER){
								if(lexer->getTokenContent() == "endif"){
									if(currentLevel == 0){
										lexer->skipToEndOfLine();
										//blockActiveStack.pop_back();
										break;
									}
									-- currentLevel;
								}else if(lexer->getTokenContent() == "else"){
									if(currentLevel == 0){
										globalScope->blockActiveStack.push_back(true);
										lexer->skipToEndOfLine();
										break;
									}
								}else if(lexer->getTokenContent() == "elif"){
									if(currentLevel == 0){
										
										PreprocessingExpr expr(globalScope);
										lexer->skipToEndOfLine(expr.content);
										if(!expr.evaluate())
										{					
											reportError("invalid expression found");
											exit(__LINE__);
										}
										if(expr.value != 0){
											globalScope->blockActiveStack.push_back(true);
											break;
										}
									}
								}else if(lexer->getTokenContent() == "if" ||
									lexer->getTokenContent() == "ifdef" ||
									lexer->getTokenContent() == "ifndef")
									++currentLevel;
							}
						}else
							continue;
					}
					token = lexer->getToken();
				}
				if(token == Lexer::END){
					reportError("'endif' preprocessing insrtuction expected");
					exit(__LINE__);
				}
			}
		}else if(name == "endif"){
			if(globalScope->blockActiveStack.size() < 2){
				reportError("invalid pre processing instruction found");
				exit(__LINE__);
			}
			globalScope->blockActiveStack.pop_back();
		}else if(name == "else" || name == "elif"){
			lexer->skipToEndOfLine();
			int token = lexer->getToken();
				int currentLevel = 0;
			while(token != Lexer::END){
				if(token == Lexer::NEWLINE){
					token = lexer->getToken();
					if(token == '#'){
						token = lexer->getToken();
						if(token == Lexer::IDENTIFIER){
							if(lexer->getTokenContent() == "endif"){
								if(currentLevel == 0){
									lexer->skipToEndOfLine();
									break;
								}
								--currentLevel;
							}else if(lexer->getTokenContent() == "if" ||
								lexer->getTokenContent() == "ifdef" ||
								lexer->getTokenContent() == "ifndef")
								++currentLevel;
						}else
							continue;
					}else
						continue;
				}
				token = lexer->getToken();
			}
			if(token == Lexer::END){
				reportError("'endif' preprocessing insrtuction expected", __LINE__);
			}
		}else{
			reportError("unregonized preprocessing instruction ", __LINE__);
		}
		return lexer->getToken();
	}
	reportError("invalid processing instruction found");
	exit(__LINE__);
}
bool HeaderParser::parse(){
	if(isParsed)
		return globalScope != NULL;
	isParsed = true;
	globalScope = new GlobalScope(NULL, &globalScopeAttribute);
	currentScope = globalScope;
	lexer = new Lexer(&block);
	std::vector<bool>& blockActiveStack = globalScope->blockActiveStack;
	blockActiveStack.push_back(true);
	//add(new Typedef("size_t", new PrimType(*PrimType::intType())));
	int token = getToken(Lexer::NEWLINE);
	while(token != Lexer::END){
		if(blockActiveStack.back())
			switch(token){ 
			case Lexer::IDENTIFIER:{		
					token = parseDeclaration();
					if(token != ';')
						reportError("';' expected", __LINE__);
					break;
			}
			case ';':break;//empty statement
			case Lexer::STRING:
			case Lexer::CHARACTER:
			case Lexer:: NUMBER:
			case Lexer::OPERATOR:
			default:
				reportError("constant not allowed here ", __LINE__);
			}
		token = getToken();
	}
	return globalScope != NULL;
}

int ExprParser::getToken0(){
	int token = lexer->getToken();
	switch(token){
	case Lexer::IDENTIFIER:return IDENTIFIER;
	case Lexer::NUMBER:
	case Lexer::CHARACTER:return NUMBER;
	case '(':return token;
	case ')':return token;
	case '&':return Expr::BITAND;
	case '|':return Expr::BITOR;
	case '^':return Expr::BITXOR;
	case '+':return Expr::PLUS;
	case '-':return Expr::MINUS;
	case '*':return Expr::MULTI;
	case '/':return Expr::DIV;
	case '%':return Expr::MODULA;
	case '<':return Expr::LESS;
	case '>':return Expr::GREATER;
	case '?': case ':':	case '~':
	case '!': case ',': return token;
	case Lexer::OPERATOR:{
		const string &op = lexer->getTokenContent();
			if(op.size() == 2){
				switch(op[1]){
				case '<':if(op[0] == '<') return Expr::LSHIFT;
				case '>':if(op[0] == '>') return Expr::RSHIFT;
				case '&':if(op[0] == '&') return Expr::AND;
				case '|':if(op[0] == '|') return Expr::OR;
				case '=':
					switch(op[0]){
					case '<':return Expr::LESSEQ;
					case '>':return Expr::GREATEREQ;
					case '=':return Expr::EQUAL;
					case '!':return Expr::NEQUAL;
					}
					break;
				}
			}
			cerr << "not supported operator" << endl;
			return 0;
	}
	case Lexer::STRING:
		cerr << "invalid token found" << endl;
		currentToken = 0;
		break;
	case Lexer::END:return END;
	}
	cerr << "invalid token found" << endl;
	return INVALID;
}
template<class Seq>
void deleteExprs(Seq& seq){
	typename Seq::iterator begin = seq.begin(), end = seq.end();
	while(begin != end){
		delete *begin;
		++begin;
	}
}
Expr *ExprParser::primary(bool get){
	if(get) getToken();
	switch(currentToken){
		case IDENTIFIER:{
			string name;
			name.swap(lexer->getTokenContent());
			getToken();
			if(currentToken == '('){				
				Expr* tmp = logical(true);
				vector<Expr*> exprs;
				while(currentToken == ','){
					if(tmp == NULL){
						cerr << "invalid expression" << endl;
						deleteExprs(exprs);
						return NULL;
					}
					exprs.push_back(tmp);
					tmp = logical(true);
				}
				if(tmp == NULL && !exprs.empty()){
					cerr << "invalid expression" << endl;
					deleteExprs(exprs);
					return NULL;
				}
				exprs.push_back(tmp);
				if(currentToken != ')'){
					cerr << "')' expected" << endl;
					deleteExprs(exprs);
					return NULL;
				}
				getToken();

				FunctionExpr *funcExpr = new FunctionExpr;
				funcExpr->name.swap(name);
				funcExpr->exprs.swap(exprs);
				return funcExpr;

			}
			IdentifierExpr *identExpr = new IdentifierExpr();
			identExpr->name.swap(name);
			return identExpr;
		}
		case NUMBER:{
			NumberExpr *expr = new NumberExpr();
			expr->isDouble = lexer->isDoubleValue();
			if(expr->isDouble)
				expr->number.doubleVal = lexer->getDoubleValue();
			else
				expr->number.intVal = lexer->getIntValue();
			getToken();
			return expr;
		}
		case '(':{
			Expr *expr = logical(true);
			if(currentToken != ')'){
				cerr << "')' expected" << endl;
				delete expr;
				return NULL;
			}
			getToken();
			return expr;
		}
		case '+':case Expr::PLUS:return primary(true);
		case Expr::MINUS:case '-':case '~':case '!':{
			int optr = (char)currentToken;
			return new UnaryExpr(optr, primary(true));
		}
		case END:
			return NULL;
		default:
			cerr << "unrecognized token" << endl;
	}
	return NULL;
}
Expr *ExprParser::bitwise(bool get){
	Expr* expr = primary(get);
	while(Expr::isBitwise(currentToken)){
		int optr = currentToken;
		Expr *tmp = primary(true);
		if(tmp == NULL){
			delete expr;
			return NULL;
		}
		return new BinaryExpr(optr, expr, tmp);
	}
	return expr;
}
Expr *ExprParser::arithmaticalTerm(bool get){
	Expr* expr = bitwise(get);
	for(;;)
		switch(currentToken){
		case Expr::MULTI:
		case Expr::DIV:
		case Expr::MODULA:{
			int optr = currentToken;
			Expr *tmp = bitwise(true);
			if(tmp == NULL){
				delete expr;
				return NULL;
			}
			expr = new BinaryExpr(optr, expr, tmp);
			break;
		}
		default:
			return expr;
		}
	return expr;
}
Expr *ExprParser::arithmatical(bool get){
	Expr* expr = arithmaticalTerm(get);
	for(;;)
		switch(currentToken){
		case Expr::PLUS:
		case Expr::MINUS:{
			int optr = currentToken;
			Expr *tmp = arithmaticalTerm(true);
			if(tmp == NULL){
				delete expr;
				return NULL;
			}
			expr = new BinaryExpr(optr, expr, tmp);
			break;
		}
		default:
			return expr;
		}
	return expr;
}
Expr *ExprParser::shift(bool get){
	Expr* expr = arithmatical(get);
	while(Expr::isShift(currentToken)){
		int optr = currentToken;
		Expr *tmp = arithmatical(true);
		if(tmp == NULL){
			delete expr;
			return NULL;
		}
		expr = new BinaryExpr(optr, expr, tmp);
	}
	return expr;
}
Expr *ExprParser::comparison(bool get){
	Expr* expr = shift(get);
	while(Expr::isComparison(currentToken)){
		int optr = currentToken;
		Expr *tmp = shift(true);
		if(tmp == NULL){
			delete expr;
			return NULL;
		}
		expr = new BinaryExpr(optr, expr, tmp);
	}
	return expr;
}
Expr *ExprParser::logical(bool get){
	Expr* expr = comparison(get);
	while(currentToken != END){
		int optr = currentToken;
		if(optr == '?'){
			optr = currentToken;
			Expr *l = logical(true);
			optr = currentToken;
			if(l == NULL || optr != ':'){
				cerr << "invalid expression" << endl;
				delete l;
				delete expr;
				return NULL;
			}
			Expr *r = logical(true);
			if(r == NULL){
				delete l;
				delete expr;
				return NULL;
			}
			expr = new BinaryExpr('?', expr, new BinaryExpr(':', l, r));
		}else if(Expr::isLogic(optr)){
			Expr *tmp = logical(true);
			if(tmp == NULL){
				delete expr;
				return NULL;
			}
			expr = new BinaryExpr(optr, expr, tmp);
		}else
			break;
	}
	return expr;
}
ExprParser::ExprParser(const Byte *start, const Byte *end)
	:block(start, end), lexer(NULL){
}
Expr *ExprParser::parse(){
	delete lexer;
	currentToken = END;
	lexer = new Lexer(&block);
	Expr *ret = logical(true);
	//getToken();
	if(currentToken != END){
		delete ret;
		cout << "invalid expression " << endl;
		return NULL;
	}
	return ret;
}

bool MacroExpansion::getToken(){
	//Note: since macro parameters maybe replaced by its value, 
	//when this function returns an identifier, it may not be an
	//identifier. But for macro expansion
	//and the algorithm here it's ok to take it as an identifier

	int token = Lexer::END;
Begin:
	token = lexer->getToken();
	switch(token){
	case '#':{
		BlockIterator tmp = lexer->iterator();
		int next = lexer->getToken();
		int index = -1;
		if((next == Lexer::IDENTIFIER) && (index  = findParam(lexer->getTokenContent())) != -1){
			lexer->getTokenContent() = this->paramValues[index];
			currentToken = Lexer::STRING;
			isChanged = true;
			return true;
		}
		lexer->iterator() = tmp;
		currentToken = token;
		return true;
	}
	case Lexer::COMMENT:
		cout << "comment ommited " << endl;
			//add comment here
		goto Begin;
	case Lexer::IDENTIFIER:{
		BlockIterator tmp = lexer->iterator();
		string value0;
		value0.swap(lexer->getTokenContent());
		int index0 = findParam(value0);
		int next = lexer->getToken();
		if(next == '#'){
			next = lexer->getToken();
			if(next == '#'){
				BlockIterator tmp0 = lexer->iterator();
				next = lexer->getToken();
				if(next == Lexer::IDENTIFIER){
					//IDENTIFIER ## IDENTIFIER
					int index1 = findParam(lexer->getTokenContent());
					if(index0 != -1){
						if(index1 != -1){
							string& str = lexer->getTokenContent();
							str = this->paramValues[index0];
							str += this->paramValues[index1];
						}
						else{
							string& str = lexer->getTokenContent();
							value0.swap(str);
							str = this->paramValues[index0];
							str += value0;
						}
					}
					else{
						if(index1 != -1){
							string& str = lexer->getTokenContent();
							str = value0;
							str += this->paramValues[index1];
						}
						else{
							string& str = lexer->getTokenContent();
							value0.swap(str);
							str += value0;
						}
					}
				}else{
					//IDENTIFIER ##
					lexer->iterator() = tmp0;
					lexer->getTokenContent().swap(value0);
				}						
				currentToken = Lexer::IDENTIFIER;
				isChanged = true;
				return true;
			}
		}

		lexer->iterator() = tmp;
		if(index0 != -1){
			lexer->getTokenContent() = this->paramValues[index0];
			isChanged = true;
		}
		else{
			Macro *macro = this->globalScope->findMacro(value0);
			if(macro != NULL && !this->isExpanded(value0)){
				MacroExpansion expansion(macro, globalScope, lexer, this);
				value0.clear();
				if(!expansion.expandMacro(value0))
					return false;
				isChanged = true;
				this->isChanged  = this->isChanged || expansion.isChanged;
			}
			lexer->getTokenContent().swap(value0);
		}
		currentToken = Lexer::IDENTIFIER;
		return true;
	}
	}
	currentToken = token;
	return true;
}

void Lexer::toString(int token, std::string& str){
	switch(token){
	case Lexer::COMMENT:
	case Lexer::END:
	case Lexer::BAD_TOKEN:return;
	case Lexer::OPERATOR:
	case Lexer::NUMBER:
	case Lexer::IDENTIFIER:
		str.push_back(' ');
		str += this->getTokenContent();
		str.push_back(' ');return;
	case Lexer::NEWLINE:
		str.push_back('\n');return;
	case Lexer::CHARACTER:
		str += " '";
		str += this->getTokenContent();
		str += "' ";return;
	case Lexer::STRING:
		str += " \"";
		str += this->getTokenContent();
		str += "\" ";return;
	default:
		str.push_back((char)token);
	}
}
void MacroExpansion::currentTokenToString(std::string& str){
	return lexer->toString(currentToken, str);
}
bool MacroExpansion::expandMacro(std::string& result){
	this->lexer = this->paramLexer;
	if(!macro->params.empty()){
		bool ret = getToken();
		if(!ret)
			return false;
		if(currentToken != '('){
			cerr << "macro invokation expects parameter" << endl;
			return false;
		}
		int currentLevel = 1;
		int paramIndex = 0;
		this->paramValues.push_back(std::string());
		while((ret = getToken()) && currentToken != Lexer::END && currentLevel != 0){
			if(currentToken == ',' && currentLevel == 1){
				++ paramIndex;
				if(paramIndex >= macro->params.size()){
					cerr << "too much parameters for macro expansion" << endl;
					return false;
				}
			}
			if(currentToken == '(')
				++currentLevel;
			else if(currentToken == ')')
				--currentLevel;
			if(currentLevel > 0){
				currentTokenToString(this->paramValues[paramIndex]);
			}
		}
		for(size_t i = 0;i < paramValues.size();++ i)
				Utility::trim(this->paramValues[i]);

		if(currentLevel != 0){
			cerr << "')' expected " << endl;
			return false;
		}
		if(!ret)
			return false;
	}

	Block block((Byte*)&macro->content[0], (Byte*)&macro->content[0] + macro->content.size());
	Lexer tmpLexer(&block);
	this->lexer = &tmpLexer;
	bool ret = false;
	while((ret = getToken()) && currentToken != Lexer::END)
		currentTokenToString(result);
	return ret;
}

/***************************************************************************/

std::pair<std::string, std::string> generateOutputFunctionForStruct(Scope* scope){
  std::string definitions;
  std::string declarations;
  for (auto node : scope->childNodes()){
    NamedCompoundType *type = dynamic_cast<NamedCompoundType*>(node);
    if (nullptr != type){
      std::string oneDef = "std::ostream& operator<<(std::ostream& os, const " + type->getName() + "& field)";
      declarations += oneDef;
      oneDef += "{\n";
      declarations += ";\n";
      auto defScope = type->getScope();
      auto chilren = defScope->childNodes();
      bool isFirstDeclaration = true;
      for (size_t i = 0; i < chilren.size(); ++i){
        Declaration* declaration = dynamic_cast<Declaration*>(chilren[i]);
        if (nullptr != declaration){
          Comment* comment = (i == 0 ? nullptr : dynamic_cast<Comment*>(chilren[i - 1]));
          if (nullptr == comment){
            std::cout << "No comment found for " <<type->getName() << "::" << declaration->getName() << std::endl;
          }
          if (isFirstDeclaration) oneDef += "  os << \"";
          else oneDef += "\n     << \"\\n";
          oneDef += (nullptr == comment ? declaration->getName() : comment->getContent());
          oneDef += " : \" << field.";
          oneDef += declaration->getName();
          isFirstDeclaration = false;
        }
      }
      oneDef += ";\n  return os;\n}\n";
      definitions += oneDef;
    }
  }
  return{ declarations, definitions };
}
int main(int argc, char* argv[])
{

	string str;
	bool tmp = Utility::fileToString("D:\\ThostFtdcUserApiStruct.h", str);
	assert(tmp);
	const char *start =str.c_str();
	int count = str.size();
	/*Block block((const Byte*)start, (const Byte *)start + count);
	Lexer lexer(&block);
	Position pos = lexer.getPosition();
	int token = lexer.getToken();
	while(token != Lexer::END){
		if(token == Lexer::COMMENT){
			cout << pos.row << ',' << pos.col << ':' << lexer.getTokenContent() << endl;
		}
		pos = lexer.getPosition();
		token = lexer.getToken();
	}*/
	HeaderParser headerParser((const Byte*)start, (const Byte*)start + count);
	headerParser.parse();
	string expression("-((!defined (_MSC_VER,func(a,bc,d),t) && (_MSC_VER >= 800))+(asdf)) + (a>b?fd+2:c>e?g:h)");
	//string expression("a+b*c");
	ExprParser parser((Byte*)&expression[0], (Byte*)&expression[0] + expression.size());
	Expr *expr = parser.parse();
	ExprToString toString;
	toString.visit(expr);
	//cout << toString.getString() << endl;
	//string result;
	//headerParser.getGlobalScope()->toString(result);
  auto ret = generateOutputFunctionForStruct(headerParser.getGlobalScope());
  std::cout << "DECLARATIONS: \n";
  std::cout << ret.first << std::endl;
  std::cout << "DEFINITIONS: \n";
  std::cout << ret.second << std::endl;
	cout << endl << endl;
	//cout << result << endl;
	cin.get();
	return 0;
}

