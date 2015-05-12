// HeaderParser.cpp
//

#include "stdafx.h"
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
		int declarationeEndToken;
		int varListSeparator;
		bool allowFunctionDefine;
		bool allowTypeDefine;
		bool allowVarList;
		bool allowCompondType;
		bool allowInitialization;
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
			globalScopeAttribute.declarationeEndToken = ';';
			globalScopeAttribute.allowFunctionDefine = true;
			globalScopeAttribute.allowTypeDefine = true;
			globalScopeAttribute.allowVarList = true;
			globalScopeAttribute.allowCompondType = true;
			globalScopeAttribute.allowInitialization = true;
			globalScopeAttribute.varListSeparator = ',';

			paramListAttribute.startToken =  '(';
			paramListAttribute.endToken = ')';
			paramListAttribute.declarationeEndToken = ',';
			paramListAttribute.allowFunctionDefine = false;
			paramListAttribute.allowTypeDefine = false;
			paramListAttribute.allowVarList = true;
			paramListAttribute.allowCompondType = false;
			globalScopeAttribute.varListSeparator = ',';
			
			compoundTypeScopeAttribute.startToken = '{';
			compoundTypeScopeAttribute.endToken = '}';
			compoundTypeScopeAttribute.declarationeEndToken = ';';
			compoundTypeScopeAttribute.allowFunctionDefine = false;
			compoundTypeScopeAttribute.allowTypeDefine = false;
			compoundTypeScopeAttribute.allowVarList = true;
			compoundTypeScopeAttribute.allowCompondType = true;
			globalScopeAttribute.varListSeparator = ',';
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
			ret.pop_back();
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
	if(modifiers.storage != this->storage || modifiers.storage * this->storage != 0)
		return false;
	if((modifiers.signedBit | this->signedBit) == 1
		&&((modifiers.unsignedBit | this->unsignedBit) == 1))
		return false;
	if((modifiers.shortBit | this->shortBit) == 1
		&& (modifiers.longBit1 | modifiers.longBit2 | this->longBit1 | this->longBit2) == 1)
		return false;
	return true;
}
static bool setNodeType(Type *type, Node *node){
	PointerType *pointer = dynamic_cast<PointerType*>(node);
	if(pointer != NULL){
		if(pointer->getSrcType() != NULL){
			cerr << "type already specified" << endl;
			return false;
		}
		pointer->setSrcType(type);
		return true;
	}
	Declaration *declare = dynamic_cast<Declaration*>(node);
	if(declare != NULL){
		if(declare->getType() != NULL){
			cerr << "type already specified" << endl;
			return false;
		}
		declare->setType(type);
		return true;
	}
	cerr << "unexpected node type found" << endl;
	return false;
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
Type *GlobalScope::findType(const std::string& name){
	map<string, int>::iterator pos = declares.find(name);
	if(pos == declares.end())
		return NULL;
	return dynamic_cast<Type*>(nodes[pos->second]);
}
Node *GlobalScope::findDeclare(const std::string& name){
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
	str.swap(tmp);
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
void Modifiers::toString(std::string& str){
	string& tmp = str;
	if(isAuto()){
		tmp += Keyword::autoKeyword().name();
		tmp.push_back(' ');
	}
	if(isStatic()){
		tmp += Keyword::staticKeyword().name();
		tmp.push_back(' ');
	}
	if(isRegister()){
		tmp += Keyword::registerKeyword().name();
		tmp.push_back(' ');
	}
	if(isTypedef()){
		tmp += Keyword::typedefKeyword().name();
		tmp.push_back(' ');
	}
	if(isExtern()){
		tmp += Keyword::externKeyword().name();
		tmp.push_back(' ');
	}


	if(isConst()){
		tmp += Keyword::constKeyword().name();
		tmp.push_back(' ');
	}
	if(isVolatile()){
		tmp += Keyword::volatileKeyword().name();
		tmp.push_back(' ');
	}
	if(isLong1()){
		tmp += Keyword::longKeyword().name();
		tmp.push_back(' ');
	}
	if(isLong2()){
		tmp += Keyword::longKeyword().name();
		tmp.push_back(' ');
	}
	if(isShort()){
		tmp += Keyword::shortKeyword().name();
		tmp.push_back(' ');
	}
	if(isSigned()){
		tmp += Keyword::signedKeyword().name();
		tmp.push_back(' ');
	}
	if(isUnsigned()){
		tmp += Keyword::unsignedKeyword().name();
		tmp.push_back(' ');
	}
}
void PrimType::toString(std::string& str)const {
	string tmp;
	tmp += keyword.name();
	tmp.push_back(' ');
	tmp += str;
	tmp.swap(str);
}
void PointerType::toString(std::string& str)const {
	string tmp;
	tmp.push_back('*');
	tmp += str;
	if(srcType != NULL){
		srcType->toString(tmp);
	}
	str.swap(tmp);
}
void ArrayType::toString(std::string& str)const{
	str += "[]";
	if(srcType != NULL){
		srcType->toString(str);
	}
}

void ParameterList::toString(std::string& str)const {
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
	tmp.push_back('(');
	if(params != NULL)
		params->toString(tmp);
	else
		tmp += "void";
	tmp.push_back(')');
	retType->toString(tmp);
	str.swap(tmp);
}
void Typedef::toString(std::string& str) const{
	string tmp;
	if(!tmp.empty())
		tmp.push_back(' ');
	tmp += name;
	tmp.push_back(' ');
	tmp += str;
	str.swap(tmp);
}
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

bool GlobalScope::add(Node *node){
	Macro *macro = dynamic_cast<Macro*>(node);
	if(macro != NULL)
	{
		map<string, int>::iterator pos = macros.find(macro->name);
		if(pos == macros.end()){
			macros[macro->name] = nodes.size();
			nodes.push_back(node);
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
		nodes.push_back(node);
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
		MacroExpansion macroExpansion(macro, this->globalScope, this->lexer);
		string result;
		if(!macroExpansion.expandMacro(result)){
			this->reportError("invalid macro invokation");
			return false;
		}
		Utility::trim(result);
		if(!macroExpansion.changed() && !macro->content.empty())
			blocks.push_back(Block((const Byte*)macro->content.c_str(), (const Byte*)macro->content.c_str() + macro->content.size(), 
				(int)Block::NOT_UPDATE_POS | (int)Block::IS_MACRO_EXPANSION));
		else if(!result.empty()){
			expandedMacros.push_back(std::string());
			string& str = expandedMacros.back();
			str.swap(result);
			blocks.push_back(Block((const Byte*)str.c_str(), (const Byte*)str.c_str() + str.size(), 
				(int)Block::NOT_UPDATE_POS | (int)Block::IS_MACRO_EXPANSION));
		}
		lexer->iterator().split(&blocks.back());
			return true;
	}
	return false;
}
int HeaderParser::parseMacro(){
		int token = lexer->getToken();
		if(token != Lexer::IDENTIFIER){
				reportError("macro definition expects an identifier");
				exit(__LINE__);
		}
		if(globalScope->findMacro(lexer->getTokenContent()) != NULL){
				reportError("macro redefined");
				exit(__LINE__);
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
						reportError("identifier expected, invalid macro definition");
						exit(__LINE__);
					}
					p->params.push_back(string());
					p->params.back().swap(lexer->getTokenContent());
					if(p->params.end() != ::find(p->params.begin(), p->params.end(), lexer->getTokenContent())){
						reportError("duplicate parameter name  in macro");
						exit(__LINE__);
					}
					token = lexer->getToken();
					if(token == ','){
						token = lexer->getToken();
						if(token != Lexer::IDENTIFIER){
							reportError("invalid macro definition");
							exit(__LINE__);
						}
					}
				}
				if(token == Lexer::NEWLINE){
					reportError("')' expected, invalid macro definition");
					exit(__LINE__);
				}
				lexer->skipToEndOfLine(p->content);
		}else{
				reportError("unexpected macro definition");
				exit(__LINE__);
		}
		globalScope->add(p);
		return getToken();
}
void HeaderParser::includeFile(const char *fileName){
	static string baseDir = "D:\\C++\\Source\\Include\\";

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
int HeaderParser::skipBraces(char ch){
	vector<char> str;
	str.push_back(ch);
	int token = getToken();
	while(token != Lexer::END){
		switch(token){
			case '(':case '[':case '{':
				str.push_back(token);break;
			case ')':case ']':case '}':{
				char expected = Utility::getMatchedBrace(str.back());
				if(expected != token){
					reportError("unmatched right brace found");
					exit(__LINE__);
				}
				str.pop_back();
				break;
			}
		}
		if(str.empty())
			break;
		token = getToken();
	}
	if(token == Lexer::END){
		reportError("right brace expected");
		exit(__LINE__);
	}
	return getToken();
}

int HeaderParser::parseVariableList(int token, ModifiedType *type){
	Modifiers storage = type->getModifiers();
	Modifiers modifiers;
	ModifiedType *decType = type;
	Node *decNode = NULL;
	const Scope::ScopeAttribute *att = this->currentScope->getAttribute();
	while(token != att->declarationeEndToken && token != att->endToken && token != Lexer::END){
start:
		switch(token){
			case '*':
				if(decNode != NULL){
					reportError("'*' not allowed", __LINE__);break;
				}
				decType->setModifiers(modifiers);
				decType = new ModifiedType(new PointerType(decType));
				modifiers = Modifiers();
				token = getToken();
				break;
			case '(':{
				//could be type or declaration
				if(decNode != NULL){
					reportError("'(' not allowed", __LINE__);break;
				}
				Node *outterMost = NULL;
				Node *innerMost = NULL;
				token = parseDeclares(innerMost, outterMost);
				token = parsePostDeclarator(token, decType);
				if(!setNodeType(decType, outterMost)){
					reportError("error occured", __LINE__);break;
				}
				FunctionType *funcType = NULL;
				Declaration *declare = dynamic_cast<Declaration*>(innerMost);
				if(declare != NULL)
					funcType = dynamic_cast<FunctionType *>(declare->getType());
				else {
					ModifiedType *tmpType = dynamic_cast<ModifiedType*>(innerMost);
					if(tmpType != NULL)
						decType = tmpType;
					funcType = dynamic_cast<FunctionType *>(innerMost);
				}
				if(funcType != NULL){
					if(storage.isStorageSet() && !storage.isStatic() && !storage.isExtern()){
						reportError("storage type not allowed for function definition", __LINE__);break;
					}
					if(declare == NULL){
						reportError("function defintion requires a function name ", __LINE__);break;
					}
					decNode = new FunctionDeclaration(declare->getName(), funcType);
				}else if(declare != NULL){
					if(storage.isTypedef()){
						decNode = new Typedef(declare->getName(), declare->getType());
					}else{
						decNode = new VariableDeclaration(declare->getName(), declare->getType());
					}
				}
				if(token == '='){
					if(storage.isTypedef()){
						reportError("could not initialize a typedef", __LINE__);break;
					}
					if(!att->allowInitialization){
						reportError("could not initialize a member or formal parameter", __LINE__);break;
					}
					token = getToken();
					while(token != att->endToken && token != att->varListSeparator &&
						token != att->declarationeEndToken && token != Lexer::END)
						token =  getToken();
					}
				if(funcType && token != att->declarationeEndToken){
					reportError("syntax error", __LINE__);
				}
				break;
			}
			case Lexer::IDENTIFIER:{
				if(decNode != NULL){
					reportError("identifier not allowed", __LINE__);break;
				}
				string name;
				name.swap(lexer->getTokenContent());
				Keyword keyword(name.c_str());
				switch(keyword.category().getIndex()){
					case KeywordCategory::MODIFIER:{
						Modifier modifier(keyword);
						if(modifier.category() != ModifierCategory::cv()){
							reportError("syntax error, invalid identifier found", __LINE__);break;
						}
						if(!modifiers.add(Modifier(keyword))){
							reportError("incompatible modifier", __LINE__);break;
						}
						if(!modifiers.check(decType)){
							reportError("modifiers incompatible with type", __LINE__);break;
						}
						token = getToken();
						break;
					}
					case KeywordCategory::INVALID:{	
						if(currentScope->findDeclare(name) != NULL){
							reportError("duplicate declaration of identifier", __LINE__);break;
						}
						if(decType != type){
							decType->setModifiers(modifiers);
							modifiers = Modifiers();
						}			
						token = getToken();
						if(token == ':'){
							cerr << "bit type here" << endl;
						}else{
							token = parsePostDeclarator(token, decType);
						}
						FunctionType *funcType = dynamic_cast<FunctionType*>(decType);						
						if(funcType != NULL){
							if(storage.isStorageSet() && !storage.isStatic() && !storage.isExtern()){
								reportError("storage type not allowed for function definition", __LINE__);break;
							}
							decNode = new FunctionDeclaration(name, funcType);							
						}else{
							if(storage.isTypedef()){
								decNode = new Typedef(name, decType);
							}else{
								decNode = new VariableDeclaration(name, decType);
							}
						}
						if(token == '='){
							if(storage.isTypedef()){
								reportError("could not initialize a typedef", __LINE__);break;
							}
							if(att->allowInitialization){
								reportError("initialization not allowed", __LINE__);break;
							}
							token = getToken();
							while(token != att->declarationeEndToken && 
								token != att->varListSeparator && token != att->endToken && token != Lexer::END)
								token =  getToken();
						}
						if(funcType && (att->allowFunctionDefine || token != att->declarationeEndToken)){
								reportError("syntax error", __LINE__);break;
						}
						break;
					}
					default:
						reportError("syntax error, invalid token found", __LINE__);break;
				}
			}
			break;
			case '[':
				token = parsePostDeclarator(token, decType);
				break;
			default:
				if(token == att->varListSeparator){
					if(!att->allowVarList){
						reportError("only one variable allowed", __LINE__);break;
					}
					if(decNode != NULL){
						currentScope->add(decNode);
						decNode = NULL;
					}else{
						currentScope->add(decType);
					}
					decType = type;
					modifiers = Modifiers();
					token = getToken();
					goto start;
				}
				reportError("syntax error", __LINE__);
		}
	}
	if(decNode != NULL)
		currentScope->add(decNode); 
	else{
		currentScope->add(decType);
	}
	return token;
	return 0;
}

int HeaderParser::parseDeclares(Node *&innerMost, Node *&outterMost){
	int token = getToken();
	Modifiers modifiers;
	ModifiedType *tmpType = new ModifiedType(NULL);
	while(token != Lexer::END && token != ')'){
		switch(token){
			case '*':{
				tmpType->setType(new PointerType(tmpType));
				tmpType->setModifiers(modifiers);
				modifiers = Modifiers();
				if(outterMost == NULL)
					outterMost = tmpType;
				innerMost = tmpType;
				token = getToken();
				break;
			}
			case '(':{
				Node *outterTemp = NULL;
				token = parseDeclares(innerMost, outterTemp);
				token = parsePostDeclarator(token, tmpType);
				if(outterMost == NULL)
					outterMost = outterTemp;
				if(!setNodeType(tmpType, outterTemp)){
					reportError("error occured", __LINE__);break;
				}
				if(token != ')'){
					reportError("')' expected", __LINE__);break;
				}
				break;
			}
			case '[':{
				if(tmpType->getType() == NULL){
					reportError("unexpected token '['", __LINE__);break;
				}
				//void (*(*func)(int, int)[])(int, int) is invalid;
				if(dynamic_cast<const FunctionType*>(tmpType->getType()) != NULL){
					reportError("array of function type is not allowned", __LINE__);break;
				}
				tmpType->setType(new ArrayType(tmpType->getType()));
				token = skipBraces('[');
				break;
			}
			case Lexer::IDENTIFIER:{
				string name;
				name.swap(lexer->getTokenContent());
				Keyword keyword(name.c_str());
				switch(keyword.category().getIndex()){
					case KeywordCategory::MODIFIER:{
						Modifier modifier(keyword);
						if(modifier.category() != ModifierCategory::cv()){
							reportError("syntax error, modifier not allowed", __LINE__);break;
						}
						if(!modifiers.add(Modifier(keyword))){
							reportError("incompatible modifier at line ", __LINE__);
						}
						if(tmpType == NULL){
							reportError("modifier not allowed", __LINE__);break;
						}
						if(!modifiers.check(tmpType)){
							reportError("modifiers incompatible with type", __LINE__);break;
						}
						tmpType->setModifiers(modifiers);
						token = getToken();
						break;
					}
					case KeywordCategory::INVALID:{	
						if(currentScope->findDeclare(name) != NULL){
							reportError("duplicate declaration of identifier ", __LINE__);
							break;
						}
						token = getToken();
						token = parsePostDeclarator(token, tmpType);
						FunctionType *funcType = dynamic_cast<FunctionType*>(tmpType);
						Node *declare = new Declaration(name, tmpType);
						if(token != ')'){
							reportError("')' expected", __LINE__);
							break;
						}
						if(outterMost == NULL)
							outterMost = declare;
						innerMost = declare;
						break;
					}
					default:
						reportError("syntax error, identifier not allowed", __LINE__);
				}
			}
			break;
			default:break;
		}
	}
	return getToken();
}

int HeaderParser::parsePostDeclarator(int token, ModifiedType *&type){
	if(token == '('){
		ParameterList *node = new ParameterList(currentScope, &paramListAttribute);
		currentScope = node;
		token = getToken();
		while(token == Lexer::IDENTIFIER){
				token = parseDeclaration();
		}
		currentScope = node->getParent();
		type = new ModifiedType(new FunctionType(type, node));
		return token;
	}
	while(token == '['){
		type->setType(new ArrayType(type->getType()));
		token = skipBraces('[');
	}
	return token;
}
int HeaderParser::parseDefinition(Type **type){
	int token = getToken();
	CompoundType *compoundType = NULL;
	ScopeAttribute *att = currentScope->getAttribute();
	while(token != Lexer::END){
		if(token == Lexer::IDENTIFIER){
			Node *node = findDeclare(lexer->getTokenContent());
			if(node == NULL){
				compoundType = new NamedCompoundType(lexer->getTokenContent());
				add(compoundType);
			}else{
				compoundType = dynamic_cast<CompoundType*>(node);
				if(compoundType == NULL){
					reportError("identifier redefined", __LINE__);break;
				}
			}
			token = getToken();
		}
		if(token == '{'){
			if(att->allowCompondType){
				reportError("struct definition not allowed", __LINE__);break;
			}
			if(compoundType == NULL)
				compoundType = new AnonymousCompoundType();
			if(compoundType->isDefined()){
				reportError("identifier redefined", __LINE__);break;
			}
			Scope *tmpNode = new ;
			currentNode = compoundType;
			token = getToken(lexer);
			while(token == Lexer::IDENTIFIER){
				token = parseDeclaration(lexer);
			}
			currentNode = tmpNode;
			if(token != '}'){
				cerr << "'}' expected at line " << lexer.getLineNo() << endl;
				exit(__LINE__);
			}
			token = getToken(lexer);
			compoundType->setDefined(true);
		}
	//}
	*type = compoundType;
	if(*type == NULL){
		cerr << "type required at line " << lexer.getLineNo() << endl;
		exit(__LINE__);
	}
	return token;
}
int HeaderParser::parseDeclaration(){ 
	int token = Lexer::IDENTIFIER;
	Modifiers modifiers;
	ModifiedType *type = NULL;
	const string& name = lexer->getTokenContent();
	bool hasVarList = false;
	const Scope::ScopeAttribute *att = this->currentScope->getAttribute();
	while(token != ';' && token != att->declarationeEndToken
		&& token != att->endToken && token != Lexer::END){
		bool isVarList = true;
		if(token == Lexer::IDENTIFIER){
				Keyword keyword(name.c_str());
				isVarList = false;
				switch(keyword.category().getIndex()){
					case KeywordCategory::OPTER:
					case KeywordCategory::CONTROL: reportError("keyword  unexpected", __LINE__);break;
					case KeywordCategory::MODIFIER:{
						if(keyword == Keyword::autoKeyword() || keyword == Keyword::registerKeyword() || !modifiers.add(Modifier(keyword))){
							reportError("modifier not allowed", __LINE__);break;
						}
						break;
					}
					case KeywordCategory::PRIM:
						if(type != NULL){
							reportError("multiple type specified", __LINE__);
							delete type;break;
						}
						{
							const PrimType * primType = PrimType::get(keyword);
							assert(primType != NULL);
							if(!modifiers.check(primType)){
								reportError("modifiers incompatible with type", __LINE__);break;
							}
							type = new ModifiedType(primType);
						}
						break;
					case KeywordCategory::METATYPE:{
						Type *tmpType = NULL;
						token = parseDefinition(&tmpType);
						assert(tmpType != NULL);
						type = new ModifiedType(tmpType);
						continue;
					}
					case KeywordCategory::INVALID:{
						bool isType = false;
						if(type == NULL){
							const Type *tmpType = currentScope->findType(name);
							isType = (tmpType != NULL);
							if(!isType){
								if(modifiers.hasIntModifiers()){
									tmpType = PrimType::intType();
									assert(tmpType != NULL);
								}
							}
							if(tmpType == NULL){
									reportError("unexpected indentifier, type required", __LINE__);break;
							}
						}
						isVarList = !isType;
					}
				}
		}
		
		if(isVarList){
			if(type == NULL && modifiers.hasIntModifiers()){
					const PrimType * primType = PrimType::intType();
					assert(primType != NULL);
					type = new ModifiedType(primType);
			}
			if(!type->setModifiers(modifiers)){
				reportError("modifiers incompatible with type", __LINE__);
			}
			hasVarList = true;
			token = parseVariableList(token, type);
			modifiers = Modifiers();
			type = NULL;
			continue;
		}
		token = getToken();
	}
	if(!hasVarList){
		currentScope->add(type);
	}
	return token;
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
				reportError("undef expects an identifier");
				exit(__LINE__);
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
					reportError("invalid expression found");
					exit(__LINE__);
				}
				isBlockActive = (expr.value != 0);
			}else{
				int token = lexer->getToken();
				string tmp;				
				lexer->skipSpace();
				lexer->skipToEndOfLine(tmp);
				if(token != Lexer::IDENTIFIER){
					reportError("invalid pre processing instruction found");
					exit(__LINE__);
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
						}
					}
				}
				token = lexer->getToken();
			}
			if(token == Lexer::END){
				reportError("'endif' preprocessing insrtuction expected");
				exit(__LINE__);
			}
		}else{
			reportError("unregonized preprocessing instruction ");
			exit(__LINE__);
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
						reportError("';' expected");
					break;
			}
			case ';':break;//empty statement
			case Lexer::STRING:
			case Lexer::CHARACTER:
			case Lexer:: NUMBER:
			case Lexer::OPERATOR:
			default:
				reportError("constant not allowed here ");
				exit(__LINE__);
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
	//when this function returns and identifier token, it may not be an
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
	
			/* COMMENT = 50000,
		  = 20002, STRING = 2000*/
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

int _tmain(int argc, _TCHAR* argv[])
{
	string str;
	bool tmp = Utility::fileToString("d:\\jpeglib.h", str);
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
	cout << toString.getString() << endl;
	cin.get();
	return 0;
}

