#include <string>
#include <vector>
#include <map>
#include <cassert>
#include <algorithm>
#include <functional>
#include <streambuf>
#include <fstream>
#include <sstream>
#include <list>
#include <iostream>
namespace tmlg{
	typedef unsigned char Byte;
	
	/****************************************/
	/*This class is used to identify keyword cateogries*/
	class Keyword;
	class KeywordCategory{
		friend class Keyword;
		Byte index;
		KeywordCategory();
		KeywordCategory(Byte index) : index(index){}
	public:
		enum { INVALID = 0, MODIFIER, CONTROL, PRIM, METATYPE, 
			OPTER};
		operator Byte(){
			return index;
		}
		bool isValid(){
			return index != 0;
		}
		Byte getIndex(){
			return index;
		}
		static KeywordCategory modifier(){
			return KeywordCategory(MODIFIER);
		}
		static KeywordCategory control(){
			return KeywordCategory(CONTROL);
		}
		static KeywordCategory prim(){
			return KeywordCategory(PRIM);
		}
		static KeywordCategory metatype(){
			return KeywordCategory(METATYPE);
		}
		static KeywordCategory opter(){
			return KeywordCategory(OPTER);
		}
		static KeywordCategory invalid(){
			return KeywordCategory(INVALID);
		}
	};
	/****************************************/
	/*This class is used to identify keyword modifier categories*/
	class ModifierCategory{
		friend class Modifier;
		Byte index;
		ModifierCategory(Byte index) : index(index){}
	public:
		enum{INVALID = 0, STORAGE, CV , SIGN, LENGTH, };
		ModifierCategory() : index(0){}
		operator Byte(){
			return index;
		}
		bool isValid(){
			return index != 0;
		}
		Byte getIndex(){
			return index;
		}
		static ModifierCategory storage(){return ModifierCategory(STORAGE);}
		static ModifierCategory cv(){return ModifierCategory(CV);}
		static ModifierCategory sign(){return ModifierCategory(SIGN);}
		static ModifierCategory length(){return ModifierCategory(LENGTH);}
		static ModifierCategory invalid(){return ModifierCategory(INVALID);}
	};
	/****************************************/
	/*This class is a wrapper class for the keywords of c language*/
	class Keyword{
		Byte index;
	public:
		/*The order of these enum defitions should be changed with care*/
		enum {INVALID = 0, 
			/*modifiers*/
			/*storage modifiers*/
			AUTO, EXTERN, REGISTER, STATIC, TYPEDEF, //5
			/*const-volatile modifiers*/
			VOLATILE, CONST, //7
			/*length modifiers*/
			LONG, SHORT, //9
			/*sign modifiers*/
			SIGNED, UNSIGNED, //11
		/*primitive tyeps*/
		CHAR, VOID, DOUBLE, INT, FLOAT, //16
		/*control keywords*/
		BREAK, CASE, CONTINUE, DEFAULT, DO, ELSE, FOR, GOTO, 
		IF, WHILE, RETURN, SWITCH, //28
		/*operator keywords*/
		SIZEOF, //29
		/*metatype key words*/
		ENUM, STRUCT, UNION, //32
		KEYWORD_COUNT,
		};
		Keyword(Byte index) : index(index){}
		Keyword(const char *name=NULL);
		const char *name() const{
			if(index > UNION)
				return keywords[INVALID];
			return keywords[index];
		}
		bool isValid() const{
			return index >= AUTO && index <= UNION;
		}
		Byte getIndex(){
			return index;
		}
		ModifierCategory modiferCatetory() const{
			if(index >= AUTO && index <= TYPEDEF)
				return ModifierCategory::storage();
			if(index >= VOLATILE && index <= CONST)
				return ModifierCategory::cv();
			if(index >= LONG && index <= SHORT)
				return ModifierCategory::length();
			if(index >= SIGNED && index <= UNSIGNED)
				return ModifierCategory::length();
			return ModifierCategory::invalid();
		}
		KeywordCategory category() const{
			if(index >= AUTO && index <= UNSIGNED)
				return KeywordCategory::modifier();
			if(index >= CHAR && index <= FLOAT)
				return KeywordCategory::prim();
			if(index >= BREAK && index <= SWITCH)
				return KeywordCategory::control();
			if(index >= SIZEOF && index <= SIZEOF)
				return KeywordCategory::opter();
			if(index >= ENUM && index <= UNION)
				return KeywordCategory::metatype();
			return KeywordCategory::invalid();
		}
		bool operator<(Keyword w) const{
			return index < w.index;
		}
		bool operator ==(Keyword w) const{
			return index == w.index;
		}
		static Keyword autoKeyword(){return Keyword(AUTO);}
		static Keyword externKeyword(){return Keyword(EXTERN);}
		static Keyword registerKeyword(){return Keyword(REGISTER);}
		static Keyword staticKeyword(){return Keyword(STATIC);}
		static Keyword typedefKeyword(){return Keyword(TYPEDEF);}

		static Keyword constKeyword(){return Keyword(CONST);}
		static Keyword volatileKeyword(){return Keyword(VOLATILE);}
		
		static Keyword longKeyword(){return Keyword(LONG);}
		static Keyword shortKeyword(){return Keyword(SHORT);}

		static Keyword signedKeyword(){return Keyword(SIGNED);}
		static Keyword unsignedKeyword(){return Keyword(UNSIGNED);}


		static Keyword charKeyword(){return Keyword(CHAR);}
		static Keyword voidKeyword(){return Keyword(VOID);}
		static Keyword doubleKeyword(){return Keyword(DOUBLE);}
		static Keyword intKeyword(){return Keyword(INT);}
		static Keyword floatKeyword(){return Keyword(FLOAT);}

		static Keyword enumKeyword(){return Keyword(ENUM);}
		static Keyword structKeyword(){return Keyword(STRUCT);}
		static Keyword unionKeyword(){return Keyword(UNION);}
		static const char *keywords[];
	};
	/****************************************/
	class Modifier{
		Keyword keyWord;
	public:
		Modifier(Keyword keyWord) : keyWord(keyWord){
		}
		bool isValid(){
			return keyWord.modiferCatetory().isValid();
		}
		ModifierCategory category(){
			return keyWord.modiferCatetory();
		}
		Keyword getKeyword() const{
			return keyWord;
		}
		const char *name() const{
			return getKeyword().name();
		}
		bool operator<(Modifier m) const{
			return keyWord < m.keyWord;
		}
	};
	/***************************************/
	class Type;
	class Modifiers{
	private:
		unsigned int storage : 3;
		unsigned int constBit : 1;
		unsigned int volatileBit : 1;
		unsigned int longBit1 : 1;
		unsigned int longBit2 : 1;
		unsigned int shortBit : 1;
		unsigned int signedBit : 1;
		unsigned int unsignedBit : 1;
	public:
		Modifiers(){
			memset(this, 0, sizeof(*this));
		}
		bool isAuto()const{return storage == 1;}
		bool isStatic()const{return storage == 2;}
		bool isRegister()const{return storage == 3;}
		bool isTypedef()const{return storage == 4;}
		bool isExtern()const{return storage == 5;}
		bool isStorageSet()const{return storage != 0;}

		bool isConst()const{return constBit != 0;}
		bool isVolatile()const{return volatileBit != 0;}
		bool isLong1()const{return longBit1 != 0;}
		bool isLong2()const{return longBit2 != 0;}
		bool isShort()const{return shortBit != 0;}
		bool isSigned()const{return signedBit != 0;}
		bool isUnsigned()const{return unsignedBit != 0;}
		bool hasIntModifiers(){
			return isLong1() || isLong2() || isShort()
					|| isSigned() || isUnsigned();
		}
		void clearStorage(){
			storage = 0;
		}
		void clearCV(){
			constBit = volatileBit = 0;
		}
		void clearSign(){
			signedBit = unsignedBit = 0;
		}
		void clearLength(){
			longBit1 = longBit2 = shortBit = 0;
		}
		bool add(Modifier modifier, bool duplicateRet = true);
		bool check(Modifier modifier);
		bool check(Modifiers modifiers);
		bool check(const Type *type);
		void toString(std::string& str);
	};
	/****************************************/
	class Node{
	public:
		virtual~Node()=0;
		virtual void toString(std::string& str)const=0;
	};

	/****************************************/
	class Type : public Node{
	};
	/****************************************/
	
	class Scope : public Node{
	public:
		class ScopeAttribute;
		Scope(Scope *parent = NULL, ScopeAttribute *att = NULL) : parent(parent), att(att){}
		virtual Type *findType(const std::string& name) = 0;
		virtual Node *findDeclare(const std::string& name) = 0;
		const ScopeAttribute *getAttribute() const {
			return att;
		}
		virtual bool add(Node *node) = 0;
		virtual void toString(std::string& str) const{
		}
		Scope *getParent(){
			return parent;
		}
	protected:
		Scope *parent;
		ScopeAttribute *att;
		std::vector<Node*> nodes;
	};
	/****************************************/
	class ModifiedType : public Type{
	protected:
		const Type *srcType;
		Modifiers modifiers;
	public:
		ModifiedType(const Type *srcType) : srcType(srcType){
			if(dynamic_cast<const ModifiedType*>(srcType) != NULL)
				assert(false);
		}
		bool setModifiers(Modifiers modifiers){
			bool ret = modifiers.check(srcType);
			if(ret)
				this->modifiers = modifiers;
			return ret;
		}
		const Type* getUnModifiedType() const{
			const ModifiedType* tmp = dynamic_cast<const ModifiedType*>(srcType);
			if(tmp == NULL)
				return srcType;
			return tmp->getUnModifiedType();
		}
		Modifiers& getModifiers(){
			return modifiers;
		}
		const Type *getType() const{
			return srcType;
		}
		void setType(Type *type){
			srcType = type;
		}
		Modifiers getModifiers() const{
			return modifiers;
		}
		virtual void toString(std::string& str) const{
		}
	};
	/****************************************/
	class Declaration : public Node{
	protected:
		Modifiers modifiers;
		std::string name;
		Type *type;
	public:
		Declaration(const std::string& name, 
			Type *type) : name(name), type(type){
		}
		const std::string& getName() const{
			return name;
		}
		Type *getType() const{
			return type;
		}
		void setType(Type *type){
			this->type = type;
		}
		virtual void toString(std::string& str)const;
	};

	class Macro : public Node{
	public:
		Macro(){}
		std::string name;
		std::vector<std::string> params;
		std::string content;
		virtual void toString(std::string& str)const;
	};

	/****************************************/
	class Comment : public Node{
		std::string content;
	public:
		Comment(){
		}
		Comment(const std::string& conent) : content(content){
		}
		const std::string& getContent() const{
			return content;
		}
		std::string& getContent(){
			return content;
		}
		virtual void toString(std::string& str)const;
	};


	/****************************************/
	class PrimType : public Type{
		Keyword keyword;
		PrimType();
		PrimType(Keyword keyword) : keyword(keyword){
		}
	public:
		bool equals(const PrimType *type) const{
			return keyword == type->keyword;
		}
		virtual void toString(std::string& str)const;
		static const PrimType *get(Keyword keyword){
			if(keyword == intType()->keyword)
				return intType();
			if(keyword == charType()->keyword)
				return charType();
			if(keyword == doubleType()->keyword)
				return doubleType();
			if(keyword == floatType()->keyword)
				return floatType();
			if(keyword == voidType()->keyword)
				return voidType();
			return NULL;
		}
		static const PrimType* intType();
		static const PrimType* charType();
		static const PrimType* doubleType();
		static const PrimType* floatType();
		static const PrimType* voidType();
	};

	/****************************************/
	class PointerType : public Type{
	protected:
		const Type *srcType;
	public:
		PointerType(const Type *srcType) : srcType(srcType){
		}
		void setSrcType(const Type *type){
			srcType = type;
		}
		const Type *getSrcType()const{
			return srcType;
		}
		virtual void toString(std::string& str)const;
	};
	/****************************************/
	class ArrayType : public PointerType{
		std::string arraySize;
	public:
		ArrayType(const Type *srcType) : PointerType(srcType){
		}
		virtual void toString(std::string& str)const;
	};
	/****************************************/
	class Typedef : public Type{
		std::string name;
		const Type *srcType;
	public:
		Typedef(const std::string& name, Type *srcType)
			: srcType(srcType), name(name){}
		const Type *getSrc() const{
			return srcType;
		}
		std::string& getName(){
			return name;
		}
		const std::string& getName()const{
			return name;
		}
		const Type *getBaseType()const{
			const Typedef *tmp = dynamic_cast<const Typedef*>(srcType);
			if(tmp == NULL)
				return srcType;
			return tmp->getBaseType();
		}
		virtual void toString(std::string& str)const;
	};
	
	/****************************************/
	class FunctionType;
	class ParameterList : public Scope{
		friend class FunctionType;
	public:
		ParameterList(Scope *parent, ScopeAttribute *att) : Scope(parent, att){
		}
		const std::vector<Node*>& getParams(){
			return nodes;
		}
		virtual Type *findType(const std::string& name){
			return NULL;
		}
		virtual Node *findDeclare(const std::string& name){
			return NULL;
		}
		virtual bool add(Node *node){
			return false;
		}
		virtual void toString(std::string& str)const;
	};
	/****************************************/
	class FunctionType : public Type{
		Type *retType;
		ParameterList* params;
	public:
		FunctionType(Type *retType, ParameterList* params)
			: retType(retType), params(params){
		}
		Type *getRetType(){
			return retType;
		}
		ParameterList *getParamList(){
			return params;
		}
		virtual void toString(std::string& str)const;
	};

	/****************************************/
	class VariableDeclaration : public Declaration{
		bool isInitialized;
	public:
		VariableDeclaration(const std::string& name, 
			Type *type) : Declaration(name, type){
		}
	};

	/****************************************/
	class FunctionDeclaration : public Declaration{
	public:
		FunctionDeclaration(const std::string& name, 
			FunctionType *type) : Declaration(name, type){
		}
	};
	/****************************************/	
	class GlobalScope : public Scope{
		std::map<std::string, int> macros;
		std::map<std::string, int> declares;
		std::map<std::string, int> metaTypes;
	public:
		GlobalScope(Scope *parent, ScopeAttribute *att) : Scope(parent, att){}
		virtual Type *findType(const std::string& name);
		virtual Node *findDeclare(const std::string& name);
		Macro *findMacro(const std::string& macro);
		std::pair<Macro*,int> findMacroIndex(const std::string& macro);
		bool removeMacro(const std::string& macro);
		bool add(Node *node);
		std::vector<bool> blockActiveStack;
	};
	/****************************************/
	class CompoundTypeScope : public Scope{
	};
	/****************************************/
	class CompoundType : public Type{
		CompoundTypeScope *scope;
	public:
		CompoundType() : scope(NULL){}
		bool isDefined(){
			return scope != NULL;
		}
	};
	class BitsType : public Type{
		int bitLen;
	public:
	};
	class AnonymousCompoundType : public CompoundType{
	public:
		virtual void toString(std::string& str);
	};
	class NamedCompoundType : public CompoundType{
		std::string name;
	public:
		NamedCompoundType(const std::string& name) : name(name){
		}
		const std::string& getName(){
			return name;
		}
		virtual void toString(std::string& str);
	};

	/****************************************/
	struct Position{
		int row;
		int col;
		bool operator<(Position pos) const{
			if(row < pos.row)
				return true;
			if(row > pos.row)
				return false;
			return col < pos.col;
		}
		bool operator==(Position pos)const{
			return pos.row == row && pos.col == col;
		}
		Position(int row = 1, int col = 0) : row(row), col(col){
		}
	};
	class BlockIterator;
	class Block{
		friend class BlockIterator;
		const Byte *start;
		const Byte *end;
		Block *next;
		Block *prev;
		int flag;
		Position savedPos;
		Position startPos;
	public:
		enum {NOT_UPDATE_POS = 1, NOT_RESET_POS = 2, NOT_RECOVER_POS = 4, IS_MACRO_EXPANSION = 8};
		Block(const Byte *start, const Byte *end, int flag = 0,
			Block *next = NULL, Block *prev = NULL) : 
		start(start), end(end), next(next), prev(prev), flag(flag){
			assert(start < end);
		}
		void insert(Block *block){
			block->next = next;
			block->prev = this;
			if(block->next != NULL)
				block->next->prev = block;
			this->next = block;
		}
		bool isUpdatePos(){
			return (flag & NOT_UPDATE_POS) == 0;
		}
		bool isResetPos(){
			return (flag & NOT_RESET_POS) == 0;
		}
		bool isRecoverPos(){
			return (flag & NOT_RECOVER_POS) == 0;
		}
		bool isMacroExpansion(){
			return (flag & IS_MACRO_EXPANSION) != 0;
		}
		void setSavedPos(Position pos){
			savedPos = pos;
		}
		void setStartPosition(Position pos){
			startPos = pos;
		}
		Position getSavedPosition(){
			return savedPos;
		}
		Position getStartPosition(){
			return startPos;
		}
	};
	class BlockIterator{
		Position pos;
		Block *block;
		const Byte *current;
	public:
		BlockIterator(Block *block, const Byte *current, Position pos)
			: pos(pos), block(block), current(current){
		}
		BlockIterator(Block *block)
			: pos(1, 0), block(block), current(block->start){
		}
		int operator*() const{
			if(current >= block->end){
				assert(block->next == NULL);
				return -1;
			}
			return *current;
		}
		Block *getBlock() const{
			return block;
		}
		const Byte *getCurrent() const{
			return current;
		}
		Position getPosition() const{
			return pos;
		}
		BlockIterator& operator++(){
			++ current;
			if(current < block->end)
			{
				if(block->isUpdatePos())
				{
					if(*current == '\n'){
						++pos.row;
						pos.col = 0;
					}else
						++pos.col;
				}
			}
			else{
				if(block->next == NULL){
					std::cerr << "attempt to access a null block "
						<< "at line " << __LINE__ << ", file " << __FILE__ << std::endl;
					//assert("attempt to access a null block" == NULL);
				}
				else
					while(block->next != NULL && current >= block->end){
						if(block->isRecoverPos()){
							this->pos = block->getSavedPosition();
						}
						if(block->next->isRecoverPos()){
							block->next->setSavedPos(this->pos);
						}
						block = block->next;
						if(block->isResetPos()){
							Position tmp = pos;
							this->pos = block->getStartPosition();							
						}
						current = block->start;
					}
			}
			return *this;
		}
		void split(Block *block){
			if(current >= this->block->end){
				this->block->insert(block);
			}else if(current == this->block->start){
				if(this->block->prev != NULL){	
					this->block->prev->insert(block);
				}else{
					block->prev = NULL;
					block->next = this->block;
					this->block->prev = block;
				}
			}else{
				Block *newBlock = new Block(current, this->block->end, this->block->flag, this->block->next, block);
				this->block->next = block;
				block->next = newBlock;
				block->prev = this->block;
				current = block->start;		
				newBlock->setStartPosition(pos);		
			}			
			if(this->block->isRecoverPos())
				pos = this->block->getSavedPosition();
			if(block->isResetPos()){
				Position tmp = pos;
				pos = block->getSavedPosition();
				block->setSavedPos(tmp);
			}
			current = block->start;
			this->block = block;
		}
		BlockIterator operator++(int){
			BlockIterator tmp(*this);
			this->operator++();
			return tmp;
		}
	};
	class Utility{
	public:
		static std::string& ltrim(std::string& s){
			s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(::isspace))));
			return s;
		}
		static std::string& rtrim(std::string& s){
			s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(::isspace))).base(), s.end());
			return s;
		}
		static std::string& trim(std::string& s){
			return ltrim(rtrim(s));
		}
		static bool isOctDigit(Byte ch){
			return ch >= '0' && ch <= '7';
		}
		static bool isHexDigit(Byte ch){
			return (ch >= '0' && ch <= '7') 
				|| (ch >= 'a' && ch <= 'f') 
				|| (ch >= 'A' && ch <= 'F');
		}
		static int oct2dec(Byte ch){
			assert(isOctDigit(ch));
			return ch - '0';
		}
		static int hex2dec(Byte ch){
			assert(isHexDigit(ch));
			if(ch >= '0' && ch <= '9')
				return ch - '0';
			if(ch >= 'a' && ch <= 'f')
				return ch - 'a' + 10;
			return ch - 'A' + 10;
		}
		static char getMatchedBrace(char ch){
			switch(ch){
				case '(':return ')';
				case ')':return '(';
				case '[':return ']';
				case ']':return '[';
				case '{':return '}';
				case '}':return '{';
			}
			return ch;
		}
		static Type *getType(Node *node){
			Type *type = dynamic_cast<Type*>(node);
			if(type == NULL){
				Declaration * declare = dynamic_cast<Declaration *>(node);
				if(declare != NULL)
					type = declare->getType();
			}
			return type;
		}
		static bool isValidIdentifier(const char *name){
			if(*name != '_' && !::isalpha(*name))
				return false;
			++ name;
			while(*name){
				if(*name != '_' && !::isalnum(*name))
					return false;
				++ name;
			}
			return true;
		}
		static bool fileToString(const char *file, std::string& content){
			std::ifstream is(file);
			if(!is)
				return false;
			content.assign((std::istreambuf_iterator<char>(is)),
                 std::istreambuf_iterator<char>());
			return true;
		}
	};
	class Lexer{
	protected:
		Block *currentBlock;
		BlockIterator itor;
		std::string curToken;
		int intValue;
		double doubleValue;
		bool _isDoubleValue;
		void readCharSequence(int delimiter, std::string& ret);
		void readMultiLineComment(std::string& ret);
		void readIdentifier(std::string& ret);
		bool readNum(BlockIterator&, std::string& ret);
	public:
		enum{
			IDENTIFIER = 10000, CONSTANT = 20000, OPERATOR = 30000, NEWLINE = 40000, COMMENT = 50000,
			NUMBER = 20001, CHARACTER = 20002, STRING = 20003, BAD_TOKEN = 0, END = EOF
		};
		Lexer(Block *data) : currentBlock(data), itor(data), 
			intValue(0), doubleValue(0), _isDoubleValue(false){}
		int next(BlockIterator& itor, bool escape = true);
		int next(bool escape = true){
			return next(itor, escape);
		}
		int getToken();
		std::string& getTokenContent(){
			return curToken;
		}
		Block *getCurrentBlock(){
			return currentBlock;
		}
		bool isDoubleValue(){
			return _isDoubleValue;
		}
		double getDoubleValue(){
			return doubleValue;
		}
		int getIntValue(){
			return intValue;
		}
		int getLineNo(){
			return itor.getPosition().row;
		}
		Position getPosition(){
			return itor.getPosition();
		}
		BlockIterator& iterator(){
			return itor;
		}
		int peek();
		void skipSpace();
		void skipToEndOfLine(std::string& ret);
		void skipToEndOfLine(){
			skipToEndOfLine(curToken);
		}
		void toString(int token, std::string& str);
	};

	class HeaderParser{
		std::vector<std::pair<Position, Position> > nodeBlocks;
		std::map<std::string, std::string> fileContents;
		std::list<Block> blocks;
		Scope *currentScope;
		GlobalScope *globalScope;
		Block block;
		Lexer *lexer;
		bool isParsed;
		std::list<std::string> expandedMacros;
		HeaderParser& operator=(const HeaderParser&);
		HeaderParser(const HeaderParser&);
		bool expandMacro(Macro *macro);
		int getToken();
		int getToken(int token);
		int preprocessInstruction();
		int parseDeclares(Node *&innterMost, Node *&outterMost);
		int parsePostDeclarator(int token, ModifiedType *&type);
		int parseDeclaration();		
		int parseDefinition(Type **type);
		int parseVariableList(int token, ModifiedType *type);
		int parseMacro();
		int skipBraces(char ch);
		void includeFile(const char *fileName);
		void reportError(const char *msg, int exitCode = 0);
	public:
		HeaderParser(const Byte *start, const Byte *end);
		bool parse();
		~HeaderParser(){
			delete lexer;
			delete globalScope;
		}
	};
	/*
		Macro expansion is more complex than we think it is.
		Take the following code as example:
		#define DD 5+EE
		#define CC DD+4
		#define BB 3+CC+DD
		#define EE 2+AA+EE+CC
		#define AA CC+BB+1
		int a = AA;

		We could use a directed graph to represent the relationships between macros.
		Each vertex in the graph represents a macro and there is an edge from one macro to
		another, if the one references another. Direct and indirect self-references cause loops
		and circuits in the directed graph, which make the expansion of a macro more difficult.
		You must avoid infinite expansion.

		Take the above code as example. What would the line of code "int a = AA;"
		produce? Especially you have to take many different compilers into consideration.
		I have tested on MSVC++. 
		The expansion of "AA" in the expression "int a = AA;" should very likely be:
			"5 + 2 + AA + EE + CC + 4 + 3 + 5 + 2 + AA + EE + CC + 4 + 5 + 2 + AA + EE + DD + 4 + 1".

		I have tested this using the following code:
		int AA=1;
		int CC=2;
		int DD=3;
		int EE=4;
		int expaned = 
			5 + 2 + AA + EE + CC + 4 + 3 + 5 + 2 + AA + EE + CC + 4 + 5 + 2 + AA + EE + DD + 4 + 1;
		#define DD 5+EE
		#define CC DD+4
		#define BB 3+CC+DD
		#define EE 2+AA+EE+CC
		#define AA CC+BB+1
		
		int a = AA;
		The procedure I worked out is as follows:
                                                  EE--->2+AA+EE+CC                            CC--->DD+4
                                                         /        \                                  /       \
              EE--->2+AA+EE+CC                    DD--->5+EE--->5+2+AA+EE+CC           EE--->2+AA+EE+CC--->2+AA+EE+DD+4
                    /        \                          /                \                     /                \
           DD--->5+EE--->5+2+AA+EE+CC           CC--->DD+4 --->  5+2+AA+EE+CC+4  DD--->5+EE    --->    5+2+AA+EE+DD+4
                 /                \                     /                    \              /                        \
          CC--->DD+4 ---> 5+2+AA+EE+CC+4  BB--->3+CC+DD      --->      3+5+2+AA+EE+CC+4+DD     --->    3+5+2+AA+EE+CC+4+5+2+AA+EE+DD+4
                /                    \              /                                                                \
        AA---->CC+BB+1    --->      5+2+AA+EE+CC+4+BB+1                --->            5+2+AA+EE+CC+4+3+5+2+AA+EE+CC+4+5+2+AA+EE+DD+4+1

		If you output the values of the variables "expaned" and "a" you should see both of them are "59"
		under MSVC++.

		Another thing to notice is that, for some compilers such as MSVC++, the following macro:
		#define A() abd

		should be invoked as:
		A()

		not
		A

		But here the macro should be invoked as:
		A

		not 
		A()

	*/
	class MacroExpansion{
		GlobalScope *globalScope;
		Macro *macro;
		MacroExpansion* parent;
		Lexer *lexer;
		Lexer *paramLexer;
		int currentToken;
		std::vector<std::string> paramValues;
		bool isChanged;
		bool isExpanded(const std::string& str){
			if(macro->name == str)
				return true;
			if(parent == NULL)
				return false;
			return parent->isExpanded(str);
		}
		int findParam(const std::string& str){
			for(int i = 0;i < macro->params.size();++ i)
				if(str == macro->params[i])
					return i;
			return -1;
		}
		void currentTokenToString(std::string& str);
		bool getToken();
	public:
		MacroExpansion(Macro *macro, GlobalScope *globalScope, Lexer *paramLexer, MacroExpansion* parent = NULL)
			: macro(macro), globalScope(globalScope), paramLexer(paramLexer), parent(parent), lexer(NULL), currentToken(Lexer::END), isChanged(false){
		}
		bool changed(){
			return isChanged;
		}
		bool expandMacro(std::string& result);		
	};
	struct Expr{
		enum{UNARY =0, BITWISE = (1 << 16), ARITH = (2 << 16), 
			SHIFT = (3 << 16), COMP = (4 << 16), LOGICAL = (5 << 16),};
		enum{
			PLUS = (ARITH + '+'), MINUS = (ARITH + '-'),
			MULTI = (ARITH + '*'), DIV = (ARITH + '/'),
			MODULA = (ARITH + '%'),
		};
		enum{
			LSHIFT = (SHIFT + ('<' << 8) + '<'),
			RSHIFT = (SHIFT + ('>' << 8) + '>'),
		};
		enum{
			LESS = (COMP + ('<' << 8)), GREATER = (COMP + ('>' << 8)), 
			LESSEQ  = (LESS + '='), GREATEREQ  = (GREATER + '='),
			EQUAL = (SHIFT + ('=' << 8) + '='),
			NEQUAL = (SHIFT + ('!' << 8) + '='),
		};
		enum{
			AND = (LOGICAL + ('&' << 8) + '&'),
			OR = (LOGICAL + ('|' << 8) + '|'),
		};
		enum{
			BITAND = (BITWISE + '&'), BITOR = (BITWISE + '|'),
			BITXOR = (BITWISE + '^'),
		};
		/*LOGICAL < COMPARISION < SHIFT < ARITH < BITWISE*/
		static int isPlusOrMinus(int op1){
			if(op1 == PLUS || op1 == MINUS)
				return 1;
			return 0;
		}
		static bool precede(int op1, int op2){
			if(op1 == '?' || op1 == ':'){
				return false;
			}
			if(op2 == '?' || op2 == ':'){
				return true;
			}
			int tmp1 = (op1 & 0xFFFF0000);
			int tmp2 = (op2 & 0xFFFF0000);
			if(tmp1 != tmp2)
				return tmp1 < tmp2;
			if(tmp1 == ARITH){
				return isPlusOrMinus(op1) < isPlusOrMinus(op2);
			}
			return false;
		}
		static bool isComparison(int token){
			return (token & 0xFFFF0000) == Expr::COMP;
		}
		static bool isLogic(int token){
			return (token & 0xFFFF0000) == Expr::LOGICAL;
		}
		static bool isShift(int token){
			return (token & 0xFFFF0000) == Expr::SHIFT;
		}
		static bool isArithmatical(int token){
			return (token & 0xFFFF0000) == Expr::ARITH;
		}
		static bool isBitwise(int token){
			return (token & 0xFFFF0000) == Expr::BITWISE;
		}
		static void toString(int op, std::string& str){
			char tmp = (char)(op >> 8);
			if(tmp != 0)
				str.push_back(tmp);
			tmp = (char)op;
			str.push_back(tmp);
		}
		virtual ~Expr(){
		}
	};
	struct NumberExpr : public Expr{
		union{
			double doubleVal;
			int intVal;
		} number;
		bool isDouble;
	};
	struct FunctionExpr : public Expr{
		std::string name;
		std::vector<Expr*> exprs;
	};
	struct CastExpr : public Expr{
		Keyword type;
		Modifier modifer;
		Expr *expr;
	};
	struct UnaryExpr : public Expr{
		int optr;
		Expr *expr;
		UnaryExpr(int optr, Expr *expr) : optr(optr), expr(expr){
		}
	};
	struct IdentifierExpr : public Expr{
		std::string name;
	};
	struct BinaryExpr : public Expr{
		int optr;
		Expr *left;
		Expr *right;
		BinaryExpr(int optr = 0, Expr * l = NULL, Expr * r = NULL) 
			: optr(optr), left(l), right(r){
		}
		~BinaryExpr(){
			delete left;
			delete right;
		}
	private:
		BinaryExpr& operator=(const BinaryExpr&);
		BinaryExpr(const BinaryExpr&);
	};
	class ExprVisitor{
	public:
		virtual bool visit(Expr *expr){
			NumberExpr * num = dynamic_cast<NumberExpr*>(expr);
			if(num != NULL)
				return visit(num);
			FunctionExpr * func = dynamic_cast<FunctionExpr*>(expr);
			if(func != NULL)
				return visit(func);
			IdentifierExpr * ident = dynamic_cast<IdentifierExpr*>(expr);
			if(ident != NULL)
				return visit(ident);
			BinaryExpr * bin = dynamic_cast<BinaryExpr*>(expr);
			if(bin != NULL)
				return visit(bin);
			UnaryExpr *unary = dynamic_cast<UnaryExpr*>(expr);
			if(unary != NULL)
				return visit(unary);
			return false;
		}
		virtual bool visit(NumberExpr *expr) = 0;
		virtual bool visit(FunctionExpr *expr) = 0;
		virtual bool visit(IdentifierExpr *expr) = 0;
		virtual bool visit(BinaryExpr *expr) = 0;
		virtual bool visit(UnaryExpr *expr) = 0;
	};
	class ExprToString : public ExprVisitor{
		std::string str;
		int lastOp;
	public:
		ExprToString() : lastOp(0){}
		const std::string& getString(){
			return str;
		}
		virtual bool visit(Expr *expr){
			return ExprVisitor::visit(expr);
		}
		virtual bool visit(NumberExpr *expr){
			std::stringstream os;
			if(expr->isDouble)
				os << expr->number.doubleVal;
			else
				os << expr->number.intVal;
			str = os.str();
			return true;
		}
		virtual bool visit(FunctionExpr *expr){
			str = expr->name;
			str.push_back('(');
			//if(expr->expr != NULL){
				for(size_t i = 0;i < expr->exprs.size();++ i){
					ExprToString tmp;
					tmp.visit(expr->exprs[i]);
					if(i != 0)
						str.push_back(',');
					str += tmp.str;
				}
			//}
			str.push_back(')');
			return true;			
		}
		virtual bool visit(IdentifierExpr *expr){
			str = expr->name;
			return true;
		}
		virtual bool visit(BinaryExpr *expr){
			ExprToString left;
			ExprToString right;
			if(expr->left != NULL)
				left.visit(expr->left);
			if(expr->right != NULL)
				right.visit(expr->right);
			if(Expr::precede(expr->optr, left.lastOp)){
				str.push_back('(');
				str += left.str;
				str.push_back(')');
			}else
				str += left.str;
			Expr::toString(expr->optr, str);
			if(Expr::precede(expr->optr, right.lastOp)){
				str.push_back('(');
				str += right.str;
				str.push_back(')');
			}else
				str += right.str;
			lastOp = expr->optr;
			return true;
		}
		virtual bool visit(UnaryExpr *expr){
			Expr::toString(expr->optr, str);
			if(expr->expr != NULL){
				ExprToString tmp;
				tmp.visit(expr->expr);
				if(Expr::precede(expr->optr, tmp.lastOp)){
					str.push_back('(');
					str += tmp.str;
					str.push_back(')');
				}else
					str += tmp.str;
			}	
			lastOp = expr->optr;
			return true;	
		}
	};
	/*
	  primary :
			NUMBER
			IDDENTIFIER
			UNARYOP(+,-,~!)	primary
			(expression)
			IDDENTIFIER(expression,expression,...)		Note: this resembles function invokation that only take one parmeter
			(TYPENAME)primary		Note: This is not supported
	  bitwise :
			bitwise OPERATOR(&,|,^) primary
			primary
	  arithmetical :
			arithmetical OPERATOR(+,-,*,/,%) bitwise
			bitwise
	  shift :
			shift OPERATOR(<<,>>) arithmetical
			arithmetical
	  comparision :
			comparision OPERATOR(>,<,>=,<=,==,!=) shift
			shift
	  logical :
			logical OPERATOR(&&,||) comparison
			comparison ? logical : logical
			comparision
	  expression :
			logical

	*/
	class ExprParser{
		int currentToken;
		Block block;
		Lexer *lexer;
		ExprParser& operator=(const ExprParser&);
		ExprParser(const ExprParser&);
		int getToken0();
		void getToken(){
			currentToken = getToken0();
		}
		enum{IDENTIFIER = Lexer::IDENTIFIER, NUMBER = Lexer::NUMBER, INVALID=0, END=Lexer::END};
		Expr *primary(bool get);
		Expr *bitwise(bool get);
		Expr *arithmatical(bool get);
		Expr *arithmaticalTerm(bool get);
		Expr *shift(bool get);
		Expr *comparison(bool get);
		Expr *logical(bool get);
	public:
		ExprParser(const Byte *start, const Byte *end);
		Expr *parse();
		~ExprParser(){
			delete lexer;
		}
	};

}
