#include <vector>
#include <cctype>
#include <iostream>
#include <unordered_map>
#include <cassert>
#include <list>
#include <set>
using namespace std;

template<class T>
class Trait {
public:
	static bool isTerminal(const T& val);
	static bool isNonTerminal(const T& val);
	static bool isStart(const T& val);
	static bool isValid(const T& val);
};
enum ProductionType {
	INVALID_PRODUCTION,
	TYPE0_PRODUCTION,
	TYPE1_PRODUCTION,
	TYPE2_PRODUCTION,
	TYPE3_PRODUCTION,
};
ostream& operator<<(ostream& os, ProductionType type) {
	switch (type) {
	case INVALID_PRODUCTION:
		os << "invalid";
		return os;
	case TYPE0_PRODUCTION:
		os << "type 0";
		return os;
	case TYPE1_PRODUCTION:
		os << "type 1";
		return os;
	case TYPE2_PRODUCTION:
		os << "type 2";
		return os;
	case TYPE3_PRODUCTION:
		os << "type 3";
		return os;
	}
	return os;
}
bool operator<(ProductionType type1, ProductionType type2) {
	return (int) type1 > (int) type2;
}
bool operator>(ProductionType type1, ProductionType type2) {
	return (int) type1 < (int) type2;
}
template<class T>
struct Production {
	vector<T> l;
	vector<T> r;
public:
	Production() {

	}
	Production(const vector<T>& l, const vector<T>& r) :
			l(l), r(r) {
	}
	const vector<T>& left() const {
		return l;
	}
	const vector<T>& right() const {
		return r;
	}
	bool isValid() const {
		for (size_t i = 0; i < l.size(); ++i)
			if (Trait<T>::isNonTerminal(l[i]))
				return true;
		return false;
	}
	bool isStart() const {
		return l.size() == 1 && Trait<T>::isStart(l[0]);
	}
	ProductionType getType() const {
		if (l.size() == 1) {
			if (Trait<T>::isStart(l[0])) {
				if (r.empty()
						|| (r.size() == 2 && Trait<T>::isTerminal(r[0])
								&& Trait<T>::isNonTerminal(r[1])))
					return TYPE3_PRODUCTION;
				return TYPE2_PRODUCTION;
			} else if (Trait<T>::isNonTerminal(l[0])) {
				return TYPE2_PRODUCTION;
			}
			return INVALID_PRODUCTION;
		}
		if (l.empty())
			return INVALID_PRODUCTION;
		if (r.empty() || r.size() < l.size() - 1) {
			return isValid() ? TYPE0_PRODUCTION : INVALID_PRODUCTION;
		}
		//production of the form, lAr--->lr
		if (r.size() == l.size() - 1) {
			size_t i = 0;
			size_t j = 0;
			while (j < r.size() && l[i] == r[j]) {
				++i;
				++j;
			}
			if (j == r.size())
				return TYPE1_PRODUCTION;
			if (!Trait<T>::isNonTerminal(l[i]))
				return isValid() ? TYPE0_PRODUCTION : INVALID_PRODUCTION;
			++i;
			while (j < r.size() && l[i] == r[j]) {
				++i;
				++j;
			}
			return j == r.size() ? TYPE1_PRODUCTION : TYPE0_PRODUCTION;
		}
		size_t startL = 0, endL = l.size() - 1;
		size_t startR = 0, endR = r.size() - 1;
		while (startL < endL) {
			while (startL < endL && Trait<T>::isTerminal(l[startL])
					&& l[startL] == r[startR]) {
				++startL;
				++startR;
			}
			if (Trait<T>::isTerminal(l[startL]))
				return isValid() ? TYPE0_PRODUCTION : INVALID_PRODUCTION;
			size_t tmpL = endL;
			size_t tmpR = endR;
			while (tmpL > startL && l[tmpL] == r[tmpR]) {
				--tmpL;
				--tmpR;
			}
			if (tmpL == startL)
				return TYPE1_PRODUCTION;
			if (l[startL] != r[startR])
				return TYPE0_PRODUCTION;
			++startL;
			++startR;
		}
		if (Trait<T>::isNonTerminal(l[startL]))
			return TYPE1_PRODUCTION;
		return TYPE0_PRODUCTION;
	}
	friend ostream& operator<<(ostream& os, const Production& p) {
		for (size_t i = 0; i < p.l.size(); ++i)
			os << p.l[i];
		os << "\t--->\t";
		for (size_t i = 0; i < p.r.size(); ++i)
			os << p.r[i];
		return os;
	}
	bool operator==(const Production& p) {
		return l == p.l && r == p.r;
	}
	void swap(Production& p) {
		l.swap(p.l);
		r.swap(p.r);
	}
};
class Symbol {
	char ch;
public:
	Symbol(char ch) {
		this->ch = 0;
		if (::isalnum(ch))
			this->ch = ch;
	}
	operator char() const {
		return ch;
	}
	bool operator==(Symbol val) const {
		return ch == val.ch;
	}
	friend ostream& operator<<(ostream& os, Symbol val) {
		os << (char) val;
		return os;
	}
	struct HashKey {
		int operator()(Symbol symbol) const {
			return (char) symbol;
		}
	};
};
class Digit {
	char ch;
public:
	Digit(char ch) {
		this->ch = 0;
		if (::isalnum(ch) || ch == '-' || ch == '+')
			this->ch = ch;
	}
	operator char() const {
		return ch;
	}
	bool operator==(Digit val) const {
		return ch == val.ch;
	}
	friend ostream& operator<<(ostream& os, Digit val) {
		os << (char) val;
		return os;
	}
	struct HashKey {
		int operator()(Digit symbol) const {
			return (char) symbol;
		}
	};
};
template<>
class Trait<Digit> {
public:
	static bool isTerminal(const Digit& val) {
		return ::isdigit((char) val) || (char) (val) == '-'
				|| (char) (val) == '+';
	}
	static bool isNonTerminal(const Digit& val) {
		return ::isalpha((char) val);
	}
	static bool isStart(const Digit& val) {
		return (char) val == 'S';
	}
	static bool isValid(const Digit& val) {
		return val != '\0';
	}
};
template<>
class Trait<Symbol> {
public:
	static bool isTerminal(const Symbol& val) {
		return ::islower((char) val) || ::isdigit((char) val);
	}
	static bool isNonTerminal(const Symbol& val) {
		return ::isupper((char) val);
	}
	static bool isStart(const Symbol& val) {
		return (char) val == 'S';
	}
	static bool isValid(const Symbol& val) {
		return val != '\0';
	}
};
template<class T>
struct GeneratedString {
	unsigned int steps;
	vector<T> val;
	GeneratedString() :
			steps(0) {

	}
	GeneratedString(unsigned int steps, const string& val) :
			steps(steps), val(val) {

	}
};
template<class T>
struct StringArrayComparator {
	const vector<GeneratedString<T> >& generated;
	StringArrayComparator(const vector<GeneratedString<T> >& generated) :
			generated(generated) {
	}
	bool operator()(unsigned int i, unsigned int j) const {
		return generated[i].val < generated[j].val;
	}
};
template<class T, class Traits = Trait<T> >
class PhraseStructureGrammer {
	vector<Production<T> > productions;
public:
	bool add(const T *l, const T *r) {
		vector<T> left;
		vector<T> right;
		while (Traits::isValid(*l)) {
			left.push_back(T(*l));
			++l;
		}
		while (Traits::isValid(*r)) {
			right.push_back(T(*r));
			++r;
		}

		Production<T> p(left, right);
		if (p.isValid()) {
			productions.push_back(Production<T>());
			productions.back().swap(p);
			return true;
		}
		return false;
	}
	void generateString(unsigned int maxStep,
			vector<GeneratedString<T> >& generated) const {
		if (generated.empty())
			for (size_t i = 0; i < productions.size(); ++i) {
				Production<T> p = productions[i];
				if (p.isStart()) {
					generated.push_back(GeneratedString<T>());

					GeneratedString<T>& str = generated.back();
					str.steps = 1;
					str.val.assign(p.right().begin(), p.right().end());
				}
			}
		std::unordered_map<T, vector<pair<unsigned int, unsigned int> >,
				typename T::HashKey> productionMap(256);
		for (size_t i = 0; i < productions.size(); ++i) {
			const vector<T>& p = productions[i].left();
			for (size_t j = 0; j < p.size(); ++j)
				if (Traits::isNonTerminal(p[j])) {
					productionMap[p[j]].push_back(
							make_pair(static_cast<unsigned int>(i),
									static_cast<unsigned int>(j)));
				}
		}

		StringArrayComparator<T> comparator = StringArrayComparator<T>(
				generated);
		std::set<unsigned int, StringArrayComparator<T> > insertedStrings(
				comparator);
		for (size_t i = 0; i < generated.size(); ++i) {
			const GeneratedString<T>& str = generated[i];
			unsigned int newStep = str.steps + 1;
			if (str.steps < maxStep) {
				std::set<pair<pair<unsigned int, unsigned int>, unsigned int> > usedRules;
				const vector<T> val = str.val;
				for (size_t j = 0; j < val.size(); ++j)
					if (Traits::isNonTerminal(val[j])) {
						vector<pair<unsigned int, unsigned int> >& useful =
								productionMap[val[j]];
						for (size_t k = 0; k < useful.size(); ++k) {
							pair<unsigned int, unsigned int>& one = useful[k];
							const vector<T>& p = productions[one.first].left();
							if (one.second <= j) {
								unsigned int start = j - one.second;
								if (start + p.size() <= val.size()) {
									pair<pair<unsigned int, unsigned int>,
											unsigned int> tmpOne;
									tmpOne.first.first = start;
									tmpOne.first.second = start + p.size();
									tmpOne.second = one.first;
									if (usedRules.find(tmpOne)
											== usedRules.end()) {
										usedRules.insert(tmpOne);
										unsigned int s = 0;
										while (s < p.size()
												&& val[s + start] == p[s])
											++s;
										if (s >= p.size()) {
											unsigned int lastIndex =
													generated.size();
											generated.push_back(
													GeneratedString<T>());

											GeneratedString<T>& last =
													generated.back();
											const vector<T>& r =
													productions[one.first].right();
											last.steps = newStep;
											for (unsigned int t = 0; t < start;
													++t)
												last.val.push_back(val[t]);
											for (unsigned int t = 0;
													t < r.size(); ++t)
												last.val.push_back((char) r[t]);
											for (unsigned int t = start
													+ p.size(); t < val.size();
													++t)
												last.val.push_back(val[t]);
											if (insertedStrings.find(lastIndex)
													!= insertedStrings.end()) {
												generated.pop_back();
											} else {
												insertedStrings.insert(
														lastIndex);
												//cout << newStep << ':' << last.val << '\t' << val << '\t' << one.first << ',' << one.second << endl;
											}
										}
									}
								}
							}
						}
					}
			}
		}
	}
	ProductionType getType() const {
		if (productions.empty())
			return INVALID_PRODUCTION;

		ProductionType type = productions[0].getType();
		for (size_t i = 1; i < productions.size(); ++i) {
			ProductionType t = productions[i].getType();
			if (t > type)
				type = t;
		}
		return type;
	}
	friend ostream& operator<<(ostream& os,
			const PhraseStructureGrammer& grammer) {

		for (size_t i = 0; i < grammer.productions.size(); ++i) {
			os << grammer.productions[i] << std::endl;
		}
		return os;
	}
};
int test101(int argc, char *argv[]) {
	PhraseStructureGrammer<Symbol> grammer;
	const char *rules[] = {
	/*"S","DCBA",
	 "B","0B",
	 "B", "1B",
	 "B", "0",
	 "B", "1",
	 "CA", "E",
	 "DE", "",

	 "1A", "FI1",
	 "0A", "GI0",
	 "0F", "F0",
	 "0G", "G0",
	 "1F", "F1",
	 "1G", "G1",
	 "CF", "FC",
	 "CG", "GC",
	 "DF", "D1H",
	 "DG", "D0H",
	 "H0", "0H",
	 "H1", "1H",
	 "HC", "CH",
	 "HI", "A",

	 "0E", "E0",
	 "1E", "E1",*/
	/*"S", "SS",
	 "S", "100",
	 "S", "010",
	 "S", "001",
	 "S", "",
	 "1S", "S1",
	 "0S", "S0",*/
	/*"S", "B0AE",
	 "AE", "C",
	 "BC", "",
	 "0C", "C0",
	 "AE", "DE",
	 "0D", "D00",
	 "BD", "BF",
	 "F0", "0F",
	 "FE", "AE",*/
	/*"S", "aBc",
	 "aBc", "adBfc",
	 "eBf", "eaBf",
	 "B", "BB",
	 "B", "f",
	 "B", "a",
	 "B", "c",*/
	/*"S", "0S1",
	 "S", "",*/
	"S", "0S", "S", "1S", "S", "",

	};
	bool ret = true;
	unsigned int count = (sizeof(rules) / sizeof(*rules)) / 2;
	for (unsigned int i = 0; i < count; ++i) {
		ret = grammer.add((Symbol*) rules[2 * i], (Symbol*) rules[2 * i + 1]);
		assert(ret);
	}
	cout << grammer << std::endl;

	vector<GeneratedString<Symbol> > generated;
	//generated.push_back(GeneratedString(1, "DC00110111011A"));
	grammer.generateString(10, generated);
	cout << grammer.getType() << endl;
	for (size_t i = 0; i < generated.size(); ++i) {
		size_t j = 0;
		for (;
				j < generated[i].val.size()
						&& Trait<Symbol>::isTerminal(generated[i].val[j]); ++j)
			;

		if (j >= generated[i].val.size()) {
			for (j = 0; j < generated[i].val.size(); ++j)
				cout << generated[i].val[j];
			cout << endl;
		}
	}
	return 0;
}

int test102(int argc, char *argv[]) {
	PhraseStructureGrammer<Digit> grammer;
	const char *rules[] = { "S", "AI",/**/
	"A", "+", /**/"A", "-", /**/"I", "D",
	/**/"I", "DI", /**/"D", "0", /**/"D", "1", /**/"D", "2",
	/**/"D", "3", /**/"D", "4", /**/"D", "5",
	/**/"D", "6", /**/"D", "7", /**/"D", "8",
	/**/"D", "9",};
	bool ret = true;
	unsigned int count = (sizeof(rules) / sizeof(*rules)) / 2;
	for (unsigned int i = 0; i < count; ++i) {
		ret = grammer.add((Digit*) rules[2 * i], (Digit*) rules[2 * i + 1]);
		assert(ret);
	}
	cout << grammer << std::endl;

	vector<GeneratedString<Digit> > generated;
	grammer.generateString(20, generated);
	cout << grammer.getType() << endl;
	for (size_t i = 0; i < generated.size(); ++i) {
		size_t j = 0;
		for (;
				j < generated[i].val.size()
						&& Trait<Digit>::isTerminal(generated[i].val[j]); ++j)
			;

		if (j >= generated[i].val.size()) {
			for (j = 0; j < generated[i].val.size(); ++j)
				cout << generated[i].val[j];
			cout << endl;
		}
	}
	return 0;
}

