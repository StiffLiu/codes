#include <cstdlib>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <vector>
#include <ctime>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <iterator>
#include <stack>
using namespace std;
#define ARSIZE(ar) ((sizeof (ar)) / (sizeof (*ar)))

/*
 * generate a random string with a given length
 */
string randString(int num) {
	string tmp;
	for (int i = 0; i < num; ++i) {
		char ch = 'a';
		if (rand() % 2 == 0)
			ch = 'A';
		ch = +rand() % 26;
		tmp.push_back(ch);
	}
	return tmp;
}

/*
 * Compare the performance of hash map(unordered_map) and
 * sorted map(map).
 * After running serveral test, I found that the performance of hash map
 * is several times better.
 * But in cases when we want the elements to be ordered, we may prefer a sorted map.
 *
 * Note that hash map become part of standard c++ only after c++11.
 * You may need to set some compiler flags, so that this program could be be compiled.
 * You can add the following compiler options under g++:
 * 		-std=gnu++0x
 *
 */
int compareMaps(int argc, char *argv[]) {
	const int numOfIntegers = 10000;
	typedef int Item;
	std::unordered_map<Item, int> hashMap(1);
	std::map<Item, int> binMap;
	vector<Item> numbers;

	srand(time(0));
	for (int i = 0; i < numOfIntegers; ++i)
		numbers.push_back(rand());//randString(10));

	const int iteration = 1000;
	clock_t start = clock();
	for (int i = 0; i < iteration; ++i)
		for (int j = 0; j < numOfIntegers; ++j)
			hashMap[numbers[j]] = j;
	clock_t end = clock();
	cout << "bucket count " << hashMap.bucket_count() << ", size "
			<< hashMap.size() << ", hash map time : " << (end - start) << endl;
	start = clock();
	for (int i = 0; i < iteration; ++i)
		for (int j = 0; j < numOfIntegers; ++j)
			binMap[numbers[j]] = j;
	end = clock();
	cout << "size " << binMap.size() << ", bin map time : " << (end - start)
			<< endl;
	return 0;
}

/*
 * I was considering to have some special hash key generator for
 * an integer. But I didn't came out some ideas, except using the
 * value of an integer directly.
 */
struct HashKeyCalculator {
	int operator()(const int& val) const {
		return val;
	}
};
template<class State, class Input, class HashKeyCalculator, class Output = Input>
class MooreMachine {
	struct HashKey {
		HashKeyCalculator calculator;
		HashKey(const HashKeyCalculator& calculator) :
				calculator(calculator) {
		}
		int operator()(const pair<State, Input>& one) const {
			return calculator(one.first) * 31 + calculator(one.second);
		}
	};
	typedef std::unordered_map<pair<State, Input>, pair<State, Output>, HashKey> Functions;
	Functions functions;
	HashKeyCalculator calculator;
public:
	MooreMachine(const HashKeyCalculator& calculator = HashKeyCalculator()) :
			functions(1, HashKey(calculator)) {

	}
	bool addFunction(const State& state, const Input& input,
			const State& nextState, const Output& output) {
		pair<State, Input> one = make_pair(state, input);
		typename Functions::iterator pos = functions.find(one);
		if (pos == functions.end()) {
			pair<State, Output>& val = functions[one];
			val.first = nextState;
			val.second = output;
			return true;
		}
		return false;
	}
	bool getFunction(const State& state, const Input& input, State& nextState,
			Output& output) const {
		pair<State, Input> one = make_pair(state, input);
		typename Functions::const_iterator pos = functions.find(one);
		if (pos == functions.end()) {
			return false;
		}
		nextState = pos->second.first;
		output = pos->second.second;
		return true;
	}
	template<class Seq, class ForwardIterator>
	bool getOutput(const State& initialState, ForwardIterator& begin,
			ForwardIterator end, State& endState, Seq& output) const {
		endState = initialState;
		while (begin != end) {
			typename Functions::const_iterator pos = functions.find(
					make_pair(endState, *begin));
			if (pos == functions.end())
				return false;
			endState = pos->second.first;
			output.push_back(pos->second.second);
			++begin;
		}
		return true;
	}
};
template<class State, class Input, class HashKeyCalculator, class Output = Input>
class MealyMachine {
	struct HashKey {
		HashKeyCalculator calculator;
		HashKey(const HashKeyCalculator& calculator) :
				calculator(calculator) {
		}
		int operator()(const pair<State, Input>& one) const {
			return calculator(one.first) * 31 + calculator(one.second);
		}
	};
	typedef std::unordered_map<pair<State, Input>, State, HashKey> Functions;
	typedef std::unordered_map<State, Output, HashKeyCalculator> Outputs;
	HashKeyCalculator calculator;
	Functions functions;
	Outputs outputs;
public:
	MealyMachine(const HashKeyCalculator& calculator = HashKeyCalculator()) :
			functions(1, HashKey(calculator)), outputs(1, calculator) {

	}
	bool addFunction(const State& state, const Input& input,
			const State& nextState) {
		pair<State, Input> one = make_pair(state, input);
		typename Functions::iterator pos = functions.find(one);
		if (pos == functions.end()) {
			functions[one] = nextState;
			return true;
		}
		return false;
	}
	bool getFunction(const State& state, const Input& input,
			State& nextState) const {
		pair<State, Input> one = make_pair(state, input);
		typename Functions::const_iterator pos = functions.find(one);
		if (pos == functions.end()) {
			return false;
		}
		nextState = pos->second;
		return true;
	}
	bool addOutput(const State& state, const Output& output) {
		typename Outputs::iterator pos = outputs.find(state);
		if (pos == outputs.end()) {
			outputs[state] = output;
			return true;
		}
		return false;
	}
	bool getOutput(const State& state, Output& output) const {
		typename Outputs::const_iterator pos = outputs.find(state);
		if (pos == outputs.end()) {
			return false;
		}
		output = pos->second;
		return true;
	}
	template<class Seq, class ForwardIterator>
	bool getOutput(const State& initialState, ForwardIterator& begin,
			ForwardIterator end, State& endState, Seq& output) const {
		endState = initialState;
		while (begin != end) {
			typename Functions::const_iterator pos = functions.find(
					make_pair(endState, *begin));
			typename Outputs::const_iterator pos1 = outputs.find(endState);
			if (pos == functions.end() || pos1 == outputs.end())
				return false;
			endState = pos->second;
			output.push_back(pos1->second);
			++begin;
		}
		return true;
	}
};

/*
 * This class represent a finite state automaton.
 */
template<class State, class Input, class HashKeyCalculator>
class FAutomaton {
public:
	typedef std::unordered_set<State, HashKeyCalculator> StateSet;
	/*
	 * Set the start state
	 */
	void setAsStart(const State& state) {
		startState = state;
	}

	/*
	 * Set the final state
	 */
	void addFinal(const State& state) {
		finalStates.insert(state);
	}

	/*
	 * Get the final states
	 */
	const StateSet& getFinals() const {
		return finalStates;
	}
protected:

	/*
	 * a hash key generator for a pair
	 */
	struct HashKey {
		HashKeyCalculator calculator;
		/*HashKey(const HashKeyCalculator& calculator) :
		 calculator(calculator) {
		 }*/
		int operator()(const pair<State, Input>& one) const {
			return calculator(int(one.first)) * 31 + calculator(one.second);
		}
	};
	StateSet finalStates;
	State startState;

	/*
	 * output the finite state machine to stream.
	 */
	friend ostream& operator<<(ostream& os, const FAutomaton& fa) {
		os << "start state : " << fa.startState << endl;
		typename StateSet::const_iterator start =
				fa.finalStates.begin(), term =
				fa.finalStates.end();
		os << "final states : ";

		bool isFirst = true;
		while (start != term) {
			if (!isFirst)
				os << ',';
			os << *start;
			isFirst = false;
			++start;
		}
		return os;
	}
};

/*
 * A deterministic finite state automaton.
 */
template<class State, class Input, class HashKeyCalculator>
class DFAutomaton: public FAutomaton<State, Input, HashKeyCalculator> {
protected:
	typedef FAutomaton<State, Input, HashKeyCalculator> SuperClass;
	typedef std::unordered_map<pair<State, Input>, State,
			typename SuperClass::HashKey> Transitions;
	Transitions transitions;
public:
	void setTransition(const State& state, const Input& input,
			const State& next) {
		transitions[make_pair(state, input)] = next;
	}
	const State* getTransition(const State& state, const Input& input) const {
		typename Transitions::const_iterator pos = transitions.find(state);
		if (pos == transitions.end())
			return NULL;
		return &pos->second;
	}
	template<class ForwardIterator>
	bool recognize(ForwardIterator& begin, ForwardIterator end) const {
		State currentState = SuperClass::startState;
		while (begin != end) {
			typename Transitions::const_iterator pos = transitions.find(
					make_pair(currentState, *begin));
			if (pos == transitions.end())
				return false;
			currentState = pos->second;
			++begin;
		}
		return SuperClass::finalStates.find(currentState)
				!= SuperClass::finalStates.end();
	}
	friend ostream& operator<<(ostream& os, const DFAutomaton& dfa) {
		os << (const SuperClass&)dfa;
		os << endl;

		typename Transitions::const_iterator begin = dfa.transitions.begin(),
				end = dfa.transitions.end();
		while (begin != end) {
			os << begin->first.first << '-' << begin->first.second << ':'
					<< begin->second << endl;
			++begin;
		}
		return os;
	}
};
template<class T1, class T2, class T3>
struct Tuple {
	T1 first;
	T2 second;
	T3 third;
	Tuple(const T1& t1, const T2& t2, const T3& t3) :
			first(t1), second(t2), third(t3) {
	}
};
template<class State, class Input, class HashKeyCalculator>
class NDFAutomaton: public FAutomaton<State, Input, HashKeyCalculator> {
protected:
	typedef FAutomaton<State, Input, HashKeyCalculator> SuperClass;
	typedef std::unordered_map<Input, typename SuperClass::StateSet,
			HashKeyCalculator> InputStateSet;
	typedef std::unordered_map<State, InputStateSet, HashKeyCalculator> Transitions;
	Transitions transitions;
	const InputStateSet *getInputStateSet(const State& state) const {
		typename Transitions::const_iterator pos = transitions.find(state);
		if (pos == transitions.end())
			return NULL;
		return &pos->second;
	}
public:
	typedef typename SuperClass::StateSet StateSet;
	void setTransition(const State& state, const Input& input,
			const State& next) {
		transitions[state][input].insert(next);
	}
	const typename SuperClass::StateSet* getTransition(const State& state,
			const Input& input) const {
		typename Transitions::const_iterator pos = transitions.find(state);
		if (pos == transitions.end())
			return NULL;
		typename InputStateSet::const_iterator pos1 = pos->second.find(input);
		if (pos1 == pos->second.end())
			return NULL;
		return &pos1->second;
	}
	template<class ForwardIterator>
	bool recognize(ForwardIterator& begin, ForwardIterator end) {
		return recognize(begin, end, (vector<State>*) NULL);
	}
	template<class ForwardIterator, class Seq>
	bool recognize(ForwardIterator& begin, ForwardIterator end,
			Seq * seq) const {
		if (begin == end)
			return SuperClass::finalStates.find(SuperClass::startState)
					!= SuperClass::finalStates.end();
		typedef typename StateSet::const_iterator SetItor;
		typedef typename Transitions::const_iterator TransItor;
		typedef Tuple<SetItor, SetItor, ForwardIterator> Element;
		stack<Element> s;
		const StateSet *next = getTransition(SuperClass::startState, *begin);
		if (next == NULL)
			return false;

		s.push(Element(next->begin(), next->end(), begin));
		if (seq != NULL)
			seq->push_back(SuperClass::startState);
		while (!s.empty()) {
			Element& tuple = s.top();
			begin = tuple.third;
			if (tuple.first == tuple.second) {
				s.pop();
				if (seq != NULL)
					seq->pop_back();
			} else {
				State state = *tuple.first;
				++tuple.first;
				++begin;
				if (begin == end) {
					if (SuperClass::finalStates.find(state)
							!= SuperClass::finalStates.end())
						return true;
				} else {
					next = getTransition(state, *begin);
					if (next != NULL) {
						s.push(Element(next->begin(), next->end(), begin));
						if (seq != NULL) {
							seq->push_back(state);
						}
					}
				}
			}
		}
		return begin == end;
	}
	static ostream& output(ostream& os, const State& state, const InputStateSet& set){
		typename InputStateSet::const_iterator start = set.begin(), end = set.end();
		while(start != end){
			typename StateSet::const_iterator s = start->second.begin(), t = start->second.end();
			os << state << '-' << start->first << ':';

			bool isFirst = true;
			while(s != t){
				if(!isFirst)
					os << ',';
				os << *s;
				isFirst = false;
				++s;
			}
			os << endl;
			++start;
		}
		return os;
	}
	friend ostream& operator<<(ostream& os, const NDFAutomaton& dfa) {
			os << (const SuperClass&)dfa;
			os << endl;

			typename Transitions::const_iterator begin = dfa.transitions.begin(),
					end = dfa.transitions.end();
			while (begin != end) {
				output(os, begin->first, begin->second);
				++begin;
			}
			return os;
		}
	class CDFAutomaton;
	struct CompoundState {
	private:
		CDFAutomaton *automaton;
		unsigned int index;
	public:
		CompoundState() :
				automaton(0), index(0) {
		}
		CompoundState(CDFAutomaton *automaton, unsigned int index) :
				automaton(automaton), index(index) {
		}
		void set(CDFAutomaton *automaton, unsigned int index) {
			this->automaton = automaton;
			this->index = index;
		}
		bool operator==(const CompoundState& another) const {
			return index == another.index;
		}
		operator int() const {
			return index;
		}
		ostream& output(ostream &os) const;
		friend ostream& operator<<(ostream& os, const CompoundState& state) {
			return state.output(os);
		}
	};
	class CDFAutomaton: public DFAutomaton<CompoundState, Input,
			HashKeyCalculator> {
		friend class NDFAutomaton<State, Input, HashKeyCalculator> ;
		friend class CompoundHashKey;
	protected:
		vector<StateSet> states;
	};
	void toDFAutomaton(CDFAutomaton& result) const {
		struct TempState {
			const vector<StateSet> *states;
			int hashKey;
			size_t index;
			bool isFinal;
			TempState() :
					states(NULL), index(0), hashKey(0), isFinal(false) {
			}
			TempState(const vector<StateSet>* states, size_t index, int hashKey,
					bool isFinal) :
					states(states), hashKey(hashKey), index(index), isFinal(
							isFinal) {
			}
			bool operator==(const TempState& another) const {
				if (another.states != states)
					return false;
				if (states == NULL)
					return true;
				return (*states)[index] == (*another.states)[another.index];
			}
			static size_t hash(const TempState& state) {
				return state.hashKey;
			}
		};

		typedef size_t (*HashFunc)(const TempState& state);
		HashKeyCalculator calculator;
		vector<StateSet>& states = result.states;
		std::unordered_set<TempState, HashFunc> processedStates(1,
				TempState::hash);
		states.push_back(StateSet());
		states.back().insert(SuperClass::startState);

		const StateSet& finalStates = SuperClass::finalStates;
		processedStates.insert(
				TempState(&states, 0, calculator(SuperClass::startState),
						finalStates.find(SuperClass::startState)
								!= finalStates.end()));
		if(finalStates.find(SuperClass::startState)
								!= finalStates.end())
			result.addFinal(CompoundState(&result, 0));
		result.setAsStart(CompoundState(&result, 0));
		for (size_t i = 0; i < states.size(); ++i) {
			StateSet& one = states[i];
			InputStateSet next;
			typename StateSet::const_iterator begin = one.begin(), end =
					one.end();

			while (begin != end) {
				const InputStateSet *stateSet = getInputStateSet(*begin);
				if (stateSet != NULL) {
					typename InputStateSet::const_iterator start =
							stateSet->begin(), term = stateSet->end();
					while (start != term) {
						next[start->first].insert(start->second.begin(),
								start->second.end());
						++start;
					}
				}
				++begin;
			}
			typename InputStateSet::iterator start = next.begin(), term =
					next.end();
			states.push_back(StateSet());
			while (start != term) {
				/*
				 * Since element in an unordered_set could be accessed in random order.
				 * We should guarantee that two unordered_sets contains the same elements must have
				 * the the same hash key.
				 * Thus a simple way is used here
				 */
				int hashKey = 0;
				bool isFinal = false;
				typename StateSet::const_iterator s = start->second.begin(), t =
						start->second.end();
				while (s != t) {
					hashKey += calculator(*s);
					isFinal = (isFinal || (finalStates.find(*s) != finalStates.end()));
					++s;
				}

				states.back().swap(start->second);
				TempState tmpState(&states, states.size() - 1, hashKey,
						isFinal);
				typename std::unordered_set<TempState, HashFunc>::const_iterator pos =
						processedStates.find(tmpState);
				if (pos == processedStates.end()) {
					processedStates.insert(tmpState);
					result.setTransition(CompoundState(&result, i),
							start->first,
							CompoundState(&result, tmpState.index));
					if (tmpState.isFinal)
						result.addFinal(CompoundState(&result, tmpState.index));
					states.push_back(StateSet());
				} else {
					result.setTransition(CompoundState(&result, i),
							start->first, CompoundState(&result, pos->index));
				}
				++start;
			}
			states.pop_back();
		}
	}
	void concat(const NDFAutomaton& another, NDFAutomaton& output){

	}
	void unite(const NDFAutomaton& another, NDFAutomaton& output){

	}
	void kleeneClosure(NDFAutomaton& output){

	}
};
template<class State, class Input, class HashKeyCalculator>
ostream& NDFAutomaton<State, Input, HashKeyCalculator>::CompoundState::output(
		ostream& os) const{
	typename StateSet::const_iterator begin =
			automaton->states[index].begin(), end =
			automaton->states[index].end();
	os << '{';

	bool isFirst = true;
	while (begin != end) {
		if (!isFirst)
			os << ',';
		os << *begin;
		isFirst = false;
		++begin;
	}
	os << '}';
	return os;
}
int test20(int argc, char *argv[]) {
	typedef char State;
	typedef bool Input;
	typedef bool Output;
	MooreMachine<State, Input, HashKeyCalculator, Output> mooreMachine;
	State nextState;
	Output output;
	vector<Output> result;
	bool start[] = { 1, 0, 1, 0, 1, 1 };
	bool *begin = start, *end = start + ARSIZE(start);
	/*
	 * An example moore machine
	 *
	 * ------------------------
	 * 	      |  f    |   g   |
	 *        |-------|--------
	 * State  | 0 | 1 | 0 | 1 |
	 * ------------------------
	 *    a   | b | d | 1 | 0 |
	 *    b   | b | c | 1 | 1 |
	 *    c   | d | e | 0 | 0 |
	 *    d   | b | a | 0 | 0 |
	 *    e   | d | e | 0 | 0 |
	 * ------------------------
	 */
	mooreMachine.addFunction('a', false, 'b', true);
	mooreMachine.addFunction('a', true, 'd', false);

	mooreMachine.addFunction('b', false, 'b', true);
	mooreMachine.addFunction('b', true, 'c', true);

	mooreMachine.addFunction('c', false, 'd', false);
	mooreMachine.addFunction('c', true, 'e', false);

	mooreMachine.addFunction('d', false, 'b', false);
	mooreMachine.addFunction('d', true, 'a', false);

	mooreMachine.addFunction('e', false, 'd', false);
	mooreMachine.addFunction('e', true, 'e', false);

	mooreMachine.getFunction('a', false, nextState, output);
	assert(nextState == 'b' && output == true);
	mooreMachine.getOutput('a', begin, end, nextState, result);
	cout << "moore machine input : ";
	copy(start, end, ostream_iterator<bool>(cout, ""));
	cout << endl;
	cout << "moore machine output : ";
	copy(result.begin(), result.end(), ostream_iterator<bool>(cout, ""));
	cout << endl;

	MealyMachine<State, Input, HashKeyCalculator, Output> mealyMachine;
	/*
	 * An example mealy machine
	 *
	 * --------------------
	 * 	      |  f    | g |
	 *        |-------|----
	 * State  | 0 | 1 | 0 |
	 * --------------------
	 *    a   | a | c | 0 |
	 *    b   | d | a | 1 |
	 *    c   | c | b | 1 |
	 *    d   | c | a | 1 |
	 * --------------------
	 */
	mealyMachine.addFunction('a', false, 'a');
	mealyMachine.addFunction('a', true, 'c');

	mealyMachine.addFunction('b', false, 'd');
	mealyMachine.addFunction('b', true, 'a');

	mealyMachine.addFunction('c', false, 'c');
	mealyMachine.addFunction('c', true, 'b');

	mealyMachine.addFunction('d', false, 'c');
	mealyMachine.addFunction('d', true, 'a');

	mealyMachine.addOutput('a', false);
	mealyMachine.addOutput('b', true);
	mealyMachine.addOutput('c', true);
	mealyMachine.addOutput('d', true);

	output = true;
	nextState = 'b';
	mealyMachine.getOutput('a', output);
	mealyMachine.getFunction('a', false, nextState);
	assert(nextState == 'a');
	assert(output == false);
	result.clear();
	begin = start;
	mealyMachine.getOutput('a', begin, start + ARSIZE(start), nextState,
			result);
	cout << "mealy machine input : ";
	copy(start, end, ostream_iterator<bool>(cout, ""));
	cout << endl;
	cout << "mealy machine output : ";
	copy(result.begin(), result.end(), ostream_iterator<bool>(cout, ""));
	cout << endl;
	return 0;
}
class IntState{
	int value;
public:
	IntState(int val = (1 << 31)){
		this->value = val;
	}
	operator int() const{
		return value;
	}
	friend ostream& operator<<(ostream& os, const IntState& state){
		os << 's' << state.value;
		return os;
	}
};
int test25(int argc, char *argv[]) {
	typedef IntState State;
	typedef bool Input;
	typedef DFAutomaton<State, Input, HashKeyCalculator> DFA;
	typedef NDFAutomaton<State, Input, HashKeyCalculator> NDFA;

	DFA dfa;
	NDFA ndfa;
	NDFA::CDFAutomaton cdfa;
	const int numOfBits = 20;
	vector<bool> vals;
	/*bool vals[]={0,1,1,0,1,1, 1, 0,0,0};
	 bool *begin = vals, *end = vals + ARSIZE(vals);*/

	srand(time(0));
	for (int i = 0; i < numOfBits; ++i)
		vals.push_back(rand() % 2 == 0);
	vector<bool>::iterator begin = vals.begin(), end = vals.end();
	/*
	 * A Deterministic Fnite-State Automaton Recognizing the Set of
	 * Bit Strings Containing an Odd Number of 1s and Ending with
	 * at Least Two 0s
	 */
	dfa.setAsStart(0);
	dfa.setTransition(0, false, 1);
	dfa.setTransition(0, true, 3);
	dfa.setTransition(1, false, 2);
	dfa.setTransition(1, true, 3);
	dfa.setTransition(2, true, 3);
	dfa.setTransition(2, false, 2);
	dfa.setTransition(3, false, 4);
	dfa.setTransition(3, true, 0);
	dfa.setTransition(4, true, 0);
	dfa.setTransition(4, false, 5);
	dfa.setTransition(5, false, 5);
	dfa.setTransition(5, true, 0);
	dfa.addFinal(5);
	cout << dfa << endl;
	cout << '\'';
	copy(begin, end, ostream_iterator<bool>(cout, ""));
	cout << '\'';
	cout << (dfa.recognize(begin, end) ? " is " : " isn't ")
			<< "recognized by the Deterministic Finite-state Automaton."
			<< endl;
	/*ndfa.setTransition(0, true, 3);
	ndfa.setTransition(0, false, 0);
	ndfa.setTransition(0, false, 1);
	ndfa.setTransition(1, false, 0);
	ndfa.setTransition(1, true, 1);
	ndfa.setTransition(1, true, 3);
	ndfa.setTransition(2, true, 0);
	ndfa.setTransition(2, true, 2);
	ndfa.setTransition(3, false, 0);
	ndfa.setTransition(3, false, 1);
	ndfa.setTransition(3, false, 2);
	ndfa.setTransition(3, true, 1);
	ndfa.setAsStart(0);
	ndfa.addFinal(3);*/
	ndfa.setTransition(0, false, 0);
	ndfa.setTransition(0, false, 2);
	ndfa.setTransition(0, true, 1);
	ndfa.setTransition(1, false, 3);
	ndfa.setTransition(1, true, 4);
	ndfa.setTransition(2, true, 4);
	ndfa.setTransition(3, false, 3);
	ndfa.setTransition(4, false, 3);
	ndfa.setTransition(4, true, 3);
	ndfa.setAsStart(0);
	ndfa.addFinal(0);
	ndfa.addFinal(4);

	NDFA tmp, result;
	ndfa.toDFAutomaton(cdfa);
	cout << ndfa << endl;
	cout << cdfa << endl;
	/*bool vals1[] =
	 { 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, };
	 {0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1,};
	 bool *begin1 = vals1, *end1 = vals1 + ARSIZE(vals1);*/
	cout << '\'';
	begin = vals.begin();
	copy(begin, end, ostream_iterator<bool>(cout, ""));
	cout << '\'';
	vector<State> states;
	bool isRecognized = ndfa.recognize(begin, end, &states);
	cout << (isRecognized ? " is " : " isn't ")
			<< "recognized by the Non-Deterministic Finite-state Automaton."
			<< endl;
	begin = vals.begin();
	cout << cdfa.recognize(begin, end) << endl;
	/*if(isRecognized)
	 copy(states.begin(), states.end(), ostream_iterator<char>(cout, ""));*/
	return 0;
}

