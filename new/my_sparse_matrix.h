#ifndef MY_LIB_SPARSE_MATRIX_H
#define MY_LIB_SPARSE_MATRIX_H
#include <exception>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <cassert>

namespace my_lib{

/**
 * stl symbol tables, like std::map, std::unordered_map, don't have the same
 * interface as my symbol table interface.
 * Need this class to adapt stl symbol table interface to our symbol table iterface
 */

template<class StdST>
class StdSTAdapter{
	StdST stdST;
public:
	template<class Base>
	class Iterator : public Base{
	public:
		Iterator(const Base& base) : Base(base){
		}
		const typename StdST::key_type& key() const {
			return (Base::operator->)()->first;
		}
		const typename StdST::mapped_type& value() const {
			return (Base::operator->)()->second;
		}
		const Iterator& operator*() const{
			return *this;
		}
	};
	StdSTAdapter(const StdST& stdST = StdST()) : stdST(stdST){
	}

	template<class K>
	const typename StdST::mapped_type* get(const K& k) const {
		auto pos = stdST.find(k);
		if(pos == stdST.end())
			return nullptr;
		return &pos->second; 
	}

	template<class K>
	void remove(const K& k){
		stdST.erase(k);
	}

	template<class K>
	bool contains(const K& k) const {
		return stdST.find(k) != stdST.end();
	}

	template<class K, class V>
	void put(const K& k, const V& v){
		stdST[k] = v;
	}

	unsigned int size() const {
		return stdST.size();
	}

	Iterator<typename StdST::const_iterator> begin() const {
		return Iterator<typename StdST::const_iterator>(stdST.begin());
	}

	Iterator<typename StdST::const_iterator> end() const {
		return Iterator<typename StdST::const_iterator>(stdST.end());
	}

};

/*
 * A sparse array is an array where most elements have the same value
 * except some few elements.
 * The frequency of a value is the number of elements that have this value.
 * Define the value that has the highest frequency as the default value.
 * A non-default elmenent is an element whose value is not the default value.
 * This implemtation uses a symbol table to map the index of a non-default elements to it's value
 * Note that the default value can be arbitrary.
 */
template<class T, class ST>
class SparseArray{
	ST st;
	unsigned int n;
	T defaultValue;
#define MY_LIB_SPARSE_VECTOR_CHECK_AGAINST(theOther)\
	if(n != theOther.dim()){\
		std::stringstream os;\
		os << "the dimensions of the two vectors must be the same. One is " << n\
		   << ", the other is " << theOther.dim() << ", line " << __LINE__ << " file " << __FILE__;\
		throw std::invalid_argument(os.str());\
	}

#define MY_LIB_SPARSE_VECTOR_OPERATE_AGAINST(theOther, Statement1)\
	MY_LIB_SPARSE_VECTOR_CHECK_AGAINST(theOther)\
	auto& theOtherNonDefaults = theOther.getNonDefaults();\
	{\
		auto begin = theOtherNonDefaults.begin();\
		auto end = theOtherNonDefaults.end();\
		while(begin != end){\
			Statement1;\
			++begin;\
		}\
	}\

#define MY_LIB_SPARSE_VECTOR_OPERATOR_ASSIGN(op) \
template<class U, class V> \
SparseArray& operator op##=(const SparseArray<U, V>& theOther) throw(std::invalid_argument) {\
	MY_LIB_SPARSE_VECTOR_OPERATE_AGAINST(theOther, \
		set(begin.key(), (*this)[begin.key()] op begin.value()))\
	auto& theOtherDefault = theOther.getDefaultValue();\
	auto begin = st.begin(), end = st.end();\
	defaultValue op##= theOtherDefault;\
	std::vector<unsigned int> keys;\
	while(begin != end){\
		if(!theOtherNonDefaults.contains(begin.key()))\
			keys.push_back(begin.key());\
		++begin;\
	}\
	for(auto key : keys)\
		set(key, (*this)[key] op theOtherDefault);\
	return *this;\
}

#define MY_LIB_SPARSE_VECTOR_OPERATOR(op) \
template<class U>\
SparseArray operator op(const U& theOther) const throw(std::invalid_argument) {\
	SparseArray vec(*this);\
	vec op##= theOther;\
	return vec;\
}

#define MY_LIB_SPARSE_VECTOR_SCALAR_OPERATOR_ASSIGN(op) \
template<class U>\
SparseArray& operator op##=(const U& val){\
	auto begin = st.begin(), end = st.end();\
	while(begin != end){\
		set(begin.key(), begin.value() op val);\
		++begin;\
	}\
	defaultValue op##= val;\
	return *this;\
}

public:
	SparseArray(unsigned int n, const T& defaultValue = T()) 
		: n(n), defaultValue(defaultValue){
	}

	T operator[](unsigned int index) const{
		auto value = st.get(index);
		if(value == nullptr)
			return defaultValue;
		return *value;
	}

	bool set(unsigned int index, const T val = T()){
		if(index >= n)
			return false;
		if(val == defaultValue){
			st.remove(index);
		}else{
			st.put(index, val);
		}
		return true;
	}

	unsigned int nonDefaultCount() const {
		return st.size();
	}

	const ST& getNonDefaults() const {
		return st;
	}

	const T& getDefaultValue() const {
		return defaultValue;
	}

	unsigned int dim() const {
		return n;
	}

	template<class U, class V>
	T dot(const SparseArray<U, V>& theOther) throw(std::invalid_argument){
		unsigned int count = 0;
		T sum = T();
		MY_LIB_SPARSE_VECTOR_OPERATE_AGAINST(theOther, \
			{sum += (*this)[begin.key()] * begin.value(); ++ count;})\
		auto& theOtherDefault = theOther.getDefaultValue();
		auto begin = st.begin(), end = st.end();
		while(begin != end){
			if(!theOtherNonDefaults.contains(begin.key())) {
				sum += begin.value() * theOtherDefault;
				++ count;
			}
			++begin;
		}
		sum += (n - count) * defaultValue * theOther.getDefaultValue();
		return sum;
	}

	MY_LIB_SPARSE_VECTOR_OPERATOR_ASSIGN(+)
	MY_LIB_SPARSE_VECTOR_OPERATOR_ASSIGN(-)
	MY_LIB_SPARSE_VECTOR_OPERATOR_ASSIGN(*)
	MY_LIB_SPARSE_VECTOR_OPERATOR_ASSIGN(/)
	MY_LIB_SPARSE_VECTOR_SCALAR_OPERATOR_ASSIGN(+)
	MY_LIB_SPARSE_VECTOR_SCALAR_OPERATOR_ASSIGN(-)
	MY_LIB_SPARSE_VECTOR_SCALAR_OPERATOR_ASSIGN(*)
	MY_LIB_SPARSE_VECTOR_SCALAR_OPERATOR_ASSIGN(/)
	MY_LIB_SPARSE_VECTOR_OPERATOR(+)
	MY_LIB_SPARSE_VECTOR_OPERATOR(-)
	MY_LIB_SPARSE_VECTOR_OPERATOR(*)
	MY_LIB_SPARSE_VECTOR_OPERATOR(/)

	template<class Stream>
	friend Stream& operator<<(Stream& os, const SparseArray& array){
		os << "[";
		for(unsigned int i = 0;i < array.dim();++ i){
			if(i > 0)
				os << '\t';
			os << array[i];
		}
		os << "]";
		return os;
	}

	/*template<class T>
	void optimize(){
		if(st.size() > n / 2){
			T counts;
			for(auto pair: st){
				auto count = st.get(pair.value());
				counts.put(pair.value(), (count == nullptr ? 1 : *count + 1));
			}
			for(auto pair: counts){
				
			}
		}
	}*/
};


template<class T, class ST>
using SparseVector = SparseArray<T, ST>;


/*
 * Let T be the type of value.
 * A sparse matrix is a matrix where most elements have the value T()
 * A non-default elmenent is an element whose value is not T()
 * This implemtation uses a symbol table to map the index of a non-default elements to it's value.
 * The index of an element is calculated as i * col + j, where the element is at ith row and jth column,
 * 	and col is the number of columns.
 *
 * It's hard to write an efficient implementation without making some assumptions, especially for matrix multiplication:
 * Let value be arbitrary, this class has the following assumptions
 * 	1. T() * value = T()
 * 	2. T() + value = value
 * 	3. value - T() = value
 * 	4. T() / value = value
 */
template<class T, class ST>
class SparseMatrix{
	ST st;
	typedef std::map<unsigned int, std::set<unsigned int> > IndexMap;
	//row number to non-zero column indices map.
	IndexMap rowIndices;
	//column number to non-zero row indices map.
	IndexMap colIndices;
	unsigned int r, c;
public:
#define MY_LIB_SPARSE_MATRIX_CHECK_AGAINST(theOther)\
	if(r != theOther.row()){\
		std::stringstream os;\
		os << "the rows of the two matrices must be the same. One is " << r\
		   << ", the other is " << theOther.row() << ", line " << __LINE__ << " file " << __FILE__;\
		throw std::invalid_argument(os.str());\
	}\
	if(c != theOther.col()){\
		std::stringstream os;\
		os << "the columns of the two matrices must be the same. One is " << c\
		   << ", the other is " << theOther.col() << ", line " << __LINE__ << " file " << __FILE__;\
		throw std::invalid_argument(os.str());\
	}

#define MY_LIB_SPARSE_MATRIX_CHECK_AGAINST0(theOther)\
	if(c != theOther.row()){\
		std::stringstream os;\
		os << "columns of the first matrix must equal rows of the second matrix. columns of first matrix is " << c\
		   << ", rows of the second matrix is " << theOther.row() << ", line " << __LINE__ << " file " << __FILE__;\
		throw std::invalid_argument(os.str());\
	}


	SparseMatrix(unsigned int row, unsigned int col) : c(col), r(row){
	}

	unsigned int row() const {
		return r;
	}

	unsigned int col() const {
		return c;
	}

	double ratio() const {
		return st.size() / (double)(r * c);
	}

	T get(unsigned int i, unsigned int j) const{
		auto value = st.get(i * c + j);
		if(value == nullptr)
			return T();
		return *value;
	}

	bool set(unsigned int i, unsigned int j, const T& val = T()){
		if(i >= r || j >= c)
			return false;
		if(val == T()){
			if(st.get(i * c + j) == nullptr)
				return true;
			auto pos = rowIndices.find(i);
			assert(pos != rowIndices.end());
			assert(!pos->second.empty());
			pos->second.erase(i);
			if(pos->second.empty())
				rowIndices.erase(pos);

			pos = colIndices.find(j);
			assert(pos != colIndices.end());
			assert(!pos->second.empty());
			pos->second.erase(j);
			if(pos->second.empty())
				colIndices.erase(pos);

			st.remove(i * c + j);
		}else{
			st.put(i * c + j, val);
			rowIndices[i].insert(j);
			colIndices[j].insert(i);
		}
		return true;
	}

	unsigned int nonDefaultCount() const {
		return st.size();
	}

	const ST& getNonDefaults() const {
		return st;
	}

#define MY_LIB_SPARSE_MATRIX_OPERATOR_ASSIGN(op) \
template<class U>\
SparseMatrix& operator op##=(const U& theOther) throw(std::invalid_argument) {\
	MY_LIB_SPARSE_MATRIX_CHECK_AGAINST(theOther)\
	for(auto e : theOther.st){\
		auto val = st.get(e.key());\
		set(e.key() / c, e.key() % c, (val == nullptr ? T() : *val) op e.value());\
	}\
	return *this;\
}

#define MY_LIB_SPARSE_MATRIX_OPERATOR(op) \
template<class U>\
SparseMatrix operator op(const U& theOther) const throw(std::invalid_argument) {\
	SparseMatrix matrix(*this);\
	matrix op##= theOther;\
	return matrix;\
}

	MY_LIB_SPARSE_MATRIX_OPERATOR(+)
	MY_LIB_SPARSE_MATRIX_OPERATOR(-)
	MY_LIB_SPARSE_MATRIX_OPERATOR(/)
	MY_LIB_SPARSE_MATRIX_OPERATOR_ASSIGN(+)
	MY_LIB_SPARSE_MATRIX_OPERATOR_ASSIGN(-)
	MY_LIB_SPARSE_MATRIX_OPERATOR_ASSIGN(/)

	template<class U>
	SparseMatrix operator*(const U& theOther) const throw(std::invalid_argument) {
		MY_LIB_SPARSE_MATRIX_CHECK_AGAINST0(theOther)

		SparseMatrix result(row(), theOther.col());
		for(auto i : rowIndices)
			for(auto j : theOther.colIndices){
				//Note:
				//	T sum = T();
				//is different from
				//	T sum;
				//for primitive types
				//
				//The later one does not initliaze a value for primitive type.
				T sum = T();
				auto& indices = (j.second.size() < i.second.size() ? j.second : i.second);
				for(auto k : indices)
			  		sum = sum + theOther.get(k, j.first) * get(i.first, k);
				result.set(i.first, j.first, sum);
				assert(result.get(i.first, j.first) == sum);
			}

		return result;
	}
	template<class Stream>
	friend Stream& operator<<(Stream& os, const SparseMatrix& matrix){
		for(unsigned int i = 0;i < matrix.row();++ i){
			if(i > 0)
				os << '\n';
			for(unsigned int j = 0;j < matrix.col();++ j){
				if(j > 0)
					os << '\t';
				os << matrix.get(i, j);
			}
		}
		return os;
	}
};

}
#endif //MY_LIB_SPARSE_MATRIX_H
