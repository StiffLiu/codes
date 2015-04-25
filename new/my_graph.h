#ifndef MY_LIB_GRAPH_H
#define MY_LIB_GRAPH_H
#include "my_math.h"
#include <map>
#include <set>
#include <cassert>
#include <sstream>
#include <exception>
#include <vector>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <stack>
#include <algorithm>
namespace my_lib{

template<class K, class V>
using DefaultStdMap = std::map<K, V>;
template<class K>
using DefaultStdSet = std::set<K>;

template<template<class U, class V> class Map = DefaultStdMap, template<class U> class List = DefaultStdSet >
class AdjListT{
public:
	typedef List<unsigned int> Vertices;
	typedef Map<unsigned int, Vertices> Edges;

	AdjListT(unsigned int vertexCount) : vertexCount(vertexCount){
	}

	unsigned int getVertexCount() const {
		return vertexCount;
	}

	void reset(unsigned int vertexCount){
		this->vertexCount = vertexCount;
		edges.clear();
	}

	const Edges& getEdges() const {
		return edges;
	}

	template<class Stream>
	friend Stream& operator<<(Stream& stream, const AdjListT& graph){
		for (auto& vertices : graph.getEdges()){
			stream << vertices.first;
			stream << "\t:";
			for (auto vertex : vertices.second)
				stream << vertex << '\t';
			stream << '\n';
		}
		return stream;
	}

protected:
	bool isAdj(unsigned int i, unsigned int j) const {
		auto pos = edges.find(i);
		if(pos != edges.end())
			return pos->second.find(j) != pos->second.end();
		return false;
	}

	void addEdgeNoCheck(unsigned int i, unsigned int j) {
		edges[i].insert(j);
	}

	void check(unsigned int i) const throw(std::invalid_argument) {
		if(i >= vertexCount){
			std::stringstream os;
			os << "vertex " << i << " out of range, max is " << vertexCount << ", line " << __LINE__ << " file " << __FILE__;
			throw std::invalid_argument(os.str());
		}
	}

	void check(unsigned int i, unsigned int j) const throw(std::invalid_argument){
		if(i == j){
			std::stringstream os;
			os << "self loop for vertex " << i << " not allowed, line " << __LINE__ << " file " << __FILE__;
			throw std::invalid_argument(os.str());
		}
	}

	void addEdgeCheck(unsigned int i, unsigned int j) throw(std::invalid_argument) {
		check(i);
		check(j);
		check(i, j);
		addEdgeNoCheck(i, j);
	}

	void addEdge(unsigned int i, unsigned int j) throw(std::invalid_argument) {
		addEdgeCheck(i, j);
	}

	unsigned int degree(unsigned int i) const {
		auto pos = edges.find(i);
		if(pos != edges.end())
			return pos->second.size();
		return 0;
	}

	const Vertices* adj(unsigned int i) const {
		auto pos = edges.find(i);
		if(pos != edges.end())
			return &pos->second;
		return nullptr;
	}

	void extend(unsigned int count){
		vertexCount += count;
	}

	Edges edges;
	unsigned int vertexCount;
};

template<template<class U, class V> class Map = DefaultStdMap, template<class U> class List = DefaultStdSet >
class UndirAdjListT : public AdjListT<Map, List>{
	typedef AdjListT<Map, List> Super;
public:
	UndirAdjListT(unsigned int vertexCount = 0) : Super(vertexCount){
	}
	using Super::isAdj;

	void addEdge(unsigned int i, unsigned int j){
		Super::check(i);
		Super::check(j);
		Super::check(i, j);
		Super::addEdgeNoCheck(i, j);
		Super::addEdgeNoCheck(j, i);
	}

	using Super::degree;
	using Super::adj;
	using Super::extend;
};

template<template<class U, class V> class Map = DefaultStdMap, template<class U> class List = DefaultStdSet >
class DirAdjListT : public AdjListT<Map, List>{
	typedef AdjListT<Map, List> Super;
public:
	DirAdjListT(unsigned int vertexCount = 0) : Super(vertexCount){
	}
	using Super::isAdj;
	using Super::addEdge;
	using Super::extend;

	unsigned int outDegree(unsigned int i) const {
		return Super::degree(i);
	}

	const typename Super::Vertices* outAdj(unsigned int i) const {
		return Super::adj(i);
	}

	bool isMap() const {
		if(Super::getVertexCount() != Super::edges.size())
			return false;
		for(auto& v : Super::edges)
			if(v.second.size() != 1)
				return false;
		return true;
	}

	DirAdjListT toReverse() const{
		DirAdjListT reverse(Super::getVertexCount());
		for (auto& vertex : Super::getEdges()){
			for (auto v : vertex.second){
				reverse.addEdge(v, vertex.first);
			}
		}
		return reverse;
	}
};

template<template<class U, class V> class Map = DefaultStdMap, template<class U> class List = DefaultStdSet >
class DblDirAdjListT : protected DirAdjListT<Map, List>{
	typedef DirAdjListT<Map, List> Super;
	DirAdjListT<Map, List> reverse;
public:
	DblDirAdjListT(unsigned int vertexCount = 0) : Super(vertexCount), reverse(vertexCount){
	}
	using Super::isAdj;
	using Super::getVertexCount;

	void addEdge(unsigned int i, unsigned int j){
		Super::addEdge(i, j);
		reverse.addEdge(j, i);
	}
	
	void extend(unsigned int count){
		Super::extend(count);
		reverse.extend(count);
	}

	using Super::outDegree;
	using Super::outAdj;
	using Super::isMap;

	unsigned int inDegree(unsigned int i) const {
		return reverse.outDegree(i);
	}
	
	auto inAdj(unsigned int i) const -> decltype(reverse.outAdj(i)) {
		return reverse.outAdj(i);
	}

	void reset(unsigned int vertexCount){
		Super::reset(vertexCount);
		reverse.reset(vertexCount);
	}

	auto getReverse() const -> decltype(reverse){
		return reverse;
	}

};

using UndirAdjList = UndirAdjListT<>;
typedef DirAdjListT<> DirAdjList;
typedef DblDirAdjListT<> DblDirAdjList;


/*
 * Graph stored in text string.
 * First number is the number of vertices.
 * Second number is the number of edges
 * For the rest numbers, every two represet an edge of the graph.
 * */
template<class Stream>
class TxtStreamGraphReader{
	Stream& stream;
public:
	TxtStreamGraphReader(Stream& stream) : stream(stream){
	}

	operator bool() const {
		return stream;
	}

	template<class Graph>
	const TxtStreamGraphReader& operator>>(Graph& graph) const throw(std::exception){
		int vertexCount;
		if(stream >> vertexCount){
			graph.reset(vertexCount);
			int edgeCount;
			if(stream >> edgeCount){
				for(int i = 0;i < edgeCount;++ i){
					int v1, v2;
					if(stream >> v1 >> v2){
						graph.addEdge(v1, v2);
					}else{
						std::stringstream os;
						os << "error reading " << i << "th edge " << ", line " << __LINE__ << " file " << __FILE__;
						throw std::invalid_argument(os.str());
					}

				}
			}
		}
		return *this;
	}
};

class TxtFileGraphReader{
	std::string filePath;
	mutable bool result = false;
public:
	TxtFileGraphReader(const char *filePath) : filePath(filePath){
	}

	TxtFileGraphReader(const std::string& filePath) : filePath(filePath){
	}

	operator bool() const {
		return result;
	}

	template<class Graph>
	const TxtFileGraphReader& operator>>(Graph& graph) const throw(std::exception) {
		std::ifstream in(filePath.c_str());
		TxtStreamGraphReader<std::ifstream> reader(in);
		reader >> graph;
		result = bool(in);
		return *this;
	}
	
};

class TxtFileStreamGraphReader{
	std::ifstream stream;
	TxtStreamGraphReader<std::ifstream> reader;
public:
	TxtFileStreamGraphReader(const char *filePath) : stream(filePath), reader(stream){
	}
	template<class Graph>
	const TxtFileStreamGraphReader& operator>>(Graph& graph) const throw(std::exception){
		reader >> graph;
		return *this;
	}
};

template<class Derived>
struct UndirSearchStrategy{
	template<class Graph, class CC>
	void operator()(const Graph& graph, CC& cc, unsigned int& componentCount) const {
		(*this)(graph, cc, componentCount, (Derived&)*this);
	}

	template<class Graph, class CC, class Strategy>
	void operator()(const Graph& graph, CC& cc, unsigned int& componentCount, Strategy& strategy) const {
		for(auto vertices : graph.getEdges()){
			if(cc(vertices.first)){
				//cc(vertext, parent of vertex, root of the tree where the vertex lives in)
				cc(vertices.first, vertices.first, vertices.first);
				strategy(vertices.first, vertices.first, graph, cc);
				++componentCount;
			}
		}
	}

};

/*
 * Depth first search algorithm for undirected graph.
 * */
struct UndirDfsStrategy : public UndirSearchStrategy<UndirDfsStrategy> {
	using UndirSearchStrategy::operator();

	template<class Graph, class CC>
	unsigned int operator()(unsigned int source, unsigned int root, const Graph& graph, CC& cc) const {
		return (*this)(source, root, graph, cc, *this);
	}

	template<class Graph, class CC, class Strategy>
	unsigned int operator()(unsigned int source, unsigned int root, const Graph& graph, CC& cc, const Strategy& strategy) const {
		auto adj = graph.adj(source);
		unsigned int depth = 0;
		if(adj != nullptr)
			for(auto v : *adj)
				//cc(source, v) should return true if v has not been vistited already
				if(cc(source, v)){
					//cc(vertext, parent of vertex, root of the tree where the vertex lives in)
					cc(v, source, root);

					auto newDepth = 1 + strategy(v, root, graph, cc);
					if(newDepth > depth)
						depth = newDepth;
				}
		return depth;
	}
};

/*
 * Breadth first search algorithm for undirected graph.
 * */
struct UndirBfsStrategy : public UndirSearchStrategy<UndirBfsStrategy>{
	using UndirSearchStrategy::operator();
	template<class Graph, class CC>
	unsigned int operator()(unsigned int source, unsigned int root, const Graph& graph, CC& cc){
		std::queue<unsigned int> parents;
		std::queue<unsigned int> children;
		parents.push(source);

		unsigned int depth = 0;
		while(!parents.empty()){
			while(!parents.empty()){
				unsigned int current = parents.front();
				parents.pop();

				auto adj = graph.adj(current);
				if(adj != nullptr)
					for(auto v : *adj)
						//cc(current, v) should return true if v has not been vistited already
						if(cc(current, v)){
							//cc(vertext, parent of vertex, root of the tree where the vertex lives in)
							cc(v, current, root);
							children.push(v);
						}

			}
			children.swap(parents);
			++depth;
		}
		return depth - 1;
	}

};

template<class Strategy>
struct DirSearchStrategyAdapter{
	Strategy strategy;
	DirSearchStrategyAdapter(const Strategy& strategy = Strategy()) : strategy(strategy){
	}
	template<class Graph>
	struct GraphAdapter{
		Graph& graph;
		GraphAdapter(Graph& graph) : graph(graph){
		}
		auto adj(unsigned int v) const -> decltype(graph.outAdj(v)) {
			return graph.outAdj(v);
		}
		auto getEdges() const -> decltype(graph.getEdges()) {
			return graph.getEdges();
		}
	};

	template<class Graph, class CC>
	void operator()(const Graph& graph, CC& cc, unsigned int& componentCount) const {
		GraphAdapter<const Graph> adapter(graph);
		strategy(adapter, cc, componentCount);
	}

	template<class Graph, class CC>
	unsigned int operator()(unsigned int source, unsigned int root, const Graph& graph, CC& cc) const {
		GraphAdapter<const Graph> adapter(graph);
		strategy(source, root, adapter, cc);
	}
};

using DirDfsStrategy = DirSearchStrategyAdapter<UndirDfsStrategy>;
using DirBfsStrategy = DirSearchStrategyAdapter<UndirBfsStrategy>;

/*
 * Calculating connnected components for undirected graph.
 * */
class ConnectedComponent{
	//<vertex, <parent of vertex, root of the dfs tree the vertex resides in> >
	//cc stands for connected component
	unsigned int componentCount = 0;
	class CC : public std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> >{
	public:
		//true if the v hasn't been visited
		bool operator()(unsigned int v){
			return find(v) == end();
		}
		bool operator()(unsigned int source, unsigned int v){
			return find(v) == end();
		}
		void operator()(unsigned int v, unsigned int p, unsigned int r){
			auto& parents = (*this)[v];
			parents.first = p;
			parents.second = r;
		}

	} cc;
public:
	template<class Graph, class Strategy>
	ConnectedComponent(const Graph& graph, Strategy strategy){
		strategy(graph, cc, componentCount);
	}

	/*
	 * Find a path between two vertices, returns true if the two vertices is connected, 
	 * else returns false.
	 */
	bool pathTo(unsigned int v1, unsigned int v2, 
		std::vector<unsigned int>* path = nullptr) const {
		if(v1 == v2)
			return true;
		auto pos1 = cc.find(v1), pos2 = cc.find(v2);
		if(pos1 != cc.end() && pos2 != cc.end()
			&& pos1->second.second == pos2->second.second){
			if(path != nullptr){
				std::vector<unsigned int> half;
				while(pos1->first != pos1->second.first){
					path->push_back(pos1->first);
					pos1 = cc.find(pos1->second.first);
				}
				path->push_back(pos1->first);

				while(pos2->first != pos2->second.first){
					half.push_back(pos2->first);
					pos2 = cc.find(pos2->second.first);
				}
				half.push_back(pos2->first);

				auto rbegin2 = half.rbegin(), rbegin1 = path->rbegin();
					while(rbegin2 != half.rend() && rbegin1 != path->rend()
						&& *rbegin2 == *rbegin1){
						++ rbegin1;
						++ rbegin2;
					}
					assert(path->empty() || rbegin1 != path->rbegin());

					auto count = path->rend() - rbegin1;
					if(count > 0)
						path->erase(path->begin() + count + 1, path->end());
				while(rbegin2 != half.rend()){
					path->push_back(*rbegin2);
					++rbegin2;
				}
			}
			return true;
		}
		return false;
	}	

	/*
	 * returns the number of the connected component of the graph.
	 * */
	unsigned int getComponentCount(){
		return componentCount;
	}

	template<class T>
	void getComponentSize(T& componentSize){
		for (auto& vertex : cc){
			++ componentSize[vertex.second.second];
		}
	}

	template<class T>
	void getComponentVertices(T& componentVertices){
		for (auto& vertex : cc){
			componentVertices[vertex.second.second].insert(vertex.first);
		}
	}
};

/*
 * Undirected graph properties.
 * */
class GraphProperties{
	struct VertexProperty{
		unsigned int eccentricity;
		unsigned int componentName;
		unsigned int shortestCycle;
	};
	std::unordered_map<unsigned int, VertexProperty> verticesProperties;
	unsigned int center;
	unsigned int peripheral;
	unsigned int girth;
	class CC : public std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> > {
	public:
		//Assuming that this will be called for every edge,
		//if not then shortestCycle won't be calculated correctly.
		//And assuming that when the vertex is found,
		//that indicate the existence of a cycle.
		bool operator()(unsigned int source, unsigned int v){
			auto pos = find(v);
			if(pos == end())
				return true;
			if(pos->second.first != v && pos->second.second + 1 < shortestCycle){
				//a back edge is found.
				//source should have been visited.
				auto pos1 = find(source);
				assert(pos1 != end());
				unsigned int len = pos1->second.second + pos->second.second + 1;
				if(len < shortestCycle)
					shortestCycle = len;
			}
			return false;
		}
		void operator()(unsigned int v, unsigned int p, unsigned int r){
			auto& info = (*this)[v];
			info.first = p;
			if(v != p){
				//p should have been visited.
				assert(find(p) != end());
				info.second += (*this)[p].second;
				info.second += 1;
			}
		}
		unsigned int shortestCycle = -1;
	};

public:
	/*
	 *Cost is O(V*E), where V is the number of vertices, E is the number of edges.
	 */
	template<class Graph>
	GraphProperties(const Graph& graph){
		UndirBfsStrategy bfs;
		unsigned int componentCount = 0;
		for(auto vertices : graph.getEdges()){
			CC tmp;
			auto& p = verticesProperties[vertices.first];
			auto& dat = tmp[vertices.first];
			dat.first = vertices.first;
			dat.second = 0;

			p.eccentricity = bfs(vertices.first, vertices.first, graph, tmp);
			p.shortestCycle = tmp.shortestCycle;
			if(p.componentName == 0){
				++componentCount;
				for(auto v : tmp)
					verticesProperties[v.first].componentName = componentCount;
			}
		}
		center = -1;
		peripheral = -1;
		girth = -1;
		for(auto& v : verticesProperties)
			if(center == -1){
				center = v.first;
				peripheral = v.first;
				girth = v.first;
			}else{
				if(v.second.eccentricity < verticesProperties[center].eccentricity){
					center = v.first;
				}
				if(v.second.eccentricity > verticesProperties[peripheral].eccentricity){
					peripheral = v.first;
				}
				if(v.second.shortestCycle < verticesProperties[girth].shortestCycle){
					girth = v.first;
				}
			}
	}

	unsigned int eccentricity(unsigned int v) const {
		auto pos = verticesProperties.find(v);
		if(pos == verticesProperties.end())
			return 0;
		return pos->second.eccentricity;
	}

	unsigned int getRadius() const {
		if(center == -1)
			return 0;
		auto pos = verticesProperties.find(center);
		assert(pos != verticesProperties.end());
		return pos->second.eccentricity;
	}

	unsigned int getCenter() const {
		return center;
	}

	unsigned int getDiameter(){
		if(peripheral == -1)
			return 0;
		auto pos = verticesProperties.find(peripheral);
		assert(pos != verticesProperties.end());
		return pos->second.eccentricity;
	}

	unsigned int getPeripheral() const {
		return peripheral;
	}

	unsigned int getGirth() const {
		if(girth == -1)
			return 0;
		auto pos = verticesProperties.find(girth);
		assert(pos != verticesProperties.end());
		return pos->second.shortestCycle;
	}

};

template<class Symbol, class UnderlyingGraph>
class SymbolAdjListT{
protected:
	typedef std::vector<Symbol> SymbolArray;
	struct Key{
		union{
			unsigned int index;
			const Symbol* symbol;
		}key;
		bool isIndex = true;
		Key(unsigned int index){
			key.index = index;
		}
		Key(const Symbol& symbol) : isIndex(false){
			key.symbol = &symbol;
		}
	};
	struct HashComparator{
		const SymbolArray* symbolArray;
		HashComparator(const SymbolArray* symbolArray = nullptr) : symbolArray(symbolArray){
		}
		//hash value for the symbol
		unsigned int operator()(Key key) const{
			assert(symbolArray != nullptr);
			return std::hash<Symbol>()(key.isIndex ? (*symbolArray)[key.key.index] : *key.key.symbol);
		}

		bool operator()(Key k1, Key k2) const{
			assert(symbolArray != nullptr);
			const Symbol& s1 = (k1.isIndex ? (*symbolArray)[k1.key.index] : *k1.key.symbol);
			const Symbol& s2 = (k2.isIndex ? (*symbolArray)[k2.key.index] : *k2.key.symbol);
			return s1 == s2;
		}
	};
	struct LessComparator{
		const SymbolArray* symbolArray;
		LessComparator(const SymbolArray* symbolArray = nullptr) : symbolArray(symbolArray){
		}
		bool operator()(Key k1, Key k2) const{
			assert(symbolArray != nullptr);
			const Symbol& s1 = (k1.isIndex ? (*symbolArray)[k1.key.index] : *k1.key.symbol);
			const Symbol& s2 = (k2.isIndex ? (*symbolArray)[k2.key.index] : *k2.key.symbol);
			return s1 < s2;
		}
	};

	typedef std::unordered_set<Key, HashComparator, HashComparator> SymbolMapBase;
	struct SymbolMap{
		SymbolMapBase base;
		SymbolArray& symbolArray;
	public:
		SymbolMap(SymbolArray& symbolArray, const HashComparator& comparator)
			: symbolArray(symbolArray), base(5, comparator, comparator){
		}
		typename SymbolMapBase::value_type  operator[](const Symbol& key){
			auto pos = base.find(key);
			if (pos == base.end()){
				symbolArray.push_back(key);
				base.insert(symbolArray.size() - 1);
				return symbolArray.size() - 1;
			}
			assert(pos->isIndex);
			return pos->key.index;
		}
		const typename SymbolMapBase::value_type*  operator[](const Symbol& key) const{
			auto pos = base.find(key);
			if (pos == base.end())
				return nullptr;
			assert(pos->isIndex);
			return &*pos;
		}
	};

	SymbolArray symbolArray;
	SymbolMap symbolMap;
	UnderlyingGraph underlyingGraph;
public:
	SymbolAdjListT() :symbolMap(symbolArray, HashComparator(&symbolArray)){
		Symbol symbol;
		symbolMap[symbol];
	}
	unsigned int getVertexCount() const {
		return underlyingGraph.getVertexCount();
	}

	void clear(){
		underlyingGraph.reset(0);
		symbolArray.clear();
		symbolMap.clear();
	}

	const typename UnderlyingGraph::Edges& getEdges() const{
		return underlyingGraph.getEdges();
	}

	bool isAdj(const Symbol& s1, const Symbol& s2) const{
		auto index1 = symbolMap[s1];
		if (index1 != nullptr){
			auto index2 = symbolMap[s2];
			if (index2 != nullptr){
				return underlyingGraph.isAdj(*index1, *index2);
			}
		}
		return false;
	}

	void addEdge(const Symbol& s1, const Symbol& s2) throw(std::invalid_argument){
		if (s1 == s2){
			std::stringstream os;
			os << "self loop for vertex " << s1 << " , " << s2 << " not allowed, line " << __LINE__ << " file " << __FILE__;
			throw std::invalid_argument(os.str());
		}
		auto index1 = symbolMap[s1];
		auto index2 = symbolMap[s2];
		assert(symbolArray.size() >= underlyingGraph.getVertexCount());
		assert(index1.isIndex && index2.isIndex);
		underlyingGraph.extend(symbolArray.size() - underlyingGraph.getVertexCount());
		return underlyingGraph.addEdge(index1.key.index, index2.key.index);
	}

	const SymbolArray& symbols() const{
		return symbolArray;
	}

	unsigned int toVertex(const Symbol& symbol) const{
		auto index = symbolMap[symbol];
		if (index == nullptr)
			return -1;
		assert(index->isIndex);
		return index->key.index;
	}
};
template<class Symbol, class UnderlyingGraph>
class UndirSymbolAdjListT : public SymbolAdjListT<Symbol, UnderlyingGraph>{
	typedef SymbolAdjListT<Symbol, UnderlyingGraph> Super;
public:
	unsigned int degree(unsigned int i) const {
		return Super::underlyingGraph.degree(i);
	}

	const typename UnderlyingGraph::Vertices* adj(unsigned int i) const {
		return Super::underlyingGraph.adj(i);
	}
};
template<class Symbol, class UnderlyingGraph>
class DirSymbolGraphT : public SymbolAdjListT<Symbol, UnderlyingGraph>{
	typedef SymbolAdjListT<Symbol, UnderlyingGraph> Super;
public:
	unsigned int outDegree(unsigned int i) const {
		return Super::underlyingGraph.degree(i);
	}

	const typename UnderlyingGraph::Vertices* outAdj(unsigned int i) const {
		return Super::underlyingGraph.adj(i);
	}
};
template<class Symbol, class UnderlyingGraph>
class DblDirSymbolAdjListT : public DirSymbolGraphT<Symbol, UnderlyingGraph>{
	typedef SymbolAdjListT<Symbol, UnderlyingGraph> Super;
public:
	unsigned int inDegree(unsigned int i) const {
		return Super::underlyingGraph.inDegree(i);
	}

	const typename UnderlyingGraph::Vertices* inAdj(unsigned int i) const {
		return Super::underlyingGraph.inAdj(i);
	}
};
template<class Symbol>
using UndirSymbolAdjList = UndirSymbolAdjListT<Symbol, UndirAdjList>;
template<class Symbol>
using DirSymbolGraph = DirSymbolGraphT<Symbol, DirAdjList>;
template<class Symbol>
using DblDirSymbolAdjList = DblDirSymbolAdjListT<Symbol, DblDirAdjList>;

/*
 * Decide whether a given permutation of a directed graph is a topological order of the directed graph.
 * Algorithm:
 * 	1. get all the in-degree of all the vertices.
 * 	2. iterate through the permutation of vertices.
 * 	3. if the in-degree of the current visited vertex is not zero, then not topological order,
 * 	   else decrease the in-degree of the vertices which are out vertcies of the current visited vertex.
 * 	4 if all the vertices are visited, then it's a topological order.
 * */
template<class Graph>
/* 
 * reverse: reverse graph.
 * graph: original graph.
 * vertices: a permutation of vertices of the graph.
 * */
bool isTopologicalOrder(const Graph& reverse, const Graph& graph, const unsigned int* vertices){
	typedef std::unordered_map<unsigned int, unsigned int> VertexDegrees;
	VertexDegrees vertexDegrees;
	for (auto& v : reverse.getEdges())
		vertexDegrees[v.first] = v.second.size();
	for(unsigned int i = 0;i < graph.getVertexCount();++ i){
		auto pos = vertexDegrees.find(vertices[i]);
		if(pos != vertexDegrees.end() && pos->second != 0)
			return false;
		//Since graph is the reverse directed graph of the original graph,
		//it's out adjacency vertices are the in adjacency vertices of the original graph.
		auto adj = graph.outAdj(vertices[i]);
		if(adj)
			for(auto v : *adj){
				pos = vertexDegrees.find(v);
				assert(pos != vertexDegrees.end());
				assert(pos->second > 0);
				--pos->second;
			}
	}
	return true;
}

/*
 * Generate random undirected graphs
 * vertexCount : number of vertex in the undirected graph.
 * ratio : number of edges divided by maximum number of edges.
 * generator : a random number generator.
 */
template<class Graph, class Generator>
void randomUndirGraph(unsigned int vertexCount, double ratio, Graph& graph, Generator generator){
	graph.reset(vertexCount);
	if(ratio > 1.0)
		ratio = 1.0;
	unsigned int totalEdges = vertexCount * (vertexCount - 1) / 2;
	unsigned int numEdges = ratio * totalEdges;
	if(numEdges == 0)
		return;
	
	std::vector<unsigned int> selectedEdges;

	selectedEdges.resize(numEdges);
	randomSelect(totalEdges, numEdges, &selectedEdges[0], generator);
	for(auto e : selectedEdges){
		e += 1;
		unsigned int v1 = floorTriSqrt1(e);
		unsigned int v2 = e - v1 * (v1 - 1) / 2;
		if(v2 == 0){
			graph.addEdge(v1 - 2, v1 -1);
		}else{
			graph.addEdge(v1, v2 -1);
		}
	}
}	

/*
 * Generate random directed graph.
 * */
template<class Graph, class Generator>
void randomDirGraph(unsigned int vertexCount, double ratio, Graph& graph, Generator generator){
	graph.reset(vertexCount);
	if(ratio > 1.0)
		ratio = 1.0;
	unsigned int totalEdges = vertexCount * (vertexCount - 1);
	unsigned int numEdges = ratio * totalEdges;
	if(numEdges == 0)
		return;
	std::vector<unsigned int> selectedEdges;

	selectedEdges.resize(numEdges);
	randomSelect(totalEdges, numEdges, &selectedEdges[0], generator);
	for(auto e : selectedEdges){
		unsigned int v1 = e / (vertexCount - 1);
		unsigned int v2 = e % (vertexCount - 1);
		if(v2 >= v1)
			v2 += 1;
		graph.addEdge(v1, v2);
	}
}

/*
 * Topological sort for directed graph.
 * Tarjan's algorithm is more efficient.
 * */
class Topological{
	std::vector<unsigned int> topoOrder;
	std::vector<unsigned int> cycle;
	struct TopoSortStrategy : public UndirDfsStrategy{
		mutable std::unordered_set<unsigned int> onStack;
		mutable std::vector<unsigned int> topoOrder;
		mutable std::vector<unsigned int> cycle;
		mutable std::unordered_set<unsigned int> visited;
		mutable bool cycleFound = false;

		using UndirDfsStrategy::operator();
		template<class Graph, class CC>
		void operator()(const Graph& graph, CC& cc, unsigned int& componentCount) const {
			UndirDfsStrategy::operator()(graph, cc, componentCount, *this);
		}

		template<class Graph, class CC>
		unsigned int operator()(unsigned int source, unsigned int root, const Graph& graph, CC& cc) const {
			unsigned int depth = 0;
			if(!cycleFound){
				onStack.insert(source);
				cycle.push_back(source);
				depth = (*this)(source, root, graph, cc, *this);
				if(!cycleFound){
					onStack.erase(source);
					cycle.pop_back();
				}
			}else{
				depth = (*this)(source, root, graph, cc, *this);
			}
			topoOrder.push_back(source);
			return depth;
		}
		bool operator()(unsigned int v){
			return visited.find(v) == visited.end();
		}
		bool operator()(unsigned int source, unsigned int v){
			if(!cycleFound && onStack.find(v) != onStack.end()){
				cycleFound = true;
				for(auto pos = cycle.begin();pos != cycle.end();++ pos){
					if(*pos == v){
						cycle.erase(cycle.begin(), pos);
						break;
					}
				}
				assert(cycle[0] == v);
				onStack.clear();
			}
			return visited.find(v) == visited.end();
		}
		bool operator()(unsigned int source, unsigned int parent, unsigned int root){
			visited.insert(source);
			return true;
		}

	};
public:
	template<class Graph>
	Topological(Graph& graph){
		DirSearchStrategyAdapter<TopoSortStrategy> adapter;
		unsigned int componentCount;
		adapter(graph, adapter.strategy, componentCount);
		topoOrder.swap(adapter.strategy.topoOrder);
		std::reverse(topoOrder.begin(), topoOrder.end());
		if(adapter.strategy.cycleFound)
			cycle.swap(adapter.strategy.cycle);
	}

	auto directedCycle() const ->decltype(&cycle) {
		if(!cycle.empty())
			return &cycle;
		return nullptr;
	}

	auto order() const -> decltype(&topoOrder) {
		return &topoOrder;
	}
};

/*
 * Compute the strong connected component of a directed graph.
 * Tarjan's algorithm is more efficient.
 * */
class StrongConnectedComponent{
};

/*
 * Calculate euler cycle for each connected component of an undirected graph if exists.
 */

/*
 * Tarjan's algorithm
 * */

/*
* Compute greatest commmon descandent of a two vertices in a DAG.
*
* */

//Euler cycle for directed graph
}
#endif //MY_LIB_GRAPH_H
