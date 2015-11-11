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
#include <list>
#include <stack>
#include <algorithm>
#include <iostream>
namespace my_lib{

template<class K, class V>
using DefaultStdMap = std::map<K, V>;
template<class K>
using DefaultStdSet = std::set<K>;

template<template<class U, class V> class Map = DefaultStdMap, template<class U> class List = DefaultStdSet >
class AdjListT{
public:
	typedef unsigned int Vertex;
	typedef List<Vertex> Vertices;
	typedef Map<Vertex, Vertices> Edges;

	AdjListT(Vertex vertexCount) : vertexCount(vertexCount){
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
	bool isAdj(Vertex i, Vertex j) const {
		auto pos = edges.find(i);
		if(pos != edges.end())
			return pos->second.find(j) != pos->second.end();
		return false;
	}

	void addEdgeNoCheck(Vertex i, Vertex j) {
		edges[i].insert(j);
	}

	void check(Vertex i) const throw(std::invalid_argument) {
		if(i >= vertexCount){
			std::stringstream os;
			os << "vertex " << i << " out of range, max is " << vertexCount << ", line " << __LINE__ << " file " << __FILE__;
			throw std::invalid_argument(os.str());
		}
	}

	void check(Vertex i, Vertex j) const throw(std::invalid_argument){
		if(i == j){
			std::stringstream os;
			os << "self loop for vertex " << i << " not allowed, line " << __LINE__ << " file " << __FILE__;
			throw std::invalid_argument(os.str());
		}
	}

	void addEdgeCheck(Vertex i, Vertex j) throw(std::invalid_argument) {
		check(i);
		check(j);
		check(i, j);
		addEdgeNoCheck(i, j);
	}

	void addEdge(Vertex i, Vertex j) throw(std::invalid_argument) {
		addEdgeCheck(i, j);
	}

	unsigned int degree(Vertex i) const {
		auto pos = edges.find(i);
		if(pos != edges.end())
			return pos->second.size();
		return 0;
	}

	const Vertices* adj(Vertex i) const {
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
	typedef typename AdjListT<Map, List>::Vertex Vertex;
	UndirAdjListT(unsigned int vertexCount = 0) : Super(vertexCount){
	}
	using Super::isAdj;

	void addEdge(Vertex i, Vertex j){
		Super::check(i);
		Super::check(j);
		Super::check(i, j);
		Super::addEdgeNoCheck(i, j);
		Super::addEdgeNoCheck(j, i);
	}

	template<class Edge>
	void addEdge(const Edge& edge){
		addEdge(edge.getV1(), edge.getV2());
	}

	using Super::degree;
	using Super::adj;
	using Super::extend;
};

template<template<class U, class V> class Map = DefaultStdMap, template<class U> class List = DefaultStdSet >
class DirAdjListT : public AdjListT<Map, List>{
	typedef AdjListT<Map, List> Super;
public:
	typedef typename AdjListT<Map, List>::Vertex Vertex;
	DirAdjListT(unsigned int vertexCount = 0) : Super(vertexCount){
	}
	using Super::isAdj;
	using Super::addEdge;
	using Super::extend;

	template<class Edge>
	void addEdge(const Edge& edge){
		addEdge(edge.getV1(), edge.getV2());
	}

	unsigned int outDegree(Vertex i) const {
		return Super::degree(i);
	}

	const typename Super::Vertices* outAdj(Vertex i) const {
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
	typedef typename DirAdjListT<Map,List>::Vertex Vertex;
	DblDirAdjListT(unsigned int vertexCount = 0) : Super(vertexCount), reverse(vertexCount){
	}
	using Super::isAdj;
	using Super::getVertexCount;

	void addEdge(Vertex i, Vertex j){
		Super::addEdge(i, j);
		reverse.addEdge(j, i);
	}

	template<class Edge>
	void addEdge(const Edge& edge){
		addEdge(edge.getV1(), edge.getV2());
	}
	
	void extend(unsigned int count){
		Super::extend(count);
		reverse.extend(count);
	}

	using Super::outDegree;
	using Super::outAdj;
	using Super::isMap;

	unsigned int inDegree(Vertex i) const {
		return reverse.outDegree(i);
	}
	
	auto inAdj(Vertex i) const -> decltype(reverse.outAdj(i)) {
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
template<class Stream, class ReadOneEdge>
class TxtStreamGraphReaderT{
public:
	TxtStreamGraphReaderT(Stream& stream) : stream(stream){
	}

	operator bool() const {
		return stream;
	}

	template<class Graph>
	const TxtStreamGraphReaderT& operator>>(Graph& graph) const throw(std::exception){
		return read(graph);
	}

protected:
	template<class Graph>
	const TxtStreamGraphReaderT& read(Graph& graph) const{
		int vertexCount;
		if(stream >> vertexCount){
			graph.reset(vertexCount);
			int edgeCount;
			if(stream >> edgeCount){
				ReadOneEdge readOneEdge;
				for(int i = 0;i < edgeCount;++ i){
					if(!readOneEdge(stream, graph)){
						std::stringstream os;
						os << "error reading " << i << "th edge " << ", line " << __LINE__ << " file " << __FILE__;
						throw std::invalid_argument(os.str());
					}
				}
			}
		}
		return *this;
	}

	Stream& stream;
};

struct ReadOneEdge{
	template<class Stream, class Graph>
	bool operator()(Stream& stream, Graph& graph) const {
		int v1, v2;
		if(stream >> v1 >> v2){
			graph.addEdge(v1, v2);
			return true;
		}
		return false;
	}
};

template<class Reader>
class TxtFileGraphReaderT{
	std::string filePath;
	mutable bool result = false;
public:
	TxtFileGraphReaderT(const char *filePath) : filePath(filePath){
	}

	TxtFileGraphReaderT(const std::string& filePath) : filePath(filePath){
	}

	operator bool() const {
		return result;
	}

	template<class Graph>
	const TxtFileGraphReaderT& operator>>(Graph& graph) const throw(std::exception) {
		std::ifstream in(filePath.c_str());
		Reader reader(in);
		reader >> graph;
		result = bool(in);
		return *this;
	}
	
};
using TxtFileGraphReader = TxtFileGraphReaderT<TxtStreamGraphReaderT<std::ifstream, ReadOneEdge> >;

template<class Reader>
class TxtFileStreamGraphReaderT{
	std::ifstream stream;
	Reader reader;
public:
	TxtFileStreamGraphReaderT(const char *filePath) : stream(filePath), reader(stream){
	}
	template<class Graph>
	const TxtFileStreamGraphReaderT& operator>>(Graph& graph) const throw(std::exception){
		reader >> graph;
		return *this;
	}
};
using TxtFileStreamGraphReader= TxtFileStreamGraphReaderT<TxtStreamGraphReaderT<std::ifstream, ReadOneEdge> >;

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

	template<class Vertex, class Graph, class CC>
	unsigned int operator()(Vertex source, Vertex root, const Graph& graph, CC& cc) const {
		return (*this)(source, root, graph, cc, *this);
	}

	template<class Vertex, class Graph, class CC, class Strategy>
	unsigned int operator()(Vertex source, Vertex root, const Graph& graph, CC& cc, const Strategy& strategy) const {
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
	template<class Vertex, class Graph, class CC>
	unsigned int operator()(Vertex source, Vertex root, const Graph& graph, CC& cc){
		std::queue<Vertex> parents;
		std::queue<Vertex> children;
		parents.push(source);

		unsigned int depth = 0;
		while(!parents.empty()){
			while(!parents.empty()){
				auto current = parents.front();
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
		typedef typename Graph::Vertex Vertex;
		auto adj(Vertex v) const -> decltype(graph.outAdj(v)) {
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
	unsigned int operator()(typename Graph::Vertex source, typename Graph::Vertex root, const Graph& graph, CC& cc) const {
		GraphAdapter<const Graph> adapter(graph);
		return strategy(source, root, adapter, cc);
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
		if(center == static_cast<unsigned int>(-1))
			return 0;
		auto pos = verticesProperties.find(center);
		assert(pos != verticesProperties.end());
		return pos->second.eccentricity;
	}

	unsigned int getCenter() const {
		return center;
	}

	unsigned int getDiameter(){
		if(peripheral == static_cast<unsigned int>(-1))
			return 0;
		auto pos = verticesProperties.find(peripheral);
		assert(pos != verticesProperties.end());
		return pos->second.eccentricity;
	}

	unsigned int getPeripheral() const {
		return peripheral;
	}

	unsigned int getGirth() const {
		if(girth == static_cast<unsigned int>(-1))
			return 0;
		auto pos = verticesProperties.find(girth);
		assert(pos != verticesProperties.end());
		return pos->second.shortestCycle;
	}

};

template<class Symbol, class UnderlyingGraph>
class SymbolAdjListBaseT{
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
	SymbolAdjListBaseT(const UnderlyingGraph& underlyingGraph = UnderlyingGraph())
		:symbolMap(symbolArray, HashComparator(&symbolArray)), underlyingGraph(underlyingGraph){
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

	const UnderlyingGraph& getUnderlyingGraph() const {
		return underlyingGraph;
	}
protected:
	UnderlyingGraph& getUnderlyingGraph() {
		return underlyingGraph;
	}
};
template<class Symbol, class UnderlyingGraph>
class SymbolAdjListT : public SymbolAdjListBaseT<Symbol, UnderlyingGraph>{
	typedef SymbolAdjListBaseT<Symbol, UnderlyingGraph> Super;
public:	
	void addEdge(const Symbol& s1, const Symbol& s2) throw(std::invalid_argument){
		if (s1 == s2){
			std::stringstream os;
			os << "self loop for vertex " << s1 << " , " << s2 << " not allowed, line " << __LINE__ << " file " << __FILE__;
			throw std::invalid_argument(os.str());
		}
		auto index1 = Super::symbolMap[s1];
		auto index2 = Super::symbolMap[s2];
		assert(Super::symbolArray.size() >= Super::getUnderlyingGraph().getVertexCount());
		assert(index1.isIndex && index2.isIndex);
		Super::getUnderlyingGraph().extend(Super::symbolArray.size() - Super::getUnderlyingGraph().getVertexCount());
		return Super::getUnderlyingGraph().addEdge(index1.key.index, index2.key.index);
	}
};
template<class Base>
class UndirSymbolAdjListBaseT : public Base{
public:
	using Base::getUnderlyingGraph;
	unsigned int degree(unsigned int i) const {
		return getUnderlyingGraph().degree(i);
	}

	auto adj(unsigned int i) const ->decltype(this->getUnderlyingGraph().adj(i)) {
		return getUnderlyingGraph().adj(i);
	}
};

template<class Base>
class DirSymbolAdjListBaseT : public Base{
public:
	using Base::getUnderlyingGraph;
	unsigned int outDegree(unsigned int i) const {
		return getUnderlyingGraph().degree(i);
	}

	auto outAdj(unsigned int i) const ->decltype(getUnderlyingGraph().adj(i)) {
		return getUnderlyingGraph().adj(i);
	}
};

template<class Super>
class DblDirSymbolAdjListBaseT : public Super{
public:
	using Super::getUnderlyingGraph;
	unsigned int inDegree(unsigned int i) const {
		return getUnderlyingGraph().inDegree(i);
	}

	auto inAdj(unsigned int i) const ->decltype(getUnderlyingGraph().inAdj(i)) {
		return getUnderlyingGraph().inAdj(i);
	}
};

template<class Symbol, class UnderlyingGraph>
using UndirSymbolAdjListT = UndirSymbolAdjListBaseT < SymbolAdjListT<Symbol, UnderlyingGraph> >;
template<class Symbol>
using UndirSymbolAdjList = UndirSymbolAdjListT<Symbol, UndirAdjList>;

template<class Symbol, class UnderlyingGraph>
using DirSymbolAdjListT = DirSymbolAdjListBaseT < SymbolAdjListT<Symbol, UnderlyingGraph> >;
template<class Symbol>
using DirSymbolAdjList = DirSymbolAdjListT<Symbol, DirAdjList>;

template<class Symbol, class UnderlyingGraph>
using DblDirSymbolAdjListT = DblDirSymbolAdjListBaseT <DirSymbolAdjListT<Symbol, UnderlyingGraph> >;
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
bool isTopologicalOrder(const Graph& reverse, const Graph& graph, const unsigned int* vertices, unsigned int n){
	typedef std::unordered_map<unsigned int, unsigned int> VertexDegrees;
	VertexDegrees vertexDegrees;
	for (auto& v : reverse.getEdges())
		vertexDegrees[v.first] = v.second.size();
	for(unsigned int i = 0;i < n;++ i){
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
		unsigned int v1 = floorTriSqrt1(e);
		unsigned int v2 = e - v1 * (v1 - 1) / 2;
		graph.addEdge(v2, v1);
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

template<class Graph, class Generator>
void randomDAG(unsigned int vertexCount, double ratio, Graph& graph, Generator generator){
	struct Edges{
		std::unordered_map<unsigned int, unsigned int> vertexMap;
		std::vector<unsigned int> vertices;
		std::vector<std::pair<unsigned int, unsigned int> > edges;
		void reset(unsigned int){
		}
		void addEdge(unsigned int v1, unsigned int v2){
			if(v1 <= v2)
				std::swap(v1, v2);

			auto pos = vertexMap.find(v1);
			if(pos == vertexMap.end()){
				auto& v = vertexMap[v1];
			        v = vertices.size();
				vertices.push_back(v1);
				v1 = v;
			}else{
				v1 = pos->second;
			}
			pos = vertexMap.find(v2);
			if(pos == vertexMap.end()){
				auto& v = vertexMap[v2];
			        v = vertices.size();
				vertices.push_back(v2);
				v2 = v;
			}else{
				v2 = pos->second;
			}
			edges.push_back({v1, v2});
		}
	}edges;
	randomUndirGraph(vertexCount, ratio, edges, generator);
	for(auto i = 0;i < edges.vertices.size(); ++ i){
		unsigned int index = i + generator() % (edges.vertices.size() - i);
		std::swap(edges.vertices[index], edges.vertices[i]);
	}
	graph.reset(vertexCount);
	for(auto& edge : edges.edges)
		graph.addEdge(edges.vertices[edge.first], edges.vertices[edge.second]);
}

/*
 * Topological sort.
 * The definition of topological order is given as follows:
 * 1. For directed acyclic graph:
 * 	a permutation of vertices with the property that if there's a directed path from v1 to v2
 * 	then v2 must be after v1 in the permutation.
 * 2. For directed graph contains cycles:
 * 	in this case, there's no topological order actually.
 * 	However the following property of the pertmutation of vertices calculated maybe usefull for some cases
 * 		Let Scc1 and Scc2 be two strong connected component of the graph, 
 * 		if there're edges from Scc1 to Scc2, then there exits a vertex in Scc1 that comes before all the vertices of Scc2 in the permutation.
 * */
class Topological{
protected:
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

		typedef unsigned int Vertex;
		template<class Graph, class CC>
		unsigned int operator()(Vertex source, Vertex root, const Graph& graph, CC& cc) const {
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
	Topological(const Graph& graph, bool reverse = true){
		DirSearchStrategyAdapter<TopoSortStrategy> adapter;
		unsigned int componentCount;
		adapter(graph, adapter.strategy, componentCount);
		topoOrder.swap(adapter.strategy.topoOrder);
		if(reverse)
			std::reverse(topoOrder.begin(), topoOrder.end());
		if(adapter.strategy.cycleFound)
			cycle.swap(adapter.strategy.cycle);
	}

	/*
	 * Returns a directed cycle in the directed graph.
	 * nullptr if no such a cycle exists.
	 * */
	auto directedCycle() const ->decltype(&cycle) {
		if(!cycle.empty())
			return &cycle;
		return nullptr;
	}

	/*
	 * Returns the topological order of the directed graph.
	 * Only meaningful for directed acyclic graph.
	 * 
	 * For directed graph that contains cycles, this order 
	 * can also be used for, say calculating strong connected component.
	 * */
	auto order() const -> decltype(&topoOrder) {
		return &topoOrder;
	}
};

/*
 * Compute the strong connected component of a directed graph using Kosaraju's algorithm.
 * Need to compute the topological order of the reverse graph.
 * I was expecting using the reverse topological order of the original graph had the same effect
 * as the topological of 
 * Tarjan's algorithm is more efficient.
 * */
class StrongConnectedComponent{
	struct SCC : public std::unordered_map<unsigned int, std::pair<unsigned int, unsigned int> >, 
		     public Topological{
		unsigned int componentCount = 0;
		template<class Graph>
		SCC(const Graph& graph) : Topological(graph.toReverse()){
			DirDfsStrategy strategy;
			for(auto vertex : topoOrder){
				if((*this)(vertex)){
					//cc(vertext, parent of vertex, root of the tree where the vertex lives in)
					(*this)(vertex, vertex, vertex);
					strategy(vertex, vertex, graph, *this);
					++componentCount;
				}
			}
			std::reverse(topoOrder.begin(), topoOrder.end());
		}
		bool operator()(unsigned int v) const {
			return find(v) == end();
		}
		bool operator()(unsigned int s, unsigned int v) const {
			return find(v) == end();
		}
		void operator()(unsigned int s, unsigned int p, unsigned int r){
			auto& dat = (*this)[s];
		        dat.first = p;
			dat.second = r;
		}
	}scc;
public:
	template<class Graph>
	StrongConnectedComponent(const Graph& graph) : scc(graph){
	}

	unsigned int getComponentCount(){
		return scc.componentCount;
	}

	template<class T>
	void getSCC(T& t){
		for(auto& kvp : scc){
			t[kvp.second.second].insert(kvp.first);
		}
	}
	auto directedCycle() const ->decltype(scc.directedCycle()) {
		return scc.directedCycle();
	}

	/*
	 * Topological order for directed acyclic graph.
	 * */
	auto order() const -> decltype(scc.order()) {
		return scc.order();
	}
};

/*
 * Calculate euler cycle for each connected component of an undirected graph if exists.
 */

/*
 * Tarjan's algorithm
 * */
class TarjanStrategy{
	typedef std::vector<std::pair<unsigned int, unsigned int> > Components;
	Components components;
public:
	template<class Graph>
	TarjanStrategy(const Graph& graph){
		struct Data{
			//<vertex, <low link, is on stack> >
			std::unordered_map<unsigned int, std::pair<unsigned int, bool> >  links;
			std::vector<unsigned int> stack;
		};
		struct Dfs{
			const Graph& graph;
			Data& data;
			unsigned int& nextComponent;
			unsigned int& index;
			Components& components;
			Dfs(const Graph& graph, Data& data, Components& components, unsigned int& nextComponent, unsigned int& index)
				: graph(graph), data(data), nextComponent(nextComponent), components(components), index(index){
			}
			unsigned int operator()(unsigned int , unsigned int v){
				auto link = &data.links[v];
				auto stackIndex = data.stack.size();
				auto vertexIndex = index;
				link->first = index;
				link->second = true;
				data.stack.push_back(v);
				++index;

				auto adj = graph.outAdj(v);
				if(adj != nullptr)
					for(auto a : *adj){
						auto pos = data.links.find(a);
						if(pos == data.links.end()){
							unsigned int lowLink = (*this)(v, a);
							//in case a rehash happens,
							//the original pointer becomes invalid.
							link = &data.links[v];
							if(lowLink < link->first) link->first = lowLink;
						}else if(pos->second.second && pos->second.first < link->first){
							link->first = pos->second.first;
						}
					}

				//in case a rehash happens, the original pointer becomes invalid.
				link = &data.links[v];
				if(link->first == vertexIndex){
					auto pos = data.stack.begin() + stackIndex, end = data.stack.end();
					for(auto begin = pos;begin != end;++begin){
						components.push_back({*begin, nextComponent});
						data.links[*begin].second = false;
					}
					data.stack.erase(pos, end);
					++ nextComponent;
				}
				return link->first;
			}
		};
		Data data;
		unsigned int nextComponent = 0, index = 0;
		for (auto& vertices : graph.getEdges()){
			Dfs dfs(graph, data, components, nextComponent, index);
			if(data.links.find(vertices.first) == data.links.end())
				dfs(vertices.first, vertices.first);
			assert(data.stack.empty());
		}
	}
	template<class T>
	void getSCC(T& t){
		for(auto& c : components){
			t[c.second].insert(c.first);
		}
	}
};

/*
* Compute least commmon ancestor of two vertices in a DAG.
* */
class LeastCommonAncestor{
	std::unordered_map<unsigned int, unsigned int> reachability;
	std::unordered_set<unsigned int> lcas;
	unsigned int v1, v2;
	template<class Graph>
	unsigned int visit(const Graph& graph, unsigned int s){
		auto r = reachability[s] = (s == v1 ? 1 : (s == v2 ? 2 : 0));

		auto adj = graph.outAdj(s);
		auto isLeast = true;
		if(adj != nullptr){
			for(auto v : *adj){
				auto pos = reachability.find(v);
				if(pos == reachability.end()){
					auto c = visit(graph, v);
					isLeast = isLeast && (c != 3);
					r |= c;
				}else{
					isLeast = isLeast && (pos->second != 3);
					r |= pos->second;
				}
			}
		}
		reachability[s] = r;
		if(r == 3 && isLeast)
			lcas.insert(s);
		return r;
	}
public:
	template<class Graph>
	LeastCommonAncestor(const Graph& graph, unsigned int v1, unsigned int v2) : v1(v1), v2(v2){
		if(v1 == v2){
			lcas.insert(v1);
			return;
		}
		for(auto& vertex : graph.getEdges()){
			if(reachability.find(vertex.first) == reachability.end()){
				visit(graph, vertex.first);
			}
		}
	}
	const std::unordered_set<unsigned int>& lca() const {
		return lcas;
	}
};

/*
 * Algorithm to calculate if a directed acyclic graph contains a Hamiltonian path
 * return true if contains.
 * throws std::invalid_argument if the given graph is not a directed acyclic graph.
 * */
template<class DAG>
bool dagContainsHamiltonianPath(const DAG& dag){
	Topological topological(dag);
	if(topological.directedCycle() != nullptr)
		throw std::invalid_argument("given graph is not a directed acyclic graph");
	auto order = topological.order();
	for(unsigned int i = 0;i + 1< order->size();++ i)
		if(!dag.isAdj((*order)[i], (*order)[i + 1]))
			return false;
	return true;
}
/*
 * return true if the given directed graph contains Euler cycle
 */
template<class Graph>
bool containsEulerianCycle(const Graph& graph){
	std::unordered_map<unsigned int, int> degDiff;
	for(auto& vertex : graph.getEdges()){
		degDiff[vertex.first] -= vertex.second.size();
		for(auto& v : vertex.second){
			++ degDiff[v];
		}
	}
	for(auto& kvp : degDiff)
		if(kvp.second != 0){
			return false;
		}
	return true;
}
/*
 * returns true if the strong connected components of the graph contains Euler cycle
 * And calculates the Euler cycle of each strong connected component.
 */
template<class Graph, class Path>
bool directedEulerianCycle(const Graph& graph, Path& path){
	typedef decltype(graph.getEdges().begin()->second.begin()) AdjRange;
	std::vector<std::pair<unsigned int, unsigned int> > paths;
	std::unordered_map<unsigned int, std::pair<AdjRange, AdjRange> > remainingEdges;
	std::unordered_map<unsigned int, unsigned int> vertexToLink;
	for(auto& vertex : graph.getEdges()){
		assert(!vertex.second.empty());
		if(!vertex.second.empty()){
			remainingEdges[vertex.first] = {vertex.second.begin(), vertex.second.end()};
			vertexToLink[vertex.first] = -1;
		}
	}

	//start link of every eurlerian cycle
	std::vector<unsigned int> begins;
	while(!vertexToLink.empty()){
		auto startVertex = *vertexToLink.begin();
		unsigned int startLink = paths.size();

		auto current = remainingEdges.find(startVertex.first);
		//loop invariant : there's unvisited edges.
		while(current != remainingEdges.end() && current->second.first != current->second.second){
			unsigned int curLink = paths.size();
			if(curLink != startLink)
				paths[curLink - 1].second = curLink;
			paths.push_back(std::make_pair(current->first, (unsigned int)-1));
			vertexToLink[current->first] = curLink;

			auto nextVertex = *current->second.first;
			++current->second.first;
			if(current->second.first == current->second.second)
				vertexToLink.erase(current->first);
			current = remainingEdges.find(nextVertex);
		}
		
		//1. if we end at a vertex with out degree zero, no Eurlerian cycle.
		//2. we end at a vertex that's not the start vertex, no Eulerian cycle.
		if(current == remainingEdges.end() || current->first != startVertex.first)
			return false;
		vertexToLink.erase(startVertex.first);
		//Note that since we're not assming self-loops, so a cycle must contains at least 2 vertices.
		assert(paths.size() > startLink && paths.size() - startLink > 1);
		if(startVertex.second == -1){
			begins.push_back(startLink);
		}else{
			unsigned int next = paths[startVertex.second].second;
			paths[startVertex.second].second = paths[startLink].second;
			paths.back().second = startLink;
			paths[startLink].second = next;
		}
	}
	for(auto begin : begins){
		while(begin != -1){
			path.push_back(paths[begin].first);
			begin = paths[begin].second;
		}
		path.push_back(-1);
	}

	return true;
}

template<class Vertex>
class EdgeBase{
protected:
	Vertex v1, v2;
public:
	EdgeBase(Vertex v1, Vertex v2) : v1(v1), v2(v2){
	}

	bool operator==(const EdgeBase& e) const {
		return v1 == e.v1 && v2 == e.v2;
	}

	bool operator<(const EdgeBase& e) const {
		if(v1 == e.v1)
			return v2 < e.v2;
		return v1 < e.v1;
	}

	Vertex getV1() const {
		return v1;
	}

	Vertex getV2() const {
		return v2;
	}

	struct Hash{
		size_t operator()(const EdgeBase& e) const {
			return e.v1 * 31 + e.v2;
		}
	};

	template<class Stream>
	friend Stream& operator<<(Stream& os, const EdgeBase& e){
		os << '(' << e.getV1() << ',' << e.getV2() << ')';
	}
};

template<class Vertex>
class UndirEdge : public EdgeBase<Vertex>{
	typedef EdgeBase<Vertex> Super;
public:
	UndirEdge(Vertex v1, Vertex v2) : Super(v2 < v1 ? Super(v2, v1) : Super(v1, v2)){
	}
};

template<class Vertex>
class DirEdge : public EdgeBase<Vertex>{
public:
	DirEdge(Vertex v1, Vertex v2) : EdgeBase<Vertex>(v1, v2){
	}
};

template<class Edge, class Weight = double>
class WeightedEdgeT : public Edge{
	Weight weight;
public:
	template<class Vertex>
	WeightedEdgeT(Vertex v1, Vertex v2, Weight w) : Edge(v1, v2), weight(w){
	}

	Weight getWeight() const {
		return weight;
	}

	template<class Stream>
	friend Stream& operator<<(Stream& os, const WeightedEdgeT& e){
		return os << '(' << e.getV1() << ',' << e.getV2()  << ',' << e.getWeight() << ')';
	}
};

template<class Vertex, class Weight = double>
using WeightedUndirEdgeT = WeightedEdgeT<UndirEdge<Vertex>, Weight>;

template<class Vertex, class Weight = double>
using WeightedDirEdgeT = WeightedEdgeT<DirEdge<Vertex>, Weight>;

template<class Vertex>
using WeightedUndirEdge = WeightedUndirEdgeT<Vertex>;

template<class Vertex>
using WeightedDirEdge = WeightedDirEdgeT<Vertex>;

template<class UnderlyingAdjList, class Weight, template<class V, class T> class WeightedEdgeT>
class WeightedAdjListT : protected UnderlyingAdjList{
	typedef UnderlyingAdjList Super;
	Weight minimum;
public:
	typedef typename UnderlyingAdjList::Vertex Vertex;
	typedef WeightedEdgeT<Vertex, Weight> WeightedEdge;
	//typedef std::unordered_set<WeightedEdge, EdgeBase::Hash> WeightedEdgeSet;
	typedef std::set<WeightedEdge> WeightedEdgeSet;
	WeightedAdjListT(unsigned int vertexCount = 0, Weight minimum = Weight()) : Super(vertexCount), minimum(minimum){
	}

	void addEdge(Vertex i, Vertex j, Weight weight){
		if (weight < minimum){
			std::stringstream os;
			os << "weight " << weight << " is smaller than the minimum " << minimum << ", line " << __LINE__ << " file " << __FILE__;
			throw std::invalid_argument(os.str());
		}
		Super::addEdge(i, j);
		weightedEdges.insert({ i, j, weight });
	}

	template<class Edge>
	void addEdge(const Edge& edge){
		addEdge(edge.getV1(), edge.getV2(), edge.getWeight());
	}

	void reset(unsigned int vertexCount){
		Super::reset(vertexCount);
		weightedEdges.clear();
	}

	const WeightedEdgeSet& getWeightedEdges() const {
		return weightedEdges;
	}

	const WeightedEdge* getWeightedEdge(Vertex i, Vertex j) const {
		auto pos = weightedEdges.find(WeightedEdge(i, j, Weight()));
		if (pos == weightedEdges.end())
			return nullptr;
		return &*pos;
	}

	Weight getWeight(Vertex i, Vertex j) const {
		auto weightedEdge = getWeightedEdge(i, j);
		return nullptr == weightedEdge ? Weight() : weightedEdge->getWeight();
	}

	template<class Stream>
	friend Stream& operator<<(Stream& stream, const WeightedAdjListT& graph){
		for (auto& edge : graph.getWeightedEdges()){
			stream << edge.getV1() << ' ' << edge.getV2() << "\t\t: " << edge.getWeight() << "\n";
		}
		return stream;
	}

	using Super::isAdj;
	using Super::getVertexCount;
	using Super::getEdges;
	using Super::extend;
protected:
	WeightedEdgeSet weightedEdges;
};

template<class Weight = double, template<class U, class V> class Map = DefaultStdMap, template<class U> class List = DefaultStdSet >
class WeightedUndirAdjListT : public WeightedAdjListT<UndirAdjListT<Map, List>, Weight, WeightedUndirEdgeT>{
	typedef WeightedAdjListT<UndirAdjListT<Map, List>, Weight, WeightedUndirEdgeT> Super;
	typedef UndirAdjListT<Map, List> GrandFather;
public:
	WeightedUndirAdjListT(unsigned int vertexCount = 0, Weight minimum = Weight()) : Super(vertexCount, minimum){
	}
	using Super::degree;
	using GrandFather::adj;
};

using WeightedUndirAdjList = WeightedUndirAdjListT<>;

template<class Weight = double, template<class U, class V> class Map = DefaultStdMap, template<class U> class List = DefaultStdSet >
class WeightedDirAdjListT : public WeightedAdjListT<DirAdjListT<Map, List>, Weight, WeightedDirEdgeT>{
	typedef WeightedAdjListT<DirAdjListT<Map, List>, Weight, WeightedDirEdgeT> Super;
	typedef DirAdjListT<Map, List> GrandFather;
public:
	WeightedDirAdjListT(unsigned int vertexCount = 0, Weight minimum = Weight()) : Super(vertexCount, minimum){
	}
	using GrandFather::outAdj;
	using GrandFather::outDegree;
	using GrandFather::isMap;
};

using WeightedDirAdjList = WeightedDirAdjListT<>;
template<class Symbol, class UnderlyingGraph>
class WeightedSymbolAdjListT : public SymbolAdjListBaseT<Symbol, UnderlyingGraph>{
	typedef SymbolAdjListBaseT<Symbol, UnderlyingGraph> Super;
public:
	template<class Weight>
	void addEdge(const Symbol& s1, const Symbol& s2, const Weight& weight) throw(std::invalid_argument){
		if (s1 == s2){
			std::stringstream os;
			os << "self loop for vertex " << s1 << " , " << s2 << " not allowed, line " << __LINE__ << " file " << __FILE__;
			throw std::invalid_argument(os.str());
		}
		auto index1 = Super::symbolMap[s1];
		auto index2 = Super::symbolMap[s2];
		assert(Super::symbolArray.size() >= Super::getUnderlyingGraph().getVertexCount());
		assert(index1.isIndex && index2.isIndex);
		Super::getUnderlyingGraph().extend(Super::symbolArray.size() - Super::getUnderlyingGraph().getVertexCount());
		return Super::getUnderlyingGraph().addEdge(index1.key.index, index2.key.index, weight);
	}
};

template<class Symbol, class UnderlyingGraph>
using WeightedUndirSymbolAdjListT = UndirSymbolAdjListBaseT < WeightedSymbolAdjListT<Symbol, UnderlyingGraph> >;
template<class Symbol>
using WeightedUndirSymbolAdjList = WeightedUndirSymbolAdjListT<Symbol, WeightedUndirAdjList>;

template<class Symbol, class UnderlyingGraph>
using WeightedDirSymbolAdjListT = DirSymbolAdjListBaseT < WeightedSymbolAdjListT<Symbol, UnderlyingGraph> >;
template<class Symbol>
using WeightedDirSymbolAdjList = WeightedDirSymbolAdjListT<Symbol, WeightedDirAdjList>;

/*
template<class Symbol, class UnderlyingGraph>
using WeightedDblDirSymbolAdjListT = DblDirSymbolAdjListBaseT <DirSymbolAdjListT<Symbol, UnderlyingGraph> >;
template<class Symbol>
using WeightedDblDirSymbolAdjList = WeightedDblDirSymbolAdjListT<Symbol, WeightedDblDirAdjList>;*/

struct ReadOneWeightedEdge{
	template<class Stream, class Graph>
	bool operator()(Stream& stream, Graph& graph) const {
		int v1, v2;
		double weight;
		if(stream >> v1 >> v2 >> weight){
			graph.addEdge(v1, v2, weight);
			return true;
		}
		return false;
	}
};

using TxtFileWeightedGraphReader = TxtFileGraphReaderT<TxtStreamGraphReaderT<std::ifstream, ReadOneWeightedEdge> >;

using TxtFileStreamWeightedGraphReader= TxtFileStreamGraphReaderT<TxtStreamGraphReaderT<std::ifstream, ReadOneWeightedEdge> >;

/*
 * Prim's algorithm for calculating minimum spanning tree of weighted undirected graph.
 * Time complexity : O(E*log(E))
 * 		     E is the number of edges
 * */
class PrimMST{
protected:
	//vertex to parent map.
	std::unordered_map<unsigned int, unsigned int> mst;
	struct CompareWeight{
		template<class T>
		bool operator()(const T *e1, const T *e2) const{
			return e1->getWeight() > e2->getWeight();
		}
	};

	template<class Graph, template<class Edge> class PriorityQueue>
	void calculate(const Graph& graph){
		for (auto& vertex : graph.getEdges()){
			if (mst.find(vertex.first) == mst.end()){
				mst[vertex.first] = vertex.first;
				//heap ordered by edge weight.
				typedef const typename Graph::WeightedEdge* WEP;
				//std::priority_queue<WEP, std::vector<WEP>, CompareWeight> pq;
				PriorityQueue<WEP> pq;
				for (auto v : vertex.second){
					auto weightedEdge = graph.getWeightedEdge(vertex.first, v);
					assert(weightedEdge != nullptr);
					//std::cout << "pushing " << *weightedEdge << std::endl;
					pq.push(v, weightedEdge);
				}
				while (!pq.empty()){
					//pick up the edge with minimum weight.
					auto e = pq.pop();
					auto pos = mst.find(e->getV1());
					auto vertexNotInTree = e->getV1();
					auto vertexInTree = e->getV2();
					//std::cout << "using " << *e << std::endl;
					//each edge in the heap must have at one vertex already in the tree.
					//so if e->v1 is not in the tree, e->v2 must be in the tree.
					assert(pos != mst.end() || mst.find(e->getV2()) != mst.end());
					if (pos != mst.end()){
						//e->v1 already in the tree

						//if e->v2 already in the tree
						//this edge is ineligible.
						if (mst.find(e->getV2()) != mst.end())
							continue;

						//e->v2 is not in the tree
						vertexNotInTree = e->getV2();
						vertexInTree = e->getV1();
					}
					auto adj = graph.adj(vertexNotInTree);
					mst[vertexNotInTree] = vertexInTree;
					if (adj != nullptr)
						for (auto v : *adj)
							if (mst.find(v) == mst.end()){
								auto weightedEdge = graph.getWeightedEdge(vertexNotInTree, v);
								assert(weightedEdge != nullptr);
								//std::cout << "pushing " << *weightedEdge << std::endl;
								pq.push(v, weightedEdge);
							}
				}
			}
		}
	}

	template<class Edge>
	class PriorityQueueAdapter{
		std::priority_queue<Edge, std::vector<Edge>, CompareWeight> pq;
	public:
		Edge pop(){
			auto edge = pq.top();
			pq.pop();
			return edge;
		}
		template<class Vertex>
		void push(Vertex v, Edge edge){
			pq.push(edge);
		}
		bool empty(){
			return pq.empty();
		}
	};

	PrimMST(){
	}
public:
	template<class Graph>
	PrimMST(const Graph& graph){
		calculate<Graph, PriorityQueueAdapter>(graph);
	}

	template<class Graph, class T>
	void getEdges(const Graph& graph, T& edges){
		for(auto& e: mst){
			if(e.first != e.second){
				auto edge = graph.getWeightedEdge(e.first, e.second);
				if(edge != nullptr)
					edges.push_back(edge);
			}
		}
	}

	template<class Graph>
	auto getMSTWeight(const Graph& graph) const -> decltype(graph.getWeightedEdge(0, 0)->getWeight()) {
		typedef decltype(getMSTWeight(graph)) Weight;
		Weight sum = Weight();
		for(auto& e: mst){
			if(e.first != e.second){
				auto edge = graph.getWeightedEdge(e.first, e.second);
				if(edge != nullptr)
					sum = sum + edge->getWeight();
			}
		}
		return sum;
	}
	
	bool operator==(const PrimMST& mst){
		return this->mst == mst.mst;
	}
};

//a specialized indexed heap for weighted edge.
template<class WeightedEdge>
struct IndexedMinWeightedEdgeHeap{
	typedef std::unordered_map<unsigned int, std::pair<WeightedEdge, unsigned int> > Maps;
	typedef std::vector<unsigned int> Vertices;
	static const int d = 2;

	//Not a real heap push operation. :-)
	//For our special case, only a smaller weighted edge could be eligible.
	//So do nothing when not this case.
	bool push(unsigned int vertex, WeightedEdge edge){
		auto pos = maps.find(vertex);
		if (pos == maps.end()){
			maps[vertex] = { edge, vertices.size() };
			vertices.push_back(vertex);
			promote(vertices.size() - 1);
			return true;
		}
		else if (edge->getWeight() < pos->second.first->getWeight()){
			pos->second.first = edge;
			promote(pos->second.second);
			return true;
		}
		return false;
	}
	bool get(unsigned int vertex, WeightedEdge& edge){
		auto pos = maps.find(vertex);
		if (pos == maps.end())
			return false;
		edge = pos->second.first;
		return true;
	}
	WeightedEdge pop(){
		auto minVertex = maps.find(vertices[0]);
		assert(minVertex != maps.end());
		WeightedEdge min = minVertex->second.first;
		maps.erase(minVertex);
		if (vertices.size() > 1){
			auto s = vertices.size() - 1;
			auto vertex = maps.find(vertices[s]);
			assert(vertex != maps.end());
			vertex->second.second = 0;
			vertices[0] = vertices[s];
			vertices.pop_back();
			heapify(0);
		}
		else{
			vertices.pop_back();
		}
		return min;
	}

	std::pair<unsigned int, WeightedEdge> popWithVertex(){
		auto vertex = vertices[0];
		return {vertex, pop()};
	}

	bool empty(){
		return vertices.empty();
	}
	void promote(size_t index) {
		if (index > 0){
			size_t parent = (index - 1) / d;
			auto parentVertex = maps.find(vertices[parent]);
			auto vertex = maps.find(vertices[index]);
			assert(parentVertex != maps.end());
			assert(vertex != maps.end());
			while (vertex->second.first->getWeight() < parentVertex->second.first->getWeight()){
				std::swap(vertex->second.second, parentVertex->second.second);
				std::swap(vertices[index], vertices[parent]);
				index = parent;
				if (index == 0)
					break;
				parent = (index - 1) / d;
				//vertex = parentVertex;
				parentVertex = maps.find(vertices[parent]);
				assert(parentVertex != maps.end());
			}
		}
	}
	void heapify(size_t index){
		auto vertex = maps.find(vertices[index]);
		auto s = vertices.size();
		assert(vertex != maps.end());
		while (index < s){
			size_t smallest = index;
			auto smallestVertex = vertex;
			if (d != 2){
				size_t start = d * index + 1;
				size_t end = std::min(s, start + d);
				for (size_t i = start; i < end; ++i){
					auto childVertex = maps.find(vertices[i]);
					assert(childVertex != maps.end());
					if (childVertex->second.first->getWeight() < smallestVertex->second.first->getWeight()){
						smallest = i;
						smallestVertex = childVertex;
					}
				}
			}else{
				//avoiding use for loop for binary heap, so that it would be a little faster
				size_t l = d * index + 1;
				size_t r = l + 1;
				if (l < s){
					auto childVertex = maps.find(vertices[l]);
					assert(childVertex != maps.end());
					if (childVertex->second.first->getWeight() < smallestVertex->second.first->getWeight()){
						smallest = l;
						smallestVertex = childVertex;
					}
				}
				if (r < s){
					auto childVertex = maps.find(vertices[r]);
					assert(childVertex != maps.end());
					if (childVertex->second.first->getWeight() < smallestVertex->second.first->getWeight()){
						smallest = r;
						smallestVertex = childVertex;
					}
				}
			}
			if (smallest == index)
				break;
			std::swap(vertex->second.second, smallestVertex->second.second);
			std::swap(vertices[index], vertices[smallest]);
			index = smallest;
			//vertex = smallestVertex;
		}
	}
private:
	Maps maps;
	Vertices vertices;
};
/*
 * Prim's algorithm eager version for calculating minimum spanning tree of weighted undirected graph.
 * Time complexity : O(E*log(V))
 * 		     E is the number of edges
 * 		     V is the number of vertices
 * */
class EagerPrimMST : public PrimMST{
public:
	template<class Graph>
	EagerPrimMST(const Graph& graph){
		calculate<Graph, IndexedMinWeightedEdgeHeap>(graph);
	}
};

class ShortestPathBase{
protected:
	//vertex to parent map.
	std::unordered_map<unsigned int, unsigned int> sp;
public:
	template<class Path>
	bool getPath(unsigned int vertex, Path *path){
		auto pos = sp.find(vertex);
		if (pos == sp.end())
			return false;
		if (path != nullptr){
			while (pos != sp.end()){
				path->push_back(pos->second);
				pos = sp.find(pos->second);
			}
		}
		std::reverse(path->begin(), path->end());
		path->push_back(vertex);
		return true;
	}
};
/*
 * Dijkstra's algorithm for computing shortest path of a weighted directed graph with non-negative weights. 
 * When edge's have negative weights, the behaviour of this algorithm is undefined.
 * Edge weights must be non-negative
 * Time complexity: O(E + V*log(V))
 * 		    With  E is the number of edge, V is the number of vertices.
 * */
class DijkstraShortestPath : public ShortestPathBase{
	unsigned int src;
public:
	template<class Graph>
	DijkstraShortestPath(const Graph& graph, unsigned int src){
		struct ShortestDist{
			unsigned int p;
			double weight;
			ShortestDist(unsigned int p = -1, double weight = 0) : weight(weight), p(p){
			}
			ShortestDist* operator->(){
				return this;
			}
			double getWeight(){
				return weight;
			}
		};
		IndexedMinWeightedEdgeHeap<ShortestDist> pq;
		pq.push(src, ShortestDist(src, 0));

		//the number of elements in the priority queue is at most V
		//and with one element poped from the priority queue each loop
		//so the number of loops is at most V
		while (!pq.empty()){
			
			//pop from the priority queue takes O(log(V)) time
			auto shortestDist = pq.popWithVertex();
			auto adj = graph.outAdj(shortestDist.first);
			sp[shortestDist.first] = shortestDist.second.p;

			if(adj)
				for(auto v : *adj)
					//every edge will be visited at most once in this inner loop
					//so the overall total cost of this inner loop is O(E)
					if(sp.find(v) == sp.end()){
						auto weightedEdge = graph.getWeightedEdge(shortestDist.first, v);
						assert(weightedEdge != nullptr);
						ShortestDist newShortestDist(shortestDist.first, shortestDist.second.getWeight() + weightedEdge->getWeight());
						pq.push(v, newShortestDist);
					}
		}
		sp.erase(src);
		this->src = src;
	}
};

/*
 * Shortest path in directed acyclic graph.
 * Edge weights can be arbitrary
 * Time complexity: O(E)
 * 		    E is the number of edge.
 * */
class DAGShortestPath : public ShortestPathBase{
public:
	template<class Graph>
	DAGShortestPath(const Graph& graph) throw(std::invalid_argument) {
		Topological topological(graph);
		if(topological.directedCycle() != nullptr)
			throw std::invalid_argument("given graph is not a directed acyclic graph");
		auto order = topological.order();
		std::unordered_map<unsigned int, double> dists;
		for(auto src : *order){
			auto adj = graph.outAdj(src);
			if(adj)
				for(auto v : *adj){
					auto weightedEdge = graph.getWeightedEdge(src, v);
					assert(weightedEdge != nullptr);
					auto pos = dists.find(v);
					auto newDist = dists[src] + weightedEdge->getWeight();
					if(pos == dists.end() || newDist < pos->second){
						dists[v] = newDist;
						sp[v] = src;
					}
				}
		}
	}
};

/*
 * BellmanFord algorithm to compute shortest path in a directed graph.
 * Edge weights can be arbitrary.
 * Worst case time complexity: O(V*E) 
 * General case time complexity: O(E + V)
 * 			       V is the number of vertices
 * 			       E is the number of edges
 */
class BellmanFordShortestPath : public ShortestPathBase{
	unsigned int src;
	std::vector<unsigned int> cycle;
public:
	template<class Graph>
	BellmanFordShortestPath(const Graph& graph, unsigned int src){
		std::unordered_map<unsigned int, double> dists;
		std::unordered_set<unsigned int> inQueue;
		std::priority_queue<unsigned int> vertices;
		/*shortest path tree used to calculate if a cycle exists*/
		class EdgeRemovableAdjList : public DirAdjList{
		public:
			EdgeRemovableAdjList(unsigned int vertexCount) : DirAdjList(vertexCount){
			}
			bool removeEdge(unsigned int v1, unsigned int v2){
				auto edges = getEdges();
				auto pos = edges.find(v1);
				if(pos == edges.end())
					return false;
				auto pos1 = pos->second.find(v2);
				if(pos1 == pos->second.end())
					return  false;
				pos->second.erase(pos1);
				return true;
			}
		}spt(graph.getVertexCount());
		vertices.push(src);

		unsigned int count = 0;
		while(!vertices.empty()){
			unsigned int src = vertices.top();
			auto parent = sp.find(src);
			vertices.pop();
			inQueue.erase(src);
			if (parent != sp.end() && inQueue.find(parent->second) != inQueue.end())
				continue;
			auto adj = graph.outAdj(src);


			if(adj)
				for(auto v : *adj){
					auto weightedEdge = graph.getWeightedEdge(src, v);
					assert(weightedEdge != nullptr);
					auto pos = dists.find(v);
					auto newDist = dists[src] + weightedEdge->getWeight();
					if(pos == dists.end()){
						dists[v] = newDist;
						spt.addEdge(src, v);
						sp[v] = src;
						assert(inQueue.find(v) == inQueue.end());
						inQueue.insert(v);
						vertices.push(v);
					}else if(newDist < pos->second){
						auto& parent = sp[v];
						bool isRemoved = spt.removeEdge(parent, v); 
						assert(isRemoved);
						dists[v] = newDist;
						spt.addEdge(src, v);
						parent = src;
						if(inQueue.find(v) == inQueue.end()){
							inQueue.insert(v);
							vertices.push(v);
						}
					}
				}
			if(++count % spt.getVertexCount() == 0){
				Topological topological(spt);
				auto cycle = topological.directedCycle();
				if(cycle != nullptr){
					this->cycle = *cycle;
					break;
				}
			}
		}
		this->src = src;
	}
	bool hasNegativeCycle(){
		return !cycle.empty();
	}
	
	auto negativeCycle() const ->decltype((cycle)){
		return cycle;
	}
};

template<class Graph, class WeightGenerator>
struct RandomWeightedGraphGenerationAdapter{
	Graph& graph;
	WeightGenerator& generator;
	RandomWeightedGraphGenerationAdapter(Graph& graph, WeightGenerator& generator) : graph(graph), generator(generator){
	}
	void reset(unsigned int vertexCount){
		graph.reset(vertexCount);
	}
	void addEdge(unsigned int v1, unsigned int v2){
		graph.addEdge(v1, v2, generator());
	}
};

template<class Graph, class Generator, class WeightGenerator>
void randomWeightedUndirGraph(unsigned int vertexCount, double ratio, Graph& graph,
	       			Generator generator, WeightGenerator weightGenerator){
	RandomWeightedGraphGenerationAdapter<Graph, WeightGenerator>  rwgga(graph, weightGenerator);
	randomUndirGraph(vertexCount, ratio, rwgga, generator);
}

template<class Graph, class Generator, class WeightGenerator>
void randomWeightedDirGraph(unsigned int vertexCount, double ratio, Graph& graph,
	       			Generator generator, WeightGenerator weightGenerator){
	RandomWeightedGraphGenerationAdapter<Graph, WeightGenerator>  rwgga(graph, weightGenerator);
	randomDirGraph(vertexCount, ratio, rwgga, generator);
}

template<class Graph, class Generator, class WeightGenerator>
void randomWeightedDAG(unsigned int vertexCount, double ratio, Graph& graph,
	       			Generator generator, WeightGenerator weightGenerator){
	RandomWeightedGraphGenerationAdapter<Graph, WeightGenerator>  rwgga(graph, weightGenerator);
	randomDAG(vertexCount, ratio, rwgga, generator);
}

template<class UnderlyingGraph, class Coordinate>
class EuclideanGraphT : public UnderlyingGraph{
	std::unordered_map<unsigned int, Coordinate> coordinates;
public:
	EuclideanGraphT(unsigned int vertexCount = 0) : UnderlyingGraph(vertexCount){
	}

	const Coordinate* getCoordinate(unsigned int v) const {
		auto pos = coordinates.find(v);
		if(pos == coordinates.end())
			return nullptr;
		return &pos->second;
	}

	auto getCoordinates() const -> decltype((coordinates)){
		return coordinates;
	}

	void setCoordinate(unsigned int v, const Coordinate& coordinate) {
		UnderlyingGraph::check(v);
		coordinates[v] = coordinate;
	}

	void clearCoordinates(){
		coordinates.clear();
	}
};

namespace{
/*
 * this class is expected only to be used in thhis namespace.
 * No-way to make something private to a namespace?
 */
template<class DirGraph, class Vertex>
struct private_IsDAG{
	private_IsDAG(const DirGraph& graph) : graph(graph){
		for(const auto& vertices : graph.getEdges()) for(const auto& v : vertices.second) ++inDegrees[v];
	}
	/*
	 * DAG : directed acyclic graph.
	 * A DAG must contain at least one vertex whose in-degree is zero.
	 * The Directed Graph generated by removing one vertex from a DAG is also a DAG.
	 * Assuming that the number of vertices of a DAG is finite, then by iteratively removing 
	 * vertices whose in-dgree are zero in a DAG and the subsequently generated DAG, all
	 * the vertices will be removed. 
	 *
	 * Time complexity of thi function is O(E), space complexity is O(V)
	 * where E is the number of edges and V is the number of vertices.
	 */
	operator bool() const {
		std::vector<Vertex> toRemove;
		for(const auto& vertices : graph.getEdges()) 
			if(inDegrees.find(vertices.first) == inDegrees.end())
				toRemove.push_back(vertices.first);

		while(!toRemove.empty()){
			for(const auto& v : toRemove){
				auto outAdj = graph.outAdj(v);
				if(decltype(outAdj)() != outAdj)
					for(auto v : *outAdj) --inDegrees[v];
			}
			toRemove.clear();
			for(const auto& v : inDegrees) if(v.second == 0) toRemove.push_back(v.first);
			for(const auto& v : toRemove) inDegrees.erase(v);
	 	}
	 	return inDegrees.empty();
 	}
protected:
 	mutable std::map<Vertex, unsigned int> inDegrees; const DirGraph& graph;
};

}

template<class DirGraph>
bool isDAG(const DirGraph& graph){
	//Relying on typedef in class, instead of the following complex expression:
	//  typedef typename std::remove_cv<decltype(graph.getEdges().begin()->first)>::type Vertex;
	return private_IsDAG<DirGraph, typename DirGraph::Vertex>(graph);
}

/*
 * A ST-network is a DAG with only one vertex whose in-degree is zero, called the source, 
 * and only one vertex whose out-degree is zero, called the sink.
 * This function checks to see is the given graph is a ST-network.
 * */
template<class DirGraph, class Vertex>
bool isSTNetWork(const DirGraph& graph, Vertex &source, Vertex& sink){
	struct IsSTNetWork : public private_IsDAG<DirGraph, Vertex>{
		IsSTNetWork(const DirGraph& graph) : private_IsDAG<DirGraph, Vertex>(graph){}
		bool getSrcAndSink(Vertex& source, Vertex &sink){
			const auto& inDegrees = private_IsDAG<DirGraph, Vertex>::inDegrees;
			const auto& edges = private_IsDAG<DirGraph, Vertex>::graph.getEdges();
		 	if(inDegrees.size() != edges.size()) return false;
			auto sourcePos = edges.end();
			{
				auto end = sourcePos;
				for(auto begin = edges.begin();begin != end;++begin)
					if(inDegrees.find(begin->first) == inDegrees.end()){
						if(sourcePos != end) return false;
						sourcePos = begin;
					}
				if(sourcePos == end) return false;
			}
			
			auto sinkPos =  inDegrees.end();
			{
				auto end = sinkPos;
				for(auto begin = inDegrees.begin();begin != end;++begin)
					if(edges.find(begin->first) == edges.end()){
						if(sinkPos != end) return false;
						sinkPos = begin;
					}
				if(sinkPos == end) return false;
			}
			source = sourcePos->first;
			sink = sinkPos->first;
			return true;
		}
	}is(graph);
	return is.getSrcAndSink(source, sink) && is;
} 

/*
 * Find the max flow for a ST-network.
 * Worst case time complexity is O(E*E*V)
 * The {@code flow[dest][src]} is the weight(or flow) choosen
 * for the edge (src, dest).
 * */
template<class WeightedDirGraph, class STFlow>
void STNetworkMaxFlow(const WeightedDirGraph& graph, STFlow& flow){
	// We're assuming a ST-network.
	typedef typename WeightedDirGraph::Vertex Vertex;
	Vertex source, sink;
	if(!isSTNetWork(graph, source, sink)) throw std::invalid_argument("not a ST-network, this algorithm works only for ST-network");

	if(source == sink) return;

	//std::cout << "source vertex is : " << source << ", sink vertex is : " << sink << std::endl;

	typedef typename std::remove_reference<typename std::remove_cv<decltype(flow[sink][source])>::type>::type Weight;

	struct AugmentPathFinder{
		STFlow& flow;
		const WeightedDirGraph& graph;

		AugmentPathFinder(const WeightedDirGraph& graph, STFlow& flow) : graph(graph), flow(flow){
		}

		struct AugmentPath{
			STFlow& flow;
			const WeightedDirGraph& graph;
			std::vector<Vertex> path;
			AugmentPath(const WeightedDirGraph& graph, STFlow& flow) : graph(graph), flow(flow){
			}
			bool adjust(){
				if(path.empty()) return false;
				Weight minWeight = Weight();
				//iterate all the forward edge and backward edge to find the 
				//minimum capacity.
				for(size_t i = 1;i < path.size();++ i){
					auto weightedEdge = graph.getWeightedEdge(path[i], path[i - 1]);
					if(nullptr != weightedEdge){
						//check if a forward edge.
						auto weight = weightedEdge->getWeight();
						auto backEdgesItor = flow.find(path[i - 1]);
						if(backEdgesItor != flow.end()){
							const auto& backEdges = backEdgesItor->second;
							auto backEdgeItor = backEdges.find(path[i]);
							if(backEdgeItor != backEdges.end()){
								weight = weight - backEdgeItor->second;
							}
						}
						if(minWeight == Weight()) minWeight = weight;
						else minWeight = std::min(minWeight, weight);
					}else{
						//check if a backward edge.
						auto backEdgesItor = flow.find(path[i]);
						if(backEdgesItor != flow.end()){
							const auto& backEdges = backEdgesItor->second;
							auto backEdgeItor = backEdges.find(path[i - 1]);
							if(backEdgeItor != backEdges.end()){
								if(minWeight == Weight()) backEdgeItor->second;
								else minWeight = std::min(minWeight, backEdgeItor->second);
							}
						}
					}
				}
				
				//adjust the capacity for each edge.
				if(Weight() < minWeight){
					for(size_t i = 1;i < path.size();++ i){
						auto backEdgesItor = flow.find(path[i]);
						if(backEdgesItor != flow.end()){
							auto& backEdges = backEdgesItor->second;
							auto backEdgeItor = backEdges.find(path[i - 1]);
							if(backEdgeItor != backEdges.end()){
								if(!(minWeight < backEdgeItor->second))
									backEdges.erase(path[i - 1]);
								else
									backEdgeItor->second = backEdgeItor->second - minWeight;
							}else{
								assert(nullptr != graph.getWeightedEdge(path[i], path[i - 1]));
								//std::cout << "add " << minWeight << " to (" << path[i - 1] << "," << path[i] << ")" << std::endl;
								flow[path[i - 1]][path[i]] += minWeight;
								assert(!(graph.getWeight(path[i], path[i - 1]) < flow[path[i - 1]][path[i]]));
							}
						}else{
							assert(nullptr != graph.getWeightedEdge(path[i], path[i - 1]));
							//std::cout << "add " << minWeight << " to (" << path[i - 1] << "," << path[i] << ")" << std::endl;
							flow[path[i - 1]][path[i]] += minWeight;
							assert(!(graph.getWeight(path[i], path[i - 1]) < flow[path[i - 1]][path[i]]));
						}
					}

					//clear vertices that doesn't have backward edges.
					std::vector<Vertex> toBeRemoved;
					for(const auto& backEdges : flow)	if(backEdges.second.empty()) toBeRemoved.push_back(backEdges.first);
					for(auto& v : toBeRemoved) flow.erase(v);
					//std::cout << "min weight is : " << minWeight << std::endl;
					return true;
				}
				return false;
			}
		};

		// an edge is a forward edge if it still has remaining capacity.
		// precondition is that (src, dest) is an edge in the ST-network.
		bool isForwardEdge(Vertex src, Vertex dest){
			auto backEdgesItor = flow.find(dest);
			if(backEdgesItor == flow.end()) return true;
			
			const auto& backEdges = backEdgesItor->second;
			auto backEdgeItor = backEdges.find(src);
			if(backEdgeItor == backEdges.end()) return true;
			//std::cout << "weight is " << graph.getWeight(src, dest) << std::endl;
			return backEdgeItor->second < graph.getWeight(src, dest);
		}

		AugmentPath construct(const std::map<Vertex, Vertex>& childToParent,
			Vertex source, Vertex sink){
			AugmentPath path{graph, flow};
			auto pos = childToParent.find(sink);
			while(pos != childToParent.end()){
				path.path.push_back(pos->first);
				if(pos->second == source) break;
				pos = childToParent.find(pos->second);
			}
			if(pos == childToParent.end()){
				assert(false);
				path.path.clear();
			}else{
				path.path.push_back(pos->second);
			}
			//for(auto& v: path.path) std::cout << v << ' ';
			//std::cout << std::endl;
			return path;
		}

		AugmentPath findOne(Vertex source, Vertex sink){
			std::map<Vertex, Vertex> childToParent;
			std::queue<unsigned int> parents;
			std::queue<unsigned int> children;
			parents.push(source);

			while(!parents.empty()){
				while(!parents.empty()){
					auto current = parents.front();

					parents.pop();

					//visiting all unvisited forward edges.
					auto adj = graph.outAdj(current);
					if(adj != nullptr){
						for(auto v : *adj)
							if(childToParent.find(v) == childToParent.end() && isForwardEdge(current, v)){
								//std::cout << "adding forward edge " << current << ',' << v << std::endl;
								childToParent[v] = current;
								if(v == sink ) return construct(childToParent, source, sink);
								children.push(v);
							}
					}

					//visiting all unvisited back edges.
					auto backEdgesItor = flow.find(current);
					if(backEdgesItor != flow.end()){
						for(const auto& backEdge : backEdgesItor->second)
							if(childToParent.find(backEdge.first) == childToParent.end()){
								//std::cout << "adding backward edge " << current << ',' << backEdge.first << std::endl;
								childToParent[backEdge.first] = current;
								if(backEdge.first == sink) return construct(childToParent, source, sink);
								children.push(backEdge.first);
							}
					}
				}

				parents.swap(children);
			}
			return {graph, flow};
		}
	};

	while(AugmentPathFinder(graph, flow).findOne(source, sink).adjust());
}
}
#endif //MY_LIB_GRAPH_H
