#include "my_test_base.h"
#include "my_graph.h"
#include <iostream>
#include <random>


struct PairToPointTraits{
	template<class Coord>
	auto x(const Coord& coord) const -> decltype(coord.first){
		return coord.first;
	}

	template<class Coord>
	auto y(const Coord& coord) const -> decltype(coord.second){
		return coord.second;
	}
};

struct RandomGraphGenerator{
	unsigned int vertexCount = 30;
	double ratio = 0.3;
	RandomGraphGenerator(unsigned int vertexCount) : vertexCount(vertexCount){
	}
	template<class Graph>
	void operator()(Graph& graph){
		randomUndirGraph(vertexCount, ratio, graph, rand);
	}
};

struct RandomGraphGenerator1{
	unsigned int vertexCount = 30;
	double ratio = 0.3;
	RandomGraphGenerator1(unsigned int vertexCount) : vertexCount(vertexCount){
	}
	template<class Graph>
	void operator()(Graph& graph) const {
		std::vector<double> weights;
		weights.resize(vertexCount * (vertexCount - 1) / 2);
		for(unsigned int i = 1;i < vertexCount;++ i){
			unsigned int base = i * (i - 1) / 2;
			auto coord1 = graph.getCoordinate(i);
			if(coord1 == nullptr){
				assert(false);
				continue;
			}
			for(unsigned int j = 0;j < i;++ j){
				unsigned int index =  base + j;
				auto coord2 = graph.getCoordinate(j);
				if(coord2 == nullptr){
					assert(false);
					continue;
				}

				double deltaX = coord1->first - coord2->first;
				double deltaY = coord1->second - coord2->second;
				weights[index] = 1.0 / std::sqrt(deltaX * deltaX + deltaY * deltaY);
			}
		}
		unsigned int edge = static_cast<unsigned int>(vertexCount * std::log(vertexCount) * ratio);
		unsigned int count = 0;

		std::random_device rd;
		std::mt19937 gen(rd());
		std::discrete_distribution<unsigned int> dd(weights.begin(), weights.end());

		clock_t start = clock();
		for(unsigned int i = 0;count < edge;++ i){
			unsigned int e = dd(gen);
			unsigned int v1 = my_lib::floorTriSqrt1(e);
			unsigned int v2 = e - v1 * (v1 - 1) / 2;
			if(!graph.isAdj(v1, v2)){
				++ count;
				graph.addEdge(v2, v1);
			}
			if(clock() - start > 1000 * CLOCKS_PER_SEC)
				break;
		}
	}
};

struct RandomGraphGenerator2{
	unsigned int vertexCount = 30;
	double ratio = 0.2;
	double startRatio = 0.1;
	RandomGraphGenerator2(unsigned int vertexCount, double ratio = 0.2, double startRatio = 0.1) : vertexCount(vertexCount), ratio(ratio), startRatio(startRatio){
	}

	typedef my_lib::WeightedUndirEdge<unsigned int> WeightedEdge;
	template<class Graph>
	void operator()(Graph& graph) const {
		using namespace my_lib;
		using namespace std;
		struct CompareWeight{
			bool operator()(const WeightedEdge& e1, const WeightedEdge& e2){
				return e1.getWeight() < e2.getWeight();
			}
		};

		std::priority_queue<WeightedEdge, std::vector<WeightedEdge>, CompareWeight> edges;
		unsigned int totalEdge = vertexCount * vertexCount / 2;
		unsigned int startEdgeCount = totalEdge * startRatio;
		unsigned int edgeCount = totalEdge * ratio + startEdgeCount;
		for(auto& one : graph.getCoordinates())
			for(auto& another : graph.getCoordinates())
				if(one.first < another.first){
					auto deltaX = one.second.first - another.second.first;
					auto deltaY = one.second.second - another.second.second;
					edges.push(WeightedEdge(one.first, another.first,
						deltaX * deltaX + deltaY * deltaY));
					if(edges.size() > edgeCount)
						edges.pop();
				}
		while(edgeCount > startEdgeCount){
			graph.addEdge(edges.top());
			edges.pop();
			--edgeCount;
		}
	}
};

struct RandomCoordinateGenerator{
	unsigned int vertexCount = 30;
	RandomCoordinateGenerator(unsigned int vertexCount) : vertexCount(vertexCount){
	}
	template<class Graph>
	void operator()(Graph& graph) const {
		for(unsigned int i = 0;i < vertexCount;++ i){
			graph.setCoordinate(i, {rand() / (double)RAND_MAX, rand() / (double)RAND_MAX});
		}
	}
};

struct RandomCoordinateGenerator1{
	unsigned int vertexCount = 30;
	double half_pi = std::asin(1.0);
	RandomCoordinateGenerator1(unsigned int vertexCount) : vertexCount(vertexCount){
	}
	template<class Graph>
	void generate(Graph& graph, double centX, double centY, double radius, unsigned int base, unsigned int count) const {
		double ratio = 1.0;
		if(count <= 10){
			radius /= 2;
			for(unsigned int i = base;i < base + count;++ i){
				double angle = rand() / (double)RAND_MAX * 4 * half_pi;
				double newRadius = rand() / (double)RAND_MAX * ratio * radius;
				double x = newRadius * std::cos(angle), y = newRadius * std::sin(angle);
				graph.setCoordinate(i, {centX + x, centY + y});
			}
			return;
		}
		const unsigned int numRegion = 4;
		const unsigned int numPoint = count / numRegion;
		double x1, y1, x2, y2, x3, y3, x4, y4;
		double angle = rand() / (double)RAND_MAX * half_pi;
		radius /= 2;
		x1 = radius * std::cos(angle), y1 = radius * std::sin(angle);
		angle += half_pi;
		x2 = radius * std::cos(angle), y2 = radius * std::sin(angle);
		angle += half_pi;
		x3 = radius * std::cos(angle), y3 = radius * std::sin(angle);
		angle += half_pi;
		x4 = radius * std::cos(angle), y4 = radius * std::sin(angle);

		radius *= ratio;
		generate(graph, centX + x1, centY + y1, radius, base, numPoint);
		generate(graph, centX + x2, centY + y2, radius, base + numPoint, numPoint);
		generate(graph, centX + x3, centY + y3, radius, base + 2 * numPoint, numPoint);
		generate(graph, centX + x4, centY + y4, radius, base + 3 * numPoint, count - 3 * numPoint);

	}
	template<class Graph>
	void operator()(Graph& graph) const {
		generate(graph, 0, 0, 2, 0, vertexCount);
	}
};

struct RandomUndirGraphPloter : public my_lib::UndirGraphPlot{
	unsigned int vertexCount = 100;
	double startRatio = 0.1, ratio = 0.2;
	mutable std::string text;
	mutable std::unordered_map<unsigned int, unsigned int> pointToVertex;
	
	void keyboard(unsigned char key, int x, int y) override {
		if(key == 'c'){
			std::cin >> vertexCount >> ratio >> startRatio;
			return;
		}
		if(key != 'n')
			return;
		using namespace my_lib;
		using namespace std;
		typedef std::pair<double, double> Coordinate;
		EuclideanGraphT<WeightedUndirAdjList, Coordinate> euclideanGraph;
		RandomCoordinateGenerator coordGenerator(vertexCount);
		RandomGraphGenerator2 generator(vertexCount, ratio, startRatio);

		euclideanGraph.reset(vertexCount);
		coordGenerator(euclideanGraph);
		generator(euclideanGraph);
		init(euclideanGraph, cout, PairToPointTraits());
		pointToVertex.clear();

		EagerPrimMST prim(euclideanGraph);
		for(auto& i : vertexToPoint) pointToVertex[i.second] = i.first;
		typedef UndirEdge<unsigned int> UndirEdge;
		typedef std::map<UndirEdge, int> EdgeToIndex;
		EdgeToIndex edgeToIndex;
		for(size_t i = 0;i < edges.size();i += 2) edgeToIndex[UndirEdge(pointToVertex[edges[i]], pointToVertex[edges[i + 1]])] = i / 2;
		edgeColors.clear();
		edgeColors.resize(edges.size() / 2 * 3);
		//default to red edge
		std::cout << "----------graph------------------" << std::endl;
		std::cout << euclideanGraph << std::endl;
		std::cout << "-----------mst-------------------" << std::endl;
		for(size_t i = 0;i < edgeColors.size();i += 3) edgeColors[i] = 1.0;
		struct MapAdapter{
			EdgeToIndex& edgeToIndex; 
			std::vector<double>& edgeColor;
			MapAdapter(EdgeToIndex& edgeToIndex, std::vector<double>& edgeColor) : edgeToIndex(edgeToIndex), edgeColor(edgeColor){
			}
			typedef my_lib::WeightedUndirEdge<unsigned int> WeightedEdge;
			void push_back(const WeightedEdge *edge) const {
				auto pos = edgeToIndex.find(UndirEdge(edge->getV1(), edge->getV2()));
				if(pos != edgeToIndex.end()){
					//mst edges colored yellow
					edgeColor[pos->second * 3] = 0.0;
					std::cout << "selected ";
				}
				std::cout << *edge << std::endl;
			}
		}edgeSet(edgeToIndex, edgeColors);
		prim.getEdges(euclideanGraph, edgeSet); 
	}

	const char *getString(unsigned int i) const override{
		std::ostringstream os;
		os << pointToVertex[i];
		text = os.str();
		return text.c_str();
	}
};
int visualizeGraph(int argc, char *argv[]){
	using namespace my_lib;
	RandomUndirGraphPloter plot;
	plot.keyboard('n', 0, 0);
	plot.run(argc, argv);
	return 0;
}
