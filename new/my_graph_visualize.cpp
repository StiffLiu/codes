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
	void keyboard(unsigned char key, int x, int y) override {
		if(key != 'n')
			return;
		using namespace my_lib;
		using namespace std;
		typedef std::pair<double, double> Coordinate;
		const unsigned int vertexCount = 100;
		EuclideanGraphT<UndirAdjList, Coordinate> euclideanGraph;
		RandomCoordinateGenerator1 coordGenerator(vertexCount);
		RandomGraphGenerator1 generator(vertexCount);

		euclideanGraph.reset(vertexCount);
		coordGenerator(euclideanGraph);
		generator(euclideanGraph);
		init(euclideanGraph, cout, PairToPointTraits());
		//redisplay();
	}
};
int visualizeGraph(int argc, char *argv[]){
	using namespace my_lib;
	RandomUndirGraphPloter plot;
	plot.keyboard('n', 0, 0);
	plot.run(argc, argv);
	return 0;
}
