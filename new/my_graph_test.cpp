#include "my_graph.h"
#include "my_math.h"
#include <iostream>
#include <map>
#include <set>
#include <ctime>

using namespace std;
void outputPath(const vector<unsigned int>& path){
	if(!path.empty()){
		cout << '\t';
		for(auto v : path)
			cout << v << ' ';
		cout << endl;
	}
}

template<class Graph>
void out_graph_property(const Graph& graph){
	using namespace my_lib;
	GraphProperties property(graph);
	cout << "center : " << property.getCenter() << endl;
	cout << "radius : " << property.getRadius() << endl;
	cout << "peripheral : " << property.getPeripheral() << endl;
	cout << "diameter : " << property.getDiameter() << endl;
	cout << "girth : " << property.getGirth() << endl;
	cout << "eccentricity 0 : " << property.eccentricity(0) << endl;
}
int graph_test(int argc, char *argv[]){
	using namespace my_lib;
	UndirAdjList ual(10);
	DirAdjList dal(20);
	DblDirAdjList ddal(30);
	ddal.getReverse();
	//const char *fileName = "/backup/documents/books/algos4/algs4.cs.princeton.edu/41undirected/tinyCG.txt";
	//const char *fileName = "/backup/documents/books/algos4/algs4.cs.princeton.edu/41undirected/mediumG.txt";
	const char *fileName = "E:\\data\\tinyG.txt";
	TxtFileGraphReader reader(fileName);
	reader >> ual >> dal >> ddal;
	cout << "undirected graph : \n" << ual << endl;
	//cout << "directed graph : \n" << dal << endl;

	ConnectedComponent dfs(ual, UndirDfsStrategy());
	vector<unsigned int> path;
	cout << "connected component count : " << dfs.getComponentCount() << endl;
	cout << "v13 is connected to v231 : " << dfs.pathTo(13, 231, &path) << endl;
	outputPath(path);
	path.clear();
	cout << "v207 is connected to v221 : " << dfs.pathTo(207, 221, &path) << endl;
	outputPath(path);

	out_graph_property(ual);
	cin.get();
	return 0;
}

#include <iostream>
#include <fstream>
#include <sstream>

namespace{
	typedef my_lib::UndirSymbolAdjList<std::string> KevinBacon;
	struct KevinBaconStrategy{
		template<class CC>
		void operator()(const KevinBacon& graph, CC& cc, unsigned int& componentCount){
			int vertex = graph.toVertex("Bacon, Kevin");
			assert(vertex >= 0);
			my_lib::UndirBfsStrategy bfs;
			bfs(vertex, vertex, graph, cc);
		}
	};
}
int symbol_graph_test(int argc, char *argv[]){
	using namespace my_lib;
	using namespace std;

	ifstream in("E:\\data\\movies.txt");
	KevinBacon kevinBacon;
	std::string line;
	int maxLen = 0;
	std::set<std::string> actors;
	while (getline(in, line)){
		std::string movie;
		std::istringstream is(line);
		if (getline(is, movie, '/')){
			std::string actor;
			while (getline(is, actor, '/')){
				kevinBacon.addEdge(movie, actor);
				if (actor.size() > maxLen)
					maxLen = actor.size();
				actors.insert(actor);
			}
		}
	}

	ConnectedComponent keveinBaconComponent(kevinBacon, KevinBaconStrategy());
	int kevinBaconVertex = kevinBacon.toVertex("Bacon, Kevin");
	bool shdOutput = false;
	std::map < int, int > kevinBaconHistogram;
	assert(kevinBaconVertex >= 0);
	for (auto& symbol : kevinBacon.symbols()){
		std::vector<unsigned int> path;
		int anotherVertex = kevinBacon.toVertex(symbol);
		if (anotherVertex == kevinBaconVertex)
			continue;
		if (keveinBaconComponent.pathTo(kevinBaconVertex, anotherVertex, &path)){
			if (path.size() % 2 != 0){
				if (shdOutput){
					assert(actors.find(symbol) != actors.end());
					cout.width(maxLen + 2);
					cout << symbol << ":\t" << path.size() / 2 << endl;
				}
				++kevinBaconHistogram[path.size() / 2];
			}
		}else{
			if (actors.find(symbol) != actors.end())
				++ kevinBaconHistogram[0];
		}
	}
	for (auto& kvp : kevinBaconHistogram){
		cout.width(10);
		cout << kvp.first << ":" << kvp.second << endl;
	}
	out_graph_property(kevinBacon);
	cin.get();
	return 0;
}

int test_random_select(int argc, char *argv[]){
	const unsigned int m = 30;
	const unsigned int n = 10;

	using namespace my_lib;
	unsigned int result[n];
	randomSelect(m, n, result, rand);
	for(unsigned int i = 0;i < n;++ i)
		std::cout << result[i] << ' ';
	std::cout << std::endl;

	unsigned int count[m];
	for (unsigned int i = 0; i < m; ++i)
		count[i] = 0;

	for (unsigned int i = 0; i < 1000000; ++i){
		randomSelect(m, n, result, rand);
		for (unsigned int j = 0; j < n; ++j)
			++count[result[j]];
	}
	std::cout << "==========" << std::endl;
	for (unsigned int i = 0; i < m; ++i)
		std::cout << i << '\t' << count[i] << std::endl;
	std::cin.get();
	return 0;
}

int test_random_undir_graph(int argc, char *argv[]){
	using namespace my_lib;
	UndirAdjList ual(10);
	srand(time(0));
	randomUndirGraph(30, 0.2, ual, rand);
	cout << ual << endl;
	out_graph_property(ual);

	ConnectedComponent dfs(ual, UndirDfsStrategy());
	cout << "connected component count : " << dfs.getComponentCount() << endl;
	return 0;
}

int test_random_dir_graph(int agrc, char *argv[]){
	using namespace my_lib;
	DirAdjList dal(10);
	srand(time(0));
	while(true){
		randomDirGraph(30, 0.2, dal, rand);
		cout << dal << endl;
		Topological topo(dal);
		auto dirCycle = topo.directedCycle();
		if(dirCycle != nullptr){
			std::cout << "directed cycle found : " << std::endl;;
			outputPath(*dirCycle);
		}else{
			std::cout << "is adg, topological order is : " << std::endl;
			assert(topo.order()->size() > 0);
			outputPath(*topo.order());
			assert(isTopologicalOrder(dal.toReverse(), dal, &(*topo.order())[0]));
			break;
		}
	}
	cin.get();
	return 0;
}

int test_math(int argc, char *argv[]){
	std::cout << "99 : " << my_lib::floorSqrt(99) << std::endl;
	std::cout << "6 : " << my_lib::floorTriSqrt0(6) << std::endl;
	std::cout << "15 : " << my_lib::floorTriSqrt0(15) << std::endl;
	std::cout << "99 : " << my_lib::floorTriSqrt1(99) << std::endl;
	return 0;
}

int main(int argc, char *argv[]){
	return test_random_select(argc, argv);
}
