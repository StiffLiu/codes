#include "my_graph.h"
#include "my_math.h"
#include <iostream>
#include <iterator>
#include <algorithm>
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

template<class Graph, class Stream>
void readAdjList(Graph& graph, Stream& stream){
	std::string line;
	while(getline(stream, line)){
		std::istringstream is(line);
		unsigned int v1;
		if(is >> v1){
			unsigned int v2;
			if(is >> v2){
				if(v1 + 1 > graph.getVertexCount())
					graph.extend(graph.getVertexCount() - v1 - 1);
				do{
					if(v2 + 1 > graph.getVertexCount())
						graph.extend(graph.getVertexCount() - v2 - 1);
					//std::cout << "(v1, v2) = (" << v1 << "," << v2 << ")\n";
					graph.addEdge(v1, v2);
				}while(is >> v2);
			}
		}
	}
}
template<class Graph>
void readAdjList(Graph& graph, const char *file){
	graph.reset(0);
	ifstream in(file);
	readAdjList(graph, in);
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

template<class SCC>
void outputSCC(const SCC& scc){
	for(auto& d : scc){
		std::cout << d.first << "\t: ";
		for(auto& v : d.second)
			std::cout << v << '\t';
		std::cout << std::endl;
	}
}

bool getCmd(std::string& cmd){
	std::string newCmd;
	if(!getline(cin, newCmd) || newCmd == "q"){
		cmd = newCmd;
		return false;
	}
	if(!newCmd.empty())
		cmd = newCmd;
	return true;
}
int test_random_dir_graph(int agrc, char *argv[]){
	using namespace my_lib;
	DirAdjList dal(10);
	srand(time(0));

	unsigned int vertexCount = 15;
	double ratio = 0.2;
	std::string cmd = "dag";
	while(true){
		if(cmd != "dag")
			randomDirGraph(vertexCount, ratio, dal, rand);
		else
			randomDAG(vertexCount, ratio, dal, rand);
		//readAdjList(dal, "issue1.txt");
		cout << dal;

		StrongConnectedComponent scc(dal);
		//typedef std::unordered_map<unsigned int, std::unordered_set<unsigned int> > SCC;
		typedef std::map<unsigned int, std::set<unsigned int> > SCC;
		SCC cc;
		scc.getSCC(cc);

		std::cout << "strong connected components : \n";
		std::cout << "===============" << std::endl;
		outputSCC(cc);
		std::cout << "===============" << std::endl;
		TarjanStrategy tarjan(dal);
		cc.clear();
		tarjan.getSCC(cc);
		outputSCC(cc);
		std::cout << "===============\n" << std::endl;

		unsigned int v1 = rand() % vertexCount, v2 = rand() % vertexCount;
		LeastCommonAncestor lca(dal, v1, v2);
		std::cout << "lca of " << v1 << " and " << v2 << " : ";
		std::copy(lca.lca().begin(), lca.lca().end(), std::ostream_iterator<unsigned int>(cout, " "));
		std::cout << std::endl;

		auto dirCycle = scc.directedCycle();
		std::cout << "order is : " << std::endl; 
		outputPath(*scc.order());
		if(dirCycle != nullptr){
			std::cout << "directed cycle found : " << std::endl;;
			outputPath(*dirCycle);
			std::cout << (containsEulerianCycle(dal) ? "" : "not ") << "contains Eulerian cycle." << std::endl;
		}else{
			std::cout << "is dag."<< std::endl;
			assert(isTopologicalOrder(dal.toReverse(), dal, &(*scc.order())[0], scc.order()->size()));
			std::cout << (dagContainsHamiltonianPath(dal) ? "" : "not ") << "contains Hamiltonian path." << std::endl;
		}
		if(!getCmd(cmd))
			break;
	}
	return 0;
}

int test_eulerian(int argc, char *argv[]){
	using namespace my_lib;
	DirAdjList dal(10), dal1(10), dal2(10);
	readAdjList(dal, "directed_eulerian_graph1.txt");
	readAdjList(dal1, "directed_eulerian_graph2.txt");
	readAdjList(dal2, "hamilton_dag.txt");
	assert(containsEulerianCycle(dal));
	assert(containsEulerianCycle(dal1));
	assert(dagContainsHamiltonianPath(dal2));
	Topological topological(dal2);
	std::cout << "topological order is : ";
	outputPath(*topological.order());

	std::vector<unsigned int> path;
	assert(directedEulerianCycle(dal, path));
	outputPath(path);
	path.clear();
	assert(directedEulerianCycle(dal1, path));
	outputPath(path);
	return 0;
}

int test_math(int argc, char *argv[]){
	std::cout << "99 : " << my_lib::floorSqrt(99) << std::endl;
	std::cout << "6 : " << my_lib::floorTriSqrt0(6) << std::endl;
	std::cout << "15 : " << my_lib::floorTriSqrt0(15) << std::endl;
	std::cout << "99 : " << my_lib::floorTriSqrt1(99) << std::endl;
	return 0;
}

#include <iostream>
#include <iterator>
#include <string>
#include <regex>

int test_regex(int argc, char *argv[]){
	std::string s = "Some people, when confronted with a problem, think "
		"\"I know, I'll use regular expressions.\" "
		"Now they have two problems.";

	std::regex self_regex("REGULAR EXPRESSIONS",
		std::regex_constants::ECMAScript | std::regex_constants::icase);
	if (std::regex_search(s, self_regex)) {
		std::cout << "Text contains the phrase 'regular expressions'\n";
	}

	std::regex word_regex("(\\S+)");
	auto words_begin =
		std::sregex_iterator(s.begin(), s.end(), word_regex);
	auto words_end = std::sregex_iterator();

	std::cout << "Found "
		<< std::distance(words_begin, words_end)
		<< " words\n";

	const int N = 6;
	std::cout << "Words longer than " << N << " characters:\n";
	for (std::sregex_iterator i = words_begin; i != words_end; ++i) {
		std::smatch match = *i;
		std::string match_str = match.str();
		if (match_str.size() > N) {
			std::cout << "  " << match_str << '\n';
		}
	}
	std::regex long_word_regex("(\\w{7,})");
	std::string new_s = std::regex_replace(s, long_word_regex, "[$&]");
	std::cout << new_s << '\n';
	std::cin.get();
	return 0;
}

int main(int argc, char *argv[]){
	return test_random_dir_graph(argc, argv);
}
