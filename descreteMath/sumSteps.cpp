#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <cassert>
using namespace std;
unsigned int totalSteps(unsigned int m, unsigned int n,
	unsigned int i, unsigned int j){
	if(m <= 0 || n <= 0)
		return 0;
	if(i <= 0 || i > m)
		return 0;
	if(j <= 0 || j > n)
		return 0;
	return i * j * (i + j - 2) / 2 +
		i * (n - j + 1) * (i + n - j + 1 - 2) / 2 +
		(m - i + 1) * j * (m - i + 1 + j - 2) / 2 +
		(m - i + 1) * (n - j + 1) * (m - i + 1 + n - j + 1 - 2) / 2 -
		i * (i - 1) / 2 -
		j * (j - 1) / 2 -
		(n - j + 1) * (n - j) / 2 -
		(m - i + 1) * (m - i) / 2;
}
void geryCodeSequence(unsigned int numBit){
	unsigned int total = (1 << numBit);
	for(unsigned int i = 0;i < total;++ i){
		unsigned int sum = 0;
		unsigned int lastNum = 0;
		unsigned int pow = 1 << (numBit - 1);
		for(unsigned int j = 0;j < numBit;++ j){
			unsigned int number = (i >> (numBit - 1 - j));
			if((lastNum + number) % 2 == 1){
				cout << '1';
				sum += pow;
			}else{
				cout << '0';
			}
			lastNum = number;
			pow >>= 1;
		}
		cout << ':' << sum << endl;
	}
}
class MultiGraph{
public:
	typedef map<int, unsigned int> AdjVertices;
	typedef map<int, AdjVertices> AdjList;
	AdjVertices* getToVertices(int i){
		AdjList::iterator pos = adjList.find(i);
		if(pos == adjList.end())
			return NULL;
		return &pos->second;
	}
	unsigned int edgeCount(int i, int j){
		AdjVertices *adjVertices = getToVertices(i);
		if(adjVertices == NULL)
			return 0;
		AdjVertices::iterator pos = adjVertices->find(j);
		if(pos == adjVertices->end())
			return 0;
		return pos->second;
	}
	unsigned int inDegree(int i){
		AdjList::iterator begin = adjList.begin(), end = adjList.end();
		int degree = 0;
		while(begin != end){
			AdjVertices::iterator pos = begin->second.find(i);
			if(pos != begin->second.end())
				degree += pos->second;
			++begin;
		}
		return degree;
	}
	unsigned int outDegree(int i){
		AdjVertices *adjVertices = getToVertices(i);
		if(adjVertices == NULL)
			return 0;
		AdjVertices::iterator begin = adjVertices->begin(), end = adjVertices->end();
		unsigned degree = 0;
		while(begin != end){
			degree += begin->second;
			++begin;
		}
		return degree;
	}
	void addEdge(int i, int j){
		++adjList[i][j];
	}
	bool removeEdge(int i, int j, bool all = false){
		AdjList::iterator posAdjVertices = adjList.find(i);
		if(posAdjVertices == adjList.end())
			return false;
		AdjVertices::iterator posVertex = posAdjVertices->second.find(j);
		if(posVertex == posAdjVertices->second.end())
			return false;
		if(posVertex->second == 1 || all){
			posAdjVertices->second.erase(posVertex);
			if(posAdjVertices->second.empty()){
				posAdjVertices = adjList.find(i);
				assert(posAdjVertices != adjList.end());
				adjList.erase(posAdjVertices);
			}
			return true;
		}
		-- posVertex->second;
		return true;			
	}
	bool removeVertice(int i){
		AdjList::iterator posAdjVertices = adjList.find(i);
		if(posAdjVertices == adjList.end())
			return false;
		vector<int> all(posAdjVertices->second.size(), 0);
		AdjVertices::iterator begin = posAdjVertices->second.begin(), end = posAdjVertices->second.end();
		int j = 0;
		while(begin != end){
			all[j] = begin->first;
			++ j;
			++ begin;
		}
		for(int k = 0;k < all.size();++ k)
			removeEdge(all[k], i, true);
		adjList.erase(posAdjVertices);
		return true;
	}
	bool isEmpty(){
		return adjList.empty();
	}
	void output(ostream& os) const{
		AdjList::const_iterator begin = adjList.begin(), end = adjList.end();
		while(begin != end){
			AdjVertices::const_iterator start = begin->second.begin(), term = begin->second.end();
			cout << begin->first << '\t';
			while(start != term){
				cout << start->first << ':' << start->second << ' ';
				++start;
			}
			cout << endl;
			++begin;
		}
	}
protected:
	AdjList adjList;
};
class UndirectedMultiGraph : private MultiGraph{
	unsigned int vertexCount;
public:
	UndirectedMultiGraph(unsigned int vertexCount) : vertexCount(vertexCount){
	}
	AdjVertices* getAdjVertices(int i){
		if(i >= vertexCount)
			return NULL;
		return MultiGraph::getToVertices(i);
	}
	bool getOneConnectedComponent(vector<int>& vertices){
		if(isEmpty())
			return false;
		int vertex = adjList.begin()->first;
		set<int> processed;
		vertices.push_back(vertex);
		for(int i = 0;i < vertices.size();++ i){
			AdjVertices* adjVertcies = getAdjVertices(vertices[i]);
			if(adjVertcies != NULL){
				int start = vertices.size();
				AdjVertices::const_iterator begin = adjVertcies->begin(), end = adjVertcies->end();
				while(begin != end){
					if(processed.find(begin->first) == processed.end()){
						processed.insert(begin->first);
						vertices.push_back(begin->first);
					}
					++begin;
				}
			}
			MultiGraph::removeVertice(vertices[i]);
		}
		return true;
	}
	void addEdge(int i, int j){
		if(i >= vertexCount || j >= vertexCount)
			return;
		MultiGraph::addEdge(i, j);
		MultiGraph::addEdge(j, i);
	}
	unsigned int edgeCount(int i, int j){
		if(i >= vertexCount || j >= vertexCount)
			return 0;
		return MultiGraph::edgeCount(i, j);
	}
	bool removeEdge(int i, int j){
		if(i >= vertexCount || j >= vertexCount)
			return false;
		bool ret1 = MultiGraph::removeEdge(i, j);
		bool ret2 = MultiGraph::removeEdge(j, i);
		assert(ret1 == ret2);
		return ret1;
	}
	unsigned int degree(int i){
		if(i >= vertexCount)
			return 0;
		return MultiGraph::outDegree(i);
	}
	unsigned int edgeCount(){
		unsigned int sum = 0;
		AdjList::const_iterator begin = adjList.begin(), end = adjList.end();
		while(begin != end){
			sum += degree(begin->first);
			++begin;
		}
		return sum / 2;
	}
	bool hasEulerCircuit(){
		AdjList::const_iterator begin = adjList.begin(), end = adjList.end();
		while(begin != end){
			if(degree(begin->first) % 2 != 0)
				return false;
			++begin;
		}
		return true;
	}
	bool findEulerCircuit(int *path){
		map<int, unsigned int> degrees;
		unsigned int edgeCount = 0;
		{
			AdjList::const_iterator begin = adjList.begin(), end = adjList.end();
			while(begin != end){
				int deg = degree(begin->first);
				degrees[begin->first] = deg;
				edgeCount += deg;
				++begin;
			}
			edgeCount /= 2;
		}
		if(edgeCount == 0)
			return true;
		unsigned int i = 1, j = edgeCount - 1;
		while(!degrees.empty()){
			assert(j > i);
			map<int, unsigned int>::iterator posi = degrees.begin();
			map<int, unsigned int>::iterator posj = degrees.begin();
			path[0] = posi->first;
			while(posi->second > 0 && posj->second > 0){
				if(posi->second == 1 || posj->second != 1){
					int verti0 = path[i - 1];
					AdjVertices* adjVerti =  getAdjVertices(verti0);
					if(adjVerti == NULL){
							//cout << "***********+++++++" << endl;
							//cout << *this << endl;
							//cout << "***********+++++++" << endl;
						return false;
					}
					int verti = adjVerti->begin()->first;
					if(posi->second == 1)
						degrees.erase(posi);
					else
						-- posi->second;
					removeEdge(verti0, verti);
					cerr << "edge " << path[i - 1] << ',' << verti << " added " << endl;
					path[i] = verti;++ i;
					posi = degrees.find(verti);
					if(posi == degrees.end()){
						if(posj->second != 0){
							//cout << "+++++++" << endl;
							//cout << *this << endl;
							//cout << "+++++++" << endl;
							return false;
						}
					}
					-- posi->second;
				}
				if(posj->second > 0 && (posj->second == 1 || posi->second != 1)){
					int vertj0 = path[(j + 1) % edgeCount];
					AdjVertices* adjVertj =  getAdjVertices(vertj0);
					if(adjVertj == NULL){
							//cout << "-------*******+++++++" << endl;
							//cout << *this << endl;
							//cout << "-------*******+++++++" << endl;
						return false;
					}
					assert(adjVertj != NULL && !adjVertj->empty());
					int vertj = adjVertj->begin()->first;
					if(posj->second == 1)
						degrees.erase(posj);
					else
						-- posj->second;
					removeEdge(vertj0, vertj);
					cerr << "edge " << vertj0 << ',' << vertj << " added " << endl;
					path[j] = vertj;-- j;
					posj = degrees.find(vertj);
					if(posj == degrees.end()){
						if(posi->second != 0){
							//cout << "+++-----------++++" << endl;
							//cout << *this << endl;
							//cout << "+++-----------++++" << endl;
						}
					}
					-- posj->second;
				}
			}
			map<int, unsigned int>::iterator start = degrees.begin(), end = degrees.end();
			while(start != end){
				map<int, unsigned int>::iterator pos = start;
				++ start;
				if(pos->second == 0)
					degrees.erase(pos);
				/*else
					cout << pos->first << ',' << pos->second << endl;*/
			}
			if(!degrees.empty()){
				
				return false;
			}
		}
		return true;
	}
	bool isEmpty(){
		return MultiGraph::isEmpty();
	}
	friend ostream& operator<<(ostream& os, const UndirectedMultiGraph& graph){
		graph.output(os);
		return os;
	}
};
int test(int argc, char *argv[]){
	UndirectedMultiGraph graph(10);
	int pathes[] = {
		1, 2, 	1, 3, 	1, 5, 	1, 6, 
		2, 3, 	2, 5, 	2, 6, 
		3, 4, 	3, 5, 	3, 7, 	3, 8,
 		4, 6, 	4, 7, 	4, 8, 
		5, 6, 
		6, 7,	6, 8,
		7, 9, 
		8, 9, 
	};
	const int pathCount = (sizeof pathes) / (sizeof *pathes) / 2;
	vector<int> path(pathCount, 0);
	for(int i = 0;i < pathCount;++ i)
		graph.addEdge(pathes[2 * i], pathes[2 * i + 1]);
	cout << graph << endl;
	//cout << "number of edges: " << graph.edgeCount() << endl;
	cout << "degrees: " << endl;
	for(int i = 1;i <=10;++ i)
		cout << graph.degree(i) << endl;
	//for(int i = 0;i < pathCount;++ i)
	//	graph.removeEdge(pathes[2 * i], pathes[2 * i + 1]);
	cout << "euler circuit" << endl;
	if(graph.findEulerCircuit(&path[0]))
		for(int i = 0;i < path.size(); ++ i)
			cout << path[i] << ' ';
	else
		cout << "don't have a euler circuit" << endl;
	cout << endl;
	return 0;
}
int test0(int argc, char *argv[]){
	unsigned int m = 5, n = 5;
	for(unsigned int i = 1;i <= m;++ i)
		for(unsigned int j = 1;j <= n;++ j)
			std::cout << '(' << i << ',' << j << ')' << totalSteps(m, n, i, j) << std::endl;
	geryCodeSequence(5);
}
int test1(int argc, char *argv[]){
	int pathes[] = {
		1, 2, 	1, 3,	2, 3, 
		4, 5, 
		6, 7, 	7, 8,  
	};
	UndirectedMultiGraph graph(10);
	const int pathCount = (sizeof pathes) / (sizeof *pathes) / 2;
	vector<int> vertices;
	for(int i = 0;i < pathCount;++ i)
		graph.addEdge(pathes[2 * i], pathes[2 * i + 1]);
	cout << graph << endl;
	while(graph.getOneConnectedComponent(vertices)){
		cout << graph << endl;
		cout << "vertices : ";
		for(int i = 0;i < vertices.size();++ i)
			cout << vertices[i] << ' ';
		cout << endl;
		vertices.clear();
	}
}
int main(int argc, char *argv[]){
	test1(argc, argv);
}
