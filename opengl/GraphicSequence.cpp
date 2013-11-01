#include <GL/gl.h>
#include <GL/glut.h>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <map>
#include <set>
using namespace std;

class Graph{
public:
	virtual int vertexCount() = 0;
	virtual bool isAdjacent(int i, int j) = 0;
	virtual ~Graph()=0;
	virtual bool isBipartite(set<int>& V1, set<int>& V2) = 0;
};
Graph::~Graph(){
}
class SimpleGraph : public Graph{
	unsigned int *data;
	int count;
	SimpleGraph(const SimpleGraph&);
	SimpleGraph& operator=(const SimpleGraph);
public:
	SimpleGraph(int cnt = 15){
		reset(cnt);
	}
	void reset(int cnt){
		delete []data;
		data = NULL;
		int total = cnt * cnt;
		const int bitOfInt = (sizeof (int)) * 8;
		int num = total / bitOfInt;
		if((total % bitOfInt) != 0)
			num += 1;
		data = new unsigned int[num];
		for(int i = 0;i < num;++ i)
			data[i] = 0;
		count = cnt;
	}
	void widthFirst(vector<int>& vertices){
		std::set<int> tmp;
		for(int i = 0;i < count;++ i)
			tmp.insert(i);
		while(!tmp.empty()){
			std::set<int>::iterator first = tmp.begin();
			int startIndex = vertices.size();
			vertices.push_back(*first);
			tmp.erase(first);
			while(startIndex < vertices.size()){
				int v = vertices[startIndex];
				std::set<int>::iterator begin = tmp.begin(), end = tmp.end();
				while(begin != end){
					std::set<int>::iterator pos = begin;
					++begin;
					if(isAdjacent(v, *pos)){
						vertices.push_back(*pos);					
						tmp.erase(pos);
					}	
				}
				++startIndex;
			}
		}
	}
	bool isBipartite(set<int>& V1, set<int>& V2){
		if(count == 0)
			return true;
		vector<char> ret(count, 0);
		vector<int> vertices;
		widthFirst(vertices);
		for(int i = 0;i < vertices.size();++ i)
			cout << vertices[i] << ' ';
		cout << endl;
		for(int i = 0;i < count;++ i){
			int v = vertices[i];
			if(ret[v] == 0)
				ret[v] = 1;

			for(int j = 0;j < count;++ j){
				if(isAdjacent(v, j)){					
					if(ret[j] == 0){
						ret[j] = (ret[v] == 1 ? 2 : 1);
					}else if(ret[j] == ret[v]){
						cout << v << ',' << j << endl;
						return false;
					}
				}	
			}
		}

		vector<char>::iterator begin = ret.begin(), end = ret.end();
		for(int i = 0;i < count;++ i){
			assert(ret[i] != 0);
			if(ret[i] == 1)
				V1.insert(i);
			else
				V2.insert(i);
		}
		return true;
	}
	virtual int vertexCount(){
		return count;
	}
	virtual bool isAdjacent(int i, int j){
		assert(i < count && j < count);
		int index = i * count + j;
		const int bitOfInt = (sizeof (int)) * 8;
		return (data[index / bitOfInt] & (1 << (index % bitOfInt))) != 0;
	}
	void set(int i, int j, bool adjacent){
		assert(i < count && j < count);
		int index = i * count + j;
		const int bitOfInt = (sizeof (int)) * 8;
		if(adjacent)
			data[index / bitOfInt] |= (1 << (index % bitOfInt));
		else
			data[index / bitOfInt] &= ~(1 << (index % bitOfInt));
	}
	~SimpleGraph(){
		delete []data;
	}
}instance;
struct GraphData{
	double len;
	int ptCount;
	vector<double> points;
	vector<string> labels;
	vector<int> lines;
	vector<float> colors;
	vector<int> colorIndex;
	Graph *graph;
	string msg;
	GraphData() : graph(&instance), len(2.0){}
	void init(){
		len = 2.0;
		updatePoints();
		resetColor();
		colors.push_back(1.0); colors.push_back(0.0); colors.push_back(0.0);
		colors.push_back(1.0); colors.push_back(1.0); colors.push_back(0.0);
		colors.push_back(0.0); colors.push_back(1.0); colors.push_back(1.0);
	}
	void updateEdges(){
		lines.clear();
		for(int i = 0;i < ptCount;++ i){
			for(int j = i + 1;j < ptCount;++ j)
				if(graph->isAdjacent(i, j)){
					lines.push_back(i);
					lines.push_back(j);
				}
		}	
	}
	void resetColor(){
		colorIndex.clear();
		for(int i = 0;i < ptCount;++ i){
			colorIndex.push_back(0);
		}
	}
	void updatePoints(){
		double radius = 0;
		const double pi2 = asin(1.0) * 4;
		double step = 0;
		points.clear();
		labels.clear();
		radius = len * 0.9;
		ptCount = graph->vertexCount();
		step = pi2 / ptCount;
		for(int i = 0;i < ptCount;++ i){
			double angle = i * step;
			ostringstream os;
			points.push_back(radius * sin(angle));
			points.push_back(radius * cos(angle));
			os << 'v' << i;
			labels.push_back(os.str());
		}
		updateEdges();
	}
}data;
class Demo{
	struct Vertex{
		int degree;
		int index;
		Vertex(int d, int i) : degree(d), index(i){}
		bool operator<(const Vertex& v)const{
			if(degree == v.degree)
				return false;
			if(degree == 0)
				return false;
			if(v.degree == 0)
				return true;
			return degree < v.degree;
		}
	};
	vector<Vertex> degrees;
	vector<Vertex> savedDegrees;
	bool done;
public:
	Demo() : done(true){}
	void updateDesc(){
		int vertexCount= degrees.size();
		data.msg.clear();
		ostringstream os;
		for(int i = 0;i < vertexCount;++ i){
			if(i != 0)
				os << ", ";
			os << degrees[i].degree;
		}
		data.msg = os.str();
	}
	bool next(){
		if(done)
			return false;
		size_t i = 0;
		for(;i < degrees.size();++ i)
			if(degrees[i].degree != 0)
				break;
		if(i < degrees.size() - 1){
			int degree = degrees[i].degree;
			if(degree > degrees.size() - i - 1){
				degrees = savedDegrees;
				updateDesc();
				data.msg += " is not a graphic sequence";
				done = true;
				return false;
			}
			degrees[i].degree = 0;
			for(int j = 1;j <= degree;++ j){
				instance.set(degrees[i].index, degrees[i + j].index, true);
				instance.set(degrees[i + j].index, degrees[i].index, true);
				if(degrees[i + j].degree == 0){
					degrees = savedDegrees;
					updateDesc();
					data.msg += " is not a graphic sequence";
					done = true;
					return false;
				}
				-- degrees[i + j].degree;
			}
			sort(degrees.begin(), degrees.end());
			reverse(degrees.begin(), degrees.end());
			data.updateEdges();
			updateDesc();
			return true;
		}
		degrees = savedDegrees;
		updateDesc();
		data.msg += " is a graphic sequence";
		done = true;
		return false;
	}
	void reset(){
		int vertexCount = 10 + rand() % 10;
		done = false;
		degrees.clear();
		for(int i = 0;i < vertexCount;++ i){			
			degrees.push_back(Vertex(rand() % (vertexCount - 1), i));
		}
		sort(degrees.begin(), degrees.end());
		reverse(degrees.begin(), degrees.end());
		savedDegrees = degrees;
		updateDesc();
		instance.reset(vertexCount);
		data.init();
	}
	void testBipartie(){
		set<int> V1, V2;
		updateDesc();
		if(data.graph->isBipartite(V1, V2)){
			data.msg = " is bipartitie";
			set<int>::iterator begin = V1.begin(), end = V1.end();
			while(begin != end){
				data.colorIndex[*begin] = 1;
				++begin;
			}
			begin = V2.begin(), end = V2.end();
			while(begin != end){
				data.colorIndex[*begin] = 2;
				++begin;
			}
		}else
			data.msg = " is not bipartite";
	}
}demo;
static void renderFunction()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glPointSize(5.0);
    glBegin(GL_POINTS);
	for(int i = 0;i < data.ptCount;++ i){
    		glColor3fv(&data.colors[3 * data.colorIndex[i]]);
		glVertex2dv(&data.points[2 * i]);
	}
    glEnd();
    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINES);
	for(int i = 0;i < data.lines.size();i += 2){
		glVertex2dv(&data.points[2 * data.lines[i]]);
		glVertex2dv(&data.points[2 * data.lines[i + 1]]);
	}
    glEnd();
    glMatrixMode(GL_MODELVIEW);
    glPointSize(1.0);
    glColor3f(0.0, 1.0, 0.0);
    for(int i = 0;i < data.ptCount;++ i){
	int index = 2 * i;

	glPushMatrix();
	glTranslated(data.points[index], data.points[index + 1], 0);
	glScaled( 0.001, 0.001, 0.001 );
	for (size_t j = 0; j < data.labels[i].size(); j ++ ){
		glutStrokeCharacter(GLUT_STROKE_ROMAN  , data.labels[i][j]);
	}
	glPopMatrix();
    }
    if(!data.msg.empty()){
	glPushMatrix();
	glTranslated(- data.len * 0.9 + 0.01, - data.len * 0.97, 0);
	glScaled( 0.0007, 0.0007, 0.0007 );
	for (size_t j = 0; j < data.msg.size(); j ++ ){
		glutStrokeCharacter(GLUT_STROKE_ROMAN  , data.msg[j]);
	}
	glPopMatrix();
    }
    glutSwapBuffers();
}
void init(){
	double len = 0;
	len = data.len;
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-len, len, -len, len, 0, 4);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void initGraphicSequence(){
	demo.reset();
	data.init();
}
void initTestbipartite(){
	int vertexCount;
	int maxEdgeCount;
	int count = 0;
	demo.reset();
	vertexCount = instance.vertexCount();
	count = vertexCount * vertexCount;
	maxEdgeCount = vertexCount * (vertexCount - 1) / 2;
	int edgeCount = rand() % maxEdgeCount;
	for(int i = 0;i < edgeCount;++ i){
		int tmp = rand() % count;
		int v1 = tmp / vertexCount, v2 = tmp % vertexCount;
		/*while((v1 = tmp / vertexCount) == (v2 = tmp % vertexCount))*/
		if(v1 == v2)
			continue;
		instance.set(v1, v2, true);
		instance.set(v2, v1, true);
	}
	data.init();
	demo.testBipartie();
}
class IterateGraph{
	unsigned int edgeCount;
	unsigned int vertexCount;
	unsigned int current;
public:
	IterateGraph(unsigned int vertexCount) : vertexCount(vertexCount){
		edgeCount = vertexCount * (vertexCount - 1) / 2;
		current = 0;
	}
	void reset(){
		instance.reset(vertexCount);
		data.init();
		current = 0;
	}
	bool next(){
		if(current >= (1 << edgeCount))
			return false;
		instance.reset(vertexCount);
		unsigned int num = current;
		unsigned int m = 0, n = 0;
		unsigned int c = vertexCount;
		for(unsigned int i = 0;i < edgeCount;++ i){
			if(num % 2 == 1){
				instance.set(m, vertexCount - 1 - n, true);
				instance.set(vertexCount - 1 - n, m, true);
			}
			++ n;
			if(n >= c){
				++ m;
				-- c;
				n = 0;
			}
			num >>= 1;
		}
		data.init();
		++ current;
		return true;
	}
}graphIterator(4);
void initGraphIterator(){
	graphIterator.reset();
}
void (*initFunc[])()={
	initGraphicSequence, initTestbipartite, initGraphIterator, 
};
static void kbd(unsigned char ch, int x, int y);
static void kbd1(unsigned char ch, int x, int y){
	switch(ch){
		case 'r' : demo.reset();break;
		case 'n' : data.resetColor();demo.next();break;
		default:
			kbd(ch, x, y);break;
	}
	glutPostRedisplay();
}
static void kbd2(unsigned char ch, int x, int y){
	switch(ch){
		case 'r' : initTestbipartite();break;
		default:
			kbd(ch, x, y);break;
	}
	glutPostRedisplay();
}
static void kbd3(unsigned char ch, int x, int y){
	switch(ch){
		case 'n' : data.resetColor();graphIterator.next();break;
		default:
			kbd(ch, x, y);break;
	}
	glutPostRedisplay();
}
void (*kbdFunc[])(unsigned char ch, int x, int y)={
	kbd1, kbd2, kbd3, 
};
#define SizeOfArray(a) ((sizeof (a)) / (sizeof (*a)))
static int current = 0;
void kbd(unsigned char ch, int x, int y){
	if(ch == 's'){
		++current;
		current %= SizeOfArray(initFunc);
		initFunc[current]();
		glutKeyboardFunc(kbdFunc[current]);
	}
}
int main(int argc, char *argv[]){
 	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE);
	glutInitWindowSize(500,500);
	glutInitWindowPosition(100,100);
	glutCreateWindow("convex hull computation");
 	glutDisplayFunc(renderFunction);
	glutKeyboardFunc(kbd);
	init();
	glutMainLoop();
	return 0;
}
