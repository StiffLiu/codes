#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <GL/gl.h>
#include <GL/glut.h>
using namespace std;
const double tolerance = 1e-8;
/* Tell to which side of the 2d vector, which has the 2d point 
 * pStart as the starting point and 2d point pEnd as the ending point, 
 * the 2d point lies.
 * return:
 * 1 : left side
 * -1 : right side
 * 0 : on the vector
 */
template<class T>
int whichSide(T *pStart, T *pEnd, T* pt){
	T deltaX1 = pEnd[0] - pStart[0];
	T deltaY1 = pEnd[1] - pStart[1];
	T deltaX2 = pt[0] - pEnd[0];
	T deltaY2 = pt[1] - pEnd[1];
	T ret = deltaX1 * deltaY2 - deltaX2 * deltaY1;
	if(ret > tolerance)
		return 1;
	if(ret < tolerance)
		return - 1;
	return 0;

}
/* calcluate the convex hull of points with brute force
 * pts : array of 2d points
 * ptCount : number of points in the array
 * This algorithm is based on the Theorem:
 *  If we organize all the edge in clockwise order, 
 *  then all the points must lie to the right of the oriented 
 *  edge. 
*/
template<class T, class Seq>
void slowConvexHull(T* pts, int ptCount, Seq& ret, int stride = 2){
	std::map<int, int> edges;
	for(int i = 0;i < ptCount;++ i)
		for(int j = 0;j < ptCount;++ j)
			if(i != j){
				bool isValid = true;
				for(int k = 0;k < ptCount;++ k)
					if(j != k && i != k){
						int side = whichSide(pts + stride * i, pts + stride * j, pts + stride * k);
						//point "k" is to the left side of vector "(i, j)",
						//thus vector "(i, j)" is not an edge of the convex hull
						if(side == 1){
							isValid = false;
							break;
						}
						//point "k" is on the vector "(i, j)",
						//thus vector "(i, j)" may or may not be an edge of the convex hull
						if(side == 0){
							float testValue = (pts[stride * k] - pts[stride * i]) * 
								(pts[stride * k] - pts[stride * j]);
							//point "k" is outside of the vector "(i, j)",
							//thus vector "(i, j)" may or may not be an edge of the convex hull
							if(testValue > tolerance){
								isValid = false;
								break;
							}
							//there exist two points that are very close to each other
							//which may lead to the failure of this algorithm
							if(testValue > -tolerance){
								std::cerr << "degnerate case with two very near points is found\n";
							}
						}
						//ok, vector "(i, j)" may be an edge of the convex hull
						//continue the test
					}
				if(isValid){
					if(edges.find(i) != edges.end()){
						std::cerr << "the algorithm may not calculate correctly";
					}					
					edges[i] = j;
				}
			}
	if(!edges.empty()){
		std::map<int, int>::iterator begin = edges.begin(), end = edges.end();
		if(begin != end){
			ret.push_back(begin->first);
			ret.push_back(begin->second);
			std::map<int, int>::iterator pos = edges.find(ret.back());
			while(pos != end && pos->second != ret[0]){
				ret.push_back(pos->second);
				pos = edges.find(ret.back());
			}
			if(pos == end)
				std::cerr << "the algorithm may not calculate correctly";
		}
	}else{
		std::cerr << "the algorithm may not calculate correctly";
	}				
}
template<class T, class Seq>
void convexHull(T* pts, int ptCount, Seq& ret, int stride = 2){
	if(ptCount < 1)
		return;
	if(ptCount <= 2){
		ret.push_back(0);
		if(ptCount > 1)
			ret.push_back(1);
		return;
	}
	ret.push_back(0);
	ret.push_back(1);
	for(int i = 2;i < ptCount;++ i){
		ret.push_back(i);
		while(ret.size() > 2 && whichSide(pts + ret[ret.size() - 3] * stride,
			pts + ret[ret.size() - 2] * stride, pts + ret[ret.size() - 1] * stride) != -1){
			ret[ret.size() - 2] = ret[ret.size() - 1];
			ret.pop_back();
		}
	}
	assert(ret.back() == ptCount - 1);
	ret.push_back(ptCount - 2);
	for(int i = ptCount - 3;i >= 0;-- i){
		ret.push_back(i);
		while(ret.size() > 2 && whichSide(pts + ret[ret.size() - 3] * stride,
			pts + ret[ret.size() - 2] * stride, pts + ret[ret.size() - 1] * stride) != -1){
			ret[ret.size() - 2] = ret[ret.size() - 1];
			ret.pop_back();
		}
	}
	assert(ret[0] == ret.back());
}
#include <cstdlib>
#include <cmath>
double genDouble(){
	return 2.0 * rand() / (double)RAND_MAX - 1.0;
}
int generateData(float *dat, int count, float maxVal, int stride = 2){
	for(int i = 0;i < stride * count;i += stride){
		float radius = genDouble() * maxVal;
		float alpha = genDouble() * 3.14159;
		dat[i] = radius * sin(alpha);
		dat[i + 1] = radius * cos(alpha);
	}
}
template<class T, class Seq>
class Demo{
public:
	virtual void init() = 0;
	virtual bool next() = 0 ;
	virtual T* getPoints() = 0;
	virtual int getPtCount() = 0;
	virtual vector<Seq>& getRet() = 0;
	virtual double getMaxVal() = 0;
};
template<class T, class Seq>
struct ConvexHullDemo : public Demo<T, Seq>{
private:
	int i;
	int state;
	ConvexHullDemo(const ConvexHullDemo&);
	ConvexHullDemo& operator=(const ConvexHullDemo&);
public:
	virtual T* getPoints(){
		return pts;
	}
	virtual int getPtCount(){
		return ptCount;
	}
	virtual vector<Seq>& getRet(){
		return ret;
	}
	virtual double getMaxVal(){
		return maxVal;
	}
	T* pts;
	int ptCount;
	vector<Seq> ret;
	int stride;
	double maxVal;
	ConvexHullDemo(){
		pts = NULL;
	}
	void init(){
		delete []pts;
		ptCount = 1000;
		pts = NULL;
		maxVal = 10;
		i = 0;
		state = 0;
		stride = 2;
		pts = new float[2 * ptCount];
		srand(time(0));
		generateData(pts, ptCount, maxVal);
		sort((pair<float,float>*)pts, (pair<float,float>*)(pts + 2 * ptCount));
		ret.clear();
		ret.push_back(Seq());
	}
	bool next(){
		switch(state){
			case 0:{
				state = 1;
				i = 2;
				if(ptCount < 1)
					return true;
				if(ptCount <= 2){
					ret[0].push_back(0);
					if(ptCount > 1)
						ret[0].push_back(1);
					return true;
				}
				ret[0].push_back(0);
				ret[0].push_back(1);
				break;
			}
			case 1:{
				if(testBack(ret[0]))
					break;
				if(i < ptCount){
					ret[0].push_back(i);
					++ i;
				}
				else{
					state = 2;
					i = ptCount - 3;
				}
				break;			
			}
			case 2:{
				if(testBack(ret[0]))
					break;
				if(i >= 0){
					ret[0].push_back(i);
					-- i;
				}else
					state = 3;
				break;
			}
			case 3:
				return false;
		}
		return true;
	}
	~ConvexHullDemo(){
		delete []pts;
	}
private:
	bool testBack(Seq& ret){
		if(ret.size() > 2 && whichSide(pts + ret[ret.size() - 3] * stride,
			pts + ret[ret.size() - 2] * stride, pts + ret[ret.size() - 1] * stride) != -1){
			ret[ret.size() - 2] = ret[ret.size() - 1];
			ret.pop_back();
			return true;
		}
		return false;
	}
};
ConvexHullDemo<float, vector<int> > instance;
template<class T>
double angle(T *vec){
	double len = sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
	if(len < tolerance)
		return 0;
	double tmp = acos(vec[1] / len);
	assert(tmp >= 0);
	static const double pi = asin(1.0) * 2;
	if(vec[0] < tolerance)
		return 2 *pi - tmp;
	return tmp;
}
template<class T, class Seq>
void convexHull0(T* pts, int ptCount, Seq& ret){
	int idx = 0;
	for(int i = 1;i < ptCount;++ i)
		if(pts[i * 2] < pts[2 * idx] ||
			(fabs(pts[i * 2] - pts[2 * idx]) < tolerance && pts[i * 2 + 1] < pts[2 * idx + 1]))
			idx = i;
	ret.push_back(idx);
	double lastAngle = 0;
	while(ret.size() <= 1 || ret.back() != ret[0]){
		double minAngle = 10;
		double maxLen = 0;
		int index = ret.back();
		for(int i = 0;i < ptCount;++ i)
			if(i != ret.back()){
				double tmp[2] = {pts[2 * i] - pts[ret.back() * 2], pts[2 * i + 1] - pts[ret.back() * 2 + 1]};
				double tmpAngle = angle(tmp);
				if(tmpAngle > lastAngle && tmpAngle < minAngle){
					maxLen = tmp[0] * tmp[0] + tmp[1] * tmp[1];
					minAngle = tmpAngle;
					index = i;
				}else if(fabs(tmpAngle - minAngle) < tolerance){
					double tmpLen = tmp[0] * tmp[0] + tmp[1] * tmp[1];
					if(tmpLen - maxLen > tolerance){
						minAngle = tmpAngle;
						maxLen = tmpLen;
						index = i;
					}
				}
			}
		assert(index != ret.back());
		lastAngle = minAngle;
		ret.push_back(index);
		if(ret.size() > 40)
			break;
	}
}
template<class T, class Seq>
struct ConvexHullDemo0 : public Demo<T, Seq>{
private:
	int i;
	int state;
	ConvexHullDemo0(const ConvexHullDemo0&);
	ConvexHullDemo0& operator=(const ConvexHullDemo0&);
public:
	virtual T* getPoints(){
		return pts;
	}
	virtual int getPtCount(){
		return ptCount;
	}
	virtual vector<Seq>& getRet(){
		return ret;
	}
	virtual double getMaxVal(){
		return maxVal;
	}
	T* pts;
	int ptCount;
	vector<Seq> ret;
	int stride;
	double maxVal;
	ConvexHullDemo0(){
		pts = NULL;
	}
	void init(){
		delete [] pts;
		ptCount = 1000;
		pts = NULL;
		maxVal = 10;
		i = 0;
		state = 0;
		stride = 2;
		pts = new float[2 * ptCount];
		srand(time(0));
		generateData(pts, ptCount, maxVal);
		ret.clear();
		ret.push_back(Seq());
		convexHull0(pts, ptCount, ret[0]);
	}
	bool next(){
		return false;
	}
	~ConvexHullDemo0(){
		delete []pts;
	}
};
ConvexHullDemo0<float, vector<int> > instance0;
template<class T, class Seq>
struct ConvexHullDemo1 : public Demo<T, Seq>{
private:
	int i;
	int state;
	ConvexHullDemo1(const ConvexHullDemo1&);
	ConvexHullDemo1& operator=(const ConvexHullDemo1&);
	int findMin(const Seq& ret){
		int idx = 0;
		for(int i = 1;i < ret.size();++ i)
			if(pts[ret[i] * 2] < pts[2 * ret[idx]] ||
				(fabs(pts[ret[i] * 2] - pts[2 * ret[idx]]) < tolerance && pts[ret[i] * 2 + 1] < pts[2 * ret[idx] + 1]))
			idx = i;
		return idx;
	}
	int findMax(const Seq& ret){
		int idx = 0;
		for(int i = 1;i < ret.size();++ i)
			if(pts[ret[i] * 2] > pts[2 * ret[idx]] ||
				(fabs(pts[ret[i] * 2] - pts[2 * ret[idx]]) < tolerance && pts[ret[i] * 2 + 1] > pts[2 * ret[idx] + 1]))
			idx = i;
		return idx;
	}
	void merge(const Seq& ret1, const Seq& ret2, Seq& ret3){
		int minIndex1 = findMin(ret1);
		int maxIndex1 = findMax(ret1);
		int minIndex2 = findMin(ret2);
		int maxIndex2 = findMax(ret2);
		for(int i = minIndex1;(i % ret1.size()) < maxIndex1;++ i)
			ret3.push_back(ret1[i % ret1.size()]);
		ret3.push_back(ret1[maxIndex1]);
		for(int i = minIndex2;(i % ret2.size()) < maxIndex2;++ i){
			ret3.push_back(ret2[i % ret2.size()]);
			while(testBack(ret3));
		}
		ret3.push_back(ret2[maxIndex2]);
		while(testBack(ret3));
		for(int i = maxIndex2 + 1;(i % ret2.size()) != minIndex2;++ i){
			ret3.push_back(ret2[i % ret2.size()]);
			while(testBack(ret3));
		}
		ret3.push_back(ret2[minIndex2]);
		while(testBack(ret3));
		for(int i = maxIndex1;(i % ret1.size()) != minIndex1;++ i){
			ret3.push_back(ret1[i % ret1.size()]);
			while(testBack(ret3));
		}
		ret3.push_back(ret1[minIndex1]);
		while(testBack(ret3));
	}
public:
	virtual T* getPoints(){
		return pts;
	}
	virtual int getPtCount(){
		return ptCount;
	}
	virtual vector<Seq>& getRet(){
		return ret;
	}
	virtual double getMaxVal(){
		return maxVal;
	}
	T* pts;
	int ptCount;
	int stride;
	double maxVal;
	vector<Seq> ret;
	ConvexHullDemo1(){
		pts = NULL;
	}
	void init(){
		delete []pts;
		ptCount = 1000;
		pts = NULL;
		maxVal = 10;
		i = 0;
		state = 0;
		stride = 2;
		pts = new float[2 * ptCount];
		srand(time(0));
		generateData(pts, ptCount, maxVal);
		sort((pair<float,float>*)pts, (pair<float,float>*)(pts + 2 * ptCount));
		ret.clear();
		for(int i = 0;i < ptCount;i += 2){
			ret.push_back(Seq());
			ret.back().push_back(i);
			ret.back().push_back(i + 1);
			ret.back().push_back(i);
		}
	}
	bool testBack(Seq& ret){
		if(ret.size() > 2 && whichSide(pts + ret[ret.size() - 3] * stride,
			pts + ret[ret.size() - 2] * stride, pts + ret[ret.size() - 1] * stride) != -1){
			ret[ret.size() - 2] = ret[ret.size() - 1];
			ret.pop_back();
			return true;
		}
		return false;
	}
	bool next(){
		if(ret.size() > 1){
			vector<Seq > tmp;
			int count = (ret.size() / 2) * 2;
			for(int i = 0;i < count;i += 2){
				tmp.push_back(Seq());
				ret[i].pop_back();
				ret[i + 1].pop_back();
				merge(ret[i], ret[i + 1], tmp.back());
			}
			if(ret.size() % 2 == 1)
				tmp.push_back(ret.back());
			ret.swap(tmp);
			return true;
		}
		return false;
	}
	~ConvexHullDemo1(){
		delete []pts;
	}
};
ConvexHullDemo1<float, vector<int> > instance1;
#include <ctime>
#include <set>
using namespace std;
Demo<float, vector<int> > *i = &instance1;
static void renderFunction()
{
    float *pts = ::i->getPoints();
    int ptCount = ::i->getPtCount();
    vector<vector<int> >& ret = ::i->getRet();
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3f(1.0, 1.0, 1.0);
    glPointSize(2.0);
    glBegin(GL_POINTS);
	for(int i = 0;i < ptCount;++ i){
		glVertex2fv(pts + 2 * i);
	}
    glEnd();
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    for(int i = 0;i < ret.size();++ i){
    glBegin(GL_POLYGON);
	for(int j = 0;j < ret[i].size();++ j){
		glVertex2fv(pts + ret[i][j] * 2);
	}
    glEnd();
    }
    glFlush();
}
void init(){
	i->init();
	float len = i->getMaxVal() * 1.05;
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-len, len, -len, len, 0, 4);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	//glOrtho(-maxVal, maxVal, -maxVal, maxVal, 0, 4);
	/*gluLookAt(1, 1, 1, 0, 0, 0, 0, 1, 0);*/
}

void animate(int timerId){
	if(i->next()){
		glutTimerFunc(100, animate, timerId);
		glutPostRedisplay();
	}
}
void kbd(unsigned char ch, int x, int y){
	switch(ch){
	case 'r':{
		i->init();
		glutTimerFunc(100, animate, 1000);
		break;
		}
	case 'n':{
		if(i == &instance)
			i = &instance1;
		else if(i == &instance1)
			i = &instance0;
		else 
			i = &instance;
		i->init();
		glutTimerFunc(100, animate, 1000);
	}
	}
	glutPostRedisplay();
}
int main(int argc, char *argv[]){
 	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE);
	glutInitWindowSize(500,500);
	glutInitWindowPosition(100,100);
	glutCreateWindow("convex hull computation");
 	glutDisplayFunc(renderFunction);
	glutKeyboardFunc(kbd);
	glutTimerFunc(100, animate, 1000);
	init();
	glutMainLoop(); 
}
