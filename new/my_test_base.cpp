#include "my_test_base.h"
#include <algorithm>
#include <GL/glew.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include <GL/glext.h>
#include <cfloat>
#include <iostream>
#include <condition_variable>
#include <chrono>
#include <sstream>
#include <string>
#include <iostream>

namespace{
using namespace my_lib;
TwoDPlot *instance = nullptr;
double pointSize = 4.0;
thread_local double charSize = 0.0003;
void openGLInit() {
	glEnableClientState(GL_VERTEX_ARRAY);
	glClearColor(0.000, 0.110, 0.392, 0.0); 
	glColor3f(0.314, 0.314, 0.000); 
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glPointSize(pointSize);
	gluOrtho2D(-1.0, 1.0, -1.0, 1.0);
}
void drawArrays(double *points, double *colors, 
	unsigned int *indices, unsigned int n, GLenum mode){
	if(points == nullptr || n <= 0)
		return;

	//GLint savedBinding = 0;
	//glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &savedBinding);
	//glBindBuffer(GL_ARRAY_BUFFER, 0);
	glVertexPointer(2, GL_DOUBLE, 0, points);
	//For some opengl versions when GL_COLOR_ARRAY is enabled
	//glColorPointer must be used to set color pointer.
	//If we don't have color pointer, we must disable it.
	glColorPointer(3, GL_DOUBLE, 0, colors);
	if(colors != nullptr){
		glEnableClientState(GL_COLOR_ARRAY);
		glColorPointer(3, GL_DOUBLE, 0, colors);
	}else{
		glDisableClientState(GL_COLOR_ARRAY);
	}
	if(indices == nullptr)
		glDrawArrays(mode, 0, n);
	else
		glDrawElements(mode, n, GL_UNSIGNED_INT, indices); 
	//glBindBuffer(GL_ARRAY_BUFFER, savedBinding);	

}
void drawStringImpl(double x, double y, double charSize, const char *s) {
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, 1.0, 0.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	
	glPointSize(1.0);
	glLoadIdentity();
	glTranslatef(x, y, 0.0);
	glScaled(charSize, charSize, 1.0);
	while(*s){	
		glutStrokeCharacter(GLUT_STROKE_ROMAN, *s); 
		++ s;
	}
	glPopMatrix();
	glPointSize(pointSize);
}
void drawStringImpl(double x, double y, const char *s){
	drawStringImpl(x, y, charSize, s);
}
void drawStringImpl(double x, double y, const std::string& s) {
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glTranslatef(x, y, 0.0);
	//glScaled(0.0003, 0.0003, 0.0003);
	for (size_t i = 0; i < s.size(); ++i)
		glutStrokeCharacter(GLUT_STROKE_ROMAN, s[i]);
	glPopMatrix();
}
void display(){
	if(instance != nullptr)
		instance->show();
}

std::mutex m;
std::condition_variable cv;
void idleFunc(){
	std::unique_lock<std::mutex> lk(m);
	cv.wait(lk, []{return true;});
	lk.unlock();
	glutPostRedisplay();
}
void kbd(unsigned char key, int x, int y){
	if(instance != nullptr)
		instance->keyboard(key, x, y);
}
}
namespace my_lib{
/*****************************************************/
void TwoDPlot::redisplay(){
	cv.notify_one();
}
int TwoDPlot::run(int argc, char *argv[]){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(640, 480);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("my test");
	openGLInit();
	glutDisplayFunc(display);
	glutIdleFunc(idleFunc);
	instance = this;	
	glutKeyboardFunc(kbd);
	init();
	glutMainLoop();
	return 0;
}
void TwoDPlot::drawPoints(double *points, 
	double *colors, unsigned int n){
	drawArrays(points, colors, nullptr, n, GL_POINTS);
}
void TwoDPlot::drawPoints(double *points,
	unsigned int n, double *color){
	if(color != nullptr)
		glColor3d(color[0], color[1], color[2]);
	drawPoints(points, nullptr, n);
}
void TwoDPlot::drawPath(double *vertices,
	double *colors, unsigned int n){
	drawArrays(vertices, colors, nullptr, n, GL_LINE_STRIP);
}
void TwoDPlot::drawPath(double *vertices,
	unsigned int n, double *color){
	if(color != nullptr)
		glColor3d(color[0], color[1], color[2]);
	drawArrays(vertices, nullptr, nullptr, n, GL_LINE_STRIP);
}
void TwoDPlot::drawEdges(double *vertices, unsigned int* edges, 
	double *colors, unsigned int n){
	drawArrays(vertices, colors, edges, n, GL_LINES);
}
void TwoDPlot::drawEdges(double *vertices, unsigned int* edges,
	unsigned int n, double *color){
	if(color != nullptr)
		glColor3d(color[0], color[1], color[2]);
	drawArrays(vertices, nullptr, edges, n, GL_LINES);
}
void TwoDPlot::drawAxis(double minX, double maxX,
	double minY, double maxY){
	glBegin(GL_LINE_STRIP);
	glVertex2f(minX, maxY);
	glVertex2f(minX, minY);	
	glVertex2f(maxX, minY);
	glEnd();
}
void TwoDPlot::drawString(double x, double y, const char *str){
	drawStringImpl(x, y, str);
}
/*****************************************************/
void StatPlotBase::collecting(StatPlotBase *plot){
	auto& graphs = plot->graphs;
	std::vector<double> values(2 * graphs.size());
	while(instance == plot){
		if(plot->collect(&values[0])){
			std::lock_guard<std::mutex> lk(plot->m);
			for(size_t i = 0;i < graphs.size();++ i){
				std::vector<double>& points = graphs[i];
				size_t index = 2 * i;
				if(plot->maxNumPoints > 0 && points.size() / 2 > plot->maxNumPoints){
					auto pos = points.begin();
					auto end = pos + 2;
					points.erase(pos, end);
				}
				points.push_back(values[index]);
				points.push_back(values[index + 1]);
			}
			cv.notify_one();
		}
		std::chrono::milliseconds dura(plot->interval);
		std::this_thread::sleep_for(dura);
	}
}
bool StatPlotBase::getBounds(double& xMin, double& xMax, double& yMin, double& yMax) const{
	for(size_t i = 0;i < graphs.size();++ i){
		const std::vector<double>& points = graphs[i];
		for(size_t j = 0;j < points.size();j += 2){
			xMin = std::min(xMin, points[j]);
			xMax = std::max(xMax, points[j]);
			yMin = std::min(yMin, points[j + 1]);
			yMax = std::max(yMax, points[j + 1]);
		}
	}
	const double eps = 1e-3;
	if(xMin == DBL_MAX)
		xMin = 0.0;
	if(yMin == DBL_MAX)
		yMin = 0.0;
	if(xMin + eps > xMax){
		xMax += eps;
		xMin -= eps;
	}
	if(yMin + eps > yMax){
		yMax += eps;
		yMin -= eps;
	}
	return true;
}
bool StatPlotBase::getTitle(RenderInfo& renderInfo) const{
	std::ostringstream os;
	os << "min : " << renderInfo.yMin << ", max : " << renderInfo.yMax;
	renderInfo.title = os.str();
	return true;
}
void StatPlotBase::show(){
	glClear(GL_COLOR_BUFFER_BIT);
	std::lock_guard<std::mutex> lg(m);
	RenderInfo renderInfo;
	getBounds(renderInfo);
	
	double deltaX = renderInfo.xMax - renderInfo.xMin;
	double deltaY = renderInfo.yMax - renderInfo.yMin;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(renderInfo.xMin - 0.01 * deltaX, renderInfo.xMax + 0.01 * deltaX, 
		renderInfo.yMin - 0.01 * deltaY, renderInfo.yMax + 0.01 * deltaY);
	for(size_t i = 0;i < graphs.size();++ i){
		std::vector<double>& points = graphs[i];
		if(!points.empty()){
			double *color = nullptr;
			if(3 * i < colors.size())
				color = &colors[3 * i]; 
			drawPoints(&points[0], points.size() / 2, color);
		}
	}
	if(renderInfo.drawAxis)
		drawAxis(renderInfo.xMin, renderInfo.xMax, renderInfo.yMin, renderInfo.yMax);
	if(getTitle(renderInfo) && !renderInfo.title.empty()){
		::charSize = renderInfo.charSize;
		drawStringImpl(renderInfo.titleX, renderInfo.titleY, renderInfo.title.c_str());
	}
	glFlush();
}
StatPlotBase::~StatPlotBase(){
	instance = nullptr;
	collectingThread.join();
	std::cout << "exited" << std::endl;	
}
void TreePlot::show(){
	double xMin = DBL_MAX, xMax = -DBL_MAX, yMin = DBL_MAX, yMax = -DBL_MAX;
	glClear(GL_COLOR_BUFFER_BIT);
	for(size_t j = 0;j < points.size();j += 2){
		xMin = std::min(xMin, points[j]);
		xMax = std::max(xMax, points[j]);
		yMin = std::min(yMin, points[j + 1]);
		yMax = std::max(yMax, points[j + 1]);

	}

	const double eps = 1e-3;
	if(xMin == DBL_MAX)
		xMin = 0.0;
	if(yMin == DBL_MAX)
		yMin = 0.0;
	if(xMin + eps > xMax){
		xMax += eps;
		xMin -= eps;
	}
	if(yMin + eps > yMax){
		yMax += eps;
		yMin -= eps;
	}

	double deltaX = xMax - xMin;
	double deltaY = yMax - yMin;
	double margin = 0.02;
	double marginX = margin * deltaX;
	double marginY = margin * deltaY;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(xMin - marginX, xMax + marginX, 
			yMin - marginY, yMax + marginY);
	double color[3] = {1.0, 0.0, 0.0};
	if(!points.empty()){
		glPointSize(3.0);
		drawPoints(&points[0], points.size() / 2, color); 

		if(!edges.empty())
			drawEdges(&points[0], &edges[0], edges.size(), color);
		unsigned int numString = points.size() / 2;
		for(unsigned int i = 0;i < numString;++ i){
			unsigned int index = 2 * i;
			const char *str = getString(i);
			if(str != nullptr){
				drawStringImpl((points[index] - xMin + marginX) / (deltaX + 2 * marginX), 
					(points[index + 1] - yMin + marginY) / (deltaY + 2 * marginY), 0.0001, str);
			}
		}
	}
	glFlush();
}
}
