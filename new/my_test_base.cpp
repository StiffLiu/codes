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

namespace{
using namespace my_lib;
TwoDPlot *instance = nullptr;
void openGLInit() {
	glEnableClientState(GL_COLOR_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);
	glClearColor(0.000, 0.110, 0.392, 0.0); 
	glColor3f(0.314, 0.314, 0.000); 
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glPointSize(5.0);
	gluOrtho2D(-1.0, 1.0, -1.0, 1.0);
}
void drawArrays(double *points, 
	double *colors, unsigned int n, GLenum mode){
	if(points == nullptr || n <= 0)
		return;

	//GLint savedBinding = 0;
	//glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &savedBinding);
	//glBindBuffer(GL_ARRAY_BUFFER, 0);
	glVertexPointer(2, GL_DOUBLE, 0, points);
	if(colors != nullptr)
		glColorPointer(3, GL_DOUBLE, 0, colors);
	glDrawArrays(mode, 0, n);
	//glBindBuffer(GL_ARRAY_BUFFER, savedBinding);	

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
}
namespace my_lib{
/*****************************************************/
int TwoDPlot::run(int argc, char *argv[]){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(640, 480);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("my test");
	glutDisplayFunc(display);
	glutIdleFunc(idleFunc);
	openGLInit();
	instance = this;	
	init();
	glutMainLoop();
	return 0;
}
void TwoDPlot::drawPoints(double *points, 
	double *colors, unsigned int n){
	drawArrays(points, colors, n, GL_POINTS);
}
void TwoDPlot::drawPoints(double *points,
	unsigned int n, double *color){
	if(color != nullptr)
		glColor3d(color[0], color[1], color[2]);
	drawPoints(points, nullptr, n);
}
void TwoDPlot::drawPath(double *vertices,
	double *colors, unsigned int n){
	drawArrays(vertices, colors, n, GL_LINE_STRIP);
}
void TwoDPlot::drawPath(double *vertices,
	unsigned int n, double *color){
	if(color != nullptr)
		glColor3d(color[0], color[1], color[2]);
	drawArrays(vertices, nullptr, n, GL_LINE_STRIP);
}
void TwoDPlot::drawAxis(double minX, double maxX,
	double minY, double maxY){
	glBegin(GL_LINE_STRIP);
	glVertex2f(minX, maxY);
	glVertex2f(minX, minY);	
	glVertex2f(maxX, minY);
	glEnd();
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
				if(points.size() > plot->maxNumPoints){
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
void StatPlotBase::show(){
	glClear(GL_COLOR_BUFFER_BIT);
	std::lock_guard<std::mutex> lg(m);
	double xMin = DBL_MAX, xMax = -DBL_MAX, yMin = DBL_MAX, yMax = -DBL_MAX; 
	for(size_t i = 0;i < graphs.size();++ i){
		std::vector<double>& points = graphs[i];
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
	if(xMin + eps > xMax)
		xMax = xMin + eps;
	if(yMin + eps > yMax)
		yMax = yMin + eps;
	
	double deltaX = xMax - xMin;
	double deltaY = yMax - yMin;
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(xMin - 0.01 * deltaX, xMax + 0.01 * deltaX, 
		yMin - 0.01 * deltaY, yMax + 0.01 * deltaY);
	for(size_t i = 0;i < graphs.size();++ i){
		std::vector<double>& points = graphs[i];
		if(!points.empty()){
			double *color = nullptr;
			if(3 * i < colors.size())
				color = &colors[3 * i]; 
			drawPoints(&points[0], points.size() / 2, color);
		}
	}
	drawAxis(xMin, xMax, yMin, yMax);
	glFlush();
}
StatPlotBase::~StatPlotBase(){
	instance = nullptr;
	collectingThread.join();
	std::cout << "exited" << std::endl;	
}
}
