#include "my_test_base.h"
#include <algorithm>
#include <GL/glut.h>
#include <float>
namespace{
thread_local TwoDPlot *instance = nullptr;
void openGLInit() {
	glEnableClientState(GL_COLOR_ARRAY);
	glEnableClientState(GL_VERTEX_ARRAY);
	glClearColor(0.000, 0.110, 0.392, 0.0); 
	glColor3f(0.314, 0.314, 0.000); 
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glPointSize(2.0);
	gluOrtho2D(-1.0, 1.0, -1.0, 1.0);
}
void drawArrays(double *points, 
	double *colors, unsigned int n, GLenum mode){
	if(points == nullptr || n <= 0)
		return;

	GLint savedBinding = 0;
	glGet(GL_ARRAY_BUFFER_BINDING, &savedBinding);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glVertexPointer(2, GL_DOUBLE, 0, points);
	glColorPointer(3, GL_DOUBLE, colors);
	glDrawArrays(mode, 0, n);
	glBindBuffer(GL_ARRAY_BUFFER, savedBinding);	

}
void display(){
	if(instance != nullptr)
		instance->draw();
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
	openGLInit();
	instance = this;	
	init();
	glutMainLoop();
	instance = nullptr;
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
	drawArrays(points, colors, n, GL_LINE_STRIP);
}
void TwoDPlot::drawPath(double *vertices,
	unsigned int n, double *color){
	if(color != nullptr)
		glColor3d(color[0], color[1], color[2]);
	drawArrays(points, nullptr, n, GL_LINE_STRIP);
}
void TwoDPlot::drawAxis(double minX, double maxX,
	double minY, double maxY){
	glBegin(GL_LINES);
	glVertex2f(minX, 0);
	glVertex2f(maxX, 0);	
	glVertex2f(0, minY);
	glVertex2f(0, maxY);
	glEnd();
}
/*****************************************************/
void StatPlot::collecting(StatPlot *plot){
	std::vector<double> values(2 * graphs.size());
	while(instance == plot){
		std::lock_guard<std::mutex> lg(plot->m);
		if(collector(&values[0], values.size())){
			for(size_t i = 0;size_t i < graphs.size();++ i){
				std::vector<double>& points = graphs[i];
				size_t index = 2 * i;
				points.erase(points.begin());
				points.erase(points.begin());
				points.push_back(values[index]);
				points.push_back(values[index + 1]);
			}
		}			
	}
}
void StatPlot::show(){
	glClear(GL_COLOR_BUFFER_BIT);
	std::lock_guard<std::mutex> lg(m);
	double xMin = DBL_MAX, xMax = -DBL_MAX, yMin = DBL_MAX, yMax = DBL_MIN; 
	for(size_t i = 0;size_t i < graphs.size();++ i){
		std::vector<double>& points = graphs[i];
		for(size_t j = 0;j < points.size();++ j){
			size_t index = 2 * j;
			xMin = std::min(xMin, points[index]);
			xMax = std::max(xMin, points[index]);
			yMin = std::min(yMin, points[index + 1]);
			yMax = std::max(yMax, points[index + 1]);
		}
	}
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(xMin, xMax, yMin, yMax);
	for(size_t i = 0;size_t i < graphs.size();++ i){
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
int StatPlot::run(int argc, char *argv[]){
	Super::run(argc, argv);
	collectingThread.join();
}
void StatPlot::init(){
	collectingThread = std::thread(collecting, this);
}
}
