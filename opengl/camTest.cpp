#include <iostream>
#include "Camera.h"
#include "Drawer.h"
#include "Color.h"
#include <GL/glut.h>
using namespace std;
using namespace tmlg;

static CameraF camera(Point3F(1, 1, 1), Point3F(0, 0, 0), Point3F(0, 0, 1));
float xLen = 5, zLen = 5, yLen= 5;
float len = 3.14159 / 2;
Point3F (*paramFunc)(float u, float v) = NULL;
GLfloat light_pos0[]={  1.0, 1.0, 1.0,  1.0 }; // first light over z-axis
GLfloat light_col0[]={  1.0,  0.0,  0.0,  1.0 }; // and red
GLfloat amb_color0[]={  0.3,  0.0,  0.0,  1.0 }; // even ambiently

GLfloat light_pos1[]={  -1.0, -1.0, -1.0,  1.0 }; // first light over z-axis
GLfloat light_col1[]={  1.0,  0.0,  0.0,  1.0 }; // and red
GLfloat amb_color1[]={  0.3,  0.0,  0.0,  1.0 }; // even ambiently
void display(){
	const float *eyePos = camera.getEyePos();
	const float *target = camera.getTarget();
	const float *upVector = camera.getUpVector();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(eyePos[0], eyePos[1], eyePos[2], target[0], target[1], target[2], 
			upVector[0], upVector[1], upVector[2]);
	glLightfv(GL_LIGHT0, GL_POSITION, light_pos0 ); // light 0
	glLightfv(GL_LIGHT1, GL_POSITION, light_pos1 ); // light 1
	Drawer drawer;
	ParametricSurfaceGeometryF g;
	g.startU = g.startV = -len;
	g.endU = g.endV = len;
	g.rows = 500;
	g.cols = 500;
	g.evalFunc = paramFunc;
	g.normFunc = NULL;
	//g.calculateNorm = false;
	//glColor3f(1.0, 1.0, 1.0);
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	drawer.draw(g);
	glLineWidth(2);
	glBegin(GL_LINES);
		glColor3f(1.0, 0.0, 0.0);
		glVertex3f(-2, 0, 0);
		glVertex3f(2, 0, 0);
		glColor3f(0.0, 1.0, 0.0);
		glVertex3f(0, -2, 0);
		glVertex3f(0, 2, 0);
		glColor3f(0.0, 0.0, 1.0);
		glVertex3f(0, 0, -2);
		glVertex3f(0, 0, 2);
	glEnd();
	glutSwapBuffers();
}
static float nearPlane = 0.1;
void setPerspective(float nearPlane){
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(30, 1, nearPlane, 100);
}
void reshape(int width, int height){
	glViewport(0, 0, width, height);
	setPerspective(nearPlane);
}
Point3F func1(float u, float v){
	double x = sin(u) * cos(v);
	double y = sin(u) * sin(v);
	double z = cos(u);
	double r = x * y * z;
	return Point3F(r * x, r * y, r * z);
}
Point3F func2(float u, float v){
	double x = sqrt(fabs(u)) - log(v*v+1);
	double y = u /(fabs(v) + 1);
	double z = cos(u);
	return Point3F(x, y, z);
}
Point3F func3(float u, float v){
	double x = u * u - v * v;
	double y = u - v;
	double z = u * v;
	return Point3F(x, y, z);
}
Point3F func4(float u, float v){
	double x = exp(u) - exp(fabs(v) + 1);
	double y = sin(u) - sqrt(fabs(v*v*u));
	double z = cos(u);
	return Point3F(x, y, z);
}
static void keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'w' :camera.orbit(1.0 / 180 * ConstantF::pi);break;
	case 'q' :camera.orbit(-1.0 / 180 * ConstantF::pi);break;
	case 's' :camera.zoom(0.5);break;
	case 'a' :camera.zoom(-0.5);break;
	case 'x' :camera.roll(1.0 / 180 * ConstantF::pi);break;
	case 'z' :camera.roll(-1.0 / 180 * ConstantF::pi);break;
	case 'u' :camera.movePositionY(0.1);break;
	case 'd' :camera.movePositionY(-0.1);break;
	case '+' :nearPlane += 0.1;setPerspective(nearPlane);break;
	case '-' :nearPlane -= 0.1;setPerspective(nearPlane);break;
	case 'l' :xLen += 0.1;break;
	case 'L' :
		if(xLen > 3)
			xLen -= 0.1;
		break;
	case 'y' :yLen += 0.1;break;
	case 'Y' :
		if(yLen > 3)
			yLen -= 0.1;
		break;
	case 'p' :zLen += 0.1;break;
	case 'P' :
		if(zLen > 3)
			zLen -= 0.1;
		break;
	case 'n':{
		static Point3F (*funcs[])(float u, float v)={func1, func2, func3, func4};
		static int current ;
		current ++;
		current %= (sizeof funcs / sizeof *funcs);
		paramFunc = funcs[current];
		break;
	}
	case 'i':
		len -= 0.1;break;
	case 'I':
		len += 0.1;break;
	}
	glutPostRedisplay();
}
void init(){

	glClearColor(0.1, 0.1, 0.1, 0);
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_BLEND);
	glLightfv(GL_LIGHT0, GL_POSITION, light_pos0 ); // light 0
	glLightfv(GL_LIGHT0, GL_AMBIENT, amb_color0 );
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_col0 );
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_col0 );

	glLightfv(GL_LIGHT1, GL_POSITION, light_pos1 ); // light 1
	glLightfv(GL_LIGHT1, GL_AMBIENT, amb_color1 );
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_col1 );
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_col1 );
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	paramFunc = func1;
}
int main(int argc, char *argv[]){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
	glutInitWindowSize(700, 700);
	glutInitWindowPosition(10, 10);
	glutCreateWindow("cam texture");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	init();
	glutMainLoop();
}
