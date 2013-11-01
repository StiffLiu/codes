#include <GL/glut.h>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cassert>
#include "Point3.h"
#include "Transform.h"
#include "Cylinder.h"
#include "Cone.h"
#include "Camera.h"

using namespace std;
using namespace tmlg;

class RotationAnimate{
	Point3F src;
	Point3F dest;
	double angle;
	double current;
	double step;
public:
	RotationAnimate(){
		src = Point3F(1, 1, -1);
		dest = Point3F(1, 2, 3);
		current = 0;
		step = ConstantF::pi / 180;
		angle = acos(src.dotProduct(dest) / (src.norm() * dest.norm()));
	}
	Point3F getAxis() const{
		return src.crossProduct(dest);
	}
	const Point3F& getSrc() const{
		return src;
	}
	const Point3F& getDest() const{
		return dest;
	}
	bool hasNext(){
		if(current >= angle)
			return false;
		current += step;
		return current < angle;
	}
	void reset(){
		current = 0;
	}
	TransformF next(){
		if(current > angle)
			current = angle;
		return TransformF(current, getAxis());
	}
}instance;
class Drawer{
	bool drawBotCap;
	bool drawTopCap;
	bool drawCap;
	int slice;
public:
	Drawer(){
		drawBotCap = drawTopCap = drawCap = false;
		slice = 20;
	}
	void draw(const CylinderF& c){
		float *ptBodys = new float[2 * (slice + 1) * 3];
		float *ptNorms = new float[2 * (slice + 1) * 3];
		float *ptBotCap = NULL, *ptTopCap = NULL;
		if(drawBotCap)
			ptBotCap = new float[(slice + 2) * 3];
		if(drawTopCap)
			ptTopCap = new float[(slice + 2) * 3];
		if(c.calcData(ptBodys, ptNorms, ptBotCap, ptTopCap, slice)){
			glBegin(GL_QUAD_STRIP);
			for(int i = 0;i < 2 * (slice + 1);++ i){
				glNormal3fv(ptNorms + 3 * i);
				glVertex3fv(ptBodys + 3 * i);
			}
			glEnd();
			if(ptBotCap != NULL){
				Point3F axis = c.getAxis();
				axis.normalize();
				glBegin(GL_TRIANGLE_FAN);
				for(int i = 0;i < slice + 2;++ i){
					glNormal3fv(axis);
					glVertex3fv(ptBotCap + 3 * i);
				}
				glEnd();
			}
			if(ptTopCap != NULL){
				Point3F axis = c.getAxis();
				axis.normalize();
				axis *= -1;
				glBegin(GL_TRIANGLE_FAN);
				for(int i = 0;i < slice + 2;++ i){
					glNormal3fv(axis);
					glVertex3fv(ptTopCap + 3 * i);
				}
				glEnd();
			}
		}
		delete []ptBodys;
		delete []ptNorms;
		delete []ptBotCap;
		delete []ptTopCap;
	}
	void draw(const ConeF& c){
		float *ptBodys = new float[(slice + 2) * 3];
		float *ptNorms = new float[(slice + 2) * 3];
		float *ptCaps = NULL;
		if(drawCap)
			ptCaps = new float[(slice + 2) * 3];
		if(c.calcData(ptBodys, ptNorms, ptCaps, slice)){
			glBegin(GL_TRIANGLE_FAN);
			for(int i = 0;i < slice + 2;++ i){
				glNormal3fv(ptNorms + 3 * i);
				glVertex3fv(ptBodys + 3 * i);
			}
			glEnd();
			if(ptCaps != NULL){
				Point3F axis = c.getAxis();
				axis.normalize();
				axis *= -1;
				glBegin(GL_TRIANGLE_FAN);
				for(int i = 0;i < slice + 2;++ i){
					glNormal3fv(axis);
					glVertex3fv(ptCaps + 3 * i);
				}
				glEnd();
			}
		}
		delete []ptBodys;
		delete []ptNorms;
		delete []ptCaps;
	}
};
void drawStatic(){
	Point3F srcPoint = instance.getSrc(), destPoint = instance.getDest(), axis = instance.getAxis();
	srcPoint.normalize(); destPoint.normalize();axis.normalize();
	Point3F srcPtEnd = srcPoint;
	Point3F destPtEnd = destPoint;
	srcPtEnd.multiply(1.08);
	destPtEnd.multiply(1.08);
	CylinderF axisX(Point3F(0, 0, 0), Point3F(1, 0, 0), 0.02, 0.02);
	ConeF capX(Point3F(1, 0, 0), Point3F(1.08, 0, 0), 0.03);
	CylinderF axisY(Point3F(0, 0, 0), Point3F(0, 1, 0), 0.02, 0.02);
	ConeF capY(Point3F(0, 1, 0), Point3F(0, 1.08, 0), 0.03);
	CylinderF axisZ(Point3F(0, 0, 0), Point3F(0, 0, 1), 0.02, 0.02);
	ConeF capZ(Point3F(0, 0, 1), Point3F(0, 0, 1.08), 0.03);
	CylinderF axisSrc(Point3F(0, 0, 0), srcPoint, 0.02, 0.02);
	ConeF capSrc(srcPoint, srcPtEnd, 0.03);
	CylinderF axisDest(Point3F(0, 0, 0), destPoint, 0.02, 0.02);
	ConeF capDest(destPoint, destPtEnd, 0.03);
	CylinderF rotateAxis(axis, -axis, 0.02, 0.02);
	Drawer d;
	d.draw(axisX);d.draw(capX);
	d.draw(axisY);d.draw(capY);
	d.draw(axisZ);d.draw(capZ);
	d.draw(axisSrc);d.draw(capSrc);
	d.draw(axisDest);d.draw(capDest);
	d.draw(rotateAxis);
}
void drawChanging(){
	TransformF t(instance.next());
	Point3F srcPoint = t.compute(instance.getSrc());
	srcPoint.normalize();
	Point3F srcPtEnd = srcPoint;
	srcPtEnd.multiply(1.0799);
	srcPoint *= 1.001;
	CylinderF axisX(t.compute(Point3F(0, 0, 0)), t.compute(Point3F(1, 0, 0)), 0.02, 0.02);
	ConeF capX(t.compute(Point3F(1, 0, 0)), t.compute(Point3F(1.08, 0, 0)), 0.03);
	CylinderF axisY(t.compute(Point3F(0, 0, 0)), t.compute(Point3F(0, 1, 0)), 0.02, 0.02);
	ConeF capY(t.compute(Point3F(0, 1, 0)), t.compute(Point3F(0, 1.08, 0)), 0.03);
	CylinderF axisZ(t.compute(Point3F(0, 0, 0)), t.compute(Point3F(0, 0, 1)), 0.02, 0.02);
	ConeF capZ(t.compute(Point3F(0, 0, 1)), t.compute(Point3F(0, 0, 1.08)), 0.03);
	CylinderF axisSrc(Point3F(0, 0, 0), srcPoint, 0.0199, 0.0199);
	ConeF capSrc(srcPoint, srcPtEnd, 0.0259);
	//CylinderF 
	Drawer d;
	d.draw(axisX);d.draw(capX);
	d.draw(axisY);d.draw(capY);
	d.draw(axisZ);d.draw(capZ);
	d.draw(axisSrc);d.draw(capSrc);
}
static CameraF camera(Point3F(3, 3, 3), Point3F(0, 0, 0), Point3F(0, 1, 0));
void display(){
	float eyePos[]={camera.getEyePos().x(), camera.getEyePos().y(), camera.getEyePos().z(), 0};
	const float *target = camera.getTarget();
	const float *upVector = camera.getUpVector();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(eyePos[0], eyePos[1], eyePos[2], target[0], target[1], target[2], upVector[0], upVector[1], upVector[2]);
	glLightfv(GL_LIGHT0, GL_POSITION, eyePos);

	GLfloat c1[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat c2[]= {0.0, 1.0, 0.0, 1.0};
	GLfloat c3[] = {0.0, 1.0, 0.0, 1.0};
	GLfloat c4[] = {0.5, 0.0, 1.0, 1.0};
	GLfloat c5[]= {0.5, 0.5, 1.0, 1.0};
	GLfloat c6[] = {0.5, 1.0, 0.0, 1.0};

	GLfloat mat_shininess[]={ 30.0 };
	glHint(GL_LINE_SMOOTH, GL_NICEST);
	glHint(GL_POLYGON_SMOOTH, GL_NICEST);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c1);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, c2);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, c3);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess );
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	drawStatic();
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c4);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, c5);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, c6);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess );
	drawChanging();
	glutSwapBuffers();
}
void idle(){
	static clock_t last = clock();
	if((clock() - last) < CLOCKS_PER_SEC / 50)
		return;
	last = clock();
	instance.hasNext();
	glutPostRedisplay();
}
struct Light{
	float position[4];
	float ambient[4];
	float specular[4];
	float diffuse[4];
};
void init(){	
	Light lit0={
		{-3, 0, 0, 1},
		{1, 0, 0, 1},
		{0.0, 0.5, 0, 1},
		{0.0, 0.5, 0, 1}
	};
	glClearColor(0.5, 0.5, 0.5,0.0);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_NORMALIZE);
	glEnable(GL_SMOOTH);
	glEnable(GL_BLEND);	
	glDisable(GL_COLOR_MATERIAL);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glLightfv(GL_LIGHT0, GL_AMBIENT, lit0.ambient);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lit0.specular);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lit0.diffuse);
	glLightfv(GL_LIGHT0, GL_POSITION, lit0.position);
	glEnable(GL_LIGHT0);
}
static void reshape(int w, int h)
{
	int t = min(w, h);
	glViewport(0, 0, t, t);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(30.0,1.0,1.0,100.0);
}
static void keyboard(unsigned char key, int x, int y)
{
	assert(camera.isValid());
	switch (key)
	{
	case 'w' :camera.orbit(1.0 / 180 * ConstantF::pi);break;
	case 'q' :camera.orbit(-1.0 / 180 * ConstantF::pi);break;
	case 's' :camera.zoom(0.5);break;
	case 'a' :camera.zoom(-0.5);break;
	case 'x' :camera.roll(1.0 / 180 * ConstantF::pi);break;
	case 'z' :camera.roll(-1.0 / 180 * ConstantF::pi);break;
	case 'r':instance.reset();break;
	}
	assert(camera.isValid());
	glutPostRedisplay();
}
static const int timerInterval = 20;
static const int timerID = 88;
static void timer(int i){
	glutTimerFunc(timerInterval, timer, timerID);
	instance.hasNext();
	glutPostRedisplay();
}
void menu(int val){
}
void submenu(int val){
}
//static void createMenu(void(*func)(int), char **items)
void createMenu(){
	static int mainMenu = glutCreateMenu(menu);
	static int subMenu = glutCreateMenu(submenu);
	int oldMenu = glutGetMenu();
	glutSetMenu(mainMenu);
	glutAddMenuEntry("F S", 1);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
	glutAddSubMenu("sub", subMenu);
	{
		int oldMenu = glutGetMenu();
		glutSetMenu(subMenu);
		glutAddMenuEntry("T S", 1);
		glutSetMenu(oldMenu);
	}
	glutSetMenu(oldMenu);
}
int main(int argc, char *argv[]){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(800, 800);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Draw Function");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	//glutIdleFunc(idle);
	glutKeyboardFunc(keyboard);
	glutTimerFunc(timerInterval, timer, timerID);
	createMenu();
	init();
	glutMainLoop();
	return 0;
}
