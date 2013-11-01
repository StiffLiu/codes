#include "glut.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <cassert>
#include "Cylinder.h"
#include "Camera.h"
#include "Cone.h"

using namespace std;
using namespace tmlg;
struct Light{
	float position[4];
	float ambient[4];
	float specular[4];
	float diffuse[4];
};
static Light lit0={
	{-3, 0, 0, 1},
	{1, 0, 0, 1},
	{1, 0, 0, 1},
	{1, 0, 0, 1}
},
lit1 = {
	{3, 0, 0, 1},
	{0, 1, 0, 1},
	{0, 1, 0, 1},
	{0, 1, 0, 1}
}, 
lit2={
	{0, 3, 0, 1},
	{0, 0, 1, 1},
	{0, 0, 1, 1},
	{0, 0, 1, 1}
};
static bool lightFollowCamera = false;
static void init()
{
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_NORMALIZE);
	glEnable(GL_SMOOTH);
	glEnable(GL_BLEND);
	glEnable( GL_COLOR_MATERIAL );
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClearColor(0.5, 0.5, 0.5,0.0);

	glLightfv(GL_LIGHT0, GL_AMBIENT, lit0.ambient);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lit0.specular);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lit0.diffuse);
	glLightfv(GL_LIGHT1, GL_AMBIENT, lit1.ambient);
	glLightfv(GL_LIGHT1, GL_SPECULAR, lit1.specular);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, lit1.diffuse);
	glLightfv(GL_LIGHT2, GL_AMBIENT, lit2.ambient);
	glLightfv(GL_LIGHT2, GL_SPECULAR, lit2.specular);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, lit2.diffuse);


	glLightfv(GL_LIGHT0, GL_POSITION, lit0.position);
	glLightfv(GL_LIGHT1, GL_POSITION, lit1.position);
	glEnable(GL_LIGHT0);
	//glEnable(GL_LIGHT1);
	//glEnable(GL_LIGHT2);
}

static void drawAxis()
{
	GLfloat white[] = {1.0, 1.0, 1.0, 1.0};
    GLfloat yellow[]= {1.0, 1.0, 0.0, 1.0};
	GLfloat tmpColor[] = {0.0, 1.0, 1.0, 1.0};
    GLfloat mat_shininess[]={ 30.0 };
	glHint(GL_LINE_SMOOTH, GL_NICEST);
	glHint(GL_POLYGON_SMOOTH, GL_NICEST);
	glDisable(GL_COLOR_MATERIAL);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, white );
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, yellow );
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, tmpColor );
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess );
	{
		float ptBot[3]={4, 0, 0}, ptTop[3]={0, 0, 0}, topRadius = 1, botRadius = 1.5, ptPt[3]={-1, 0, 0}, radius = 1.5;
		const int slice = 50;
		float body[2 * (slice + 1) * 3];
		float bodyNorm[2 * (slice + 1) * 3];
		float topCap[(slice + 2) * 3];
		float botCap[(slice + 2) * 3];
		CylinderF cylinderAxisX(ptBot, ptTop, topRadius, botRadius);
		ConeF conef(ptTop, ptPt, radius);
		if(cylinderAxisX.calcData(body, bodyNorm, topCap, botCap, slice))
		{
			Point3F axis = cylinderAxisX.getAxis();
			glBegin(GL_QUAD_STRIP);
			for(int i = 0;i < 2 * (slice + 1);++ i)
			{
				glNormal3fv(bodyNorm + 3 * i);
				glVertex3fv(body + 3 * i);
			}
			glEnd();
			glBegin(GL_TRIANGLE_FAN);
			for(int i = 0;i < (slice + 2);++ i)
			{
				glNormal3fv(axis);
				glVertex3fv(topCap + 3 * i);
			}
			glEnd();
			axis *= -1;
			glBegin(GL_TRIANGLE_FAN);
			for(int i = 0;i < (slice + 2);++ i)
			{
				glNormal3fv(axis);
				glVertex3fv(botCap + 3 * i);
			}
			glEnd();
		}
		if(conef.calcData(body, bodyNorm, botCap, slice))
		{
			Point3F axis = conef.getAxis();
			glBegin(GL_TRIANGLE_FAN);
			for(int i = 0;i < (slice + 2);++ i)
			{
				glNormal3fv(bodyNorm + 3 * i);
				glVertex3fv(body + 3 * i);
			}
			glEnd();
			axis *= -1;
			glBegin(GL_TRIANGLE_FAN);
			for(int i = 0;i < (slice + 2);++ i)
			{
				glNormal3fv(axis);
				glVertex3fv(botCap + 3 * i);
			}
			glEnd();
		}		
	}
	glEnable(GL_COLOR_MATERIAL);
}
static CameraF camera(Point3F(8, 8, 8), Point3F(1, 1, 1), Point3F(0, 1, 0));
static void display(void) 
{
	const float *eyePos = camera.getEyePos();
	const float *target = camera.getTarget();
	const float *upVector = camera.getUpVector();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if(!lightFollowCamera)
		glLightfv(GL_LIGHT0, GL_POSITION, lit0.position);
	else
		glLightfv(GL_LIGHT0, GL_POSITION, eyePos);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(eyePos[0], eyePos[1], eyePos[2], target[0], target[1], target[2], upVector[0], upVector[1], upVector[2]);
	drawAxis();
	glutSwapBuffers();
}
static void reshape(int w, int h)
{
	glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(30.0,1.0,1.0,100.0);
}
static void keyboard(unsigned char key, int x, int y);
static void specialKey(int key, int x, int y);
int test(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(800, 800);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Draw Function");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutSpecialFunc(specialKey);
	init();
	glutMainLoop();
	return 0;
}
void outInfo()
{
	Point3F& viewDirect = camera.getTarget() - camera.getEyePos();
	Point3F upVector = camera.getUpVector();
	Point3F ortho;
	upVector.normalize();
	upVector *= viewDirect.dotProduct(upVector);
	ortho = viewDirect - upVector;
	cout << camera << endl;
	cout << "dist : " << viewDirect.norm() << ", horizontal dist : " << ortho.norm() << endl;
}
static void keyboard(unsigned char key, int x, int y)
{
	assert(camera.isValid());
	switch (key)
	{
	case 'w' :camera.orbit(1.0 / 180 * ConstantF::pi);outInfo();break;
	case 'q' :camera.orbit(-1.0 / 180 * ConstantF::pi);outInfo();break;
	case 's' :camera.zoom(0.5);outInfo();break;
	case 'a' :camera.zoom(-0.5);outInfo();break;
	case 'x' :camera.roll(1.0 / 180 * ConstantF::pi);outInfo();break;
	case 'z' :camera.roll(-1.0 / 180 * ConstantF::pi);outInfo();break;	
	case 'i':lightFollowCamera = !lightFollowCamera;break;
	}
	assert(camera.isValid());
	glutPostRedisplay();
}
static void specialKey(int key, int x, int y)
{
	assert(camera.isValid());
	switch(key)
	{
	case GLUT_KEY_UP:camera.movePositionX(0.5);break;
	case GLUT_KEY_DOWN:camera.movePositionX(-0.5);break;
	case GLUT_KEY_LEFT:camera.movePositionY(0.5);break;
	case GLUT_KEY_RIGHT:camera.movePositionY(-0.5);break;
	}
	assert(camera.isValid());
	glutPostRedisplay();
}
int main(int argc, char *argv[]){
	extern int hsl_double_cone(int argc, char** argv);
	extern int heatDistributionBar(int argc, char** argv);
	extern int testColorRamp(int argc, char *argv[]);
	extern int testTransparency(int argc, char *argv[]);
	extern int testCheckerBoard(int argc, char *argv[]);
	extern int genData(int totalLen, const char *file);
	extern bool compareFile(const char *file1, const char *file2);
	extern int renameProjFile(const char *baseDir, const char *fileName);
	return test(argc, argv);
	
}
