#include <GL/glut.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <cassert>
#include "Cylinder.h"
#include "Camera.h"
using namespace tmlg;
using namespace std;
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
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHT2);
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
	glEnable(GL_CULL_FACE);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, white );
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, yellow );
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, tmpColor );
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess );
	{
		float ptBot[3]={0, 0, 0}, ptTop[3]={4, 0, 0}, ptPt[3]={4.2, 0, 0}, 
			color[4] = {1.0, 0.0, 0.0, 0.5}, topRadius = 0.5, botRadius = 0.1;
		int slice = 20;
		float *body = new float[2 * (slice + 1)* 3];
		float *bodyNorm = new float[2 * (slice + 1)* 3];
		float *bot = new float[(slice + 2)*3];
		float *top = new float[(slice + 2)*3];
		CylinderF cylinderAxisX(ptBot, ptTop, topRadius, botRadius);
		if(cylinderAxisX.calcData(body, bodyNorm, top, bot, slice)){
			Point3F axis = cylinderAxisX.getAxis();
			glPolygonMode(GL_FRONT, GL_LINE);
			glColor3f(0, 0, 1);
			glBegin(GL_QUAD_STRIP);
			for(int i = 0;i < 2*(slice + 1);++ i){
				glNormal3fv(bodyNorm + 3 * i);
				glVertex3fv(body + 3 * i);
			}
			glEnd();
			glColor3f(0, 1, 0);
			glPushMatrix();
			glTranslatef(ptTop[0], ptTop[1], ptTop[2]);
			glutSolidSphere(0.1, 20, 20);
			glPopMatrix();
			glBegin(GL_TRIANGLE_FAN);
			for(int i = 0;i < slice + 2;++ i){
				glNormal3fv(axis);
				glVertex3fv(top + 3 * i);
			}
			glEnd();
			glColor3f(1, 0, 0);
			glPushMatrix();
			glTranslatef(ptBot[0], ptBot[1], ptBot[2]);
			glutSolidSphere(0.1, 20, 20);
			glPopMatrix();
			glBegin(GL_TRIANGLE_FAN);
			axis *= -1;
			for(int i = 0;i < slice + 2;++ i){
				glNormal3fv(axis);
				glVertex3fv(bot + 3 * i);
			}
			glEnd();
		}
		delete [] body;
		delete [] bodyNorm;
		delete [] bot;
		delete [] top;
	}
	glEnable(GL_COLOR_MATERIAL);
}
static char ch;
static int Angle;
static CameraF camera(Point3F(8, 0, 8), Point3F(), Point3F(0, 1, 0));
static void display(void) 
{
	const float *eyePos = camera.getEyePos();
	const float *target = camera.getTarget();
	const float *upVector = camera.getUpVector();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if(!lightFollowCamera)
		glLightfv(GL_LIGHT2, GL_POSITION, lit2.position);
	else
		glLightfv(GL_LIGHT2, GL_POSITION, eyePos);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(eyePos[0], eyePos[1], eyePos[2], 
		target[0], target[1], target[2], upVector[0], upVector[1], upVector[2]);
	cout << camera << endl;
	 switch(ch)
    {
		case 'w':
			camera.orbit(Angle / 180.0 * ConstantF::pi);break;
			//glRotatef( Angle, 1.0, 0.0, 0.0); break;
		case 'q':
			glRotatef(-Angle, 1.0, 0.0, 0.0); break;
		case 's':
			glRotatef( Angle, 0.0, 1.0, 0.0); break;
		case 'a':
			glRotatef(-Angle, 0.0, 1.0, 0.0); break;
		case 'x':
			glRotatef( Angle, 0.0, 0.0, 1.0); break;
		case 'z':
			glRotatef(-Angle, 0.0, 0.0, 1.0); break;
	}
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
	init();
	glutMainLoop();
	return 0;
}

static void keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 'w' :
		case 'q' :
		case 's' :
		case 'a' :
		case 'x' :
		case 'z' :
			Angle += 2;ch = key; break;
		case 'i':
			lightFollowCamera = !lightFollowCamera;
	}
	glutPostRedisplay();
}
int main(int argc, char *argv[]){
	return test(argc, argv);
}
