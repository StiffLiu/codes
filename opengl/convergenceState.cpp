#include <glut.h>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include "Camera.h"
#include "Color.h"
#include "Constant.h"

using namespace std;
using namespace tmlg;
struct Light{
	Point3F position;
	ColorAF ambitient;
	ColorAF diffuse;
	ColorAF specular;
};
struct GlobalData{
	CameraD camera;
	Light lit0;
	Light lit1;
	int particalCount;
	Point3D *particalPos;
	Point3D *particalVelocity;
	double  sphereRadius;
	double  particalRadius;
	bool isStarted;
	double  timeInterval;
	clock_t lastTime;
	GlobalData() : camera(Point3D(6, 6 ,6), Point3D(0, 0, 0), Point3D(0, 1, 0))
	{
		lit0.position = Point3F(4, 4, 4);
		lit0.ambitient = ColorAF(0.3,  0.0,  0.0,  1.0);
		lit0.diffuse = ColorAF(1.0,  0.0,  0.0,  1.0);
		lit0.specular = ColorAF(1.0,  0.0,  0.0,  1.0);


		lit1.position = Point3F(-4, -4, -4);
		lit1.ambitient = ColorAF(0.3,  0.0,  0.0,  1.0);
		lit1.diffuse = ColorAF(1.0,  0.0,  0.0,  1.0);
		lit1.specular = ColorAF(1.0,  0.0,  0.0,  1.0);
		sphereRadius = 3.0;
		particalRadius = 0.05;
		particalPos = NULL;
		particalVelocity = NULL;
		tmpPosition = NULL;
		tmpVelocity = NULL;
		isStarted = false;
		initParticalData();
		timeInterval = 0.0005;
		lastTime = clock();
	}
	void reinit()
	{
		isStarted = false;
		initParticalData();
	}
	void next()
	{
		if(tmpPosition == NULL)
			tmpPosition = new Point3D[particalCount];
		if(tmpVelocity == NULL)
			tmpVelocity = new Point3D[particalCount];
		Point3D *tmp = particalPos;
		particalPos = tmpPosition;
		tmpPosition = tmp;
		tmp = particalVelocity;
		particalVelocity = tmpVelocity;
		tmpVelocity = tmp;

		for(int i = 0;i < particalCount;++ i)
		{
			Point3D totalForce;
			Point3D newVolocity;
			for(int j = 0;j < particalCount;++ j)
			{
				Point3D force = tmpPosition[i] - tmpPosition[j];
				double normSquare = force.normSquare();
				if(normSquare > ConstantD::tolerance)
				{
					force.multiply(1 /(normSquare * sqrt(normSquare)));
					totalForce.add(force);
				}
			}
			//totalForce.multiply(0.5);
			particalPos[i] = tmpVelocity[i];
			particalPos[i].multiply(timeInterval);
			particalPos[i].add(tmpPosition[i]);
			particalPos[i].normalize();
			particalPos[i].multiply(sphereRadius);
			newVolocity = tmpPosition[i].orth(totalForce);
			newVolocity.multiply(timeInterval);
			newVolocity.add(tmpVelocity[i]);
			particalVelocity[i] = particalPos[i].orth(newVolocity);
			double norm = particalVelocity[i].norm();
			if(norm > 10000)
				particalVelocity[i].multiply(10000 / norm);
		}
	}
	~GlobalData()
	{
		delete []particalPos;
		delete []particalVelocity;
		delete []tmpPosition;
		delete []tmpVelocity;
	}
private:
	GlobalData(const GlobalData&);
	GlobalData& operator=(const GlobalData&);
	void initParticalData()
	{
		particalCount = 50;
		srand(time(0));
		if(particalPos == NULL)
			particalPos = new Point3D[particalCount];
		if(particalVelocity == NULL)
			particalVelocity = new Point3D[particalCount];
		for(int i = 0;i < particalCount;++ i)
		{
			double phi = rand() / (float)RAND_MAX * ConstantF::pi;
			double alpha = rand() / (float)RAND_MAX * ConstantF::pi * 2;
			Point3D *pt = (particalPos + i);
			double sinPhi = sin(phi);
			pt->x() = sphereRadius * sinPhi * sin(alpha);
			pt->y() = sphereRadius * sinPhi * cos(alpha);
			pt->z() = sphereRadius * cos(phi);
		}
	}
	Point3D *tmpPosition;
	Point3D *tmpVelocity;
}globalData;
static void init()
{
	glEnable(GL_DEPTH);
	glClearColor(0, 0, 0, 1);
	glShadeModel(GL_SMOOTH);
	glLightfv(GL_LIGHT0, GL_POSITION, globalData.lit0.position);
	glLightfv(GL_LIGHT0, GL_AMBIENT, globalData.lit0.ambitient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, globalData.lit0.diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, globalData.lit0.specular);

	glLightfv(GL_LIGHT1, GL_POSITION, globalData.lit1.position);
	glLightfv(GL_LIGHT1, GL_AMBIENT, globalData.lit1.ambitient);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, globalData.lit1.diffuse);
	glLightfv(GL_LIGHT1, GL_SPECULAR, globalData.lit1.specular);
	glEnable(GL_BLEND);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}
static void display()
{
	const double *eyePos = globalData.camera.getEyePos();
	const double *target = globalData.camera.getTarget();
	const double *upVector = globalData.camera.getUpVector();
	ColorAF color(1.0, 0.5, 0.5, 0.1);
	glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(eyePos[0], eyePos[1], eyePos[2], target[0], target[1], target[2], upVector[0], upVector[1], upVector[2]);
	//glColor4f(1.0, 0.5, 0.5, 0.1);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glutSolidSphere(globalData.sphereRadius, 30, 30);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, ColorAF(1.0, 0.0, 0.0, 1.0));
	for(int i = 0;i < globalData.particalCount;++ i)
	{
		Point3D *pt = globalData.particalPos + i;
		glTranslatef(pt->x(), pt->y(), pt->z());
		glutSolidSphere(globalData.particalRadius, 10, 10);
		glTranslatef(-pt->x(), -pt->y(), -pt->z());
	}
	glutSwapBuffers();
}
static void reshape(int width, int height)
{
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, width / (float)height, 1, 100);
}
static void idle()
{
	if(globalData.isStarted)
	{
		clock_t current = clock();
		if((current - globalData.lastTime) < (CLOCKS_PER_SEC / 20))
			return;
		globalData.lastTime = current;
		for(int i = 0;i < 200;++ i)
			globalData.next();
		glutPostRedisplay();
	}
}
static void keyboard(unsigned char key, int x, int y)
{
	switch(key)
	{
	case 'r':globalData.reinit();break;
	case 's':globalData.isStarted = !globalData.isStarted;break;
	case 'w' :globalData.camera.orbit(1.0 / 180 * ConstantF::pi);break;
	case 'q' :globalData.camera.orbit(-1.0 / 180 * ConstantF::pi);break;
	}
	glutPostRedisplay();
}
int convergenceState(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE  | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(500, 500);
	glutCreateWindow("convergence state");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutIdleFunc(idle);
	glutKeyboardFunc(keyboard);
	init();
	glutMainLoop();
	return 0;
}