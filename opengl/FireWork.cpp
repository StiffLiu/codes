#include <GL/glut.h>
#include <cmath>
#include <list>
#include <iostream>
#include "Point3.h"
#include "Color.h"
#include "Camera.h"
using namespace std;
using namespace tmlg;

struct FireWork{
	float mass;
	float time;
	float gravity;
	float resistence;
	Point3F v0;
	Point3F s0;
	Point3F pos;
	ColorF color;
	FireWork() : 
		color(1.0, 1.0, 1.0), mass(1.0),
		gravity(9.8), resistence(1.0), time(0){}
	void update(){
		//z:(g/r+v0)(1-exp(-r*t))/r-g/r*t
		//x,y:v0/r*(1-exp(-r*t))
		double g = gravity;
		double r = resistence;
		double tmp = 1 - exp(-r * time);
		double ratio = g / r;
		pos.x() = v0.x() / r * tmp + s0.x(); 
		pos.y() = v0.y() / r * tmp + s0.y(); 
		pos.z() = (ratio + v0.z()) * tmp / r - ratio * time + s0.z();
		//cout << pos.z() << endl;
	}
};
struct FireWorks{
	float maxX, maxY, maxZ, gravity, resitence;
	list<FireWork> first;
	list<FireWork> second;
	list<FireWork> newest;
	CameraF camera;
	static float randFloat(){
		return rand() / (double)RAND_MAX;
	}
	FireWorks() : camera(Point3F(), Point3F(), Point3F(0, 0, 1)){
		maxX = maxY = 80;
		gravity = 9.8;
		resitence = 1;
		maxZ = (1 / exp(-3) - 1)*(gravity / resitence);
	}
	void newFireWork(){
		FireWork tmp;
		tmp.mass = 3.0;
		tmp.gravity = gravity;
		tmp.resistence = resitence;
		tmp.v0 = Point3F((2*randFloat() - 1) * maxX, (2*randFloat() - 1) * maxY, 
			(randFloat() * 0.6 + 0.4) * maxZ);
		first.push_back(tmp);		
	}
	static void draw(list<FireWork>& fws){
		list<FireWork>::iterator begin = fws.begin(), end = fws.end();
		while(begin != end){
			glColor3fv(begin->color);
			//glPointSize(begin->mass);
			glVertex3fv(begin->pos);
			++begin;
		}
	}
	void update(list<FireWork>& fws, float interval, list<FireWork>* next = NULL){
		list<FireWork>::iterator begin = fws.begin(), end = fws.end();
		while(begin != end){
			begin->time += interval;
			begin->update();
			if(begin->pos.z() < -0.1){
				list<FireWork>::iterator pos = begin;
				++begin;
				fws.erase(pos);
			}else if(next != NULL && rand() % 100 == 0){
				int count = 3 + rand() % 8;
				FireWork tmp;
				tmp.mass = 3.0;
				tmp.gravity = begin->gravity;
				tmp.resistence = begin->resistence;
				for(int i = 0;i < count; ++ i){
					tmp.v0 = Point3F((2*randFloat() - 1) * 80 / 2, (2*randFloat() - 1) * 80 / 2, 
						(2*randFloat() - 1) * maxZ / 2);
					tmp.s0 = begin->pos;
					next->push_back(tmp);
				}
				list<FireWork>::iterator pos = begin;
				++begin;
				fws.erase(pos);
			}else			
				++begin;
		}
	}
	void init(){
		float lenX = maxZ;
		float lenY = maxZ;
		float lenZ = maxZ;
		float fov = 60;
		float aspectRatio = lenX / lenY;
		float nearClip = 1;//lenX / tan(fov / 360 * ConstantF::pi);
		float farClip = nearClip + 2 * lenY;
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(fov, aspectRatio, nearClip, farClip);
		camera.setEyePos(Point3F(0, nearClip + lenY, lenZ / 2));
	}
}instance;
void display(){
	const float *eyePos = instance.camera.getEyePos();
	const float *target = instance.camera.getTarget();
	const float *upVector = instance.camera.getUpVector();
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(eyePos[0], eyePos[1], eyePos[2], target[0], target[1], target[2], 
			upVector[0], upVector[1], upVector[2]);
	glPointSize(3);
	glBegin(GL_POINTS);
	glVertex3f(0.0, 0.0, 0.0);
	FireWorks::draw(instance.first);
	FireWorks::draw(instance.second);
	FireWorks::draw(instance.newest);
	glEnd();
	glAccum( GL_MULT,   0.90 );
	glAccum( GL_ACCUM,  0.10 );
	glAccum( GL_RETURN, 1. );
	glutSwapBuffers();
}
void reshape(int w, int h){
        glViewport(0,0,(GLsizei)w,(GLsizei)h);
}
void animate(int timerId){
	float interval = 0.001 * 20;
	glutTimerFunc(10, animate, 100);
	instance.update(instance.first, interval, &instance.second);
	instance.update(instance.second, interval, &instance.newest);
	instance.update(instance.newest, interval);
	static int count = 0;
	count = (count + 1) % 100;
	if(count == 0)
		instance.newFireWork();
	glutPostRedisplay();
}
static void keyboard(unsigned char key, int x, int y){
	CameraF& camera = instance.camera;
	switch (key){
	case 'w' :camera.orbit(1.0 / 180 * ConstantF::pi);break;
	case 'q' :camera.orbit(-1.0 / 180 * ConstantF::pi);break;
	case 's' :camera.zoom(0.5);break;
	case 'a' :camera.zoom(-0.5);break;
	case 'x' :camera.roll(1.0 / 180 * ConstantF::pi);break;
	case 'z' :camera.roll(-1.0 / 180 * ConstantF::pi);break;
	case 'u' :camera.movePositionY(0.1);break;
	case 'd' :camera.movePositionY(-0.1);break;
	}
	glutPostRedisplay();
}
void init(){
	instance.init();
}
int main(int argc, char *argv[]){
	glutInit(&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_ACCUM);
	glutInitWindowSize(800, 800);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("fire works");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutTimerFunc(10, animate, 100);
	init();
	glutMainLoop();
}
