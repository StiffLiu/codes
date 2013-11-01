#include <GL/glut.h>
#include <cassert>
#include <cmath>
#include "Point3.h"
using namespace std;
using namespace tmlg;

class Spoke{
	unsigned int count;
	float radius;
	float angle;
public:
	Spoke(unsigned int n, float radius)
		: count(n), radius(radius){
		assert(n > 0);
		assert(radius > 0);
		angle=0.1;
	}
	unsigned int getCount(){
		return count;
	}
	void calculatePoints(Point3F *pts){
		if(count > 0){
			float unit = 2 * ConstantF::pi / count;
			for(unsigned int i = 0;i < count;++ i){
				unsigned int index = 3 * i;
				float ang = i * unit;
				float ang1 = ang - angle / 2;
				float ang2 = ang1 + angle;
				pts[index] = Point3F(0, 0, 0);
				pts[index + 1] = Point3F(radius * cos(ang1), radius * sin(ang1), 0);
				pts[index + 2] = Point3F(radius * cos(ang2), radius * sin(ang2), 0);
			}
		}		
	}
	void draw(){
		Point3F *pts = new Point3F[count * 3];
		calculatePoints(pts);
		glBegin(GL_TRIANGLES);
			unsigned int pointsCount = count * 3;
			for(unsigned int i = 0;i < pointsCount;++ i)
				glVertex3fv(pts[i]);
		glEnd();
		delete []pts;
	}
};
class Animation{
	Spoke spoke;
	unsigned long long frameId;
	float rotAngle;
	bool _useAccumBuf;
	bool _pause;
public:
	Animation() : spoke(8, 4.0), 
		frameId(0), rotAngle(0), _useAccumBuf(false){
	}
	void next(){
		if(!pause()){
			const float maxSpeed = 0.2;
			float t = frameId  * ConstantF::pi / 1000;
			/*if(t < ConstantF::pi){
				rotAngle = maxSpeed / 2 * (t - sin(t));
			}else
				rotAngle += maxSpeed;*/
			rotAngle = 360 * sin(t*2);
			++frameId;
		}
	}
	bool pause(){
		return _pause;
	}
	bool pause(bool val){
		_pause = val;
	}
	bool useAccumBuf(){
		return _useAccumBuf;
	}
	bool useAccumBuf(bool val){
		return _useAccumBuf = val;
	}
	void display(){
		glColor3f(1.0, 1.0, 1.0);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glRotatef(rotAngle, 0, 0, 1.0);
		spoke.draw();
		if(useAccumBuf()){
			glAccum( GL_MULT,   0.90 );
			glAccum( GL_ACCUM,  0.10 );
			glAccum( GL_RETURN, 1. );
		}
	}
}animation;
void init(){
	glClearColor(0.0, 0.0, 1.0, 0.0);
}
void reshape(int w, int h){
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(-4.5, 4.5, -4.5, 4.5);
}
void display(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	animation.display();
	glutSwapBuffers();
}
void animate(int timerId){
	glutTimerFunc(10, animate, timerId);
	animation.next();
	glutPostRedisplay();
}
void kbd(unsigned char key, int x, int y){
	switch(key){
		case 'p': animation.pause(!animation.pause());break;
		case 'a': animation.useAccumBuf(!animation.useAccumBuf());break;
	}
}
int main(int argc, char *argv[]){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_ACCUM);
	glutInitWindowSize(400, 400);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("moving spoke");

	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutTimerFunc(10, animate, 100);
	glutKeyboardFunc(kbd);
	init();
	glutMainLoop();	
}
