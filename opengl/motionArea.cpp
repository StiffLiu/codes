#include <GL/glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <set>
#include <vector>
#include <algorithm>
#include <iterator>
using namespace std;

class MotionArea{
	vector<pair<int, int> > ptAll;
	vector<pair<int, int> > ptRect;
public:
	MotionArea(){
		winWidth = winHeight = length = 600;
		const int ptCountAll = 9000;
		const int ptCountRect = ptCountAll / 9;
		const int thirdLen = length / 3;
		set<pair<int, int> > setAll;
		set<pair<int, int> > setRect;
		pause = false;
		current = 0;
		/*while(setAll.size() < ptCountAll && setRect.size() < ptCountRect){
			pair<int, int> tmp(rand() % length, rand() % length);
			if(tmp.first < thirdLen || tmp.second < thirdLen ||
				tmp.first > 2 * thirdLen || tmp.second > 2 * thirdLen){
				if(setAll.size() < ptCountAll)
					setAll.insert(tmp);
			}else{
				if(setRect.size() < ptCountRect)
					setRect.insert(pair<int, int>(tmp.first - thirdLen, tmp.second - thirdLen));
			}
		}*/
		while(setAll.size() < ptCountAll){
			pair<int, int> tmp(rand() % length, rand() % length);
			setAll.insert(tmp);
		}
		while(setRect.size() < ptCountRect){
			pair<int, int> tmp(rand() % thirdLen, rand() % thirdLen);
			setRect.insert(tmp);
		}
		copy(setAll.begin(), setAll.end(), back_inserter(ptAll));
		copy(setRect.begin(), setRect.end(), back_inserter(ptRect));
	}
	void draw(int length, vector<pair<int, int> >& pts, int displaceX, int displaceY){
		int startX = (winWidth - length) / 2;
		int startY = (winHeight - length) / 2;
		size_t count = pts.size();
		if(startX < 0)
			startX = 0;
		if(startY < 0)
			startY = 0;
		glViewport(startX, startY, length, length);
    		glMatrixMode(GL_PROJECTION);
   		glLoadIdentity();
		gluOrtho2D(0., (float)length, 0., (float)length);
		glBegin(GL_POINTS);
		for (size_t i=0; i < count; i++){
			glVertex2f((float)((pts[i].first + displaceX + length) % length), 
				(float)((pts[i].second + displaceY + length) % length));
		}
		glEnd();
	}
	void draw(){
		int tmpDisplace =  current % (length / 3);
		glPointSize(2);
		draw(length, ptAll, current, current);
		draw(length / 3, ptRect, -tmpDisplace, tmpDisplace);
		glAccum( GL_MULT,   0.90 );
		glAccum( GL_ACCUM,  0.10 );
		glAccum( GL_RETURN, 1. );
	}
	int current;
	int length;
	int winWidth;
	int winHeight;
	bool pause;
}instance;

void display( void ){
	glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	instance.draw();
	glutSwapBuffers();
}

void reshape(int w,int h){
	instance.winWidth = w;
	instance.winHeight = h;
}

void keyboard(unsigned char key, int x, int y){
	instance.pause = !instance.pause;
	glutPostRedisplay();
}

void animate(int timerId){
	glutTimerFunc(10, animate, 100);
	if(!instance.pause){
		instance.current = (instance.current + 1) % instance.length;
		glutPostRedisplay();
	}
}

int main(int argc, char** argv)
{
	glutInit(&argc,argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(instance.winWidth, instance.winHeight);
	glutInitWindowPosition(70, 70);
	glutCreateWindow("Find the shape from the motion");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutTimerFunc(10, animate, 100);
	glutKeyboardFunc(keyboard);
	glutMainLoop();
}
