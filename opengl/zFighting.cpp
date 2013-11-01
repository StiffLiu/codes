#include <GL/glut.h>
#include <iostream>
using namespace std;

void init(){
	glEnable(GL_DEPTH_TEST);
	glClearColor(0, 0, 0, 1.0);
}
double diff = 2.1341e-07;
double inc = 2.0f;
void keyboard(unsigned char ch, int x, int y){
	switch(ch){
		case 'd':
			diff /= inc;break;
		case 'i':
			diff *= inc;break;
		case 'w':
			if(inc >= 1.01)
				inc -= 0.01;
			break;
		case 's':
			if(inc <= 1.99)
				inc += 0.01;
		case 'p':
			cout << "diff:" << diff << ", inc:" << inc << endl;
			break;
	}
	glutPostRedisplay();
}
void display(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glBegin(GL_TRIANGLES);
	glColor3f(1.0, 0.0, 0.0);
	glVertex3f(1.0, 0.0, 0.0);
	glVertex3f(0.0, 1.0, 0.0);
	glVertex3f(0.0, 0.0, 1.0);

	glColor3f(0.0, 1.0, 0.0);
	glVertex3f(1 + diff, 0.0, 0.0);
	glVertex3f(0.0, 1 - diff, 0.0);
	glVertex3f(0.0, 0.0, 1 + diff);	
	glEnd();
	glutSwapBuffers();
}
void reshape(int width, int height){
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1, 1, -1, 1, 0, 2);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(1, 1, 1, 0, 0, 0, 0, 1, 0);
}

int main(int argc, char *argv[]){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(800, 800);
	glutCreateWindow("z fighting");
	init();
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutMainLoop();
}
