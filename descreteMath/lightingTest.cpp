#include <GL/glut.h>
#include <iostream>
#include <sstream>
using namespace std;
static int width = 900, height = 600;
static void init(){
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glClearColor(0, 0, 0, 1);

	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
}
static void showInfo(int xLeft, int yBottom, int winWidth, const char *info )
{
	int textHeight = 20;
	glPushMatrix();
	glViewport(xLeft, yBottom + textHeight / 3 , winWidth, textHeight);
	glLoadIdentity();
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	gluOrtho2D(0, winWidth, 0, textHeight);
	glPushMatrix();
	glTranslatef( winWidth / 3, 5, 0);
	glScalef( 0.1, 0.1, 0.1 );
	glColor4f( 1.0, 0.0, 0.0, 1.0 );
	while(*info){
		glutStrokeCharacter(GLUT_STROKE_ROMAN  , *info);
		++ info;
	}
	glPopMatrix();
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glColor4f( 0.8, 0.8, 0.8, 0.7 );
	glRectf(1, 1,  winWidth - 1, textHeight - 1);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glPopMatrix();
}
static void reshape(int width, int height){
	::width = width;
	::height = height;
}
static void setTransformation(){
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1.5, 1.5, -1.5, 1.5, 0, 4);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(1, 1, 1, 0, 0, 0, 0, 1, 0);
}
static string setSpecLightAndMat(int i){
	GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat mat_shininess[][1]= {{25.0}, {50.0}, {100.0}, {150.0}, {200.0}, {500.0}};
	GLfloat light_position[4] = {1.0, 1.0, 1.0, 0.0};
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess[i]);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLighti(GL_LIGHT0, GL_SPOT_CUTOFF, 90);
	glColorMaterial(GL_FRONT,GL_DIFFUSE);
	ostringstream os;
	os << "shiness:" << mat_shininess[i][0];
	return os.str();
}
static string setPositionalLightAndMat(int i){
	GLint constant = 1, linear = 0, quadratic = 0;
	GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat mat_shininess[][1]= {{25.0}};
	GLfloat light_position[4] = {2.0, 2.0, 2.0, 1.0};
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialf(GL_FRONT, GL_SHININESS, 25.0);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLighti(GL_LIGHT0, GL_SPOT_CUTOFF, 90);
	switch(i){
		case 0:break;
		case 1:constant = 0;linear = 1;break;
		case 2:constant = 0;quadratic = 1;break;
		case 3:constant *= 2;break;
		case 4:constant = 0;linear = 2;break;
		case 5:constant = 0;quadratic = 2;break;
	}
	glLighti(GL_LIGHT0, GL_CONSTANT_ATTENUATION, constant);
	glLighti(GL_LIGHT0, GL_LINEAR_ATTENUATION, linear);
	glLighti(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, quadratic);
	glColorMaterial(GL_FRONT,GL_DIFFUSE);
	ostringstream os;
	os << "attenuation:C-" << constant << ", L-" << linear << ", Q-" << quadratic;
	return os.str();
}
static string setSpotLight(int i){
	GLfloat mat_specular[] = {1.0, 1.0, 1.0, 1.0};
	GLint cutoff[][1]= {{5}, {10}, {20}, {35}, {55}, {80}};
	GLfloat light_position[4] = {2.0, 2.0, 2.0, 1.0};
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialf(GL_FRONT, GL_SHININESS, 25.0);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightiv(GL_LIGHT0, GL_SPOT_CUTOFF, cutoff[i]);
	glLighti(GL_LIGHT0, GL_CONSTANT_ATTENUATION, 1);
	glLighti(GL_LIGHT0, GL_LINEAR_ATTENUATION, 0);
	glLighti(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, 0);
	glColorMaterial(GL_FRONT,GL_DIFFUSE);
	ostringstream os;
	os << "cutoff:" << cutoff[i][0];
	return os.str();
}
static string (*func)(int) = setSpecLightAndMat;
void display(){
	int halfWidth = width / 3;
	int halfHeight = height / 2;
	int sideLen = (halfWidth > halfHeight ? halfHeight : halfWidth);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	for(int i = 0;i < 2;++ i)
		for(int j = 0;j < 3;++ j){
			glPushMatrix();
			glViewport(j * sideLen, i * sideLen, sideLen, sideLen);
			setTransformation();
			string ret = func(i * 3 + j);
			glColor3f(1.0, 0.0, 0.0);
			glutSolidSphere(1.0, 30, 30);
			showInfo(j * sideLen, i * sideLen, sideLen, ret.c_str());
			glPopMatrix();	
		}
	glutSwapBuffers();
}
void keyboard(unsigned char ch, int x, int y){
	switch(ch){
		case 'a':func = setSpecLightAndMat;cout << "specular" << endl;break;
		case 's':func = setPositionalLightAndMat;cout << "point" << endl;break;
		case 'd':func = setSpotLight;cout << "spot" << endl;break;
	}
	glutPostRedisplay();
}
int main(int argc, char *argv[]){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(width, height);
	glutCreateWindow("Lighting Test");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	init();
	glutMainLoop();
}
