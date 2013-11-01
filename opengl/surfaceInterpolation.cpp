#include <GL/glut.h>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cmath>
#include <iostream>
#include "Point3.h"
#include "Camera.h"
#include "Drawer.h"
#include "Color.h"
using namespace std;
using namespace tmlg;
struct CubicSurface{
	float surf[10];
	ColorAF color;
};
class SurfInterpolation{
	CameraF camera;
	vector<CubicSurface> surfs;
	CubicSurface interpolated;
	ColorAF lit_clr;
	static float randFloat(){
		double tmp = rand() / (double)RAND_MAX;
		return 2 * tmp - 1;
	}
	int current;
	int total;
	int all;
	static float * surf;
	static float calc(float x, float y){
		float x2 = x * x;
		float y2 = y * y;
		return x2 * (surf[0] * x + surf[1] * y + surf[2]) +
			y2 * (surf[3] * x + surf[4] * y + surf[5]) +
			(surf[6] * x + surf[7] * y + surf[8]) +
			surf[9];
	}
	static void calcRainbow(double value, ColorAF& c){
		float& r = c.x();
		float& g = c.y();
		float& b = c.z();
		c.a() = 1;
		if (value < 0.2)	// purple to blue ramp
			{ r=0.5*(1.0-value/0.2);g=0.0;b=0.5+(0.5*value/0.2);return;}
		if ((value >= 0.2) && (value < 0.40))	// blue to cyan ramp
			{ r= 0.0; g= (value-0.2)*5.0; b = 1.0; return; }
		if ((value >= 0.40) && (value < 0.6))	// cyan to green ramp
			{ r= 0.0; g= 1.0; b = (0.6-value)*5.0; return; }
		if ((value >= 0.6) && (value < 0.8))	// green to yellow ramp
			{ r= (value-0.6)*5.0; g= 1.0; b = 0.0; return; }
		if (value >= 0.8)	// yellow to red ramp
			{ r= 1.0; g= (1.0-value)*5.0; b= 0.0; }
		return;
	}
public:
	SurfInterpolation() : camera(Point3F(0, -1, 0), Point3F(0, 0, 0), Point3F(0, 0, 1)){
		reset();
	}
	void reset(){
		int count = 4;
		srand(time(0));
		surfs.resize(count);
		for(int i = 0;i < count;++ i){
			surfs[i].surf[i] =  randFloat() * 2;
			//cout << i << endl;
			calcRainbow(i / (double)(count - 1), surfs[i].color);
			/*cout << surfs[i].color.x() << ' ' << surfs[i].color.y() << ' ' << surfs[i].color.z() << endl;
			if(i % 2 == 0)
				surfs[i].color = ColorAF(1.0, 0.0, 0.0, 1.0);
			else
				surfs[i].color = ColorAF(0.0, 0.0, 1.0, 1.0);*/
		}
		total = 1000;
		all = count * total;
		current = 0;
		//clrA = ColorAF(1.0, 0.0, 0.0, 1.0);
		//clrB = ColorAF(0.0, 0.0, 1.0, 1.0);
		lit_clr = ColorAF(1.0, 1.0, 1.0, 1.0); // and red
	}
	void mat(ColorAF * clr){
		GLfloat shiness = 10;
	   	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, (const GLfloat*)clr );
	   	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, (const GLfloat*)clr );
	   	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, (const GLfloat*)&lit_clr );
	   	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, &shiness );
	}
	void draw(){
		const float *eyePos = camera.getEyePos();
		const float *target = camera.getTarget();
		const float *upVector = camera.getUpVector();
		int index = current / total;
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt(eyePos[0], eyePos[1], eyePos[2], target[0], target[1], target[2], 
			upVector[0], upVector[1], upVector[2]);
		SurfaceGeometryF tmp;
		Drawer drawer;
		tmp.rows = tmp.cols = 100;
		tmp.startX = tmp.startY = -1.0;
		tmp.endX = tmp.endY = 1.0;
		tmp.calculateNorm = true;
		tmp.evalFunc = calc;
		for(size_t i = 0;i < surfs.size();++ i){
			surf = surfs[i].surf;
			mat(&surfs[i].color);
			//cout << surfs[i].color.x() << ' ' << surfs[i].color.y() << ' ' << surfs[i].color.z() << endl;
			drawer.draw(tmp);
		}
		if(index + 1 < surfs.size()){
			surf = interpolated.surf;			
			mat(&interpolated.color);
			drawer.draw(tmp);
		}
	}
	void kbd(unsigned char key){
		switch (key){
			case 'w' :camera.orbit(1.0 / 180 * ConstantF::pi);break;
			case 'q' :camera.orbit(-1.0 / 180 * ConstantF::pi);break;
			case 's' :camera.zoom(0.5);break;
			case 'a' :camera.zoom(-0.5);break;
			case 'x' :camera.roll(1.0 / 180 * ConstantF::pi);break;
			case 'z' :camera.roll(-1.0 / 180 * ConstantF::pi);break;
			case 'u' :camera.movePositionY(0.1);break;
			case 'd' :camera.movePositionY(-0.1);break;
			case 'r' :reset();break;
		}
	}
	void next(){
		++current;
		int index = current / total;
		if(index + 1 >= surfs.size())
			return;
		double ratio = sin((current % total) / (double)total * ConstantF::pi / 2);
		float *surfA = surfs[index].surf;
		float *surfB = surfs[index + 1].surf;
		ColorAF clrA = surfs[index].color;
		ColorAF clrB = surfs[index + 1].color;
		double one_minus_ratio = 1 - ratio;
		for(int i = 0;i < 10;++ i){
			interpolated.surf[i] = one_minus_ratio * surfA[i] + ratio * surfB[i];
		}
		interpolated.color = ColorAF(one_minus_ratio * clrA.x() + ratio * clrB.x(),
				one_minus_ratio * clrA.y() + ratio * clrB.y(),
				one_minus_ratio * clrA.z() + ratio * clrB.z(),
				one_minus_ratio * clrA.a() + ratio * clrB.a());
	}

}instance;
float *SurfInterpolation::surf = NULL;
GLfloat light_pos0[]={  1.0, 1.0, 1.0,  1.0 }; // first light over z-axis
GLfloat light_col0[]={  1.0,  1.0,  1.0,  1.0 }; // and red
GLfloat amb_color0[]={  0.3,  0.0,  0.0,  1.0 }; // even ambiently

GLfloat light_pos1[]={  -1.0, -1.0, -1.0,  1.0 }; // first light over z-axis
GLfloat light_col1[]={  1.0,  1.0,  1.0,  1.0 }; // and red
GLfloat amb_color1[]={  0.3,  0.0,  0.0,  1.0 }; // even ambiently
void init(){
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
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
}
void display(){
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	instance.draw();
	glAccum( GL_MULT,   0.90 );
	glAccum( GL_ACCUM,  0.10 );
	glAccum( GL_RETURN, 1. );
	glutSwapBuffers();
}
void keyboard(unsigned char ch, int x, int y){
	instance.kbd(ch);
}
void reshape(int w, int h){
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(30.0f, 1.0f, 1.0f, 20.0f);
}
void animate(int timerId){
	glutTimerFunc(20, animate, timerId);
	instance.next();
	glutPostRedisplay();
}
int main(int argc, char *argv[]){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_ACCUM);
	glutInitWindowSize(400, 400);
	glutInitWindowPosition(100, 100);
	glutCreateWindow("Interpolating surface");

	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);
	glutTimerFunc(20, animate, 100);
	init();

	glutMainLoop();
}
