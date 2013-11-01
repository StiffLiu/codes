/*
	Create a view of the HSV color space for classroom use

	Source file to be used with
	Cunningham, Computer Graphics: Programming in OpenGL for Visual Communication, Prentice-Hall, 2007

	Intended for class use only
*/
#include <glut.h>
#include <cstdlib>
#include <iostream>
#include <cmath>

#define SIZE 20.0

static char ch;
typedef GLfloat point3[3];
typedef GLfloat color [4];
using namespace std;

void hsl2rgb(double h, double sl, double l, double& r, double& g, double& b)
{
	double v;

	r = l;   // default to gray
	g = l;
	b = l;
	v = (l <= 0.5) ? (l * (1.0 + sl)) : (l + sl - l * sl);
	if (v > 0)
	{
		double m;
		double sv;
		int sextant;
		double fract, vsf, mid1, mid2;

		m = l + l - v;
		sv = (v - m ) / v;
		h *= 6.0;
		sextant = (int)h;
		fract = h - sextant;
		vsf = v * sv * fract;
		mid1 = m + vsf;
		mid2 = v - vsf;
		switch (sextant)
		{
		case 0:	r = v;g = mid1;b = m;break;
		case 1:	r = mid2;g = v;b = m;break;
		case 2: r = m;g = v;b = mid1;break;
		case 3:	r = m;g = mid2;b = v;break;
		case 4:	r = mid1;g = m;b = v;break;
		case 5:	r = v;g = m;b = mid2;break;
		}
	}
}
void hsl2rgb(double h, double sl, double l, float& rVal, float& gVal, float& bVal)
{
	double r,g,b;
	hsl2rgb(h, sl, l, r, g, b);
	rVal = r;
	gVal = g;
	bVal = b;
}
void hsl2rgb(double h, double sl, double l, int& rVal, int& gVal, int& bVal)
{
	double r,g,b;
	hsl2rgb(h, sl, l, r, g, b);
	rVal = ceil(r * 255.0f);
	gVal = ceil(g * 255.0f);
	bVal = ceil(b * 255.0f);
}
void rgb2hsl (int rVal, int gVal, int bVal, double& h, double& s, double& l)
{
	double r = rVal/255.0;
	double g = gVal/255.0;
	double b = bVal/255.0;
	double v;
	double m;
	double vm;
	double r2, g2, b2;

	h = 0; // default to black
	s = 0;
	l = 0;
	v = max(r,g);
	v = max(v,b);
	m = min(r,g);
	m = min(m,b);
	l = (m + v) / 2.0;
	if (l <= 0.0)	return;
	vm = v - m;
	s = vm;
	if (s > 0.0)
		s /= (l <= 0.5) ? (v + m ) : (2.0 - v - m) ;
	else
		return;
	r2 = (v - r) / vm;
	g2 = (v - g) / vm;
	b2 = (v - b) / vm;
	if (r == v)
		h = (g == m ? 5.0 + b2 : 1.0 - g2);
	else if (g == v)
		h = (b == m ? 1.0 + r2 : 3.0 - b2);
	else
		h = (r == m ? 3.0 + g2 : 5.0 - r2);
	h /= 6.0;
}
// function prototypes
static void myinit(void);
static void display(void);
static void hslDoubleCone(void);
static void reshape(int,int);
static void keyboard(unsigned char,int,int);

static void myinit(void)
{
	glClearColor( 0.5, 0.5, 0.5, 0.0 );
	glEnable(GL_DEPTH_TEST); // allow z-buffer display
	glShadeModel(GL_SMOOTH);
}

static void display( void )
{
    typedef GLfloat point3[3];
	#define Angle 2.0

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //  NOTE that we do not take ourselves back to the original raw-world state.
    //  We let the modelview matrix hold our global view of the model, and
    //  we modify it by including the new rotation each time a key is pressed
    switch(ch)
    {
		case 'w':
			glRotatef( Angle, 1.0, 0.0, 0.0); break;
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

	//  NOW we save the state of the modelview transformation so we can restore
	//  it after the axes and cube have been drawn.
	glPushMatrix();
	hslDoubleCone();
	glPopMatrix();
	glutSwapBuffers();
 }

static void reshape(int w,int h)
{
	glViewport(0,0,(GLsizei)w,(GLsizei)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0,1.0,1.0,30.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt( 0.0,  10.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
}


static void hslDoubleCone(void)
{
#define NSTEPS 144
#define steps (float)NSTEPS
#define TWOPI 6.28318

	int i;
	double r, g, b;
	
	glBegin(GL_TRIANGLE_FAN);	//	cone of the HSV space
		glColor3f(0.0, 0.0, 0.0);
		glVertex3f(0.0, 0.0, -2.0);
		for (i=0; i<=NSTEPS + 1; i++) {
			hsl2rgb((float)(i%NSTEPS)/steps, 1.0, 0.5, r, g, b);
			glColor3f(r, g, b);
			glVertex3f(2.0*cos(TWOPI*(float)i/steps),2.0*sin(TWOPI*(float)i/steps),0.0);
		}
	glEnd();
	glBegin(GL_TRIANGLE_FAN);	//	top plane of the HSV space
		glColor3f(1.0, 1.0, 1.0);
		glVertex3f(0.0, 0.0, 2.0);
		for (i=0; i<=NSTEPS + 1; i++) {
			hsl2rgb((float)(i%NSTEPS)/steps, 1.0, 0.5, r, g, b);
			glColor3f(r, g, b);
			glVertex3f(2.0*cos(TWOPI*(float)i/steps),2.0*sin(TWOPI*(float)i/steps),0.0);
		}
	glEnd();
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
			ch = key; break;
	}
	glutPostRedisplay();
}

int hsl_double_cone(int argc, char** argv)
{
	glutInit(&argc,argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(500,500);
	glutInitWindowPosition(70,70);
	glutCreateWindow("HSV Color Space");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);

	myinit();
	glutMainLoop();
	return 0;
}