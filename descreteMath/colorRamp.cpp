#include <GL/glut.h>
#include <cmath>
static unsigned int count = 5;
static unsigned long long values[100 * 100];
static int winId = 0;
double identity(double x){
	return x;
}
static double(*func)(double) = log;
static double h = 1.0;
static double s = 1.0;
extern void hsl2rgb(double h, double sl, double l, double& r, double& g, double& b);
// n
//m
static void compute(unsigned int n, unsigned long long *ret){
	for(unsigned int i = 0;i < n;++ i)
		ret[i] = 1;
	for(unsigned int i = 0;i < n;++ i)
		ret[i * n] = 1;
	for(unsigned int i = 1;i < n;++ i){
		for(unsigned int j = 1;j < i;++ j){
			unsigned int temp = ret[i * n + j - 1] + ret[(i - j - 1) * n + j];
			ret[i * n + j] = temp;
		}
		ret[i * n + i] = 1 + ret[i * n + i - 1];
		for(unsigned int j = i + 1;j < n;++ j)
			ret[i * n + j] = ret[i * n + i];
	}
}
#include <fstream>
#include <iomanip>
static void out(std::ostream& os, unsigned int n, unsigned long long *ret){
	for(unsigned int i = 0;i < n;++ i)
		for(unsigned int j = 0;j < n;++ j)
			os  << ret[i * n + j] << (j == n -1 ? '\n' : '\t');
}
static void init(){
	glClearColor(0, 0, 0, 1);
	glEnable(GL_SHADE_MODEL);
	glMatrixMode(GL_MODELVIEW);
	glDisable(GL_POLYGON_SMOOTH);
	glShadeModel(GL_FLAT);	
	glEnable(GL_BLEND);
	glLoadIdentity();
	gluLookAt(0, 0, 1, 0, 0, 0, 0, 1, 0);
	compute(count, values);
}
static bool useRgb = false;
static void calcUniform(double value, double& r, double& g, double& b){
	hsl2rgb(h, s, value, r, g, b);
}
void calcRainbow(double value, double& r, double& g, double& b)
{
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

void calcMulti(double value, double& r, double& g, double& b)
{
	if (value < 0.2)	// purple to blue ramp, starting with black
		{ r=5.0*value*0.5*(1.0-value/0.2);
		  g=0.0;b=5.0*value*(0.5+(0.5*value/0.2));return;}
	if ((value >= 0.2) && (value < 0.40))	// blue to cyan ramp, starting with black
		{ r= 0.0; g= 5.0*(value-0.2)*(value-0.2)*5.0;
		  b = 5.0*value*1.0; return; }
	if ((value >= 0.40) && (value < 0.6))	// cyan to green ramp, starting with black
		{ r= 0.0; g= 5.0*(value-0.4)*1.0;
		  b = 5.0*(value-0.4)*(0.6-value)*5.0; return; }
	if ((value >= 0.6) && (value < 0.8))	// green to yellow ramp, starting with black
		{ r= 5.0*(value-0.6)*(value-0.6)*5.0; 
		  g= 5.0*(value-0.6)*1.0; b = 0.0; return; }
	if (value >= 0.8)	// yellow to red ramp, starting with black
		{ r= 5.0*(value-0.8)*1.0;
		  g= 5.0*(value-0.8)*(1.0-value)*5.0; b= 0.0; }
	return;
}
void calcRWB(double value, double& r, double& g, double& b)
{
	//printf("in calcRWB value is %f\n",value);
	if (value < 0.5)
	{
		r = 1.; // red
		g = 2. * value; // green
		b = 2. * value; // blue
	}
	else
	{
		r = 2. * (1. - value); // red
		g = 2. * (1. - value); // green
		b = 1.; // blue
	}
}
void calcHeat(double value, double& r, double& g, double& b)
{
	if (fabs(value-0.5)<0.02){r=g=b=0.;return;}
	if (value < 0.30) { r=value/0.3; g=0.0; b = 0.0; return; }
	if ((value >= 0.30) && (value < 0.89))
		{ r=1.0; g=(value-0.3)/0.59; b = 0.0; return; }
	if (value >= 0.89) { r=1.0; g=1.0; b=(value-0.89)/0.11; }
	return;
}
void calcEven(double value, double& r, double& g, double& b)
//	commented out line creates a green "bathtub ring" in the surface
{
	if (value < 0.33) { r=value/0.3; g=0.0; b = 0.0; return; }
	if ((value >= 0.33) && (value < 0.66))
		{ r=1.0; g=(value-0.33)/0.33; b = 0.0; return; }
	if (value >= 0.66) { r=1.0; g=1.0; b=(value-0.66)/0.33; }
	return;
}
static void (*rampFunc)(double value, double& r, double& g, double& b) = calcUniform;
static void nextRampFunc(){
	static void (*allFunc[])(double, double&, double&, double&)={
		calcUniform,calcRainbow,calcMulti,calcRWB,calcHeat,calcEven,
	};
	static const char *names[]={
		"calcUniform", "calcRainbow", "calcMulti", "calcRWB", "calcHeat", "calcEven",
	};
	static int current = -1;
	++current;current %= (sizeof allFunc / sizeof *allFunc);
	rampFunc = allFunc[current];
	glutSetWindowTitle(names[current]);
}
static void display(){
	double hasfSide = 3.5;
	double interVal = 2 * hasfSide / count;
	double max = func((double)values[count * count - 1]);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	if(useRgb)
		for(unsigned int i = 0;i < count - 1;++ i){
			glBegin(GL_QUAD_STRIP);
			for(unsigned int j = 0;j < count;++ j){
				double r, g, b;
				unsigned long long val = values[j * count + i];
				//hsl2rgb(h, s, func((double)values[j * count + i]) / max, r, g, b);
				//glColor3d(r, g, b);
				glColor4f((val % 255 / 255.0), (val >> 8) % 255 / 255.0, 
					(val >> 12) % 255 / 255.0, (val >> 16) % 255 / 255.0);
				glVertex3d(-hasfSide + j * interVal, hasfSide - i * interVal, 0);
				//hsl2rgb(h, s, func((double)values[j * count + i + 1]) / max, r, g, b);
				val = values[j * count + i + 1];
				glColor4f((val % 255 / 255.0), (val >> 8) % 255 / 255.0, 
					(val >> 12) % 255 / 255.0, (val >> 16) % 255 / 255.0);
				//glColor3d(r, g, b);
				glVertex3d(-hasfSide + j * interVal, hasfSide - (i + 1) * interVal, 0);
			}
			glEnd();
		}
	else
		for(unsigned int i = 0;i < count - 1;++ i){
			glBegin(GL_QUAD_STRIP);
			for(unsigned int j = 0;j < count;++ j){
				double r, g, b;
				rampFunc(func((double)values[j * count + i]) / max, r, g, b);
				glColor3d(r, g, b);
				glVertex3d(-hasfSide + j * interVal, hasfSide - i * interVal, 0);
				rampFunc(func((double)values[j * count + i + 1]) / max, r, g, b);
				glColor3d(r, g, b);
				glVertex3d(-hasfSide + j * interVal, hasfSide - (i + 1) * interVal, 0);
			}
			glEnd();
		}
	glutSwapBuffers();
}
static void reshape(int width, int height){
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-4, 4, -4, 4, 0, 1);
}
static void keyboardFunc(unsigned char ch, int x, int y){
	switch(ch){
		case 's':
			useRgb = !useRgb;break;
		case 'n':
			{
				static double (*allFunc[])(double) = {log, identity, sqrt, };
				static int current = 0;
				current = (current + 1) % (sizeof allFunc / sizeof *allFunc);
				func = allFunc[current];
				break;
			}
		case '+':
			if(count < 100){
				++count;
				compute(count, values);
			}
			break;
		case '-':
			if(count > 5){
				--count;
				compute(count, values);
			}
			break;
		case 'i':
			h += 0.01;if(h > 1.0) h = 1.0;break;
		case 'd':
			s += 0.01;if(s > 1.0) s = 1.0;break;
		case 'r':
			h -= 0.01;if(h < 0) h = 0;break;
		case 'f':
			s -= 0.01;if(s < 0) s = 0;break;
		case 'c':
			nextRampFunc();
			break;
	}
	glutPostRedisplay();
}
int main(int argc, char *argv[]){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	init();
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(400, 400);
	winId = glutCreateWindow("Test Color Ramp");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboardFunc);
	nextRampFunc();
	glutMainLoop();
	return 0;
}
