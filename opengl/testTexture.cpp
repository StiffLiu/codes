#include <iostream>
#include "Camera.h"
#include "Drawer.h"
#include "Color.h"
#include <GL/glut.h>
#include <Magick++.h>
#include <cassert>
using namespace std;
using namespace Magick;
using namespace tmlg;
double t = 0.1;
struct Texture{
	int textWidth;
	int textHeight;
	GLubyte *texImage;
};
Texture textures[6];
static CameraF camera(Point3F(0, 0, 8), Point3F(0, 0, 4), Point3F(0, 1, 0));
static const int width = 256;
ColorF ramp[width];
Texture texture;
float (*tmpFunc)(float x, float y);
Point3F (*normFunc)(float x, float y);
int slice = 400;
void display(){
	GLfloat texParams[4] = {0, 0, 0.5, 0.5};
	const float *eyePos = camera.getEyePos();
	const float *target = camera.getTarget();
	const float *upVector = camera.getUpVector();
	GLfloat white[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat yellow[]= {1.0, 1.0, 0.0, 1.0};
	GLfloat mat_shininess[]={ 30.0 };
	SurfaceGeometryF tmp;
	Drawer drawer;
	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	tmp.startX = tmp.startY = -4;
	tmp.endX = tmp.endY = 4;
	tmp.rows = tmp.cols = slice;
	tmp.evalFunc = tmpFunc;
	tmp.normFunc = normFunc;
	tmp.hasTexture = true;
	if(tmp.normFunc == NULL)
		tmp.calculateNorm = true;
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(eyePos[0], eyePos[1], eyePos[2], target[0], target[1], target[2], 
			upVector[0], upVector[1], upVector[2]);
	//glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);
	//glTexGenfv(GL_S, GL_OBJECT_PLANE, texParams);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, white );
        glMaterialfv(GL_BACK, GL_DIFFUSE, yellow );
        glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess );
	{
		GLuint texName;
		glGenTextures(1,&texName);                 // define texture for sixth face
		//glEnable(GL_TEXTURE_2D);
		glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
		glBindTexture(GL_TEXTURE_2D,texName);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
		glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA8,texture.textWidth,texture.textHeight,
		             0,GL_RGB,GL_UNSIGNED_BYTE,texture.texImage);
		drawer.draw(tmp);
		glDeleteTextures(1, &texName);
	}
	glutSwapBuffers();
}
void reshape(int width, int height){
	//int len = min(width, height);
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(30, 1, 1, 100);
}
float func1(float x, float y)
{
	return cos(x + y + t) + cos(x * x + y * y + t * 2);
}
Point3F normFunc1(float x, float y){
	float tmp0 = sin(x + y + t);
	float tmp1 = sin(2*x + 3 * y + t * 2);
	Point3F pt(tmp0 + 2 * x * tmp1, tmp0 + 2 * y * tmp1, 1);
	pt.normalize();
	return pt;
}
float func2(float x, float y)
{
	float tmp = x * x + y * y;
	float norm = (1 + sin(t));
	if( tmp < 0.1)
		return norm * 10;
	return norm / tmp;
}
Point3F normFunc2(float x, float y)
{
	float tmp = x * x + y * y;
	float norm = (1 + sin(t));
	if( tmp < 0.1)
		return Point3F(0, 0, 1);
	Point3F pt(2 * norm * x / tmp, 2 * norm * y / tmp, -1);
	pt.normalize();
	return pt;	
}
float func3(float x, float y){
	return rand() * 10 / (float)RAND_MAX;
}
float func4(float x, float y){
	return exp(x + y);
}
Point3F normFunc4(float x, float y)
{
	float tmp = exp(x+y);
	Point3F pt(-tmp, -tmp, -1);
	pt.normalize();
	return pt;	
}
float func5(float x, float y){
	return exp(sin(x * x + y * y + t));
}
Point3F normFunc5(float x, float y)
{
	float tmp = x * x + y * y + t;
	float tmp1 = cos(tmp);
	float tmp2 = tmp1 * exp(sin(tmp));
	Point3F pt(-2 * x * tmp2, -2 * y * tmp2, 1);
	pt.normalize();
	return pt;	
}
const char *greyImage = "hs-2011-05-d-web.jpg";
float func6(float x, float y){
	static class TmpClass{
	public:
		Texture texture;
		unsigned int mx;
		unsigned int mn;
		TmpClass(){
			Image tmp;
			GLubyte *texImage = NULL;
			mx = 0;mn = 255;
			tmp.read(greyImage);
			texture.textWidth = tmp.columns();
			texture.textHeight = tmp.rows();
			texture.texImage = new GLubyte[texture.textHeight * texture.textWidth * 3];
			texImage = texture.texImage;
			for(unsigned int i = 0;i < texture.textHeight;++ i){
				for(unsigned int j = 0, index = 3 * i * texture.textWidth;j < texture.textWidth;++ j, index += 3){
					Magick::Color c = tmp.pixelColor(j, texture.textHeight - i);
					texImage[index] = (unsigned char)c.redQuantum();
					texImage[index + 1] = (unsigned char)c.greenQuantum();
					texImage[index + 2] = (unsigned char)c.blueQuantum();
					//texImage[index + 3] = rand() / (float)RAND_MAX;
					if(texImage[index] > mx)
						mx = texImage[index];
					if(texImage[index] < mn)
						mn = texImage[index];
				}
			}
			cout << mx << ' ' << mn << endl;
			tmp.display();
		}
		~TmpClass(){
			delete [] texture.texImage;
		}
	}instance;
	int i = (x + 4) * instance.texture.textWidth / 8 + 0.5;
	int j = (y + 4) * instance.texture.textHeight / 8 + 0.5;
	if(i == instance.texture.textWidth || j == instance.texture.textHeight)
		return 0;
	int index = 3 * (j * instance.texture.textWidth + i);
	float tmp = instance.texture.texImage[index];
	assert(j <= instance.texture.textHeight);
	assert(i <= instance.texture.textWidth);
	assert(instance.texture.texImage[index] == 
		instance.texture.texImage[index + 1]);
	assert(instance.texture.texImage[index + 1] == 
		instance.texture.texImage[index + 2]);
	tmp += instance.texture.texImage[index + 1];
	tmp += instance.texture.texImage[index + 2];
	tmp /= (3.0 * (instance.mx - instance.mn));
	return tmp * 4 - 2;
}
struct FunStruct{
	float (*func)(float x, float y);
	Point3F (*normFunc)(float x, float y);
};
static FunStruct funcs[] = {
	{func1, normFunc1}, {func2, normFunc2},{func3, NULL},{func4, normFunc4},{func5, normFunc5}, {func6, NULL},
	};
static int current = 0 ;
static void keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'w' :camera.orbit(1.0 / 180 * ConstantF::pi);break;
	case 'q' :camera.orbit(-1.0 / 180 * ConstantF::pi);break;
	case 's' :camera.zoom(0.5);break;
	case 'a' :camera.zoom(-0.5);break;
	case 'x' :camera.roll(1.0 / 180 * ConstantF::pi);break;
	case 'z' :camera.roll(-1.0 / 180 * ConstantF::pi);break;
	case 'u' :camera.movePositionY(0.1);break;
	case 'd' :camera.movePositionY(-0.1);break;
	case 'n' :
		current = (current + 1) % (sizeof funcs / sizeof *funcs);
		tmpFunc = funcs[current].func;
		normFunc = funcs[current].normFunc;break;
	case 't':{
		static bool enabled = true;
		if(enabled)
			glDisable(GL_TEXTURE_2D);
		else
			glEnable(GL_TEXTURE_2D);
		enabled = !enabled;
		break;
		}
	}
	glutPostRedisplay();
}
void makeRamp(void)
{
	// Make color ramp for 1D texture: starts at 0, ends at 240, 256 steps
	for (int i=0; i<width; i++) {
		ColorF tmp((float)i/255.0, 1.0, 0.5);
		ramp[i] = tmp.toRGB();
	}
}
void init(const char *fileName){
	GLfloat light_position[]={ 0.0, 0.0, 10.0, 1.0 };
        GLfloat light_color[]   ={ 1.0, 1.0, 1.0, 1.0 };
        GLfloat ambient_color[] ={ 0.9, 0.3, 0.3, 1.0 };
        GLfloat mat_specular[]  ={ 1.0, 1.0, 1.0, 1.0 };
	int i = 1;

	glClearColor(0.1, 0.1, 0.1, 0);
	glEnable(GL_DEPTH);
	glShadeModel(GL_SMOOTH);
	glLightModeliv(GL_LIGHT_MODEL_TWO_SIDE, &i ); // two-sided lighting
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular );
        glLightfv(GL_LIGHT0, GL_POSITION, light_position );
        glLightfv(GL_LIGHT0, GL_AMBIENT, ambient_color );
        glLightfv(GL_LIGHT0, GL_SPECULAR, light_color );
        glLightfv(GL_LIGHT0, GL_DIFFUSE, light_color );
	makeRamp();
	/*glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameterf(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexImage1D(GL_TEXTURE_1D, 0, 3, width, 0, GL_RGB, GL_FLOAT, ramp[0]);
	glEnable(GL_TEXTURE_GEN_S);
	glEnable(GL_TEXTURE_1D);*/
	tmpFunc = funcs[current].func;
	normFunc = funcs[current].normFunc;
	{
		Image tmp;
		GLubyte *texImage = NULL;
		tmp.read(fileName);
		texture.textWidth = tmp.columns();
		texture.textHeight = tmp.rows();
		texture.texImage = new GLubyte[texture.textHeight * texture.textWidth * 3];
		texImage = texture.texImage;
		for(unsigned int i = 0;i < texture.textHeight;++ i){
			for(unsigned int j = 0, index = 3 * i * texture.textWidth;j < texture.textWidth;++ j, index += 3){
				Magick::Color c = tmp.pixelColor(j, texture.textHeight - i);
				texImage[index] = (unsigned char)c.redQuantum();
				texImage[index + 1] = (unsigned char)c.greenQuantum();
				texImage[index + 2] = (unsigned char)c.blueQuantum();
				//texImage[index + 3] = rand() / (float)RAND_MAX;
			}
		}
	}
        glEnable(GL_TEXTURE_2D); // allow 2D texture maps
	glEnable(GL_LIGHTING);   // so lighting models are used
        glEnable(GL_LIGHT0);     // we'll use LIGHT0
}
void idle(){
	t += 0.1;
	glutPostRedisplay();
}
int main(int argc, char *argv[]){
	const char *fileName = "/home/tmlg/Pictures/W7 Walls/Australia/3.jpg";
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
	glutInitWindowSize(700, 700);
	glutInitWindowPosition(10, 10);
	glutCreateWindow("test texture");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutIdleFunc(idle);
	if(argc >= 2)
		fileName = argv[1];
	cout << fileName << endl;
	init(fileName);
	glutMainLoop();
}
