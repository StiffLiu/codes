#include <GL/glut.h>
#include <Magick++.h>
#include <iostream>
#include "opengl/Point3.h"
#include "opengl/Camera.h"
#include <vector>
#include <string>
#include <cstdlib>
#include <cassert>
#include <ctime>
using namespace std;
using namespace Magick;
using namespace tmlg;
struct Texture{
	int textWidth;
	int textHeight;
	GLubyte *texImage;
};
vector<string> fileNames;
Texture textures[6];
void init(){
	glEnable(GL_DEPTH_TEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_NORMALIZE);
	glEnable(GL_SMOOTH);
        glEnable(GL_TEXTURE_2D); // allow 2D texture maps
	srand(time(0));
	for(int k = 0;k < fileNames.size();++ k){
		Image texture;
		GLubyte *texImage = NULL;
		texture.read(fileNames[k].c_str());
		textures[k].textWidth = texture.columns();
		textures[k].textHeight = texture.rows();
		textures[k].texImage = new GLubyte[textures[k].textHeight * textures[k].textWidth * 4];
		texImage = textures[k].texImage;
		for(unsigned int i = 0;i < textures[k].textHeight;++ i)
			for(unsigned int j = 0, index = 3 * i * textures[k].textWidth;j < textures[k].textWidth;++ j, index += 3){
				Color c = texture.pixelColor(j, textures[k].textHeight - i);
				texImage[index] = (unsigned char)c.redQuantum();
				texImage[index + 1] = (unsigned char)c.greenQuantum();
				texImage[index + 2] = (unsigned char)c.blueQuantum();
				texImage[index + 3] = rand() / (float)RAND_MAX;
			}
	}
}
GLenum textMode[] = {GL_DECAL, GL_BLEND, GL_MODULATE, GL_REPLACE};
GLenum filterMode[] = {GL_LINEAR, GL_NEAREST};
int currentFilter = 0;
const char * textModeNames[] = {"GL_DECAL", "GL_BLEND", "GL_MODULATE", "GL_REPLACE"};
const char * filterModeNames[] = {"GL_LINEAR", "GL_NEAREST"};
int current = 0;
void drawOneQuad(const Point3F& p1, const Point3F& p2, const Point3F& p3, const Point3F& p4,
	int textWidth, int textHeight, GLubyte *texImage){
	GLuint texName;
	glGenTextures(1,&texName);                 // define texture for sixth face
	glEnable(GL_TEXTURE_2D);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, textMode[current]);
	glBindTexture(GL_TEXTURE_2D,texName);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,filterMode[currentFilter]);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,filterMode[currentFilter]);
	glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA8,textWidth,textHeight,
                     0,GL_RGB,GL_UNSIGNED_BYTE,texImage);
	glBegin(GL_QUADS);
		glTexCoord2f(0.0, 0.0); 
		glVertex3fv(p1);
		glTexCoord2f(0.0, 1.0);  
		glVertex3fv(p2);
		glTexCoord2f(1.0, 1.0); 
		glVertex3fv(p3);
		glTexCoord2f(1.0, 0.0);  
		glVertex3fv(p4);
	glEnd();
	glDeleteTextures(1, &texName);
}
static CameraF camera(Point3F(0, 0, 3), Point3F(0, 0, 0), Point3F(0, 1, 0));
float  len2 = 9, dist = 3;
void display(){
	const float *eyePos = camera.getEyePos();
	const float *target = camera.getTarget();
	const float *upVector = camera.getUpVector();
	float len = len2, margin = 0.1;
	int counts = 2;
	int total = fileNames.size() / 2;
	float angle = 360.0 / total;
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(eyePos[0], eyePos[1], eyePos[2], target[0], target[1], target[2], 
			upVector[0], upVector[1], upVector[2]);
	for(int i = 0;i < total;++ i){
		glPushMatrix();
		glRotatef(i * angle, 0, 1, 0);
		for(int j = 0;j < counts;++ j){
			int index = i * counts + j;
			float zStart = j * (len + margin);
			float zEnd = zStart + len;
			drawOneQuad(Point3F(0, -len2, zStart), Point3F(0, len2, zStart), Point3F(0, len2, zEnd), Point3F(0, -len2, zEnd), 
					textures[index].textWidth, textures[index].textHeight, textures[index].texImage);
		}
		glPopMatrix();
	}
	glutSwapBuffers();
}

float ratio = 3;
void reshape(int width, int height){
	glViewport(0, 0, width, height);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//gluOrtho2D(-0.1 * textWidth, textWidth * 1.1, -0.1 * textHeight, textHeight * 1.1);
	//gluOrtho2D(0, textWidth, 0, textHeight);
	gluPerspective(60.0,ratio,1.0,1000.0);
}
static void keyboard(unsigned char key, int x, int y)
{
	assert(camera.isValid());
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
	case '+' :len2+=1;break;
	case 'i' :ratio += 0.1;break;
	case 'n' :current += 1;current %= (sizeof textMode/ sizeof *textMode);
		cout << textModeNames[current] << endl;
		break;
	case 't' :currentFilter += 1;currentFilter %= (sizeof filterMode/ sizeof *filterMode);
		cout << filterModeNames[currentFilter] << endl;
		break;
	}
	assert(camera.isValid());
	glutPostRedisplay();
}
int main(int argc, char *argv[]){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(800, 800);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Draw Function");
	glutDisplayFunc(display);
	glutKeyboardFunc(keyboard);
	glutReshapeFunc(reshape);
	string tmp;
	while(getline(cin, tmp))
		fileNames.push_back(tmp);
	init();
	glutMainLoop();
	return 0;
}
