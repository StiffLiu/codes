#include <GL/glut.h>
#include <iostream>
#include <ctime>
#include <cmath>
#include <cassert>
#include <map>
#include <set>
#include "Point3.h"
#include "Transform.h"
#include "Cylinder.h"
#include "Cone.h"
#include "Camera.h"

using namespace std;
using namespace tmlg;
class Drawer{
	bool drawBotCap;
	bool drawTopCap;
	bool drawCap;
	unsigned int slice;
public:
	Drawer(){
		drawBotCap = drawTopCap = drawCap = true;
		slice = 100;
	}
	unsigned int getSlice()
	{
		return slice;
	}
	void draw(const CylinderF& c){
		float *ptBodys = new float[2 * (slice + 1) * 3];
		float *ptNorms = new float[2 * (slice + 1) * 3];
		float *ptBotCap = NULL, *ptTopCap = NULL;
		if(drawBotCap)
			ptBotCap = new float[(slice + 2) * 3];
		if(drawTopCap)
			ptTopCap = new float[(slice + 2) * 3];
		if(c.calcData(ptBodys, ptNorms, ptBotCap, ptTopCap, slice)){
			glBegin(GL_QUAD_STRIP);
			for(int i = 0;i < 2 * (slice + 1);++ i){
				glNormal3fv(ptNorms + 3 * i);
				glVertex3fv(ptBodys + 3 * i);
			}
			glEnd();
			if(ptBotCap != NULL){
				Point3F axis = c.getAxis();
				axis.normalize();
				glBegin(GL_TRIANGLE_FAN);
				for(int i = 0;i < slice + 2;++ i){
					glNormal3fv(axis);
					glVertex3fv(ptBotCap + 3 * i);
				}
				glEnd();
			}
			if(ptTopCap != NULL){
				Point3F axis = c.getAxis();
				axis.normalize();
				axis *= -1;
				glBegin(GL_TRIANGLE_FAN);
				for(int i = 0;i < slice + 2;++ i){
					glNormal3fv(axis);
					glVertex3fv(ptTopCap + 3 * i);
				}
				glEnd();
			}
		}
		delete []ptBodys;
		delete []ptNorms;
		delete []ptBotCap;
		delete []ptTopCap;
	}
	void draw(const ConeF& c){
		float *ptBodys = new float[(slice + 2) * 3];
		float *ptNorms = new float[(slice + 2) * 3];
		float *ptCaps = NULL;
		if(drawCap)
			ptCaps = new float[(slice + 2) * 3];
		if(c.calcData(ptBodys, ptNorms, ptCaps, slice)){
			glBegin(GL_TRIANGLE_FAN);
			for(int i = 0;i < slice + 2;++ i){
				glNormal3fv(ptNorms + 3 * i);
				glVertex3fv(ptBodys + 3 * i);
			}
			glEnd();
			if(ptCaps != NULL){
				Point3F axis = c.getAxis();
				axis.normalize();
				axis *= -1;
				glBegin(GL_TRIANGLE_FAN);
				for(int i = 0;i < slice + 2;++ i){
					glNormal3fv(axis);
					glVertex3fv(ptCaps + 3 * i);
				}
				glEnd();
			}
		}
		delete []ptBodys;
		delete []ptNorms;
		delete []ptCaps;
	}
	void draw(const Point3F& minPt, const Point3F& maxPt, bool hasTexture = false)
	{
		glBegin(GL_QUADS);
		if(hasTexture){
			glNormal3f(-1, 0, 0);
			glTexCoord2f(0.0, 0.0); glVertex3fv(minPt);
			glTexCoord2f(0.0, 1.0); glVertex3f(minPt.x(), minPt.y(), maxPt.z());
			glTexCoord2f(1.0, 1.0); glVertex3f(minPt.x(), maxPt.y(), maxPt.z());
			glTexCoord2f(1.0, 0.0); glVertex3f(minPt.x(), maxPt.y(), minPt.z());
			glNormal3f(0, 1, 0);
			glTexCoord2f(0.0, 0.0); glVertex3f(minPt.x(), maxPt.y(), minPt.z());
			glTexCoord2f(0.0, 1.0); glVertex3f(minPt.x(), maxPt.y(), maxPt.z());
			glTexCoord2f(1.0, 1.0); glVertex3fv(maxPt);
			glTexCoord2f(1.0, 0.0); glVertex3f(maxPt.x(), maxPt.y(), minPt.z());
			glNormal3f(1, 0, 0);
			glTexCoord2f(0.0, 0.0); glVertex3f(maxPt.x(), maxPt.y(), minPt.z());
			glTexCoord2f(0.0, 1.0); glVertex3fv(maxPt);
			glTexCoord2f(1.0, 1.0); glVertex3f(maxPt.x(), minPt.y(), maxPt.z());
			glTexCoord2f(1.0, 0.0); glVertex3f(maxPt.x(), minPt.y(), minPt.z());
			glNormal3f(0, -1, 0);
			glTexCoord2f(0.0, 0.0); glVertex3f(maxPt.x(), minPt.y(), minPt.z());
			glTexCoord2f(0.0, 1.0); glVertex3f(maxPt.x(), minPt.y(), maxPt.z());
			glTexCoord2f(1.0, 1.0); glVertex3f(minPt.x(), minPt.y(), maxPt.z());
			glTexCoord2f(1.0, 0.0); glVertex3fv(minPt);
			glNormal3f(0, 0, 1);
			glTexCoord2f(0.0, 0.0); glVertex3fv(maxPt);
			glTexCoord2f(0.0, 1.0); glVertex3f(minPt.x(), maxPt.y(), maxPt.z());
			glTexCoord2f(1.0, 1.0); glVertex3f(minPt.x(), minPt.y(), maxPt.z());
			glTexCoord2f(1.0, 0.0); glVertex3f(maxPt.x(), minPt.y(), maxPt.z());
			glNormal3f(0, 0, -1);
			glTexCoord2f(0.0, 0.0); glVertex3fv(minPt);
			glTexCoord2f(0.0, 1.0); glVertex3f(minPt.x(), maxPt.y(), minPt.z());
			glTexCoord2f(1.0, 1.0); glVertex3f(maxPt.x(), maxPt.y(), minPt.z());
			glTexCoord2f(1.0, 0.0); glVertex3f(maxPt.x(), minPt.y(), minPt.z());
			glEnd();
		}else{
			glNormal3f(-1, 0, 0);
			glVertex3fv(minPt);
			glVertex3f(minPt.x(), minPt.y(), maxPt.z());
			glVertex3f(minPt.x(), maxPt.y(), maxPt.z());
			glVertex3f(minPt.x(), maxPt.y(), minPt.z());
			glNormal3f(0, 1, 0);
			glVertex3f(minPt.x(), maxPt.y(), minPt.z());
			glVertex3f(minPt.x(), maxPt.y(), maxPt.z());
			glVertex3fv(maxPt);
			glVertex3f(maxPt.x(), maxPt.y(), minPt.z());
			glNormal3f(1, 0, 0);
			glVertex3f(maxPt.x(), maxPt.y(), minPt.z());
			glVertex3fv(maxPt);
			glVertex3f(maxPt.x(), minPt.y(), maxPt.z());
			glVertex3f(maxPt.x(), minPt.y(), minPt.z());
			glNormal3f(0, -1, 0);
			glVertex3f(maxPt.x(), minPt.y(), minPt.z());
			glVertex3f(maxPt.x(), minPt.y(), maxPt.z());
			glVertex3f(minPt.x(), minPt.y(), maxPt.z());
			glVertex3fv(minPt);
			glNormal3f(0, 0, 1);
			glVertex3fv(maxPt);
			glVertex3f(minPt.x(), maxPt.y(), maxPt.z());
			glVertex3f(minPt.x(), minPt.y(), maxPt.z());
			glVertex3f(maxPt.x(), minPt.y(), maxPt.z());
			glNormal3f(0, 0, -1);
			glVertex3fv(minPt);
			glVertex3f(minPt.x(), maxPt.y(), minPt.z());
			glVertex3f(maxPt.x(), maxPt.y(), minPt.z());
			glVertex3f(maxPt.x(), minPt.y(), minPt.z());
			glEnd();
		}
	}
};
class SelectionBuffer{
	GLuint *buffer;
	unsigned int bufSize;
	GLint recordCount;
	SelectionBuffer(const SelectionBuffer&);
	SelectionBuffer& operator=(const SelectionBuffer&);
public:
	SelectionBuffer(unsigned int bufSize){
		buffer = new GLuint[bufSize];
		this->bufSize = bufSize;
		recordCount = -1;
	}
	unsigned int getBufSize()const{
		return bufSize;
	}
	operator GLuint *(){
		return buffer;
	}
	operator const GLuint *() const{
		return buffer;
	}
	void setRecordCount(int recordCount){
		this->recordCount = recordCount;
	}
	int getRecordCount()const {
		return recordCount;
	}
	~SelectionBuffer(){
		delete []buffer;
	}
	int getLeaves(GLuint *leaves ,int bufSize = 0) const{
		int count = recordCount;
		GLuint *current = buffer;
		int index = 0;
		while(index < bufSize && count > 0){
			GLuint cnt = *current;
			current += 3;
			if(cnt > 0){
				leaves[index] = current[cnt - 1];
				++index;
			}
			current += cnt;
			-- count;
		}
		return index;
	}
	friend ostream& operator<<(ostream& os, const SelectionBuffer& buf){
		int count = buf.recordCount;
		const GLuint *current = buf.buffer;
		while(count > 0){
			GLuint cnt = *current;
			if(cnt > 0){
				for(int i = 3;i <= cnt + 2;++ i){
					for(int j = 2;j < i;++ j)
						os << '\t';
					os << current[i];
					if(i == cnt + 2)
						os << ' ' << current[1] << '-' << current[2];
					os << endl;
				}
			}
			current += (cnt + 3);
			-- count;
		}
		return os;
	}
	GLuint getNearest()const{
		GLuint nearest = 0;
		GLuint dist = -1;
		int count = recordCount;
		const GLuint *current = buffer;
		while(count > 0){
			GLuint cnt = *current;
			if(cnt > 0 && (dist > current[1] || dist > current[2])){
				dist = min(current[1], current[2]);
				nearest = current[cnt + 2];
			}
			current += (cnt + 3);
			-- count;
		}
		return nearest;
	}
};
static GLuint lastSelPart = 0;
enum DrawMode{RENDER, SELECTION, BACK_BUFFER};
struct ColorMap{
	static map<GLuint, int> colorMap;
	static GLuint getByColor(int color){
		map<GLuint, int>::iterator begin = colorMap.begin(), end = colorMap.end();
		while(begin != end){
			if(begin->second == color)
				return begin->first;
			++begin;
		}
		return 0;
	}
};
class Car{
	float width;
	float height;
	float length;
	float rInner;
	float rOutter;
	float rCannon;
	float hCannon;
	Point3F dist;//(w, h, l)
public:
	Car(){
		width = 2.0;
		height = 1.0;
		length = 3.5;
		rInner = 1.0/24;
		rOutter = 1.0/2;
		rCannon = 0.1;
		hCannon = 1.5;
		dist = Point3F(0.2, 1.0/12, 2.0/5);
		map<GLuint, int>& colorMap = ColorMap::colorMap;
		//if(colorMap.empty()){
			colorMap[111] = 1;
			colorMap[112] = 2;
			colorMap[113] = 3;

			colorMap[121] = 4;
			colorMap[122] = 5;
			colorMap[123] = 6;


			colorMap[13] = 7;


			colorMap[141] = 8;
			colorMap[142] = 9;
			colorMap[143] = 10;
			colorMap[144] = 11;
		//}
	}
	void draw(DrawMode mode = RENDER){
		switch(mode){
			case RENDER : draw0(0);break;
			case SELECTION : draw0(1);break;
			default:draw0(2);break;
		}
	}
private:
	void draw0(int mode){
		double frontWheelAxisZ = length / 2 - dist.z(), frontWheelAxisX = width / 2 + dist.x(),
			frontWheelAxisY = dist.y();
		double frontRightWheelX = frontWheelAxisX - dist.x() / 2;
		Drawer drawer;
		GLfloat c1[] = {1.0, 1.0, 1.0, 1.0};
		GLfloat c2[]= {0.0, 1.0, 0.0, 1.0};
		GLfloat c3[] = {0.0, 1.0, 0.0, 1.0};
		GLfloat mat_shininess[]={ 30.0 };
		glHint(GL_LINE_SMOOTH, GL_NICEST);
		glHint(GL_POLYGON_SMOOTH, GL_NICEST);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c1);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, c2);
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, c3);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess );
		//front wheel;
		CylinderF frontWheelAxis(Point3F(-frontWheelAxisX, frontWheelAxisY, frontWheelAxisZ),
			Point3F(frontWheelAxisX, frontWheelAxisY, frontWheelAxisZ), rInner, rInner);
		CylinderF frontRightWheel(Point3F(frontRightWheelX, frontWheelAxisY, frontWheelAxisZ),
			Point3F(frontWheelAxisX, frontWheelAxisY, frontWheelAxisZ), rOutter, rOutter);
		CylinderF frontLeftWheel(Point3F(-frontRightWheelX, frontWheelAxisY, frontWheelAxisZ),
			Point3F(-frontWheelAxisX, frontWheelAxisY, frontWheelAxisZ), rOutter, rOutter);
		//back wheel;
		CylinderF backWheelAxis(Point3F(-frontWheelAxisX, frontWheelAxisY,-frontWheelAxisZ),
			Point3F(frontWheelAxisX, frontWheelAxisY, -frontWheelAxisZ), rInner, rInner);
		CylinderF backRightWheel(Point3F(frontRightWheelX, frontWheelAxisY, -frontWheelAxisZ),
			Point3F(frontWheelAxisX, frontWheelAxisY, -frontWheelAxisZ), rOutter, rOutter);
		CylinderF backLeftWheel(Point3F(-frontRightWheelX, frontWheelAxisY, -frontWheelAxisZ),
			Point3F(-frontWheelAxisX, frontWheelAxisY, -frontWheelAxisZ), rOutter, rOutter);
		//cannon;
		CylinderF cannon(Point3F(0, height, 0),
			Point3F(0, height + hCannon, 0), rCannon, rCannon);
		double sqrt2 = sqrt(2.0);
		CylinderF cannonTube(Point3F(0, height + hCannon, 0),
			Point3F(0, height + hCannon + hCannon / 3 / sqrt2, hCannon / 3 / sqrt2), rCannon, rCannon);
		ConeF cannonArrow(Point3F(0, height + hCannon + hCannon / 3 / sqrt2, hCannon / 3 / sqrt2),
			Point3F(0, height + hCannon + hCannon / 2 / sqrt2, hCannon / 2 / sqrt2), rCannon * 1.1);
		if(mode == 0){
#define DRAW_PART(obj, val) if(lastSelPart == val){\
			GLfloat c[] = {0.5, 0.5, 0.0, 1.0};\
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, c);\
			drawer.draw(obj);\
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, c3);\
		}else{ drawer.draw(obj);}
			//front wheel;
			DRAW_PART(frontWheelAxis, 111);
			DRAW_PART(frontRightWheel, 112);
			DRAW_PART(frontLeftWheel, 113);

			//back wheel;
			DRAW_PART(backWheelAxis, 121);
			DRAW_PART(backRightWheel, 122);
			DRAW_PART(backLeftWheel, 123);

			//body
			if(lastSelPart == 13){	
				GLfloat c[] = {0.5, 0.5, 0.0, 1.0};
				glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, c);
				drawer.draw(Point3F(-width /  2, 0, -length / 2), Point3F(width / 2, height, length / 2));
				glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, c3);
			}else{
				drawer.draw(Point3F(-width /  2, 0, -length / 2), Point3F(width / 2, height, length / 2));
			}

			//cannon;
			DRAW_PART(cannon, 141);
			DRAW_PART(cannonTube, 142);
			DRAW_PART(cannonArrow, 143);			
			if(lastSelPart == 144){	
				GLfloat c[] = {0.5, 0.5, 0.0, 1.0};
				glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, c);
			}
			glPushMatrix();
				glTranslated(0, height + hCannon, 0);
				glutSolidSphere(rCannon, drawer.getSlice(),drawer.getSlice());
			glPopMatrix();
			if(lastSelPart == 144)
				glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, c3);
		}else if(mode == 1){
			glPushName(1);
				glPushName(11);
					//front wheel;
					glPushName(111);
						drawer.draw(frontWheelAxis);
					glLoadName(112);
						drawer.draw(frontRightWheel);
					glLoadName(113);
						drawer.draw(frontLeftWheel);
					glPopName();
				glLoadName(12);
					//back wheel;
					glPushName(121);
						drawer.draw(backWheelAxis);
					glLoadName(122);
						drawer.draw(backRightWheel);
					glLoadName(123);
						drawer.draw(backLeftWheel);
					glPopName();
				glLoadName(13);
					//body
					drawer.draw(Point3F(-width /  2, 0, -length / 2), Point3F(width / 2, height, length / 2));
				glLoadName(14);
					//cannon;
					glPushName(141);
						drawer.draw(cannon);
					glLoadName(142);
						drawer.draw(cannonTube);
					glLoadName(143);
						drawer.draw(cannonArrow);
					glPushName(144);
						glPushMatrix();
							glTranslated(0, height + hCannon, 0);
							glutSolidSphere(rCannon, drawer.getSlice(),drawer.getSlice());
						glPopMatrix();
					glPopName();
				glPopName();
			glPopName();
		}else{
			//front wheel;
			GLint oldShadeModel;
			GLboolean lightingEnabled;
			int color = 0;
			map<GLuint, int>& colorMap = ColorMap::colorMap;
			glGetIntegerv(GL_SHADE_MODEL, &oldShadeModel);
			lightingEnabled = glIsEnabled(GL_LIGHTING);
			glDisable(GL_LIGHTING);
			glShadeModel(GL_FLAT);
			color = colorMap[111]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF); drawer.draw(frontWheelAxis);
			color = colorMap[112]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF); drawer.draw(frontRightWheel);
			color = colorMap[113]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF); drawer.draw(frontLeftWheel);

			//back wheel;
			color = colorMap[121]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF); drawer.draw(backWheelAxis);
			color = colorMap[122]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF); drawer.draw(backRightWheel);
			color = colorMap[123]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF); drawer.draw(backLeftWheel);
		
			//body
			color = colorMap[13]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF); drawer.draw(Point3F(-width /  2, 0, -length / 2), Point3F(width / 2, height, length / 2));

			//cannon;
			color = colorMap[141]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF); drawer.draw(cannon);
			color = colorMap[142]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF); drawer.draw(cannonTube);
			color = colorMap[143]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF); drawer.draw(cannonArrow);
			color = colorMap[144]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF); 
			glPushMatrix();
			glTranslated(0, height + hCannon, 0);
			glutSolidSphere(rCannon, drawer.getSlice(),drawer.getSlice());
			glPopMatrix();
			glPopName();
			if(oldShadeModel != GL_FLAT)
				glShadeModel(oldShadeModel);
			if(lightingEnabled)
				glEnable(GL_LIGHTING);
		}

	}
};
static const int TEX_WIDTH = 255;
static const int TEX_HEIGHT = 255;
GLubyte texImage[TEX_WIDTH][TEX_HEIGHT][3];
class Chair{
public:
	Chair(){
		map<GLuint, int>& colorMap = ColorMap::colorMap;
		//if(colorMap.empty()){
			colorMap[111] = 1;
			colorMap[112] = 2;
			colorMap[113] = 3;
			colorMap[114] = 4;
			colorMap[115] = 5;
			colorMap[121] = 6;
			colorMap[122] = 7;
			colorMap[123] = 8;
			colorMap[124] = 9;
			colorMap[125] = 10;
			colorMap[126] = 11;
			colorMap[127] = 12;
			colorMap[128] = 13;
		//}
	}
	void draw(DrawMode mode = RENDER){
		glPushMatrix();
		glRotatef(90, 1, 0, 0);
		switch(mode){
			case RENDER : draw0(0);break;
			case SELECTION : draw0(1);break;
			default:draw0(2);break;
		}
		glPopMatrix();
	}
private:
	inline void draw1(Drawer& drawer, const Point3F& p1, const Point3F& p2, GLuint name, GLfloat *oldColor){
		if(lastSelPart == name){
			GLfloat c[] = {0.5, 0.5, 0.0, 1.0};
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, c);
			drawer.draw(p1, p2, true);
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, oldColor);
		}else{ 
			drawer.draw(p1, p2, true);
		}
	}
	void draw0(int i){
		Drawer drawer;
		float halfWidth = 1.5, halfLength = 2, thickness = 0.2,
			height = 4;
		GLfloat c1[] = {1.0, 1.0, 1.0, 1.0};
		GLfloat c2[]= {1.0, 1.0, 1.0, 1.0};
		GLfloat c3[] = {1.0, 1.0, 1.0, 1.0};
		GLfloat mat_shininess[]={ 30.0 };
		Point3F p1(-halfWidth, 0, -halfLength);
		Point3F p2(halfWidth, thickness, halfLength);
		Point3F p3(halfWidth - thickness, -height, halfLength - thickness);
		Point3F p4(halfWidth - thickness, -thickness, -halfLength);
		Point3F p5(-halfWidth, -thickness, halfLength - thickness);
		Point3F p6(halfWidth, 0, halfLength);
		glHint(GL_LINE_SMOOTH, GL_NICEST);
		glHint(GL_POLYGON_SMOOTH, GL_NICEST);
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, c1);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, c2);
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, c3);
		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess );
		if(i == 0){
			GLuint texName;
			glGenTextures(1, &texName);
			glBindTexture(GL_TEXTURE_2D, texName);
			glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP);
			glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP);
			glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
			glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
			glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
			glTexImage2D(GL_TEXTURE_2D,0,GL_RGB8,TEX_WIDTH,TEX_HEIGHT,
					     0,GL_RGB,GL_UNSIGNED_BYTE,texImage);			
			glBindTexture(GL_TEXTURE_2D, texName);
			//111
			draw1(drawer, p1, p2, 111, c3);
			//112
			draw1(drawer, p4, p6, 112, c3);
			//113
			draw1(drawer, p5, p6, 113, c3);
			//114
			glPushMatrix();
			glTranslatef(thickness - 2 * halfWidth, 0, 0);
			draw1(drawer, p4, p6, 114, c3);
			glPopMatrix();
			//115
			glPushMatrix();
			glTranslatef(0, 0, thickness - 2 * halfLength);
			draw1(drawer, p5, p6, 115, c3);
			glPopMatrix();
			//121
			draw1(drawer, p3, p2, 121, c3);
			//122
			glPushMatrix();
			glTranslatef(thickness - 2 * halfWidth, 0, 0);
			draw1(drawer, p3, p2, 122, c3);
			glPopMatrix();
			//123
			glPushMatrix();
			glTranslatef(0, 0, thickness - 2 * halfLength);
			draw1(drawer, p3, p2, 123, c3);
			glPopMatrix();
			//124
			glPushMatrix();
			glTranslatef(thickness - 2 * halfWidth, 0, thickness - 2 * halfLength);
			draw1(drawer, p3, p2, 124, c3);
			glPopMatrix();
			//125
			glPushMatrix();
			glTranslatef(0, -height / 2, 0);
			draw1(drawer, p4, p6, 125, c3);
			glPopMatrix();
			//126
			glPushMatrix();
			glTranslatef(thickness - 2 * halfWidth, -height / 2, 0);
			draw1(drawer, p4, p6, 126, c3);
			glPopMatrix();
			//127
			glPushMatrix();
			glTranslatef(0, -height / 2, 0);
			draw1(drawer, p5, p6, 127, c3);
			glPopMatrix();
			//128
			glPushMatrix();
			glTranslatef(0, -height / 2, thickness - 2 * halfLength);
			draw1(drawer, p5, p6, 128, c3);
			glPopMatrix();			
			glDeleteTextures(1, &texName);
		}else if(i == 1){
			glPushName(1);
				glPushName(11);
					glPushName(111);
						//111
						drawer.draw(p1, p2);
					glLoadName(112);
						//112
						drawer.draw(p4, p6);
					glLoadName(113);
						//113
						drawer.draw(p5, p6);
					glLoadName(114);
						//114
						glPushMatrix();
						glTranslatef(thickness - 2 * halfWidth, 0, 0);
						drawer.draw(p4, p6);
						glPopMatrix();
					glLoadName(115);
						//115
						glPushMatrix();
						glTranslatef(0, 0, thickness - 2 * halfLength);
						drawer.draw(p5, p6);
						glPopMatrix();
					glPopName();
				glPushName(12);
					glPushName(121);
						//121
						drawer.draw(p3, p2);
					glLoadName(122);
						//122
						glPushMatrix();
						glTranslatef(thickness - 2 * halfWidth, 0, 0);
						drawer.draw(p3, p2);
						glPopMatrix();
					glLoadName(123);
						//123
						glPushMatrix();
						glTranslatef(0, 0, thickness - 2 * halfLength);
						drawer.draw(p3, p2);
						glPopMatrix();
					glLoadName(124);
						//124
						glPushMatrix();
						glTranslatef(thickness - 2 * halfWidth, 0, thickness - 2 * halfLength);
						drawer.draw(p3, p2);
						glPopMatrix();
					glLoadName(125);
						//125
						glPushMatrix();
						glTranslatef(0, -height / 2, 0);
						drawer.draw(p4, p6);
						glPopMatrix();
					glLoadName(126);
						//126
						glPushMatrix();
						glTranslatef(thickness - 2 * halfWidth, -height / 2, 0);
						drawer.draw(p4, p6);
						glPopMatrix();
					glLoadName(127);
						//127
						glPushMatrix();
						glTranslatef(0, -height / 2, 0);
						drawer.draw(p5, p6);
						glPopMatrix();
					glLoadName(128);
						//128
						glPushMatrix();
						glTranslatef(0, -height / 2, thickness - 2 * halfLength);
						drawer.draw(p5, p6);
						glPopMatrix();
					glPopName();
				glPopName();
			glPopName();
		}else{
			GLint oldShadeModel;
			GLboolean lightingEnabled;
			int color = 0;
			map<GLuint, int>& colorMap = ColorMap::colorMap;
			glGetIntegerv(GL_SHADE_MODEL, &oldShadeModel);
			lightingEnabled = glIsEnabled(GL_LIGHTING);
			glDisable(GL_LIGHTING);
			glShadeModel(GL_FLAT);
			color = colorMap[111]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF);
			//111			
			drawer.draw(p1, p2);
			color = colorMap[112]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF);
			//112		
			drawer.draw(p4, p6);
			color = colorMap[113]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF);
			//113	
			drawer.draw(p5, p6);
			color = colorMap[114]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF);
			//114
			glPushMatrix();
			glTranslatef(thickness - 2 * halfWidth, 0, 0);	
			drawer.draw(p4, p6);
			glPopMatrix();
			color = colorMap[115]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF);
			//115
			glPushMatrix();
			glTranslatef(0, 0, thickness - 2 * halfLength);
			drawer.draw(p5, p6);
			glPopMatrix();
			color = colorMap[121]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF);
			//121
			drawer.draw(p3, p2);
			color = colorMap[122]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF);
			//122
			glPushMatrix();
			glTranslatef(thickness - 2 * halfWidth, 0, 0);
			drawer.draw(p3, p2);
			glPopMatrix();
			color = colorMap[123]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF);
			//123
			glPushMatrix();
			glTranslatef(0, 0, thickness - 2 * halfLength);
			drawer.draw(p3, p2);
			glPopMatrix();
			color = colorMap[124]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF);
			//124
			glPushMatrix();
			glTranslatef(thickness - 2 * halfWidth, 0, thickness - 2 * halfLength);
			drawer.draw(p3, p2);
			glPopMatrix();
			color = colorMap[125]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF);
			//125
			glPushMatrix();
			glTranslatef(0, -height / 2, 0);	
			drawer.draw(p4, p6);
			glPopMatrix();
			color = colorMap[126]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF);
			//126
			glPushMatrix();
			glTranslatef(thickness - 2 * halfWidth, -height / 2, 0);	
			drawer.draw(p4, p6);
			glPopMatrix();
			color = colorMap[127]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF);
			//127
			glPushMatrix();
			glTranslatef(0, -height / 2, 0);
			drawer.draw(p5, p6);
			glPopMatrix();
			color = colorMap[128]; glColor3ub((color >> 16) & 0xFF, color >> 8 & 0xFF, color & 0xFF);
			//128
			glPushMatrix();
			glTranslatef(0, -height / 2, thickness - 2 * halfLength);
			drawer.draw(p5, p6);
			glPopMatrix();
			if(oldShadeModel != GL_FLAT)
				glShadeModel(oldShadeModel);
			if(lightingEnabled)
				glEnable(GL_LIGHTING);
		}
	}
};
map<GLuint,int> ColorMap::colorMap;
static CameraF camera(Point3F(3, 3, 3), Point3F(0, 0, 0), Point3F(0, 1, 0));
static float rotateAngle = 0, curX = 0, curZ = 0;
static int startX = -1, startY = -1;

void draw(DrawMode mode){
	const float *eyePos = camera.getEyePos();
	const float *target = camera.getTarget();
	const float *upVector = camera.getUpVector();
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(eyePos[0], eyePos[1], eyePos[2], target[0], target[1], target[2], upVector[0], upVector[1], upVector[2]);
	glTranslatef(curX, 0, curZ);	
	glRotatef(rotateAngle, 0, 1, 0);
	//Car car;
	//car.draw(mode);
	Chair chair;
	chair.draw(mode);
}
void drawDetail()
{
	if(startX >= 0 && startY >= 0){
		GLuint winWidth, winHeight;
		GLint viewPort[4];
		GLint oldScissor[4];
		GLfloat oldMatrix[16];
		glGetIntegerv(GL_VIEWPORT, viewPort);
		glGetIntegerv(GL_SCISSOR_BOX, oldScissor);
		glEnable(GL_SCISSOR_TEST);
		glViewport(0, 0, viewPort[2] / 5, viewPort[3] / 5);
		glScissor(0, 0, viewPort[2] / 5, viewPort[3] / 5);
		glClearColor(0.0, 0.0, 0.0, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		winWidth = glutGet(GLUT_WINDOW_WIDTH);
		winHeight = glutGet(GLUT_WINDOW_HEIGHT);
		glGetFloatv (GL_PROJECTION_MATRIX, oldMatrix);   
		glMatrixMode(GL_PROJECTION); 
		glLoadIdentity();
		gluPickMatrix(startX, winHeight - startY, 40, 40, viewPort);
		glMultMatrixf(oldMatrix);
		draw(RENDER);
		glMatrixMode(GL_PROJECTION); 
		glLoadIdentity();
		glMultMatrixf(oldMatrix);
		glViewport(viewPort[0], viewPort[1], viewPort[2] ,viewPort[3]);//, 0, 30, 30);
		glDisable(GL_SCISSOR_TEST);
		glClearColor(0.5, 0.5, 0.5, 1.0);
		glScissor(oldScissor[0], oldScissor[1] , oldScissor[2], oldScissor[3]);
	}
}
struct FrameRate{
	clock_t start;
	clock_t last;
	long long count;
	FrameRate() : start(clock()), last(clock()), count(0){}
}instance;
void display(){
	static class Init{
	public:
		Init(){instance.start = clock();}
	}inst;
	instance.last = clock();
	++instance.count;
	float eyePos[]={camera.getEyePos().x(), camera.getEyePos().y(), camera.getEyePos().z(), 0};
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLightfv(GL_LIGHT0, GL_POSITION, eyePos);
	draw(RENDER);
	drawDetail();
	glutSwapBuffers();
}
struct Light{
	float position[4];
	float ambient[4];
	float specular[4];
	float diffuse[4];
};
extern void hsl2rgb(double h, double sl, double l, double& r, double& g, double& b);
extern void rgb2hsl (int rVal, int gVal, int bVal, double& h, double& s, double& l);
void init(){	
	Light lit0={
		{-3, 0, 0, 1},
		{1, 1.0, 1.0, 1},
		{1.0, 1.0, 1.0, 1},
		{1.0, 1.0, 1.0, 1}
	};
	glClearColor(0.5, 0.5, 0.5,0.0);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_NORMALIZE);
	glEnable(GL_SMOOTH);
	glEnable(GL_BLEND);	
	glDisable(GL_COLOR_MATERIAL);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glLightfv(GL_LIGHT0, GL_AMBIENT, lit0.ambient);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lit0.specular);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lit0.diffuse);
	glLightfv(GL_LIGHT0, GL_POSITION, lit0.position);
	glEnable(GL_LIGHT0);
        glEnable(GL_TEXTURE_2D); // allow 2D texture maps
	const unsigned int brown = 0xA52A2A, tan = 0xD2B48C;
	double hb, s, l, ht;
	rgb2hsl((brown >> 16) & 0xFF, (brown >> 8) & 0xFF, brown & 0xFF, hb, s, l);
	rgb2hsl((tan >> 16) & 0xFF, (tan >> 8) & 0xFF, tan & 0xFF, ht, s, l);
	double interval = (ht - hb) / TEX_WIDTH;
	for(int i = 0;i < TEX_WIDTH;++ i){
		const double h = hb + i * interval;
		double r, g ,b;
		hsl2rgb(h, s, l, r, g, b);
		for(int j = 0;j < TEX_HEIGHT;++ j){
			texImage[i][j][0] = (r * 255 + 0.5);
			texImage[i][j][1] =  (g * 255 + 0.5);
			texImage[i][j][2] =  (b * 255 + 0.5);
		}
	}
}
static bool isMoving = false;
static  bool useBackBuffer = true;
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
	case 'r' :isMoving = !isMoving;break;
	case 'o' :cout << "frame rate:" << instance.count * CLOCKS_PER_SEC / 
		(double)(instance.last - instance.start) << " FPS" << endl;break;
	case 'b':useBackBuffer = !useBackBuffer;
	}
	assert(camera.isValid());
	glutPostRedisplay();
}
static const int timerInterval = 50;
static const int timerID = 88;
static const float moveStep = 0.01;
static void animate(){
	if(rand() % 2 == 0){
		switch(rand() % 10){
			case 0:curX += moveStep;break;
			case 1:curX -= moveStep;break;
			case 2:curZ += moveStep;break;
			default:curZ -= moveStep;break;
		}
	}else{
		switch(rand() % 4){
			case 0:rotateAngle += 0.1;break;
			default:rotateAngle -= 0.1;break;
		}
		if(rotateAngle > 360)
			rotateAngle -= 360;
		if(rotateAngle < -360)
			rotateAngle += 360;
	}
}
static void timer(int i){
	glutTimerFunc(timerInterval, timer, timerID);
	animate();
	glutPostRedisplay();
}
void idle(){
	animate();
	glutPostRedisplay();
}
void menu(int val){
}
void submenu(int val){
}
//static void createMenu(void(*func)(int), char **items)
void createMenu(){
	static int mainMenu = glutCreateMenu(menu);
	static int subMenu = glutCreateMenu(submenu);
	int oldMenu = glutGetMenu();
	glutSetMenu(mainMenu);
	glutAddMenuEntry("F S", 1);
	glutAttachMenu(GLUT_RIGHT_BUTTON);
	glutAddSubMenu("sub", subMenu);
	{
		int oldMenu = glutGetMenu();
		glutSetMenu(subMenu);
		glutAddMenuEntry("T S", 1);
		glutSetMenu(oldMenu);
	}
	glutSetMenu(oldMenu);
}
static int sideLen = 0;
static void reshape(int w, int h)
{
	int t = min(w, h);
	sideLen = t;
	glViewport(0, 0, t, t);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(30.0,1.0,1.0,100.0);
}
void mouse(int button, int state, int x, int y)
{
	if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		int winWidth, winHeight;
		GLint viewPort[4];
		GLfloat oldMatrix[16];
		glGetIntegerv(GL_VIEWPORT, viewPort);

		winWidth = glutGet(GLUT_WINDOW_WIDTH);
		winHeight = glutGet(GLUT_WINDOW_HEIGHT);
		glGetFloatv (GL_PROJECTION_MATRIX, oldMatrix);   
		//glMatrixMode(GL_PROJECTION); 
			//glLoadIdentity();
			//gluPickMatrix(x, winHeight - y, 1, 1, viewPort);
			//glMultMatrixf(oldMatrix);
			if(useBackBuffer){
				GLubyte backBuffer[3];
				set<GLuint> results;
				glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
				glDrawBuffer(GL_BACK);
				draw(BACK_BUFFER);
				glReadBuffer(GL_BACK);
				glReadPixels(x, winHeight - y, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, backBuffer);
				/*for(int i = 0;i < 5;++ i){
					for(int j = 0;j < 5;++ j){*/
						int index = 0;//3 * (i * sideLen + j);
						int value = backBuffer[index];value <<= 8;
						value += backBuffer[index + 1];value <<= 8;
						value += backBuffer[index + 2];
						//cout << hex << (int)backBuffer[index] << '-' << 
						//	(int)backBuffer[index + 1] << '-' << (int)backBuffer[index + 2] << ' ';
						results.insert(ColorMap::getByColor(value));
				//	}
				//	cout << endl;
				//}
				//delete [] backBuffer;*/
				if(!results.empty()){
					lastSelPart = *results.begin();
				}
			}else{
				SelectionBuffer selectBuffer(100);
				glSelectBuffer(selectBuffer.getBufSize(), selectBuffer);
				glRenderMode(GL_SELECT);
				glInitNames();
				draw(SELECTION);
				selectBuffer.setRecordCount(glRenderMode(GL_RENDER));
				cout << "-----------------------------" << endl;
				cout << selectBuffer << endl;
				lastSelPart = selectBuffer.getNearest();
			}
		glMatrixMode(GL_PROJECTION); 
			glLoadIdentity();
			glMultMatrixf(oldMatrix);
		glutPostRedisplay();
	}
}
void passiveMotion(int x, int y){
	startX = x;
	startY = y;
	glutPostRedisplay();
}

int main(int argc, char *argv[]){
	srand(time(0));
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(800, 800);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Draw Function");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	//glutIdleFunc(idle);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(mouse);
	glutPassiveMotionFunc(passiveMotion);
	//glutTimerFunc(timerInterval, timer, timerID);
	createMenu();
	init();
	glutMainLoop();
	return 0;
}

