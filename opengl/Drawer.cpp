#include <GL/glut.h>
#include "Drawer.h"
#include <cassert>
#include <iostream>
using namespace std;
using namespace tmlg;
void Drawer::draw(const CylinderGeometryF& g){
	assert(g.c != NULL);
	float *ptBodys = new float[2 * (g.slice + 1) * 3];
	float *ptNorms = new float[2 * (g.slice + 1) * 3];
	float *ptBotCap = NULL, *ptTopCap = NULL;
	if(g.drawBotCap)
		ptBotCap = new float[(g.slice + 2) * 3];
	if(g.drawTopCap)
		ptTopCap = new float[(g.slice + 2) * 3];
	if(g.c->calcData(ptBodys, ptNorms, ptBotCap, ptTopCap, g.slice)){
		glBegin(GL_QUAD_STRIP);
		for(int i = 0;i < 2 * (g.slice + 1);++ i){
			glNormal3fv(ptNorms + 3 * i);
			glVertex3fv(ptBodys + 3 * i);
		}
		glEnd();
		if(ptBotCap != NULL){
			Point3F axis = g.c->getAxis();
			axis.normalize();
			glBegin(GL_TRIANGLE_FAN);
			for(int i = 0;i < g.slice + 2;++ i){
				glNormal3fv(axis);
				glVertex3fv(ptBotCap + 3 * i);
			}
			glEnd();
		}
		if(ptTopCap != NULL){
			Point3F axis = g.c->getAxis();
			axis.normalize();
			axis *= -1;
			glBegin(GL_TRIANGLE_FAN);
			for(int i = 0;i < g.slice + 2;++ i){
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
void Drawer::draw(const ConeGeometryF& g){
	assert(g.c != NULL);
	float *ptBodys = new float[(g.slice + 2) * 3];
	float *ptNorms = new float[(g.slice + 2) * 3];
	const ConeF& c = *g.c;
	float *ptCaps = NULL;
	if(g.drawCap)
		ptCaps = new float[(g.slice + 2) * 3];
	if(c.calcData(ptBodys, ptNorms, ptCaps, g.slice)){
		glBegin(GL_TRIANGLE_FAN);
		for(int i = 0;i < g.slice + 2;++ i){
			glNormal3fv(ptNorms + 3 * i);
			glVertex3fv(ptBodys + 3 * i);
		}
		glEnd();
		if(ptCaps != NULL){
			Point3F axis = c.getAxis();
			axis.normalize();
			axis *= -1;
			glBegin(GL_TRIANGLE_FAN);
			for(int i = 0;i < g.slice + 2;++ i){
				glNormal3fv(axis);
				glVertex3fv(ptCaps + 3 * i);
			}
			glEnd();
		}
	}
}
void Drawer::draw(const BoxGeometryF& g){
	const Point3F& minPt = g.minPt;
	const Point3F& maxPt = g.maxPt;
	glBegin(GL_QUADS);
		glNormal3f(-1, 0, 0);
		if(g.hasTexture) glTexCoord2f(0.0, 0.0); glVertex3fv(minPt);
		if(g.hasTexture) glTexCoord2f(0.0, 1.0); glVertex3f(minPt.x(), minPt.y(), maxPt.z());
		if(g.hasTexture) glTexCoord2f(1.0, 1.0); glVertex3f(minPt.x(), maxPt.y(), maxPt.z());
		if(g.hasTexture) glTexCoord2f(1.0, 0.0); glVertex3f(minPt.x(), maxPt.y(), minPt.z());
		if(g.hasTexture) glNormal3f(0, 1, 0);
		if(g.hasTexture) glTexCoord2f(0.0, 0.0); glVertex3f(minPt.x(), maxPt.y(), minPt.z());
		if(g.hasTexture) glTexCoord2f(0.0, 1.0); glVertex3f(minPt.x(), maxPt.y(), maxPt.z());
		if(g.hasTexture) glTexCoord2f(1.0, 1.0); glVertex3fv(maxPt);
		if(g.hasTexture) glTexCoord2f(1.0, 0.0); glVertex3f(maxPt.x(), maxPt.y(), minPt.z());
		glNormal3f(1, 0, 0);
		if(g.hasTexture) glTexCoord2f(0.0, 0.0); glVertex3f(maxPt.x(), maxPt.y(), minPt.z());
		if(g.hasTexture) glTexCoord2f(0.0, 1.0); glVertex3fv(maxPt);
		if(g.hasTexture) glTexCoord2f(1.0, 1.0); glVertex3f(maxPt.x(), minPt.y(), maxPt.z());
		if(g.hasTexture) glTexCoord2f(1.0, 0.0); glVertex3f(maxPt.x(), minPt.y(), minPt.z());
		glNormal3f(0, -1, 0);
		if(g.hasTexture) glTexCoord2f(0.0, 0.0); glVertex3f(maxPt.x(), minPt.y(), minPt.z());
		if(g.hasTexture) glTexCoord2f(0.0, 1.0); glVertex3f(maxPt.x(), minPt.y(), maxPt.z());
		if(g.hasTexture) glTexCoord2f(1.0, 1.0); glVertex3f(minPt.x(), minPt.y(), maxPt.z());
		if(g.hasTexture) glTexCoord2f(1.0, 0.0); glVertex3fv(minPt);
		glNormal3f(0, 0, 1);
		if(g.hasTexture) glTexCoord2f(0.0, 0.0); glVertex3fv(maxPt);
		if(g.hasTexture) glTexCoord2f(0.0, 1.0); glVertex3f(minPt.x(), maxPt.y(), maxPt.z());
		if(g.hasTexture) glTexCoord2f(1.0, 1.0); glVertex3f(minPt.x(), minPt.y(), maxPt.z());
		if(g.hasTexture) glTexCoord2f(1.0, 0.0); glVertex3f(maxPt.x(), minPt.y(), maxPt.z());
		glNormal3f(0, 0, -1);
		if(g.hasTexture) glTexCoord2f(0.0, 0.0); glVertex3fv(minPt);
		if(g.hasTexture) glTexCoord2f(0.0, 1.0); glVertex3f(minPt.x(), maxPt.y(), minPt.z());
		if(g.hasTexture) glTexCoord2f(1.0, 1.0); glVertex3f(maxPt.x(), maxPt.y(), minPt.z());
		if(g.hasTexture) glTexCoord2f(1.0, 0.0); glVertex3f(maxPt.x(), minPt.y(), minPt.z());
	glEnd();
}
static void drawInternal(Point3F* points,  Point3F* norms, Point3F* texCoord,
	int rows, int cols, bool calculateNorm){
	glBegin(GL_TRIANGLES);
	for(int i = 0;i < rows;++ i){
		int index = i * (cols + 1);
		for(int j = 0;j < cols;++ j){
			int i1 = index + j;
			int i2 = i1 + cols + 1;
			if(calculateNorm){
				Point3F dir1 = points[i1 + 1];
				Point3F dir2 = points[i2];
				Point3F normal;
				dir1 -= points[i1];
				dir2 -= points[i1];
				normal = dir1.crossProduct(dir2);
				normal.normalize();
				glNormal3fv(normal);		
			}
			if(norms) 
				glNormal3fv(norms[i1]);			
			if(texCoord)
				glTexCoord2fv(texCoord[i1]);
			glVertex3fv(points[i1]);
			if(norms) 
				glNormal3fv(norms[i1 + 1]);
			if(texCoord)
				glTexCoord2fv(texCoord[i1 + 1]);
			glVertex3fv(points[i1 + 1]);
			if(norms) 
				glNormal3fv(norms[i2]);
			if(texCoord)
				glTexCoord2fv(texCoord[i2]);
			glVertex3fv(points[i2]);


			if(calculateNorm){
				Point3F dir1 = points[i2];
				Point3F dir2 = points[i1 + 1];
				Point3F normal;
				dir1 -= points[i2 + 1];
				dir2 -= points[i2 + 1];
				normal = dir1.crossProduct(dir2);
				normal.normalize();
				glNormal3fv(normal);	
			}
			if(norms) 
				glNormal3fv(norms[i1 + 1]);
			if(texCoord)
				glTexCoord2fv(texCoord[i1 + 1]);
			glVertex3fv(points[i1 + 1]);
			if(norms) 
				glNormal3fv(norms[i2 + 1]);
			if(texCoord)
				glTexCoord2fv(texCoord[i2 + 1]);
			glVertex3fv(points[i2 + 1]);			
			if(norms) 
				glNormal3fv(norms[i2]);
			if(texCoord)
				glTexCoord2fv(texCoord[i2]);
			glVertex3fv(points[i2]);
		}		
	}
	glEnd();
}
void Drawer::draw(const ParametricSurfaceGeometryF& g){
	EvaluatorF e;
	int count = (g.rows + 1) * (g.cols + 1);
	Point3F* points = new Point3F[count];
	Point3F* norms = NULL;
	Point3F* texCoord = NULL;
	if(g.normFunc != NULL && !g.calculateNorm)
		norms = new Point3F[count];
	if(g.hasTexture)
		texCoord = new Point3F[count];
	e.eval(g.startU, g.endU, g.startV, g.endV, g.rows, g.cols,
		points, g.evalFunc,
		norms, g.normFunc, texCoord);
	drawInternal(points, norms, texCoord, g.rows, g.cols, g.calculateNorm);
	delete [] norms;
	delete [] points;
	delete [] texCoord;
}
void Drawer::draw(const SurfaceGeometryF& g){
	EvaluatorF e;
	int count = (g.rows + 1) * (g.cols + 1);
	Point3F* points = new Point3F[count];
	Point3F* norms = NULL;
	Point3F* texCoord = NULL;
	if(g.normFunc != NULL && !g.calculateNorm)
		norms = new Point3F[count];
	if(g.hasTexture)
		texCoord = new Point3F[count];
	e.eval(g.startX, g.endX, g.startY, g.endY, g.rows, g.cols, points, g.evalFunc, norms, g.normFunc, texCoord);
	drawInternal(points, norms, texCoord, g.rows, g.cols, g.calculateNorm);
	delete [] norms;
	delete [] points;
	delete [] texCoord;
}

