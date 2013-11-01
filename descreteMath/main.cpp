#include "glut.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <cassert>
template<class T>
struct Constant{
	typedef T type;
	const static T pi;
	const static T tolerance;
};
typedef Constant<double> ConstantD;
typedef Constant<float> ConstantF;
const ConstantD::type ConstantD::pi = asin(1.0)*2;
const ConstantD::type ConstantD::tolerance = 1.e-10;
const ConstantF::type ConstantF::pi =  ConstantD::pi;
const ConstantF::type  ConstantF::tolerance = 1.e-6;

template<class T>
struct Point3{
	T coord[3];
	Point3(){
		coord[0] = coord[1] = coord[2] = 0;
	}
	Point3(T x, T y, T z){
		this->coord[0] = x;
		this->coord[1] = y;
		this->coord[2] = z;
	}
	Point3(const T coord[3]){
		this->coord[0] = coord[0];
		this->coord[1] = coord[1];
		this->coord[2] = coord[2];
	}
	template<class T>
	Point3& operator=(const T coord[3]){
		this->coord[0] = coord[0];
		this->coord[1] = coord[1];
		this->coord[2] = coord[2];
		return *this;
	}
	T x() const{
		return coord[0];
	}
	T y() const{
		return coord[1];
	}
	T z() const{
		return coord[2];
	}
	operator T*(){
		return coord;
	}
	operator const T*() const{
		return coord;
	}
	T dist(const Point3& pt) const{
		Point3 tmp(*this);
		return tmp.subtract(pt).norm();
	}
	T norm() const{
		return sqrt(coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2]);
	}
	T normSquare() const{
		return coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2];
	}
	void normalize(){
		double norm = this->norm();
		coord[0] /= norm;
		coord[1] /= norm;
		coord[2] /= norm;
	}
	Point3& operator*=(T val){
		return multiply(val);
	}
	Point3& multiply(T val){
		coord[0] *= val;
		coord[1] *= val;
		coord[2] *= val;
		return *this;
	}
	Point3& add(const Point3& pt){
		coord[0] += pt.coord[0];
		coord[1] += pt.coord[1];
		coord[2] += pt.coord[2];
		return *this;
	}
	Point3 proj(const Point3& vec){
		T normS = normSquare();
		if(normS <= Constant<T>::tolerance)
			return vec;
		Point3 vecProj = vec;
		vecProj *=  vec.dotProduct(vec) / normS;
		return vecProj;
	}
	Point3 orth(const Point3& vec){
		Point3 vecOrth = vec;
		vecOrth.subtract(proj(vec));
		return vecOrth;
	}
	Point3& subtract(const Point3& pt){
		coord[0] -= pt.coord[0];
		coord[1] -= pt.coord[1];
		coord[2] -= pt.coord[2];
		return *this;
	}
	static Point3 rotateBy(const Point3& vec, T angle, const Point3& axis){
		T vecNormSquare = vec.normSquare(), axisNormSquare = axis.normSquare();
		if(vecNormSquare <= Constant<T>::tolerance || axisNormSquare <= Constant<T>::tolerance || angle <= Constant<T>::tolerance)
			return vec;
		Point3 vecProj = axis;
		Point3 vecOrtho = vec;
		Point3 tmp;
		vecProj *=  vec.dotProduct(axis) / axisNormSquare;
		vecOrtho.subtract(vecProj);
		tmp = axis.crossProduct(vecOrtho);
		if(tmp.normSquare() <= Constant<T>::tolerance){
			tmp = axis;
			if(vec.dotProduct(axis) < 0)
				tmp *= -1;
			tmp.normalize();
			return tmp;
		}
		else{
			Point3 rotatedVecOrth = vecOrtho;
			T vecOrthoNorm = vecOrtho.norm();
			vecOrtho.normalize();
			tmp.normalize();
			tmp *= sin(angle) * vecOrthoNorm;
			rotatedVecOrth *= cos(angle);
			rotatedVecOrth.add(tmp);
			rotatedVecOrth.add(vecProj);
			return rotatedVecOrth;
		}
	}
	double dotProduct(const Point3& pt) const{
		return coord[0] * pt.coord[0] + coord[1] * pt.coord[1] + coord[2] * pt.coord[2];
	}
	Point3 crossProduct(const Point3& pt) const{
		return Point3(coord[1] * pt.coord[2] - coord[2] * pt.coord[1],
			coord[2] * pt.coord[0] - coord[0] * pt.coord[2],
			coord[0] * pt.coord[1] - coord[1] * pt.coord[0]);
	}
	bool operator==(const Point3& pt){
		return coord[0] == pt.coord[0] && coord[1] = pt.coord[1] && coord[2] = pt.coord[2];
	}
};
template<class T>
class Transform{
	T m[16];
public:
	Transform(){
		setIdentity();
	}
	Transform(T angle, const Point3<T>& axis){
		setIdentity();
		Point3<T> rotatedX = Point3<T>::rotateBy(Point3<T>(1, 0, 0), angle, axis);
		Point3<T> rotatedY = Point3<T>::rotateBy(Point3<T>(0, 1, 0), angle, axis);
		Point3<T> rotatedZ = Point3<T>::rotateBy(Point3<T>(0, 0, 1), angle, axis);
		m[0] = rotatedX.x();m[4] = rotatedX.y();m[8]=rotatedX.z();
		m[1] = rotatedY.x();m[5] = rotatedY.y();m[9]=rotatedY.z();
		m[2] = rotatedZ.x();m[6] = rotatedZ.y();m[10]=rotatedZ.z();
	}
	Transform(T x, T y, T z){
		setIdentity();
		m[3] = x;
		m[7] = y;
		m[11] = z;
	}
	//(*this) t
	Transform& multiply(const Transform& t){
		for(unsigned int i = 0;i < 16;i += 4){
			T tmp0 = m[i] * t.m[0] + m[i+1] * t.m[4] + m[i+2] * t.m[8] + m[i+3] * t.m[12];
			T tmp1 = m[i] * t.m[1] + m[i+1] * t.m[5] + m[i+2] * t.m[9] + m[i+3] * t.m[13];
			T tmp2 = m[i] * t.m[2] + m[i+1] * t.m[6] + m[i+2] * t.m[10] + m[i+3] * t.m[14];
			T tmp3 = m[i] * t.m[3] + m[i+1] * t.m[7] + m[i+2] * t.m[11] + m[i+3] * t.m[15];
			m[i] = tmp0;m[i+1] = tmp1;m[i+2] = tmp2;m[i+3] = tmp3;
		}
		return *this;
	}
	bool operator==(const Transform& t){
		T sum = 0;
		T tolerance = Constant<T>::tolerance * Constant<T>::tolerance;
		for(int i = 0;i < 16;++ i){
			double delta = m[i] - t.m[i];
			sum += delta * delta;
			if(sum >= tolerance)
				return false;
		}
		return true;
	}
	void setIdentity(){
		m[1] = m[2] = m[3] = 0;
		m[4] = m[6] = m[7] = 0;
		m[8] = m[9] = m[11] = 0;
		m[12] = m[13] = m[14] = 0;
		m[0] = m[5] = m[10] = m[15] = 1;
	}
	Transform(const T* m){
		for(unsigned int i = 0;i < 16;++ i)
			this->m[i] = m[i];
	}
	void transpose(){
		using std::swap;
		swap(m[1], m[4]);
		swap(m[2], m[8]);
		swap(m[3], m[12]);
		swap(m[6], m[9]);
		swap(m[7], m[13]);
		swap(m[11], m[14]);
	}
	template<class U>
	Point3<U> compute(const Point3<U>& pt) const{
		return Point3<U>(m[0] * pt.x() + m[1] * pt.y() + m[2] * pt.z() + m[3], 
			m[4] * pt.x() + m[5] * pt.y() + m[6] * pt.z() + m[7], 
			m[8] * pt.x() + m[9] * pt.y() + m[10] * pt.z() + m[11]);
	}
	friend std::ostream& operator<<(std::ostream& os, const Transform& t){
		using namespace std;
		const char *sep = "\t  ";
		os << t.m[0] << sep << t.m[1] << sep << t.m[2] << sep << t.m[3] << endl;
		os << t.m[4] << sep << t.m[5] << sep << t.m[6] << sep << t.m[7] << endl;
		os << t.m[8] << sep << t.m[9] << sep << t.m[10] << sep << t.m[11] << endl;
		os << t.m[12] << sep << t.m[13] << sep << t.m[14] << sep << t.m[15] << endl;
		return os;
	}
};
struct Light{
	float position[4];
	float ambient[4];
	float specular[4];
	float diffuse[4];
};
typedef Point3<float> Point3f;
typedef Transform<float> TransformF;
class Cylinder{
	
	Point3f ptBottom;
	Point3f ptTop;
	float color[4];
	float topRadius;
	float botRadius;
	double rotateAngle;
	Point3f rotateAxis;
	unsigned int slice;
	char hasCap;
public:
	Cylinder(const float ptB[3], const float ptT[3], const float color[4], float topRadius, float botRadius,
		unsigned int slice = 20, bool hasTopCap = true, bool hasBotCap = true){
		ptBottom = ptB;
		ptTop = ptT;
		this->topRadius = topRadius;
		this->botRadius = botRadius;
		this->slice = slice;
		this->hasCap = 0;
		if(hasTopCap) this->hasCap |= 1;
		if(hasBotCap) this->hasCap |= 2;
		
		Point3f tmp(ptTop);
		Point3f zVec(0, 1, 0);
		rotateAngle = 0;
		rotateAxis = Point3f(zVec.crossProduct(tmp.subtract(ptBottom)));
		if(rotateAxis.norm() > ConstantF::tolerance){
			rotateAngle = acos(tmp.dotProduct(zVec) / tmp.norm());
			rotateAxis.normalize();
		}else if(tmp.dotProduct(zVec) < 0){
			rotateAngle = ConstantF::pi / 2;
			rotateAxis = Point3f(0, 1, 0);
		}
		rotateAngle = rotateAngle / ConstantF::pi * 180;
		memcpy(this->color, color, 4 * sizeof(float));
	}
	void draw(){
		double height = ptBottom.dist(ptTop);
		glColor4fv(color);
		if(height <= ConstantF::tolerance){
			if(std::max(topRadius, botRadius) <= ConstantF::tolerance){
				glBegin(GL_POINTS);
				glVertex3fv(ptBottom);
				glEnd();
				return;
			}
			if(slice <= 1)
				return;
			glPushMatrix();
			applyTransform();
			if(slice == 2){
				glBegin(GL_LINE);
				glVertex3f(botRadius, 0, 0);
				glVertex3f(topRadius, 0, 0);
				glEnd();
			}else{
				double inc = 2 * ConstantF::pi / slice;
				if(fabs(topRadius - botRadius) <= ConstantF::tolerance){
					glBegin(GL_LINES);
					for(unsigned int i = 0;i <= slice;++ i)
						glVertex3f(cos(i * inc)* topRadius, 0, sin(i * inc)* topRadius);
					glEnd();
					if(hasCap){
						drawCap(0, topRadius, NULL);
					}
				}else if(hasCap){

				}
			}
			glPopMatrix();
			return;
		}
		if(slice <= 1)
			return;
		if(std::max(topRadius, botRadius) <= ConstantF::tolerance){
			glBegin(GL_LINE);
			glLineWidth(3);
			glVertex3fv(ptBottom);
			glVertex3fv(ptTop);
			glEnd();
			return;
		}
		glPushMatrix();
		applyTransform();
		if(slice == 2){
			glBegin(GL_QUADS);
			glVertex3f(botRadius, 0, 0);
			glVertex3f(topRadius, height, 0);
			glVertex3f(-topRadius, height, 0);
			glVertex3f(-botRadius, 0, 0);
			glEnd();
		}else{
				double inc = 2 * ConstantF::pi / slice;
				TransformF t(ptBottom.x(), ptBottom.y(), ptBottom.z());
				t.multiply(TransformF(rotateAngle / 180 * ConstantF::pi , rotateAxis));
				glBegin(GL_QUAD_STRIP);
				for(unsigned int i = 0;i <= slice;++ i){
					Point3f pt(cos(i * inc)* botRadius, 0, sin(i * inc)* botRadius);
					Point3f pt1(cos(i * inc)* topRadius, height, sin(i * inc)* topRadius);
					Point3f ptB = t.compute(pt);
					Point3f ptDir = pt;
					ptB.subtract(ptBottom);
					ptDir.subtract(pt1);
					Point3f normVec = ptDir.orth(ptB);
					glNormal3fv(normVec);
					glVertex3fv(pt);
					glNormal3fv(normVec);
					glVertex3fv(pt1);
				}
				glEnd();
				Point3f axis(ptTop);
				axis.subtract(ptBottom);
				axis.normalize();
				if((hasCap & 1) && topRadius > ConstantF::tolerance)
					drawCap(height, topRadius, axis);
				axis *= -1;
				if((hasCap & 2) && botRadius > ConstantF::tolerance){
					drawCap(0, botRadius, axis);
				}
		}
		glPopMatrix();
	}
private:
	void applyTransform(){
		//{
			float matrix[16];
			glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
			TransformF t(matrix);
			t.transpose();
		//}
		if(ptBottom.norm()> ConstantF::tolerance){
			glTranslatef(ptBottom.x(), ptBottom.y(), ptBottom.z());
		}
		if(rotateAngle  > ConstantF::tolerance){
			glRotatef(rotateAngle, rotateAxis.x(), rotateAxis.y(), rotateAxis.z());
		}
		{
			using namespace std;
			glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
			TransformF trans(ptBottom.x(), ptBottom.y(), ptBottom.z());
			TransformF rot(rotateAngle / 180 * ConstantF::pi , rotateAxis);
			TransformF t1(matrix);
			t.multiply(trans);
			t.multiply(rot);
			t1.transpose();
			assert(t == t1);
		}
	}
	void drawCap(float height, float radius, float *norm){
			double inc = 2 * ConstantF::pi / slice;
			glBegin(GL_TRIANGLE_FAN);
			if(norm != NULL)
				glNormal3fv(norm);
			glVertex3f(0, height, 0);
			for(unsigned int i = 0;i <= slice;++ i){
				if(norm != NULL)
					glNormal3fv(norm);
				glVertex3f(cos(i * inc)* radius, height, sin(i * inc)* radius);
			}
			glEnd();
	}
};
static Light lit0={
	{-3, 0, 0, 1},
	{1, 0, 0, 1},
	{1, 0, 0, 1},
	{1, 0, 0, 1}
},
lit1 = {
	{3, 0, 0, 1},
	{0, 1, 0, 1},
	{0, 1, 0, 1},
	{0, 1, 0, 1}
}, 
lit2={
	{0, 3, 0, 1},
	{0, 0, 1, 1},
	{0, 0, 1, 1},
	{0, 0, 1, 1}
};
static bool lightFollowCamera = false;
static void init()
{
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_NORMALIZE);
	glEnable(GL_SMOOTH);
	glEnable(GL_BLEND);
	glEnable( GL_COLOR_MATERIAL );
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClearColor(0.5, 0.5, 0.5,0.0);

	glLightfv(GL_LIGHT0, GL_AMBIENT, lit0.ambient);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lit0.specular);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lit0.diffuse);
	glLightfv(GL_LIGHT1, GL_AMBIENT, lit1.ambient);
	glLightfv(GL_LIGHT1, GL_SPECULAR, lit1.specular);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, lit1.diffuse);
	glLightfv(GL_LIGHT2, GL_AMBIENT, lit2.ambient);
	glLightfv(GL_LIGHT2, GL_SPECULAR, lit2.specular);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, lit2.diffuse);


	glLightfv(GL_LIGHT0, GL_POSITION, lit0.position);
	glLightfv(GL_LIGHT1, GL_POSITION, lit1.position);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHT2);
}

static void drawCylinder(float radius, float x, float y, float z, int slices, int stacks){

}
static void drawAxis()
{
	GLfloat white[] = {1.0, 1.0, 1.0, 1.0};
    GLfloat yellow[]= {1.0, 1.0, 0.0, 1.0};
	GLfloat tmpColor[] = {0.0, 1.0, 1.0, 1.0};
    GLfloat mat_shininess[]={ 30.0 };
	glHint(GL_LINE_SMOOTH, GL_NICEST);
	glHint(GL_POLYGON_SMOOTH, GL_NICEST);
	glDisable(GL_COLOR_MATERIAL);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, white );
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, yellow );
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, tmpColor );
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess );
	{
		float ptBot[3]={0, 0, 0}, ptTop[3]={4, 0, 0}, ptPt[3]={4.2, 0, 0}, color[4] = {1.0, 0.0, 0.0, 0.5}, topRadius = 0.03, botRadius = 0.03;
		Cylinder cylinderAxisX(ptBot, ptTop, color, topRadius, botRadius, 20, false);
		Cylinder cylinderAxisXTop(ptTop, ptPt, color, 0, botRadius * 1.5, 20, false, false);
		cylinderAxisX.draw();
		cylinderAxisXTop.draw();
	}
	{
		float ptBot[3]={0, 0, 0}, ptTop[3]={0, 4, 0}, ptPt[3]={0, 4.2, 0}, color[4] = {0.0, 1.0, 0.0, 0.5}, topRadius = 0.03, botRadius = 0.03;
		Cylinder cylinderAxisY(ptBot, ptTop, color, topRadius, botRadius, 20, false);
		Cylinder cylinderAxisYTop(ptTop, ptPt, color, 0, botRadius * 1.5, 20, false, false);
		cylinderAxisY.draw();
		cylinderAxisYTop.draw();
	}
	{
		float ptBot[3]={0, 0, 0}, ptTop[3]={0, 0, 4}, ptPt[3]={0, 0, 4.2}, color[4] = {0.0, 0.0, 1.0, 0.5}, topRadius = 0.03, botRadius = 0.03;
		Cylinder cylinderAxisZ(ptBot, ptTop, color, topRadius, botRadius, 20, false);
		Cylinder cylinderAxisZTop(ptTop, ptPt, color, 0, botRadius * 1.5, 20, false, false);
		cylinderAxisZ.draw();
		cylinderAxisZTop.draw();
	}
	glEnable(GL_COLOR_MATERIAL);
}
static char ch;
static int Angle;
static void display(void) 
{
	float eyePos[4]={8, 8, 8, 1};
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if(!lightFollowCamera)
		glLightfv(GL_LIGHT2, GL_POSITION, lit2.position);
	else
		glLightfv(GL_LIGHT2, GL_POSITION, eyePos);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(eyePos[0], eyePos[1], eyePos[2], 0, 0, 0, 0, 1, 0);
	
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
	drawAxis();
	glutSwapBuffers();
}
static void reshape(int w, int h)
{
	glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(30.0,1.0,1.0,100.0);
}
static void keyboard(unsigned char key, int x, int y);
int test(int argc, char *argv[])
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize(800, 800);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Draw Function");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	init();
	glutMainLoop();
	return 0;
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
			Angle += 2;ch = key; break;
		case 'i':
			lightFollowCamera = !lightFollowCamera;
	}
	glutPostRedisplay();
}
void arryIndexTest(){
	int tmp[10] = {1, 2, 3, 4,};
}
int main(int argc, char *argv[]){
	extern int hsl_double_cone(int argc, char** argv);
	extern int heatDistributionBar(int argc, char** argv);
	extern int testColorRamp(int argc, char *argv[]);
	extern int testTransparency(int argc, char *argv[]);
	return testTransparency(argc, argv);
}