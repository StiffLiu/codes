#include <GL/glut.h>
#include <cstdlib>
#include <iostream>
static int width = 900, height = 600;
static int startX = 600, startY = 300;
static bool mouseIn = false;
static bool drawDecomposedOrderedPlane = false;
static void init(){
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glClearColor(0, 0, 0, 1);
}
static float xyPlane[][3] = {
	{-1, -1, 0},{1, -1, 0},{1, 1, 0},{-1, 1, 0},
};
static float yzPlane[][3] = {
	{0, -1, -1},{0, 1, -1},{0, 1, 1},{0, -1, 1},
};
static float zxPlane[][3] = {
	{-1, 0, -1},{-1, 0, 1},{1, 0, 1},{1, 0, -1},
};
static float red[4]={
	1.0, 0, 0, 0.5,
};
static float green[4]={
	0, 1.0, 0, 0.5,
};
static float blue[4]={
	0.0, 0, 1.0, 0.5,
};
static float decomposedXYPlane[][3] = {
	{0, 0, 0},{1, 0, 0},{1, 1, 0},{0, 1, 0},
	{0, 0, 0},{0, 1, 0},{-1, 1, 0},{-1, 0, 0},
	{0, 0, 0},{-1, 0, 0},{-1, -1, 0},{0, -1, 0},
	{0, 0, 0},{0, -1, 0},{1, -1, 0},{1, 0, 0},
};
static float decomposedYZPlane[][3] = {
	{0, 0, 0},{0, 0, -1},{0, 1, -1},{0, 1, 0},
	{0, 0, 0},{0, 1, 0},{0, 1, 1},{0, 0, 1},
	{0, 0, 0},{0, 0, 1},{0, -1, 1},{0, -1, 0},
	{0, 0, 0},{0, -1, 0},{0, -1, -1},{0, 0, -1},
};
static float decomposedZXPlane[][3] = {
	{0, 0, 0},{1, 0, 0},{1, 0, -1},{0, 0, -1},
	{0, 0, 0},{0, 0, -1},{-1, 0, -1},{-1, 0, 0},
	{0, 0, 0},{-1, 0, 0},{-1, 0, 1},{0, 0, 1},
	{0, 0, 0},{0, 0, 1},{1, 0, 1},{1, 0, 0},
};
static float orderedDecomposedPlanes[][3]={
	{0, 0, 0},{0, 0, -1},{-1, 0, -1},{-1, 0, 0}, //blue
	{0, 0, 0},{0, 1, 0},{-1, 1, 0},{-1, 0, 0},   //red
	{0, 0, 0},{-1, 0, 0},{-1, -1, 0},{0, -1, 0}, //red
	{0, 0, 0},{0, 0, -1},{0, 1, -1},{0, 1, 0},   //green
	{0, 0, 0},{0, -1, 0},{0, -1, -1},{0, 0, -1}, //green
	{0, 0, 0},{-1, 0, 0},{-1, 0, 1},{0, 0, 1},   //blue
	{0, 0, 0},{1, 0, 0},{1, 0, -1},{0, 0, -1},   //blue
	{0, 0, 0},{0, 1, 0},{0, 1, 1},{0, 0, 1},     //green
	{0, 0, 0},{0, 0, 1},{0, -1, 1},{0, -1, 0},   //green
	{0, 0, 0},{1, 0, 0},{1, 1, 0},{0, 1, 0},     //red
	{0, 0, 0},{0, -1, 0},{1, -1, 0},{1, 0, 0},   //red
	{0, 0, 0},{0, 0, 1},{1, 0, 1},{1, 0, 0},     //blue
};
static float *colors[]={blue, red, red, green, green, blue, blue, green, green, red, red, blue};
#define ASIZE(a) (sizeof a / sizeof *a)
static void drawCube(float (*planes[])[4][3], float (*colors[])[4], int n){
	glBegin(GL_QUADS);
	for(int i = 0;i < n;++ i){
		float (*plane)[4][3] = planes[i];
		glColor4fv((*(colors + i))[0]);
		glVertex3fv((*planes[i])[0]);
		glVertex3fv((*planes[i])[1]);
		glVertex3fv((*planes[i])[2]);
		glVertex3fv((*planes[i])[3]);
	}
	glEnd();
}
static void drawDecomposedCube(){
	glBegin(GL_QUADS);
	glColor4f(1.0, 0.0, 0.0, 0.5);
	for(int i = 0;i < ASIZE(decomposedXYPlane);++i)
		glVertex3fv(decomposedXYPlane[i]);
	glColor4f(0.0, 1.0, 0.0, 0.5);
	for(int i = 0;i < ASIZE(decomposedYZPlane);++i)
		glVertex3fv(decomposedYZPlane[i]);
	glColor4f(0.0, 0.0, 1.0, 0.5);
	for(int i = 0;i < ASIZE(decomposedZXPlane);++i)
		glVertex3fv(decomposedZXPlane[i]);
	glEnd();
}
static void setTransformation(){
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-1.5, 1.5, -1.5, 1.5, 0, 4);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(1, 1, 1, 0, 0, 0, 0, 1, 0);
}
static void mousePassiveMotionFunc(int x, int y){
	int halfWidth = width / 3;
	int halfHeight = height / 2;
	int sideLen = (halfWidth > halfHeight ? halfHeight : halfWidth);
	y = height - y;
	mouseIn = (x > startX && x < startX + sideLen &&
	   y > startY && y < startY + sideLen);
	//std::cout << x << ',' << y << std::endl;
	glutPostRedisplay();		
}
static void drawByOrder(){
	int halfWidth = width / 3;
	int halfHeight = height / 2;
	int sideLen = (halfWidth > halfHeight ? halfHeight : halfWidth);
	int params[4];
	if(mouseIn)
		glClearColor(0.5, 0.5, 0.5, 1.0);
	glEnable(GL_SCISSOR_TEST);
	glGetIntegerv(GL_SCISSOR_BOX, params);
	glScissor(startX, startY , sideLen, sideLen);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPushMatrix();
	glViewport(startX, startY , sideLen, sideLen);
	setTransformation();
	glBegin(GL_QUADS);
	for(unsigned int i = 0;i < ASIZE(colors);++ i){
		glColor4fv(colors[i]);
		glVertex3fv(orderedDecomposedPlanes[4 * i]);
		glVertex3fv(orderedDecomposedPlanes[4 * i + 1]);
		glVertex3fv(orderedDecomposedPlanes[4 * i + 2]);
		glVertex3fv(orderedDecomposedPlanes[4 * i + 3]);
	}
	glEnd();
	glPopMatrix();
	glScissor(params[0], params[1] , params[2], params[3]);
	glDisable(GL_SCISSOR_TEST);
	glClearColor(0.0, 0.0, 0.0, 1.0);
}
void showInfo(int xLeft, int yBottom, int winWidth, const char *info )
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
static void display(){
	int halfWidth = width / 3;
	int halfHeight = height / 2;
	int sideLen = (halfWidth > halfHeight ? halfHeight : halfWidth);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	const char *reprents[] = {
		"XY:R-YZ:G-ZX:B", "XY:R-ZX:B-YZ:G",
		"YZ:G-XY:R-ZX:B", "YZ:G-ZX:B-XY:R",
		"ZX:B-XY:R-YZ:G", "ZX:B-YZ:G-XY:R",
	};
	/*float (*planes[3])[4][3]={
		xyPlane, yzPlane, zxPlane,
	};*/
	typedef float (*Plane[3])[4][3];
	Plane planes[]={
		{&xyPlane, &yzPlane, &zxPlane},{&xyPlane, &zxPlane, &yzPlane},
		{&yzPlane, &xyPlane, &zxPlane},{&yzPlane, &zxPlane, &xyPlane},
		{&zxPlane, &xyPlane, &yzPlane},{&zxPlane, &yzPlane, &xyPlane},
	};
	typedef float(*Colors[3])[4];
	Colors colors[]={
		{&red, &green, &blue}, {&red, &blue, &green}, 
		{&green, &red, &blue}, {&green, &blue, &red}, 
		{&blue, &red, &green}, {&blue, &green, &red}, 
	};
	if(true)
		for(int i = 0;i < 2;++ i)
			for(int j = 0;j < 3;++ j){
				glPushMatrix();
				glViewport(j * sideLen, i * sideLen, sideLen, sideLen);
				setTransformation();
				drawCube(planes[i * 3 + j], colors[i * 3 + j],3);
				showInfo(j * sideLen, i * sideLen, sideLen, reprents[i * 3 + j]);
				glPopMatrix();	
			}
	if(drawDecomposedOrderedPlane){
		drawByOrder();
	}
	glutSwapBuffers();
}
static void reshape(int width, int height){
	::width = width;
	::height = height;
	glutPostRedisplay();
}
static void keyboardFunc(unsigned char ch, int x, int y){
	switch(ch){
		case 'S':drawDecomposedOrderedPlane = !drawDecomposedOrderedPlane;break;
		case 'w':++startY;break;
		case 's':--startY;break;
		case 'a':--startX;break;
		case 'd':++startX;break;
		
		
	}
	glutPostRedisplay();
}
int main(int argc, char *argv[]){
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(width, height);
	glutCreateWindow("Test Tranparency");
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboardFunc);
	glutPassiveMotionFunc(mousePassiveMotionFunc);
	init();
	glutMainLoop();
	return 0;
}
