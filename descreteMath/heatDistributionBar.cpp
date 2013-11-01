/*
	Project example - heat flow in a thin rectangular body
	
	Assumptions are simple -- body is modeled as a set of small squares, and heat flows
	between squares based on the relative heats of the two squares and a diffusion
	function that includes a constant that models the heat conductivity of the bodies.
	
	Source file to be used with
    Cunningham, Computer Graphics: Programming in OpenGL for Visual Communication, Prentice-Hall, 2007
   
    Intended for class use only
*/
#include <glut.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define ROWS 10			//	body is assumed to be ROWSxCOLS (unitless) squares
#define COLS 30
#define AMBIENT 25.0	//	initial ambient temperature in degrees Celsius
#define HOT	50.0		//	hot temperature of heat-source cell
#define COLD 0.0		//	cold temperature of cold-sink cell
#define NHOTS 4
#define NCOLDS 5

static GLfloat angle = 0.0;
static GLfloat temps[ROWS][COLS], back[ROWS+2][COLS+2];
static GLfloat theta = 0.0, vp = 30.0;
//	set locations of fixed hot and cold spots on the bar
static int hotspots[NHOTS][2] = { {ROWS/2,0},{ROWS/2-1,0},{ROWS/2-2,0},{0,3*COLS/4} };
static int coldspots[NCOLDS][2] = { {ROWS-1,COLS/3}, {ROWS-1,1+COLS/3}, {ROWS-1,2+COLS/3}, {ROWS-1,3+COLS/3}, {ROWS-1,4+COLS/3} };

static void myinit(void);
static void cube(void);
static void display(void);
static void setColor(float);
static void reshape(int, int);
static void animate(void);
static void iterationStep(void);

static void myinit(void)
{
	int i,j;
	
	glEnable (GL_DEPTH_TEST);
	glClearColor(0.6, 0.6, 0.6, 1.0);

//	set up initial temperatures in cells
	for (i=0; i<ROWS; i++) {
		for (j=0; j < COLS; j++) {
			temps[i][j] = AMBIENT;
		}
	}
	for (i=0; i<NHOTS; i++) {
		temps[hotspots[i][0]][hotspots[i][1]]=HOT; }
	for (i=0; i<NCOLDS; i++) {
		temps[coldspots[i][0]][coldspots[i][1]]=COLD; }
}

//	Unit cube in first octant
static void cube (void)
{
	typedef GLfloat point [3];

	point v[8] = {
		{0.0, 0.0, 0.0}, {0.0, 0.0, 1.0},
		{0.0, 1.0, 0.0}, {0.0, 1.0, 1.0},
		{1.0, 0.0, 0.0}, {1.0, 0.0, 1.0},
		{1.0, 1.0, 0.0}, {1.0, 1.0, 1.0} };

	glBegin (GL_QUAD_STRIP);
		glVertex3fv(v[4]);
		glVertex3fv(v[5]);
		glVertex3fv(v[0]);
		glVertex3fv(v[1]);
		glVertex3fv(v[2]);
		glVertex3fv(v[3]);
		glVertex3fv(v[6]);
		glVertex3fv(v[7]);
	glEnd();

	glBegin (GL_QUAD_STRIP);
		glVertex3fv(v[1]);
		glVertex3fv(v[3]);
		glVertex3fv(v[5]);
		glVertex3fv(v[7]);
		glVertex3fv(v[4]);
		glVertex3fv(v[6]);
		glVertex3fv(v[0]);
		glVertex3fv(v[2]);
	glEnd();
}

static void display( void )
{
	#define SCALE 10.0
	
	int i,j;
	
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	//           eye point      center of view       up
	gluLookAt(vp, vp/2., vp/4., 0.0, 0.0, 0.0,  0.0, 0.0, 1.0);
	
	//rotate through given angle
	glPushMatrix();
	glRotatef(angle, 0.,0.,1.);
    
	//	Draw the bars
	for (i = 0; i < ROWS; i++) {
		for (j = 0; j < COLS; j++) {
			//glColor3f( temps[i][j]/HOT, 0.0, 1.0-temps[i][j]/HOT );	// hotter redder; colder bluer
			setColor(temps[i][j]);
			glPushMatrix();
			glTranslatef((float)i-(float)ROWS/2.0,(float)j-(float)COLS/2.0,0.0);
			glScalef(1.0,1.0,0.1+4.9*temps[i][j]/HOT);			// 0.1 cold, 4.0 hot
			cube();
			glPopMatrix();
			}
		}
	glutSwapBuffers();
	glPopMatrix();
}

static void setColor(float t)
{
	float r, g, b;
	
	r = t/HOT;
	g = 0.0;
	b = 1.0 - t/HOT;
	glColor3f(r, g, b);
}

static void reshape(int w,int h)
{
	glViewport(0,0,(GLsizei)w,(GLsizei)h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60.0, (float)w/(float)h, 1.0, 300.0);
}

static void animate(void)
{
	iterationStep();
	glutPostRedisplay();
}

static void iterationStep(void)
{
	int i, j, m, n;
	// this filter was changed from the value in the text to model heat moving more slowly
	float filter[3][3]={{ 0.     , 0.0125, 0.     },
						{ 0.0125 , 0.95,   0.0125 },
						{ 0.     , 0.0125, 0.     } };
	
	//	increment temperatures throughout the material
	for (i=0; i<ROWS; i++)	// first back temps up so you can recreate temps
		for (j=0; j<COLS; j++)
			back[i+1][j+1] = temps[i][j]; // note that we leave boundaries on back
	//	now fill the boundaries with what's next to them from the original temps[][]
	for (i=1; i<ROWS+2; i++) {
		back[i][0]=back[i][1];
		back[i][COLS+1]=back[i][COLS];
		}
	for (j=0; j<COLS+2; j++) {
		back[0][j] = back[1][j];
		back[ROWS+1][j]=back[ROWS][j];
		}
	for (i=0; i<ROWS; i++)	// use diffusion based on back values to compute temp
		for (j=0; j<COLS; j++) {
			temps[i][j]=0.0;
			for (m=-1; m<=1; m++)
				for (n=-1; n<=1; n++)
					temps[i][j]+=back[i+1+m][j+1+n]*filter[m+1][n+1];
		}
	for (i=0; i<NHOTS; i++) {
		temps[hotspots[i][0]][hotspots[i][1]]=HOT; }
	for (i=0; i<NCOLDS; i++) {
		temps[coldspots[i][0]][coldspots[i][1]]=COLD; }
		
	angle += 0.1; if (angle > 360.) angle -= 360.;
}

int heatDistributionBar(int argc, char** argv)
{
/* Standard GLUT initialization */
	    glutInit(&argc,argv);
	    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH); 
	    glutInitWindowSize(700,700);
	    glutInitWindowPosition(0,0);
		glutCreateWindow("Temperature in bar");
	    glutDisplayFunc(display);
	    glutReshapeFunc(reshape);
	    glutIdleFunc(animate);
		
        myinit();
	    glutMainLoop(); /* enter event loop */
		return 0;
}