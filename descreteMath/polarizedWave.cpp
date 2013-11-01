/*******************************************************************

 Source courtesy of Virginia Muncy with some modifications by the author
 
 Shows the way that a wave of light, when passed through a vertical
 polarizer, has a change in amplitude and angle where the angel is
 set back to 90 deg and the amplitude is changed to match the max
 vertical height of the original wave. The idea of polarizations is
 to eliminate all light waves that are not entering at the angle of
 polarization. This tends to dim the light but alows for interesting
 experiments in which the angle of entrance of light is significant.
 I also checked the math and science behind the model and it appears
 to be correct.
 
 I added a floor and back wall which show the projections of the
 wave for the xy and xz planes. This helps the user understand what
 is happening to the wave as it is rotated.
 
 I also have the polarized wave change from white to black as its
 amplitude is reduced. This shows the way that light is dimmed as
 it passes through a polarizer based on the angle of entry.
 
 For better display purposes I added the gluLookAt rather than
 accepting the default. I also enabled the depth buffer so that as
 the image is rotated you can see the projections on a solid surface
 by looking at the back side of the back wall of the image.
 
 The user can still control the view and wave angle, but now they
 choose between using the keyboard or a left-mose-button menu.
 
 Source file to be used with
 Cunningham, Computer Graphics: Programming in OpenGL for Visual Communication, Prentice-Hall, 2007

 Intended for class use only

*******************************************************************/
#include <glut.h>
#include 	<math.h>

GLfloat rot   =  0.0;	//the variable angle of rotation of the view
GLfloat	theta = 45.0;	//the variable angle for the initial wave
GLfloat t     =  0.0;	//the animation variable for the wave
GLfloat xx, yy;			//temp variables for later

/**********function key***********
 -sets up the keyboard interface
  for rotation the original wave
  and the view of the model.
*********************************/
static void key(unsigned char c, int x, int y)
 {
	switch(c){	//Keyboard Controls
		case 'z': rot -= 5.0;	break;	//rotates the view
		case 'x': rot += 5.0;	break;
		case 'a': theta -= 5.0;	break;	//rotates the wave
		case 's': theta += 5.0;	break;
	}
	glutPostRedisplay();
 }

/**********function menu***********
 -sets up the menu interface
  for rotation the original wave
  and the view of the model.
*********************************/
static void menu(int choice)
 {
	switch(choice){	//Mouse Controlls
		case 0:	theta += 5.0;	break;	//rotates the wave
		case 1:	theta -= 5.0;	break;
		case 2:	rot -= 5.0;	break;	//rotates the view
		case 3:	rot += 5.0;	break;
	}
	glutPostRedisplay();	
 }

/**********function init**********
 -initializes the background color
  of the window as well as setting
  up the shading and transparency
  for the program.
*********************************/
static void init(void) 
 {
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	glutCreateMenu(menu);
	glutAddMenuEntry("Rotate Wave Left",0);
	glutAddMenuEntry("Rotate Wave Right",1);
	glutAddMenuEntry("Rotate View Left",2);
	glutAddMenuEntry("Rotate View Right",3);
	glutAttachMenu(GLUT_LEFT_BUTTON);
 }

/********function display*********
 -sets the colors and draws the
  waves and the polorizing film.
  Translates them all into their
  propper positions.
*********************************/
static const double M_PI = asin(1.0)*2;
static void display(void) 
 {
	float	wavex, wavey, wavez;
	float	x, y;
	float	bright;
	int	i, angle;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-2.0,2.0,-2.0,2.0,-20.0,20.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(0.0,0.0,0.0,-0.5,-0.5,-1.0,0.0,1.0,0.0);	//default view
	
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_DEPTH);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT,GL_NICEST);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT,GL_NICEST);
    glEnable(GL_POLYGON_SMOOTH);
    glHint(GL_POLYGON_SMOOTH_HINT,GL_NICEST);

	glRotatef(rot,0,1,0);
	glLineWidth(3.0);

	//update animation parameter
	t += .01;
	
	glPushMatrix();
	// background of scene
	// translate and rotate to get right placement
	glRotatef(25., 0., 1., 0.);
	glBegin(GL_QUADS);
		glColor3f(0., .3, 0.);
		glVertex3f(4.,4.,-2.);
		glVertex3f(4.,-4.,-2.);
		glVertex3f(-4.,-4.,-2.);
		glVertex3f(-4.,4.,-2.);
	glEnd();
	glBegin(GL_TRIANGLE_FAN);
		glColor3f(0., .9, 0.);
		glVertex3f(0., -.5, -1.95);
		glColor3f(0., .3, 0);
		for (i=0; i<21; i++)
		{
			glVertex3f( 3.*cos(2.*M_PI*float(i)/20.), 3.*sin(2.*M_PI*float(i)/20.), -1.95);
		}
	glEnd();
	glPopMatrix();

	glPushMatrix();
	//back panel
	glTranslatef(-1.0,-0.5,-0.5);
	glBegin(GL_QUADS);
	glColor4f(0.5,0.7,0.7,1.0);
	glVertex3f(0.0,0.0,0.0);
	glVertex3f(2.0,0.0,0.0);
	glVertex3f(2.0,1.0,0.0);
	glVertex3f(0.0,1.0,0.0);
	glEnd();
	glPopMatrix();
	
	glPushMatrix();
	//floor panel
	glTranslatef(-1.0,-0.5,-0.5);
	glRotatef(90,1,0,0);
	glBegin(GL_QUADS);
	glColor4f(0.5,0.7,0.7,1.0);
	glVertex3f(0.0,0.0,0.0);
	glVertex3f(2.0,0.0,0.0);
	glVertex3f(2.0,1.0,0.0);
	glVertex3f(0.0,1.0,0.0);
	glEnd();
	glPopMatrix();
	
	glPushMatrix();
	//original wave
	glRotatef(theta,1,0,0);
	glBegin(GL_LINE_STRIP);
		glColor4f(1.0,1.0,1.0,1.0);
		for (wavex=-1; wavex<=0; wavex+=.005){
			wavey=sin(16.0*wavex-t)/4.0;
			if (wavex == -1.) {xx = wavex; yy = wavey;}
			glVertex3f(wavex,wavey,0.0);
		}
		glVertex3f(0., 0., 0.);
		glVertex3f(-1., 0., 0.);
		glVertex3f(xx, yy, 0.);
	glEnd();
	glPopMatrix();
	
	glPushMatrix();				//projection on back wall
	glTranslatef(0.0,0.0,-0.5);
	glBegin(GL_LINE_STRIP);
		glColor4f(0.8,0.8,0.8,1.0);
		for (wavex=-1; wavex<=0; wavex+=.005){
			wavey=sin(16.0*wavex-t)/4.0;
			wavey=wavey*sin((90.0-theta)*M_PI/180);
			glVertex3f(wavex,wavey,0.0);
		}
	glEnd();
	glPopMatrix();

	glPushMatrix();				//projection on floor
	glTranslatef(0.0,-0.5,0.0);
	glBegin(GL_LINE_STRIP);
	glColor4f(0.8,0.8,0.8,1.0);
	for(wavex=-1; wavex<=0; wavex+=.005){
		wavez=sin(16.0*wavex-t)/4.0;
		wavez=wavez*cos((90.0-theta)*M_PI/180);
		glVertex3f(wavex,0.0,wavez);
	}
	glEnd();
	glPopMatrix();

	glPushMatrix();				//polarizer
	glRotatef(90,0,1,0);
	glTranslatef(-0.5,-0.5,0.0);
	glBegin(GL_QUAD_STRIP);
	for(y=0; y<10; y+2){
		x=y/10;
		glColor4f(1.0,1.0,0.5,0.5);
		glVertex3f(x,0.0,0.0);
		glColor4f(1.0,1.0,0.5,0.5);
		glVertex3f(x,1.0,0.0);
		y++;
		x=y/10;
		glColor4f(0.5,0.5,0.5,0.5);
		glVertex3f(x,0.0,0.0);
		glColor4f(0.5,0.5,0.5,0.5);
		glVertex3f(x,1.0,0.0);		
	}
	glEnd();
	glPopMatrix();

	glPushMatrix();
	//polarized wave
	glBegin(GL_LINE_STRIP);
	angle=fabs(theta);
	if(angle%180<90){
		angle=angle%90;
		bright=angle/90.0;
		bright=1-bright;
	}
	else{
		angle=angle%90;
		bright=angle/90.0;
	}
	glColor4f(bright,bright,bright,1.0);
	for(wavex=0.; wavex<=1.; wavex+=.005){
		wavey=sin(16.0*wavex-t)/4.0;
		wavey=wavey*sin((90.0-theta)*M_PI/180);	//changes amplitude based
		if (wavex == 0.) { xx = wavex; yy = wavey; }
		glVertex3f(wavex,wavey,0.0);			//on the first wave's angle
	}
	glVertex3f(1., 0., 0.);
	glVertex3f(0., 0., 0.);
	glVertex3f(xx, yy, 0.);
	glEnd();
	glPopMatrix();
	
	glPushMatrix();				//projection on back wall
	glTranslatef(0.0,0.0,-0.5);
	glBegin(GL_LINE_STRIP);
	glColor4f(0.8,0.8,0.8,1.0);
	for(wavex=0; wavex<=1; wavex+=.005){
		wavey=sin(16.0*wavex-t)/4.0;
		wavey=wavey*sin((90.0-theta)*M_PI/180);
		glVertex3f(wavex,wavey,0.0);
	}
	glEnd();
	glPopMatrix();

	glPushMatrix();				//projection on floor
	glTranslatef(0.0,-0.5,0.0);
	glBegin(GL_LINE_STRIP);
		glColor4f(0.8,0.8,0.8,1.0);
		for(wavex=0; wavex<=1; wavex+=.005){
			glVertex3f(wavex,0.0,0.0);
		}
	glEnd();
	glPopMatrix();

	glutSwapBuffers(); 
}

/********function resahape*********
 -sets up the window so that the
  picture is not distorted if the
  window is reshaped.
**********************************/
static void reshape(int w, int h)
 {
	glViewport ( 0, 0, (GLsizei)w, (GLsizei)h);
	//glMatrixMode(GL_PROJECTION);
	//glLoadIdentity();
	//gluPerspective(60.0, 1.0, 1.0, 30.0);
	//glMatrixMode(GL_MODELVIEW);
	//glLoadIdentity();
	//gluLookAt(0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0);

	if (w >= h)			//prevents distortion in the window
		glViewport(0, 0, (GLsizei)h, (GLsizei)h);
	else
		glViewport(0, 0, (GLsizei)w, (GLsizei)w);
 }
 
static void idle(void)
	{
		glutPostRedisplay();
	}

int polarized_wave(int argc, char** argv) 
 { 
	glutInit (&argc, argv);
	glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutInitWindowSize (800, 800);
	glutInitWindowPosition (0, 0);
	glutCreateWindow ("Vertical Polarization of a Wave of Light");
	
	init();

	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(key);
	glutIdleFunc(idle);

	glutMainLoop();
	
	return 0;
  }