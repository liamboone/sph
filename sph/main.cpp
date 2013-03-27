// Include standard headers
#include <stdio.h>
#include <stdlib.h>

#include "memleak.h"

#define DEFAULT_WIDTH 640
#define DEFAULT_HEIGHT 480

#include "glew.h"
#include "../freeglut/include/GL/glut.h"
#include "../glm/glm.hpp"

#include "fluid.h"
#include "display.h"

using namespace glm;

int buttonPress;
int old_X;
int old_Y;
bool play = false;
bool singleStep = false;

void resize_cb(int, int);
void display_cb(void);
void keypress_cb(unsigned char, int, int);
void cleanup_cb(void);
void mouseDrag_cb(int, int);
void mouseClick_cb(int, int, int, int);

Display display;
Fluid theFluid;

int main(int argc, char** argv) 
{
	
	//GLUT stuff and making a window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(DEFAULT_WIDTH, DEFAULT_HEIGHT);
	glutCreateWindow("SPH");

	display.init();
	display.setFluids(&theFluid); 

	//Setup Callbacks
	glutDisplayFunc(display_cb);
	glutReshapeFunc(resize_cb);
	glutKeyboardFunc(keypress_cb);
	glutMouseFunc(mouseClick_cb);
	glutMotionFunc(mouseDrag_cb);
	glutIdleFunc(display_cb);

	//Start it off
	glutMainLoop();
	return 0;
}

void mouseClick_cb(int button, int state, int x, int y)
{
	if( state == GLUT_DOWN )
	{
		buttonPress = button;
		if( button == 3 )
		{
			display.zoomCamera( 10 );
		}
		else if( button == 4 )
		{
			display.zoomCamera( -10 );
		}
	}
	display.updateCamera();
	old_X = x;
	old_Y = y;
}

void mouseDrag_cb(int x, int y)
{
	if( buttonPress == GLUT_LEFT_BUTTON )
	{
		display.orbitCamera( x-old_X, y-old_Y );
	}
	if( buttonPress == GLUT_RIGHT_BUTTON )
	{
		display.zoomCamera( y-old_Y );
	}
	old_X = x;
	old_Y = y;
}

void keypress_cb(unsigned char key, int x, int y) {
	switch(key) 
	{
	case 'q': 
	case 'Q':
	case 27: // ascii code of esc key
		_CrtDumpMemoryLeaks();
		system( "pause" );
		exit(0);
		break;
	case '=':
		play = !play;
		break;
	case '>':
		play = true;
		break;
	case '.':
		singleStep = true;
		break;
	}
	glutPostRedisplay();
}

void display_cb() {
	//Always and only do this at the start of a frame, it wipes the slate clean
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	char title[100];
	sprintf( title, "SPH frame:%d", theFluid.frame ); 
	glutSetWindowTitle( title );
	if( play || singleStep )
	{
		singleStep = false;
		theFluid.Update(0.003, glm::vec3(0, -9.8, 0)); 
	}
	display.draw();

	glutSwapBuffers();
}

void resize_cb(int width, int height) {
	//set the viewport, more boilerplate
	glViewport(0, 0, width, height);

	display.updateViewport( width, height );

	glutPostRedisplay();
}