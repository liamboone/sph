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
#include "..\ObjCore\obj.h"
#include "..\ObjCore\objloader.h"

#include <IL/il.h>
#include <IL/ilu.h>
#include <IL/ilut.h>

using namespace glm;


int theFrameNum = 0; 
bool isRecording = false;
bool displayOn = true;

int buttonPress;
int old_X;
int old_Y;
bool play = false;
bool singleStep = false;
int screenWidth = 640;
int screenHeight = 480;

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

	//Recording - devIL
	ilInit();
    iluInit();
    ilEnable(IL_FILE_OVERWRITE);
    ilutRenderer(ILUT_OPENGL);

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
	case 'd':
		displayOn = !displayOn;
		break;
	case 'r':
		isRecording = !isRecording; 
		if (isRecording) theFrameNum = 0;
		break;
	}
	glutPostRedisplay();
}

//Taken from previous 593 assignments (JelloSim) 
void grabScreen()  
{
    unsigned int image;
    ilGenImages(1, &image);
    ilBindImage(image);

    ILenum error = ilGetError();
    assert(error == IL_NO_ERROR);

	ilTexImage(screenWidth, screenHeight, 1, 3, IL_RGB, IL_UNSIGNED_BYTE, NULL);

    error = ilGetError();
    assert(error == IL_NO_ERROR);

    unsigned char* data = ilGetData();

    error = ilGetError();
    assert(error == IL_NO_ERROR);

    for (int i=screenHeight-1; i>=0; i--) 
    {
	    glReadPixels(0,i,screenWidth,1,GL_RGB, GL_UNSIGNED_BYTE, 
		    data + (screenWidth * 3 * i));
    }

    char anim_filename[100];
    sprintf_s(anim_filename, 100, "output/%04d.png", theFrameNum++); 

    ilSave(IL_PNG, anim_filename);

    error = ilGetError();
    assert(error == IL_NO_ERROR);

    ilDeleteImages(1, &image);

    error = ilGetError();
    assert(error == IL_NO_ERROR);
}

void display_cb() {
	//Always and only do this at the start of a frame, it wipes the slate clean
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	char title[100];
	sprintf( title, "SPH frame: %d, #particles: %d", theFluid.frame, theFluid.getParticles().size() ); 
	glutSetWindowTitle( title );
	if( play || singleStep )
	{
		singleStep = false;
		theFluid.Update(0.004, glm::vec3(0, -9.8, 0)); 
		/*theFluid.Update(0.0008, glm::vec3(0, -9.8, 0)); 
		theFluid.Update(0.0008, glm::vec3(0, -9.8, 0)); 
		theFluid.Update(0.0008, glm::vec3(0, -9.8, 0)); 
		theFluid.Update(0.0008, glm::vec3(0, -9.8, 0));  */
	}
	if( displayOn )
	{
		display.draw();
	}
	else
	{
		display.march();
	}
	if (isRecording) grabScreen(); 
	glutSwapBuffers();
}

void resize_cb(int width, int height) {
	//set the viewport, more boilerplate
	screenWidth = width;
	screenHeight = height;
	glViewport(0, 0, width, height);

	display.updateViewport( width, height );

	glutPostRedisplay();
}

