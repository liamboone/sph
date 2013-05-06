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

#include <iostream>
#include <fstream>

using namespace glm;


int theFrameNum = 0; 
bool isRecording = false;
bool writeToFile = false; 
ofstream myfile;
bool displayOn = true;
bool useRaymarch = false;
bool gravity = true;
bool vorticity = false; 

int buttonPress;
int old_X;
int old_Y;
bool play = false;
bool singleStep = false;
int screenWidth = 640;
int screenHeight = 480;
int displayFlags = 0xFF;

void resize_cb(int, int);
void display_cb(void);
void keypress_cb(unsigned char, int, int);
void cleanup_cb(void);
void mouseDrag_cb(int, int);
void mouseClick_cb(int, int, int, int);

Display display;

//Fluid theFluid;
int theMenu = 0;


Fluid theFluid(vec3(-1.5f, 0.0f, -1.5f), vec3(1.5f, 5.3f, 1.5f));

//Fluid theFluid(vec3(-0.5f, 0.0f, -0.5f), vec3(0.5f, 1.0f, 0.5f));



int frame = 0;
int timebase = 0;
float fps = 0;

void onMenuCb(int value)
{
   switch (value)
   {
   case -1: exit(0);
   default: keypress_cb(value, 0, 0); break;
   }
}

int main(int argc, char** argv) 
{
	
	//GLUT stuff and making a window
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowSize(DEFAULT_WIDTH, DEFAULT_HEIGHT);
	glutCreateWindow("SPH");

	//Window for controls (middle mouse click)
	int subMenu1 = glutCreateMenu(onMenuCb);
    glutAddMenuEntry("Pause/Resume    =", '=');
    glutAddMenuEntry("Single step     .", '.');
    glutAddMenuEntry("Display on/off  d", 'd');
	glutAddMenuEntry("Record          r", 'r');

    int subMenu2 = glutCreateMenu(onMenuCb);
    glutAddMenuEntry("Color by temperature t", 't');
	glutAddMenuEntry("Color by velocity    v", 't');
    glutAddMenuEntry("Raymarch on/off      m", 'm');
	glutAddMenuEntry("Toggle fluids on/off 4", '4'); 

	int subMenu3 = glutCreateMenu(onMenuCb); 
	glutAddMenuEntry("Gravity on/off    g", 'g');

    theMenu = glutCreateMenu(onMenuCb);
    glutAddMenuEntry("Start    >", '>');
    glutAddSubMenu("Play Controls", subMenu1);
    glutAddSubMenu("Display Controls", subMenu2);
	glutAddSubMenu("Simulation Controls", subMenu3);
    glutAddMenuEntry("_________________", -1);
    glutAddMenuEntry("Exit", 27);
    glutAttachMenu(GLUT_MIDDLE_BUTTON);

	display.init();
	display.setFluids(&theFluid); 
	display.setFlags( displayFlags );

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
	//Close file if we are writing to txt 
	if (writeToFile) myfile.close(); 
	return 0;
}



void mouseClick_cb(int button, int state, int x, int y)
{
	if( state == GLUT_DOWN )
	{
		buttonPress = button;
		if( button == 3 )
		{
			display.zoomCamera( -10 );
		}
		else if( button == 4 )
		{
			display.zoomCamera( 10 );
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
	case 'v':
		displayFlags ^= DFLAG_VEL;
		display.setFlags( displayFlags );
		break;
	case 't':
		displayFlags ^= DFLAG_TEMP;
		display.setFlags( displayFlags );
		break;
	case 'r':
		isRecording = !isRecording; 
		if (isRecording) theFrameNum = 0;
		break;
	case 'm':
		useRaymarch = !useRaymarch; 
		break;
	case 'g':
		gravity = !gravity;
		break;
	case 's':
		vorticity = !vorticity; 
		gravity = !gravity; 
		break;
	case 'w':
		writeToFile = !writeToFile; 
		break;
	case '1':
	case '2':
	case '3':
	case '4':
		int ord = key - '1';
		displayFlags ^= 1 << ord;
		display.setFlags( displayFlags );
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

//Stores position for each frame for rendering in maya
void writeToText()
{

	if (myfile.is_open()) {
		myfile<<theFluid.frame;
		myfile<<"\n";

		myfile<<theFluid.getParticles().size();
		myfile<<"\n";
		for (int i = 0; i < theFluid.getParticles().size(); i++) {
			myfile<<theFluid.getParticles().at(i)->getIndex(); 
			myfile<<" "; 
			myfile<<theFluid.getParticles().at(i)->getPosition().x; 
			myfile<<" "; 
			myfile<<theFluid.getParticles().at(i)->getPosition().y;
			myfile<<" "; 
			myfile<<theFluid.getParticles().at(i)->getPosition().z; 
			myfile<<" "; 
		}
		myfile<<"\n";
	} else {		
		myfile.open ("mayaImport.txt");
		myfile<<theFluid.frame;
		myfile<<"\n";
		myfile<<theFluid.getParticles().size();
		myfile<<"\n";
		for (int i = 0; i < theFluid.getParticles().size(); i++) {
		    myfile<<theFluid.getParticles().at(i)->getIndex(); 
			myfile<<" "; 
			myfile<<theFluid.getParticles().at(i)->getPosition().x; 
			myfile<<" "; 
			myfile<<theFluid.getParticles().at(i)->getPosition().y;
			myfile<<" "; 
			myfile<<theFluid.getParticles().at(i)->getPosition().z; 
			myfile<<" "; 
		}
		myfile<<"\n";
	}
}

vec3 gravityForce( vec3 p )
{
	return vec3( 0, -9.81, 0 );
}

vec3 noForce( vec3 p )
{
	return vec3( 0 );
}

vec3 vortexForce( vec3 p )
{
	vec2 f( p.x, p.z );
	float r = glm::length( f );
	f = glm::normalize( f );
	vec3 F = 1.0f/(r+0.001f)*vec3( -f.y, 0, f.x );
	return vec3( 0, -9.8, 0 ) + F;
}

void display_cb() {
	//Always and only do this at the start of a frame, it wipes the slate clean
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	char title[100];
	
	frame++;
	int time=glutGet(GLUT_ELAPSED_TIME);

	if (time - timebase > 1000) {
		fps = frame*1000.0f/(time-timebase);
	 	timebase = time;
		frame = 0;
	}
	
	sprintf_s( title, 100, "%0.2f FPS, SPH frame: %d, #particles: %d", 
		fps,
		theFluid.frame, 
		theFluid.getParticles().size() ); 
	glutSetWindowTitle( title );

	if( play || singleStep )
	{
		force_t externalForce;
		if (vorticity )
		{
			externalForce = vortexForce; 
		} else {
			externalForce = gravity ? gravityForce : noForce;
		}
		singleStep = false;
		theFluid.Update(0.004f, externalForce);
	}
	if( displayOn )
	{
		if( useRaymarch )
			display.march();
		else
			display.draw();
	}
	if (isRecording) grabScreen(); 
	if (writeToFile) writeToText();
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

