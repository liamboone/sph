#ifndef _SPH_CAMERA_H_
#define _SPH_CAMERA_H_

// c++ libraries
#include <iostream>

#include "memleak.h"

// glm
#include "../glm/glm.hpp"
#include "../glm/gtc/matrix_transform.hpp"

using namespace glm;

class Camera
{
public:
	Camera( float fovy = 60.0f, 
			vec3 pos = vec3(5.0f, 5.0f, 5.0f), 
			vec3 target = vec3(0.0f, 0.0f, 0.0f), 
			vec3 up = vec3(0.0f, 1.0f, 0.0f), 
			float znear = 0.1f, 
			float zfar = 30.0f, 
			int width = 640, 
			int height = 480 )
	{ 
		init( fovy, pos, target, up, znear, zfar, width, height ); 
	}

	void init( float fovy, vec3 pos, vec3 target, vec3 up, float znear, float zfar, int width, int height );
	void setViewport( int w, int h );
	void setPos( vec3 newpos ) { pos = newpos; }
	
	int getWidth() { return width; }
	int getHeight() { return height; }
	vec3 getPos() { return pos; }
	mat4 getMat4();

	void zoom( int dz );
	void orbit( int dx, int dy );

private:
	vec3 pos;
	vec3 target;
	vec3 up;

	float fovy;
	float zNear;
	float zFar;

	int height;
	int width;
};

#endif