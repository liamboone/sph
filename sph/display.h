#ifndef _SPH_DISPLAY_H_
#define _SPH_DISPLAY_H_

// c++ libraries
#include <fstream>
#include <iostream>
#include <map>

#include "memleak.h"

// glew
#include "glew.h"

// glm
#include "../glm/glm.hpp"
#include "../glm/gtc/matrix_transform.hpp"

#include "camera.h"
#include "world.h"
#include "fluid.h"

#define DFLAG_VEL 0x8000
#define DFLAG_TEMP 0x4000

using namespace glm;

struct Shader
{
	unsigned int vertex;
	unsigned int fragment;
	unsigned int geometry;
	unsigned int program;
};

class Display
{
public:
	Display();
	~Display();
	char* textFileRead(const char* filename);
	void printLinkInfoLog(int prog);
	void printShaderInfoLog(int shader);
	void setFluids(Fluid *fluid); 

	void init();
	void draw();
	void march();
	void updateCamera();
	void updateViewport( int w, int h ) { camera->setViewport( w, h ); updateCamera(); }
	void zoomCamera( int dz ) { camera->zoom( dz ); updateCamera(); }
	void orbitCamera( int dx, int dy ) { camera->orbit( dx, dy ); ; updateCamera(); }
	void setFlags( int f ) { flags = f; }

private:
	std::map<int, vec3> colorMap;

	int flags;

	void initShaders();
	void loadShader( const char* vertFile, const char* fragFile, const char* geomFile, Shader & shader );

	World * world;
	Camera * camera;
	Fluid * theFluid; 

	vec3 lightPos;
	vec3 lightCol;

	//frame buffer, textures
	unsigned int fbo;
	unsigned int depthTexture;
	unsigned int fluidTexture;

	//shader stuff
	unsigned int shaderProgram;
	unsigned int shadowShaderProgram;
	unsigned int raymarchShaderProgram;

	//attributes
	unsigned int shadowPositionLocation;

	unsigned int positionLocation;
	unsigned int colorLocation;
	unsigned int normalLocation;

	//uniforms
	unsigned int u_modelMatrixLocation;
	unsigned int u_projMatrixLocation;
	unsigned int u_lightPositionLocation;
	unsigned int u_lightColorLocation;
	unsigned int u_camPositionLocation;

	unsigned int u_shadowModelMatrixLocation;
	unsigned int u_shadowProjMatrixLocation;
	unsigned int u_shadowMapLocation;
	unsigned int u_shadowBiasMatrixLocation;

	unsigned int u_resolutionLocation;
	unsigned int u_rayCamPositionLocation;
	unsigned int u_distanceMapLocation;
};

#endif