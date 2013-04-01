#pragma once
#include "..\glm\glm.hpp";
#include "..\glm/gtc/matrix_transform.hpp"
#include <vector>
#include <stdlib.h>
#include <time.h>

#include "memleak.h"

#include "particle.h"
#include "container.h"

#include "..\ObjCore\obj.h"
#include "..\ObjCore\objloader.h"

using namespace glm;

class Fluid
{
public:
	Fluid(void);
	~Fluid(void);

	const std::vector<Particle*> &getParticles();
	//************************************************************************************************
	//Initialize fluid from some point
	void addFluid(float dt);

	//Simulation functions
	void Update(float dt, glm::vec3& externalForces);
	void findNeighbors();
	void computeDensity(float dt);
	
	//Forces 
	glm::vec3 computePressure(float dt, int i);
	glm::vec3 computeViscosity(float dt, int i); 
	glm::vec3 computeSurfaceTension(float dt, int i); 
	void computeForces(float dt, glm::vec3 externalForces);

	//Position & velocity integration
	void integrate(float dt); 

	void resolveCollisions(); //TODO - for now, just pushing inside, later add in fancier collisions
	
	//Kernal functions
	float wPoly6(float r, float h); 
	glm::vec3 wPoly6Grad(glm::vec3 r, float h); 
	float wPoly6Lap(glm::vec3 r, float h); 

	float wViscosity(float r, float h);
	glm::vec3 wViscosityGrad(glm::vec3 r, float h); 
	float wViscosityLap(glm::vec3 r, float h); 

	float wSpiky(float r, float h); 
	glm::vec3 wSpikyGrad(glm::vec3 r, float h);
	float wSpikyLap(glm::vec3 r, float h); 
    
	// Reset to the initial state
    virtual void Reset();

	//Creates the particles from an imported mesh obj
	bool insideOutside(vec3 p);
	void createParticlesFromMesh(); 
	float rayTriangleIntersect(vec3 const& P0, vec3 const& V0, vec3 const& p1, vec3 const& p2, vec3 const& p3); 
	bool rayCubeIntersect(vec3 const& P0);

	//Draws the current frame 
    virtual void Draw(const glm::vec3& eyePos);

	int frame;
	Container container;

protected:
	//If we create the particles from a mesh
	obj* theMesh; 

	//Info to draw to the screen
	int numRows, numCols, numStacks; 
	unsigned int drawFlags;
	
	std::vector<Particle *> theParticles; 
	std::vector<glm::vec3> accel0; 
};

