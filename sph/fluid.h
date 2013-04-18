#pragma once
#include "..\glm\glm.hpp"
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

typedef vec3 ( *force_t )( vec3 );

class Fluid
{
public:
	Fluid(vec3 cMin, vec3 cMax);
	Fluid(void);
	~Fluid(void);

	const std::vector<Particle*> &getParticles();
	//************************************************************************************************
	//Initialize fluid from some point
	void addFluid(float dt);

	//Test functions to draw two different liquids 
	void addMultiFluid(); //Pours fluids of two densities together 
	void addLavaLamp();  //Recreates the lavalamp example from http://www.matthiasmueller.info/publications/sca05.pdf

	//Simulation functions
	void Update(float dt, force_t externalForce);
	void findNeighbors();
	void computeDensity(float dt);
	vec4 field( vec3 pos, std::map<int,vec3> colorMap );
	//Functions for multiple fluid & bubble simulation
	void computeDiffusion(float dt); 
	glm::vec3 computeArtificialBuoyancy(int i); 
	void manageAirBubbles();
	

	//Forces 
	glm::vec3 computePressure(float dt, int i);
	glm::vec3 computeViscosity(float dt, int i); 
	glm::vec3 computeSurfaceTension(float dt, int i); 
	void computeSurfaceAndInterfaceTension(int i, glm::vec3 &interfaceForce, glm::vec3 &surfaceForce);
	void computeForces(float dt, force_t externalForce);

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
	std::vector<float> findCrossing(vec3 p);
	void createParticlesFromMesh(); 
	float rayTriangleIntersect(vec3 const& P0, vec3 const& V0, vec3 const& p1, vec3 const& p2, vec3 const& p3); 
	bool footprintIntersect(vec3 const& P0);

	//Draws the current frame 
    virtual void Draw(const glm::vec3& eyePos);

	int frame;
	Container container;
	vec3 containerMax;
	vec3 containerMin;

protected:
	//If we create the particles from a mesh
	obj* theMesh; 

	//Info to draw to the screen
	int numRows, numCols, numStacks; 
	unsigned int drawFlags;
	
	std::vector<Particle *> theParticles; 
	std::vector<glm::vec3> accel0; 
	std::vector<vec3> vel0; 
};

