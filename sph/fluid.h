#pragma once
#include "..\glm\glm.hpp";
#include "..\glm/gtc/matrix_transform.hpp"
#include <vector>

#include "particle.h"
//using namespace glm;


class Fluid
{
public:
	Fluid(void);
	~Fluid(void);

	const std::vector<Particle> &getParticles();
	//************************************************************************************************
	//Initialize fluid from some point
	void addFluid(float dt);

	//Simulation functions
	virtual void Update(float dt, glm::vec3& externalForces);
	void findNeighbors();
	void computeDensity(float dt);
	void computePressure(float dt);
	void computeForces(float dt, glm::vec3 externalForces);
	void integrate(float dt); 
	//TODO - remove if we don't need later
	void computeVelocity(float dt);
	void computePosition(float dt);
	void resolveCollisions(); //TODO - for now, just pushing inside, later add in fancier collisions
	//Kernal functions
	float wPoly6(float r, float h); 
	float wViscosity(float r, float h);
	float wSpiky(float r, float h); 

    // Reset to the initial state
    virtual void Reset();

	//************************************************************************************************
	//Grid functions 
    //Set/Get our screen resolution
    virtual void SetGridSize(int cols, int rows, int stacks);
    virtual int GetGridCols() const;
    virtual int GetGridRows() const;
    virtual int GetGridStacks() const;


    // Set/Get our draw flags -> TODO - make this handle multiple types of drawing? 
    virtual void SetDrawFlags(unsigned int flags);
    virtual unsigned int GetDrawFlags() const;
	//Draws the current frame 
    virtual void Draw(const glm::vec3& eyePos);

protected:
	//Info to draw to the screen
	int numRows, numCols, numStacks; 
	unsigned int drawFlags;
	
	std::vector<Particle> theParticles; 
	
};

