#include "fluid.h"


//Globals
const double h = 1.0;
const double containerSizeX = 10.0; 
const double containerSizeY = 10.0; 
const double containerSizeZ = 10.0; 
const double k = 10.0; //TODO - make this function dependant on the temp

const double PI = 3.14159265; 

Fluid::Fluid(void)
{

	//Sets up the size of the grid that the particles fall into 
	SetGridSize(containerSizeX, containerSizeY, containerSizeZ);
}


Fluid::~Fluid(void)
{
}




//Draws the current frame 
void Fluid::Draw(const glm::vec3& eyePos)
{

}

    // Reset to the initial state
void Fluid::Reset()
{

}

//Initialize fluid from some point
void Fluid::addFluid(float dt)
{
	Particle p(1000.0, 1.0, glm::vec3(0, 10.0, 0), glm::vec3(0.5, 0, 0));
	theParticles.push_back(p);
}

//**********************************************************************************************
//The main fluid simulation 
//**********************************************************************************************

//Calls all the SPH fns
void Fluid::Update(float dt, glm::vec3& externalForces)
{
	if (theParticles.size() < 10000)
		addFluid(dt);
	findNeighbors();
	computeDensity(dt);
	computePressure(dt);
	computeForces(dt, externalForces);
	integrate(dt); 
	//computeVelocity(dt);
	//computePosition(dt);
	resolveCollisions();
}

//Finds all the the particles current neighbors and stores them in the particle's neighbors vector
void Fluid::findNeighbors()
{
	for (int i = 0; i < theParticles.size(); i++)
	{
		for (int j = 0 ; j < theParticles.size(); j++)
		{
			theParticles.at(i).clearNeighbors(); 
			//Compute the distance from particle i to j & add it if its within the kernal radius
			if (glm::distance(theParticles.at(i).getPosition(), theParticles.at(j).getPosition()) < h && i != j) {
				theParticles.at(i).addNeighbor(theParticles.at(j)); 
			}
		}
	}
}

//Finds the density at the current timestep
void Fluid::computeDensity(float dt)
{
	for (int i = 0; i < theParticles.size(); i++)
	{
		float density = 0;
		std::vector<Particle> neighbors = theParticles.at(i).getNeighbors();
		for (int j = 0; j < neighbors.size(); j++)
		{
			float r = glm::distance(theParticles.at(i).getPosition(), neighbors.at(j).getPosition()); 
			density += (neighbors.at(j).getMass() * wPoly6(r, h)); 
		}
		theParticles.at(i).setDensity(density); 
	}
}

//Computes the current pressure of the particles
//TODO - check.  The paper doesn't store this, but uses it directly in the forces (noted in particle.h as well)
void Fluid::computePressure(float dt)
{
	//Get the pressure: p = k * (currDens - restDens)

	for (int i = 0; i < theParticles.size(); i++)
	{
		float pressure = 0;
		float pi = (k* (theParticles.at(i).getDensity() - theParticles.at(i).getRestDensity())); 
		std::vector<Particle> neighbors = theParticles.at(i).getNeighbors();
		for (int j = 0; j < neighbors.size(); j++) {
			float pj = (k* (neighbors.at(j).getDensity() - neighbors.at(j).getRestDensity())); 
			float r = glm::distance(theParticles.at(i).getPosition(), neighbors.at(j).getPosition()); 
			//fPressure = - sum (mj (tempPi + tempPj) / 2 pj * gradient(W(ri - rj, h))
			pressure += (neighbors.at(j).getMass() * ((pi + pj) / (2.0 * neighbors.at(j).getDensity())) * wSpiky(r, h)); 
		}
		theParticles.at(i).setPressure(-pressure);
	}

}

void Fluid::computeForces(float dt, glm::vec3 externalForces)
{
	//TODO - add other forces for now, just add gravity
	for (int i = 0; i < theParticles.size(); i++) 
	{
		theParticles.at(i).setForce(glm::vec3(0, -9.8, 0) * theParticles.at(i).getMass()); 
	}
}

void Fluid::integrate(float dt)
{
	//Euler just in case leapfrog is wrong
	for (int i = 0; i < theParticles.size(); i++)
	{
		Particle& p = theParticles.at(i); 
		p.setVelocity(p.getVelocity() + dt * p.getForce() / p.getMass());
		p.setPostion(p.getPosition() + dt * p.getVelocity());
	}

	//************Leap frog start **************************
	//http://en.wikipedia.org/wiki/Leapfrog_integration	
	//std::vector<vec3> currAccel;
	//for (int i = 0; i < theParticles.size(); i++)
	//{
	//		Particle p = theParticles.at(i); 
	//		glm::vec3 accel = p.getForce() / p.getMass(); 
	//		currAccel.push_back(accel);
	//		glm::vec3 newPos = p.getPosition() + p.getVelocity() * dt + accel * dt * dt * (float) 0.5;
	//		p.setPostion(newPos); 
	//}

	////TODO - figure out what to add here/ if we should use different scheme
	//computeForces(dt, vec3(0, 0, 0));

	//for (int i = 0; i < theParticles.size(); i++)
	//{
	//	Particle& p = theParticles.at(i);
	//	vec3 newAccel = p.getForce()/ p.getMass();

	//	p.setVelocity(p.getVelocity() + dt * (float) 0.5 * (currAccel.at(i) + newAccel));
	//	
	//}

}

void Fluid::computeVelocity(float dt)
{
}

void Fluid::computePosition(float dt)
{

}

void Fluid::resolveCollisions()
{
	for (int i = 0; i < theParticles.size(); i++)
	{
		if (theParticles.at(i).getPosition().y < 0) {
			theParticles.at(i).setPostion(glm::vec3(theParticles.at(i).getPosition().x, 0, theParticles.at(i).getPosition().z)); 
			theParticles.at(i).setVelocity(glm::vec3(theParticles.at(i).getVelocity().x, -theParticles.at(i).getVelocity().y, theParticles.at(i).getVelocity().z));  
		}
	}
}


float Fluid::wPoly6(float r, float h)
{
	if (0 <= r && r <= h) {
		float c = 315.0 / (64.0 * PI * pow(h, 9)); 
		float w = pow(pow(h, 2) - pow(r, 2), 3); 
		return c * w; 
	} else {
		return 0.0;
	}
}

//Used for pressure calcs
float Fluid::wSpiky(float r, float h)
{
	if (0 <= r && r <= h) {
		float c = 15.0 / (PI * pow(h, 6)); 
		float w = pow (h - r, 3); 
		return c*w;
	} else {
		return 0; 
	}
}

float Fluid::wViscosity(float r, float h)
{
	if (0 <= r && r <= h) {
		float c = 15.0 / (2.0 * PI * pow(h, 3)); 
		float w = (- pow (r, 3) / (2.0 * pow(h, 3))) + (pow(r, 2) / pow(h, 2)) + h / (2.0 *r) - 1; 
		return c*w;
	} else {
		return 0; 
	}
}

//************************************************************************************************
//Grid section 
//TODO - use this for speedup
//************************************************************************************************

//Set/Get our screen resolution
void Fluid::SetGridSize(int cols, int rows, int stacks)
{
	numRows = rows;
	numCols = cols;
	numStacks = stacks; 
}

int Fluid::GetGridCols() const
{
	return numCols;
}
int Fluid::GetGridRows() const
{
	return numRows;
}

int Fluid::GetGridStacks() const
{
	return numStacks; 
}

// Set/Get our draw flags
void Fluid::SetDrawFlags(unsigned int flags)
{
	drawFlags = flags;
}

unsigned int Fluid::GetDrawFlags() const
{
	return drawFlags;
}

const std::vector<Particle>& Fluid::getParticles()
{
	return theParticles; 
}