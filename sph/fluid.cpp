#include "fluid.h"


//Globals
const double h = 0.05;
const double containerSizeX = 10.0; 
const double containerSizeY = 10.0; 
const double containerSizeZ = 10.0; 
const double k = 1.5; //TODO - make this function dependant on the temp

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
	Particle p(15, 0.002, glm::vec3(0, 2.0, 0), glm::vec3(0.1, 0.2, 0.3));
	theParticles.push_back(p);
}

//**********************************************************************************************
//The main fluid simulation 
//**********************************************************************************************

//Calls all the SPH fns
void Fluid::Update(float dt, glm::vec3& externalForces)
{
	if (theParticles.size() < 1000)
		addFluid(dt);
	findNeighbors();
	computeDensity(dt);
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
		theParticles.at(i).clearNeighbors(); 
		for (int j = 0 ; j < theParticles.size(); j++)
		{
			//Compute the distance from particle i to j & add it if its within the kernal radius
			if (glm::distance(theParticles.at(i).getPosition(), theParticles.at(j).getPosition()) <= h ) {
				theParticles.at(i).addNeighbor(&(theParticles.at(j))); 
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
		std::vector<Particle*> neighbors = theParticles.at(i).getNeighbors();
		for (int j = 0; j < neighbors.size(); j++)
		{
			float r = glm::distance(theParticles.at(i).getPosition(), neighbors.at(j)->getPosition()); 
			density += (neighbors.at(j)->getMass() * wPoly6(r, h)); 
		}
		theParticles.at(i).setDensity(density); 
	}
}

//Computes the current pressure of the particles
//TODO - check.  The paper doesn't store this, but uses it directly in the forces (noted in particle.h as well)
glm::vec3 Fluid::computePressure(float dt, int i )
{
	//Get the pressure: p = k * (currDens - restDens)
	glm::vec3 pressure(0.0);
	float pi = (k * (theParticles.at(i).getDensity() - theParticles.at(i).getRestDensity())); 
	std::vector<Particle*> neighbors = theParticles.at(i).getNeighbors();
	for (int j = 0; j < neighbors.size(); j++) {
		float pj = (k * (neighbors.at(j)->getDensity() - neighbors.at(j)->getRestDensity())); 
		vec3 r = theParticles.at(i).getPosition() - neighbors.at(j)->getPosition(); 
		//fPressure = - sum (mj (tempPi + tempPj) / 2 pj * gradient(W(ri - rj, h))
		pressure += (neighbors.at(j)->getMass() * ((pi + pj) / (2.0f * neighbors.at(j)->getDensity())) * wSpikyGrad(r, h)); //TODO - use gradient of spiky
	}
	
	return -pressure;
}

glm::vec3 Fluid::computeViscosity(float dt, int i)
{
	glm::vec3 v(0.0); 
	std::vector<Particle*> neighbors = theParticles.at(i).getNeighbors();
	for (int j = 0; j < neighbors.size(); j++)
	{
		float c = theParticles.at(i).getMass() / theParticles.at(i).getDensity();
		glm::vec3 r = theParticles.at(i).getPosition() - neighbors.at(j)->getPosition();
		glm::vec3 vel = neighbors.at(j)->getVelocity() - theParticles.at(i).getVelocity();
		v += neighbors.at(j)->getMass()/neighbors.at(j)->getDensity()*vel*wViscosityGrad( r, h );
		//v += (theParticles.at(i).getViscosity() + theParticles.at(j).getViscosity() / 2.0) *
		//	(theParticles.at(j).getMass() / theParticles.at(j).getDensity()) *
		//	(theParticles.at(j).getVelocity() - theParticles.at(i).getVelocity()) * wViscosityGrad(r, h); 
	}
	return 1.5f*v; 
}

void Fluid::computeForces(float dt, glm::vec3 externalForces)
{
	//TODO - add other forces for now, just add gravity
	for (int i = 0; i < theParticles.size(); i++) 
	{
		glm::vec3 pressureForce = computePressure(dt, i); 
		glm::vec3 viscosityForce = computeViscosity(dt, i); 
		glm::vec3 finalForce = pressureForce + viscosityForce + externalForces; 
		theParticles.at(i).setForce( finalForce ); 
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
		glm::vec3 pos = theParticles.at(i).getPosition();
		glm::vec3 vel = theParticles.at(i).getVelocity();
		bool updated = false;

		if (pos.x < -3) {
			updated = true;
			pos.x = -3;
			vel.x *= -1;
		}
		if (pos.y < 0) {
			updated = true;
			pos.y = 0;
			vel.y *= -1;  
		}
		if (pos.z < -3) {
			updated = true;
			pos.z = -3;
			vel.z *= -1; 
		}

		if (pos.x > 3) {
			updated = true;
			pos.x = 3;
			vel.x *= -1;
		}
		if (pos.y > 6) {
			updated = true;
			pos.y = 6;
			vel.y *= -1;  
		}
		if (pos.z > 3) {
			updated = true;
			pos.z = 3;
			vel.z *= -1; 
		}

		if (updated)
		{
			theParticles.at(i).setPostion(pos); 
			theParticles.at(i).setVelocity(vel);
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

glm::vec3 Fluid::wPoly6Grad(glm::vec3 r, float h)
{
	float lr = length( r );
	if (0 <= lr && lr <= h) {
		float lrs = lr*lr-h*h;
		float c = 315.0 / (64.0 * PI * pow(h, 9)); 
		float x = -6*r.x*lrs;
		float y = -6*r.y*lrs;
		float z = -6*r.z*lrs;
		vec3 w(x,y,z); 
		return c * w; 
	} else {
		return vec3( 0.0 );
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
vec3 Fluid::wSpikyGrad(vec3 r, float h)
{
	float lr = length( r );
	if (0 <= lr && lr <= h) {
		float c = 15.0 / (PI * pow(h, 6)); 
		float numerator = pow(lr - h, 2);
		float x = -3*r.x*numerator/(lr+1e-15);
		float y = -3*r.y*numerator/(lr+1e-15);
		float z = -3*r.z*numerator/(lr+1e-15);
		vec3 w = vec3( x, y, z );
		return c*w;
	} else {
		return vec3(0); 
	}
}

float Fluid::wViscosity(float r, float h)
{
	return 0.0f;
}

glm::vec3 Fluid::wViscosityGrad(glm::vec3 r, float h)
{
	return vec3(0.0f);
}

float Fluid::wViscosityLap(glm::vec3 r, float h)
{
	float lr = length( r );
	if (0 <= lr && lr <= h) {
		float c = 45.0f / (PI * pow(h, 6)); 
		float w = h - lr;
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

