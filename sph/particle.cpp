#include "particle.h"


Particle::Particle(void)
{
}

Particle::Particle(float den, float m, glm::vec3 pos, glm::vec3 vel)
{
	restDensity = den; 
	density = den; 
	mass = m; 
	position.x = pos.x; position.y = pos.y; position.z = pos.z; 
	velocity.x = vel.x; velocity.y = vel.y; velocity.z = vel.z; 
	viscosity = 1.13; 
}

Particle::~Particle(void)
{
}


glm::vec3 Particle::getPosition()
{
	return position; 
}
void Particle::setPostion(glm::vec3 pos)
{
	position = pos; 
}

glm::vec3 Particle::getVelocity()
{
	return velocity; 
}

void Particle::setVelocity(glm::vec3 vel)
{
	velocity = vel; 
}

void Particle::setViscosity(float v)
{
	viscosity = v; 
}

float Particle::getViscosity()
{
	return viscosity; 
}

std::vector<Particle*> Particle::getNeighbors()
{
	return neighbors; 
}

void Particle::setNeighbors(std::vector<Particle*> n)
{
	neighbors = n; 
}

void Particle::addNeighbor(Particle *p)
{
	neighbors.push_back(p); 
}

void Particle::addNeighbors( std::vector<Particle *> vp)
{
	neighbors.insert(neighbors.begin(), vp.begin(), vp.end());
}

void Particle::clearNeighbors()
{
	neighbors.clear(); 
}

float Particle::getMass()
{
	return mass; 
}

void Particle::setMass(int m){
	mass = m; 
}

glm::vec3 Particle::getForce()
{
	return force; 
}

void Particle::setForce(glm::vec3 f)
{
	force = f; 
}

void Particle::setDensity(float d)
{
	density = d; 
}

float Particle::getDensity()
{
	return density; 
}

void Particle::setRestDensity(float d)
{
	restDensity = d; 
}

float Particle::getRestDensity()
{
	return restDensity;
}

void Particle::setPressure(float p)
{
	pressure = p;
}

float Particle::getPressure()
{
	return pressure; 
}