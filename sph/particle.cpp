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
	mu = 50.f; 
	/*const float sigma = 5.f;*/
	k = 40.0f; //TODO - do we make this dependent on t?
	t = 10; 

	//Creating 0 cs gradient to start
	csGrad = glm::vec3(0.0); 
	cpGrad = glm::vec3(0.0); 
	cs = 1; 
	cp = 1;
	ci = 1; 
}

Particle::~Particle(void)
{
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