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
	mu = 15.f; 
	/*const float sigma = 5.f;*/
	k = 5.0f; //TODO - do we make this dependent on t?
	t = 50; 
	//Creating 0 cs gradient to start
	csGrad = glm::vec3(0.0); 
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
	mu = v; 
}

float Particle::getViscosity()
{
	return mu; 
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

void Particle::setMass(float m){
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

void  Particle::setTemp(float temp)
{
	t = temp; 
}
float  Particle::getTemp()
{
	return t; 
}

void  Particle::setCi(float _ci)
{
	ci = _ci; 
}
float  Particle::getCi()
{
	return ci; 
}

void  Particle::setCp(float _cp)
{
	cp = _cp; 
}
float  Particle::getCp()
{
	return cp;
}

void  Particle::setCs(float _cs)
{
	cs = _cs; 
}
float  Particle::getCs()
{
	return cs;
}

void  Particle::setK(float _k)
{
	k = _k;
}

float  Particle::getK()
{
	return k;
}

void Particle::setColor (glm::vec3 c)
{
	color = c; 
}

glm::vec3 Particle::getColor()
{
	return color;
}

void Particle::setCsGrad(glm::vec3 cs)
{
	csGrad = cs; 
}

glm::vec3 Particle::getCsGrad()
{
	return csGrad; 
}