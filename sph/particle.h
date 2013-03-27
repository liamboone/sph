#pragma once
//#include "vec.h"  //TODO - use vec or glm? 
#include <vector> 

#include "memleak.h"

// glm
#include "..\glm\glm.hpp";
#include "..\glm/gtc/matrix_transform.hpp"
//using namespace glm;

class Particle
{
public:
	Particle(void);
	Particle(float den, float m); 
	Particle(float den, float mass, glm::vec3 position, glm::vec3 velocity); 
	~Particle(void);

	float getMass();
	void setMass(int m); 

	glm::vec3 getForce(); 
	void setForce(glm::vec3 f); 
	
	glm::vec3 getPosition();
	void setPostion(glm::vec3 pos);
	
	glm::vec3 getVelocity();
	void setVelocity(glm::vec3 vel); 
	
	void setDensity(float d);
	float getDensity();
	
	void setRestDensity(float d);
	float getRestDensity();
	
	void setPressure(float p);
	float getPressure(); 

	void setViscosity(float v);
	float getViscosity(); 

	std::vector<Particle*> getNeighbors();
	void setNeighbors(std::vector<Particle*> n); 
	void addNeighbor(Particle *p); 
	void addNeighbors( std::vector<Particle *> vp ); 
	void clearNeighbors(); 


private:
	float restDensity;

	float density; 
	float mass; 
	float pressure; //TODO - this is not stored in the paper

	glm::vec3 position;
	glm::vec3 velocity;

	glm::vec3 force; 
	glm::vec3 color;
	float viscosity; 
	float temp; 
	std::vector<Particle*> neighbors; 
};

