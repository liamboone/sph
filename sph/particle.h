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
	void setMass(float m); 

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

	void setTemp(float temp);
	float getTemp();

	void setCi(float _ci); 
	float getCi();

	void setCp(float _cp); 
	float getCp();

    void setCs(float _cs); 
	float getCs();

	void setK(float _k);
	float getK(); 

	void setColor (glm::vec3 c);
	glm::vec3 getColor(); 

	void setCsGrad(glm::vec3 cs);
	glm::vec3 getCsGrad(); 

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

	std::vector<Particle*> neighbors; 

	//Variables for fluid-fluid interaction
	float mu; 
	float t; 
	float ci;
	float cp;
	float cs;
	glm::vec3 csGrad; 
	float k; 

};

