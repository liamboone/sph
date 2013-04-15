#pragma once
//#include "vec.h"  //TODO - use vec or glm? 
#include <vector> 

#include "memleak.h"

// glm
#include "..\glm\glm.hpp"
#include "..\glm/gtc/matrix_transform.hpp"
//using namespace glm;

class Particle
{
public:
	Particle(void);
	Particle(float den, float m); 
	Particle(float den, float mass, glm::vec3 position, glm::vec3 velocity); 
	~Particle(void);

	void setMass(float m)			{ mass = m;			} 
	void setDensity(float d)		{ density = d;		}
	void setRestDensity(float d)	{ restDensity = d;	}
	void setPressure(float p)		{ pressure = p;		}
	void setViscosity(float v)		{ mu = v;			}
	void setTemp(float temp)		{ t = temp;			}
	void setCi(float _ci)			{ ci = _ci;			}
	void setCp(float _cp)			{ cp = _cp;			}
	void setCs(float _cs)			{ cs = _cs;			}
	void setK(float _k)				{ k = _k;			}
	void setIndex(int idx)			{ index = idx;		}
	void setCsGrad(glm::vec3 cs)	{ csGrad = cs;		}
	void setCpGrad(glm::vec3 cp)	{ cpGrad = cp;		}
	void setForce(glm::vec3 f)		{ force = f;		}
	void setPostion(glm::vec3 pos)	{ position = pos;	}
	void setVelocity(glm::vec3 vel)	{ velocity = vel;	}

	float getMass()			{ return mass;			}
	float getDensity()		{ return density;		}
	float getRestDensity()	{ return restDensity;	}
	float getPressure()		{ return pressure;		}
	float getViscosity()	{ return mu;			}
	float getTemp()			{ return t;				}
	float getCi()			{ return ci;			}
	float getCp()			{ return cp;			}
	float getCs()			{ return cs;			}
	float getK()			{ return k;				}
	int getIndex()			{ return index;			}
	glm::vec3 getCsGrad()	{ return csGrad;		}
	glm::vec3 getCpGrad()	{ return cpGrad;		}
	glm::vec3 getForce()	{ return force;			}
	glm::vec3 getPosition()	{ return position;		}
	glm::vec3 getVelocity()	{ return velocity;		}

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
	int index;

	std::vector<Particle*> neighbors; 

	//Variables for fluid-fluid interaction
	float mu; 
	float t; 
	float ci;
	float cp;
	float cs;
	glm::vec3 csGrad; 
	glm::vec3 cpGrad; 
	float k;
};

