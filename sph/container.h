#ifndef _SPH_CONTAINER_H_
#define _SPH_CONTAINER_H_

#include <vector>
#include <map>

#include "memleak.h"

#include "particle.h"

// glm
#include "..\glm\glm.hpp";

using namespace glm;

struct Box
{
	int frame;
	std::vector< Particle * > particles;
	Box() : frame(-1) { }
};

class Container
{
public:
	Container( int x, int y, int z, vec3 lower, vec3 upper );
	~Container();
	Box * operator()( vec3 p );
	void clear() { grid.clear(); }

private:
	std::map<int, Box> grid;
	vec3 lBound;
	vec3 uBound;
	int width;
	int height;
	int depth;
};

#endif