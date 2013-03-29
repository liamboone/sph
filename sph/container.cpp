#include "container.h"

Container::Container( int x, int y, int z, vec3 lower, vec3 upper )
{
	lBound = lower;
	uBound = upper;
	width = x;
	height = y;
	depth = z;
}

Container::~Container() { }

Box * Container::operator()( vec3 p )
{
	vec3 P = (p - lBound)/(uBound - lBound);

	float x = P.x;
	float y = P.y;
	float z = P.z;

	int i,j,k;
	
	i = floor( width * x );
	j = floor( height * y );
	k = floor( depth * z );

	return &grid[i + j*width + k*width*height];
}