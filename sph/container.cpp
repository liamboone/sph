#include "container.h"

Container::Container( float h, vec3 lower, vec3 upper )
{
	lBound = lower;
	uBound = upper;

	span = uBound-lBound;

	width = (int)(span.x/h);
	height = (int)(span.y/h);
	depth = (int)(span.z/h);
}

Container::~Container() { }

Box * Container::operator()( vec3 p )
{
	vec3 P = (p - lBound) / span;

	float x = P.x;
	float y = P.y;
	float z = P.z;

	int i,j,k;
	
	i = (int) floor( width * x );
	j = (int) floor( height * y );
	k = (int) floor( depth * z );

	return &grid[i + j*width + k*width*height];
}