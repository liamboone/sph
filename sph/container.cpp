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
	
	if( x >= 0 && x < 1.0 )
	{
		i = floor( width * x );
	}
	else
		return NULL;

	if( y >= 0 && y < 1.0 )
	{
		j = floor( height * y );
	}
	else
		return NULL;

	if( z >= 0 && z < 1.0 )
	{
		k = floor( depth * z );
	}
	else
		return NULL;

	return &grid[i + j*width + k*width*height];
}