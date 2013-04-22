#include "fluid.h"
#include <iostream>
#include <algorithm>
#include <cmath>

//Globals
const float radius = 0.0451f;
const float h = 3*radius;
const float PI = 3.14159265f; 
const float sigma = 2.f;
const float sigmaI = 2.f;
const float restDensity = 1000.0f;
const float mass = restDensity*4.0f/3.0f*PI*radius*radius*radius;

/*
In the results section, they say: 
support radius = 0.045 m
viscosity constant = 50 Ns / m^2
k = 20 Nm / kg
omegai & omegas = 0.6
*/

//***************************************
const bool loadFromMesh = !true; 

Fluid::Fluid(vec3 cMin, vec3 cMax) : container( h, cMin, cMax)
{
	containerMin = cMin;
	containerMax = cMax;
	frame = 0;
	/*srand (time(NULL));*/

	if (loadFromMesh == true)
	{
		theMesh = new obj();
		string file = "..\\stanford_bunny\\bunny_scaled4.obj";  
		//string file = "..\\bunny.obj";
		objLoader* loader = new objLoader( file, theMesh );
		theMesh->buildVBOs();
		delete loader;

		cout<<"The mesh size: "<<(*theMesh->getFaces()).size()<<endl; 
	}
}


Fluid::~Fluid(void)
{
	if ( loadFromMesh == true) 
		delete theMesh;
}

//Draws the current frame 
void Fluid::Draw(const glm::vec3& eyePos)
{

}

    // Reset to the initial state
void Fluid::Reset()
{

}

/*
P0 = the point. 
V0 = the ray direction
p1, p2, p3 are the points on the triangle
returns NaN for no intersection
*/
float Fluid::rayTriangleIntersect(vec3 const& P0, vec3 const& V0, vec3 const& p1, vec3 const& p2, vec3 const& p3)
{
	vec3 rayDirection = V0;
	vec3 rayOrigin = P0;

	/*
	point(u,v) = (1-u-v)*p0 + u*p1 + v*p2
	where u, v >= 0
	intersect ray: rayOrigin + t * rayDirection = (1 - u - v)*p0  + u*p1 + v*p2
	*/
	float f, u, v; //cooefficient & barycentric coords

	//Calculate the triangle's edges
	vec3 e1 = p2 - p1;
	vec3 e2 = p3 - p1;

	//N cross X = d
	//using normals from triangles & ray dir, calc cooefficient
	vec3 h = cross(rayDirection, e2);
	f = 1 / dot(e1, h);

	//Calc barycentric coordinate for u
	vec3 s = rayOrigin - p1;
	u = f * (dot(s, h));

	if (u < 0.0 || u > 1.0)
		return std::numeric_limits<float>::quiet_NaN();

	//Calc the barycentric coord for v
	vec3 q = cross(s, e1);
	v = f * dot(rayDirection, q);

	if (v < 0.0 || u + v > 1.0) {
		return std::numeric_limits<float>::quiet_NaN();
	}

	//Ray intersection
	//float dotTest = dot(e2, q);
	float t = f * dot(e2, q);
	return t;
}


bool Fluid::footprintIntersect(vec3 const& P0) {

	float xMax = theMesh->xmax;
	float xMin = theMesh->xmin;
	float zMax = theMesh->zmax;
	float zMin = theMesh->zmin;

	if (P0.x > xMin && P0.x < xMax)
	{
		if (P0.z > zMin && P0.z < zMax)
		{
			return true; 
		}
	}


	return false;
}

std::vector<float> Fluid::findCrossing(vec3 p)
{
	std::vector<float> crossings;
	//Check if inside the bounding box
	if( footprintIntersect( p ) )
	{
		vector<vec4> points = *theMesh->getPoints(); 
		for (unsigned int i = 0; i < (*theMesh->getFaces()).size(); i++)
		//Then check if its actually inside the mesh
		{
			vector<int> face = (*theMesh->getFaces()).at(i); 
			glm::vec3 p1(points[face[0]].x, points[face[0]].y, points[face[0]].z);
			glm::vec3 p2(points[face[1]].x, points[face[1]].y, points[face[1]].z);
			glm::vec3 p3(points[face[2]].x, points[face[2]].y, points[face[2]].z);
			float t = rayTriangleIntersect(p, vec3(0, 1, 0), p1, p2, p3);
			if( t == t ) //NaN check
			{
				crossings.push_back(t);
			}
		}
	}
	std::sort( crossings.begin(), crossings.end() );
	return crossings; 
}
//********************************************************************************************
//Particle creation 
//********************************************************************************************
//Creates particles inside of theMesh obj
void Fluid::createParticlesFromMesh()
{
	float stepsize = 0.07f;
	for( float x = containerMin.x; x < containerMax.x; x+= stepsize )
	{
		for( float z = containerMin.z; z < containerMax.z; z+=stepsize )
		{
			vec3 castPos( x, containerMin.y, z );
			std::vector<float> crossings = findCrossing(castPos);
			if( crossings.size() <= 1 ) continue;
			std::vector<float>::iterator it = crossings.begin();
			int idx = 0;
			bool noMore = false;
			for( float y = containerMin.y; y < containerMax.y && !noMore; y+=stepsize )
			{
				if( *it > (y-containerMin.y) ) continue;
				vec3 pos(x, y, z); 
				while( *it < (y-containerMin.y) )
				{
					++it;
					idx++;
					if( it == crossings.end() )
					{
						noMore = true;
						break;
					}
				}
				--it;
				idx--;
				if( !noMore && idx % 2 == 0 )
				{
					//Make a little random
					float rx = 0.001f*((float) rand() / (RAND_MAX)); 
					float ry = 0.001f*((float) rand() / (RAND_MAX)); 
					float rz = 0.001f*((float) rand() / (RAND_MAX)); 
					pos.x += rx;
					pos.y += ry;
					pos.z += rz;
					Particle * p = new Particle(restDensity, mass, pos, glm::vec3(0.0, 0.0, 0.0));
					p->setIndex(1); 
					theParticles.push_back(p);
					Box * box = container( pos );
					if( box->frame < frame )
					{
						box->particles.clear();
						box->frame = frame;
					}
					box->particles.push_back( p );
				}
			}
		}
	}
}

//Add two different types of fluids
void Fluid::addMultiFluid()
{
	if( frame == 0 )
	{
		//Normal fluid as floor
		for( float z = containerMin.z; z < containerMax.z; z+=0.07f )
		{
			for( float y = 0; y < 0.5; y += 0.07f )
			{
				for( float x = containerMin.x; x < containerMax.x; x+=0.07f )
				{
					float rx = 0.01f*((float) rand() / (RAND_MAX)); 
					float ry = 0.01f*((float) rand() / (RAND_MAX)); 
					float rz = 0.01f*((float) rand() / (RAND_MAX)); 

					vec3 pos(x+rx, y+ry, z+rz);
					Particle * p = new Particle(restDensity, mass, pos, glm::vec3(0));
					p->setIndex(1); 
					theParticles.push_back(p);
					Box * box = container( pos );
					if( box->frame < frame )
					{
						box->particles.clear();
						box->frame = frame;
					}
					box->particles.push_back( p );
				}
			}
		}
	}
	//Heavier fluid pouring from top 
	for (float z = -0.1f; z <= 0.12f; z+=0.1f)
	{
		for (float x = -0.1f; x <= 0.12f; x+=0.1f)
		{
			float rx = 0.01f*((float) rand() / (RAND_MAX)); 
			float ry = 0.01f*((float) rand() / (RAND_MAX)); 
			float rz = 0.01f*((float) rand() / (RAND_MAX)); 

			vec3 pos(x+rx-0.5, containerMax.y+ry, z+rz-0.5);
			Particle * p = new Particle(restDensity+200, mass, pos, glm::vec3(0,-1,0));
			p->setIndex(2); 
			theParticles.push_back(p);
			Box * box = container( pos );
			if( box->frame < frame )
			{
				box->particles.clear();
				box->frame = frame;
			}
			box->particles.push_back( p );
		}
	}
	//Lighter fluid pouring from top 
	for (float z = -0.1f; z <= 0.2f; z+=0.12f)
	{
		for (float x = -0.1f; x <= 0.2f; x+=0.12f)
		{
			float rx = 0.01f*((float) rand() / (RAND_MAX)); 
			float ry = 0.01f*((float) rand() / (RAND_MAX)); 
			float rz = 0.01f*((float) rand() / (RAND_MAX)); 

			vec3 pos(x+rx+0.5, containerMax.y+ry, z+rz+0.5); 
			Particle * p = new Particle(restDensity-200, mass, pos, glm::vec3(0,-1,0));
			p->setIndex(4); 
			theParticles.push_back(p);
			Box * box = container( pos );
			if( box->frame < frame )
			{
				box->particles.clear();
				box->frame = frame;
			}
			box->particles.push_back( p );
		}
	}
}

//Adds the lava lamp from the paper & affects temperature
void Fluid::addLavaLamp()
{
	if (frame == 0) {
		//Add red fluid (this will start at the bottom & rise) 
		for( float y = 0.0f; y < 0.55f; y += 0.07f )
		{
			for( float x = containerMin.x; x < containerMax.x; x += 0.07f )
			{
				for( float z = containerMin.z; z < containerMax.z; z += 0.07f )
				{
					float rx = 0.01f*(float)rand() / RAND_MAX; 
					float ry = 0.01f*(float)rand() / RAND_MAX; 
					float rz = 0.01f*(float)rand() / RAND_MAX; 
					vec3 pos(x+rx, y+ry, z+rz);
					//TODO - mass freaks out when its 0.012 like they give in the paper!!! :(
					//Particle * p = new Particle(restDensity-200, mass, pos, glm::vec3(0));
					Particle* p = new Particle(375, mass-0.05f, pos, glm::vec3(0)); 
					theParticles.push_back(p);
					p->setIndex(2);
					p->setTemp(5.f); 
					p->setCi(-0.5); 
					p->setViscosity(50.f);
					Box * box = container( pos );
					if( box->frame < frame )
					{
						box->particles.clear();
						box->frame = frame;
					}
					box->particles.push_back( p );
				}
			}
		}
		//Add blue fluid at the top (it will sink)
		for( float y = 0.556f; y < 1.4f; y += 0.07f )
		{
			for( float x = containerMin.x + 0.15f; x <containerMax.x-0.15f; x += 0.07f )
			{
				for( float z = containerMin.z + 0.15f; z < containerMax.z-0.15f; z += 0.07f )
				{
					float rx = 0.01f*(float)rand() / RAND_MAX;
					float ry = 0.01f*(float)rand() / RAND_MAX;
					float rz = 0.01f*(float)rand() / RAND_MAX;
					vec3 pos(x+rx, y+ry, z+rz);
					//was 0.006
					//TODO - figure out if we just wanna use the defaults or if we want to change other params
					/*Particle * p = new Particle(restDensity+200, mass, pos, glm::vec3(0));*/
					Particle * p = new Particle(1200, mass+0.05f, pos, glm::vec3(0));
					theParticles.push_back(p);
					p->setIndex(1);
					p->setTemp(15.f); 
					p->setViscosity(50.f); 
					Box * box = container( pos );
					if( box->frame < frame )
					{
						box->particles.clear();
						box->frame = frame;
					}
					box->particles.push_back( p );
				}
			}
		}
	}
}

//Initialize fluid from some point
void Fluid::addFluid(float dt)
{
	//Density, mass, position, velocity (particle inputs)
	for (float z = -0.2f; z < 0.2f; z+=0.1f)
    {
		for (float y = -0.2f; y < 0.2f; y += 0.1f )
		{
			float rx = 0.01f*((float) rand() / (RAND_MAX)); 
			float ry = 0.01f*((float) rand() / (RAND_MAX)); 
			float rz = 0.01f*((float) rand() / (RAND_MAX)); 
			
			vec3 pos(rx+containerMax.x, y+ry+containerMax.y/2, z+rz-0.2);
			//Density, mass, position, velocity (particle inputs)
			Particle * p = new Particle(restDensity+300, mass, pos, glm::vec3(-3, 0, 0.0));
			p->setIndex(1); 
			theParticles.push_back(p);
			Box * box = container( pos );
			if( box->frame < frame )
			{
				box->particles.clear();
				box->frame = frame;
			}
			box->particles.push_back( p );
		}
	}
}

void Fluid::addHeavyFluid(float dt)
{
	//Density, mass, position, velocity (particle inputs)
	for (float z = -0.2f; z < 0.2f; z+=0.1f)
    {
		for (float y = -0.2f; y < 0.2f; y += 0.1f )
		{
			float rx = 0.01f*((float) rand() / (RAND_MAX)); 
			float ry = 0.01f*((float) rand() / (RAND_MAX)); 
			float rz = 0.01f*((float) rand() / (RAND_MAX)); 
			
			vec3 pos(rx-containerMax.x, y+ry+containerMax.y/2, z+rz+0.2);
			//Density, mass, position, velocity (particle inputs)
			Particle * p = new Particle(restDensity-300, mass, pos, glm::vec3(3, 0, 0.0));
			p->setIndex(2); 
			theParticles.push_back(p);
			Box * box = container( pos );
			if( box->frame < frame )
			{
				box->particles.clear();
				box->frame = frame;
			}
			box->particles.push_back( p );
		}
	}
}

//**********************************************************************************************
//The main fluid simulation 
//**********************************************************************************************

//Calls all the SPH fns
void Fluid::Update(float dt, force_t externalForce)
{
    //addLavaLamp(); 
	if (loadFromMesh == true && frame == 0) {
		createParticlesFromMesh();
	} else if (loadFromMesh == false && theParticles.size() < 3000 && frame % 4 == 0) {
		addFluid(dt);
		addHeavyFluid(dt);
		//addMultiFluid();
		if( frame == -1 )
		{
			for( float y = 0; y < 0.4; y += 0.07f )
			{
				for( float x = containerMin.x; x < containerMax.x; x += 0.07f )
				{
					for( float z = containerMin.z; z < containerMax.z; z += 0.07f )
					{
						float rx = 0.01f*(float)rand() / RAND_MAX; 
						float ry = 0.01f*(float)rand() / RAND_MAX; 
						float rz = 0.01f*(float)rand() / RAND_MAX; 
						vec3 pos(x+rx, y+ry, z+rz);
						Particle * p = new Particle(restDensity-200, mass, pos, glm::vec3(0));
						p->setIndex(2); 
						p->setCi(-0.5); 
						theParticles.push_back(p);
						Box * box = container( pos );
						if( box->frame < frame )
						{
							box->particles.clear();
							box->frame = frame;
						}
						box->particles.push_back( p );
					}
				}
			}
		}
	}
	findNeighbors();

	computeDensity(dt);
	computeForces(dt, externalForce);
	integrate(dt); 
	resolveCollisions();
	frame++;
	
	//container.clear();
	for (unsigned int i = 0; i < theParticles.size(); i++)
	{
		Particle * p = theParticles.at(i);
		Box * box = container( p->getPosition() );
		if( box == NULL )
			continue;
		if( box->frame < frame )
		{
			box->particles.clear();
			box->frame = frame;
		}
		box->particles.push_back( p );
	}
}

//Finds all the the particles current neighbors and stores them in the particle's neighbors vector
void Fluid::findNeighbors()
{
	#pragma omp parallel for
	for (int i = 0; i < (int) theParticles.size(); i++)
	{
		Particle * p = theParticles.at(i);
		p->clearNeighbors();
		//*
		vec3 offset( 0 );
		Box * visited[27];
		int j = 0;
		for ( int x = -1; x <= 1; x ++ )
		{
			offset.x = x*h;
			for ( int y = -1; y <= 1; y ++ )
			{
				offset.y = y*h;
				for ( int z = -1; z <= 1; z ++ )
				{
					offset.z = z*h;
					Box * box = container( p->getPosition() + offset );
					visited[j] = box;
					if( box != NULL )
					{
						if( box->frame == frame )
						{
							bool flag = true;
							for( int k = 0; k < j; k ++ )
							{
								flag = flag && ( visited[k] != box );
							}
							if( flag )
							{
								p->addNeighbors( box->particles );
							}
						}
					}
					j++;
				}
			}
		}
	}
}

//Finds the density at the current timestep
void Fluid::computeDensity(float dt)
{
	#pragma omp parallel for
	for (int i = 0; i < (int)theParticles.size(); i++)
	{
		float density = 0;
		vec3 pos_i = theParticles.at(i)->getPosition();
		std::vector<Particle*> neighbors = theParticles.at(i)->getNeighbors();
		
		for (unsigned int j = 0; j < neighbors.size(); j++)
		{
			float r = glm::distance(pos_i, neighbors.at(j)->getPosition()); 
			density += (neighbors.at(j)->getMass() * wPoly6(r, h));
		}
		theParticles.at(i)->setDensity(density); 

	}
} 

//Computes the current pressure of the particles
glm::vec3 Fluid::computePressure(float dt, int i )
{
	//Get the pressure: p = k * (currDens - restDens)
	glm::vec3 pressure(0.0);
	float pi = (theParticles.at(i)->getK() * (theParticles.at(i)->getDensity() - theParticles.at(i)->getRestDensity())); 
	std::vector<Particle*> neighbors = theParticles.at(i)->getNeighbors();
	for (unsigned int j = 0; j < neighbors.size(); j++) {
		float pj = (neighbors.at(j)->getK() * (neighbors.at(j)->getDensity() - neighbors.at(j)->getRestDensity())); 
		vec3 r = theParticles.at(i)->getPosition() - neighbors.at(j)->getPosition(); 
		//fPressure = - sum (mj (tempPi + tempPj) / 2 pj * gradient(W(ri - rj, h))
		pressure += neighbors.at(j)->getMass() * ((pi + pj) / (1e-15f + 2.0f * neighbors.at(j)->getDensity())) * wSpikyGrad(r, h); 
	}


	float sentinel = glm::length( pressure );
	return -pressure;
}

glm::vec3 Fluid::computeViscosity(float dt, int i)
{
	glm::vec3 v(0.0); 
	std::vector<Particle*> neighbors = theParticles.at(i)->getNeighbors();
	for (unsigned int j = 0; j < neighbors.size(); j++)
	{
		glm::vec3 r = theParticles.at(i)->getPosition() - neighbors.at(j)->getPosition();
		glm::vec3 vel = (neighbors.at(j)->getVelocity() - theParticles.at(i)->getVelocity()) / (neighbors.at(j)->getDensity() + 1e-15f);
		float muAvg = (neighbors.at(j)->getViscosity() + theParticles.at(i)->getViscosity()) / 2.0f; 
		v += muAvg * neighbors.at(j)->getMass() * vel * wViscosityLap( r, h );
	}
	return v; 
}

glm::vec3 Fluid::computeSurfaceTension(float dt, int i)
{
	glm::vec3 n(0.0); 
    float k = 0.0; 
	std::vector<Particle*> neighbors = theParticles.at(i)->getNeighbors();
	for (unsigned int j = 0; j < neighbors.size(); j++)
	{
		glm::vec3 r = theParticles.at(i)->getPosition() - neighbors.at(j)->getPosition();

		float mass = neighbors.at(j)->getMass(); 
		float denInv = 1.0f / (neighbors.at(j)->getDensity() + 1e-15f); 
		n += mass * denInv * wPoly6Grad(r, h); 
		k += mass * denInv * wPoly6Lap(r, h); 
	}

	k = -k / (glm::length(n) + 1e-15f); 
	n =  n / (glm::length(n) + 1e-15f);
	return sigma * k * n; 
}

void Fluid::computeSurfaceAndInterfaceTension(int i, glm::vec3 &interfaceForce, glm::vec3 &surfaceForce)
{
	//Interface
	glm::vec3 ni(0.0);
	float ki = 0.0f; 
	//Surface
	glm::vec3 n(0.0); 
    float k = 0.0f; 
	//Cp 
	glm::vec3 nCp(0.0); 

	std::vector<Particle*> neighbors = theParticles.at(i)->getNeighbors();
	for (unsigned int j = 0; j < neighbors.size(); j++)
	{
		glm::vec3 r = theParticles.at(i)->getPosition() - neighbors.at(j)->getPosition();
		float mass = neighbors.at(j)->getMass(); 
		//Cs should be 0 for air particles & 1 for liquids
		float denInv = theParticles.at(i)->getCs() / (neighbors.at(j)->getDensity() + 1e-15f); 
		glm::vec3 poly6Grad = wPoly6Grad(r, h);
		float poly6Lap = wPoly6Lap(r, h); 

		//Surface force - TODO - set to 0 if air (its 1 now)
		n += mass * denInv *  poly6Grad;


		//Cp stored for air particles
		nCp += mass * 1.f / (neighbors.at(j)->getDensity() + 1e-15f) * poly6Grad;
		
		//Interface force
		ni += mass * neighbors.at(j)->getCi() * denInv * poly6Grad;
		ki += mass * neighbors.at(j)->getCi() * denInv * poly6Lap; 

	}

	//The normal is the gradient of cs.  Store it to use in bubble creation
	theParticles.at(i)->setCsGrad(n); 
	theParticles.at(i)->setCpGrad(nCp); 

	k = -k / (glm::length(n) + 1e-15f); 
	n =  n / (glm::length(n) + 1e-15f);
	surfaceForce = sigma * k * n; 


	//TODO - should we be normalizing k? 
	ni = ni / (glm::length(ni) + 1e-15f); 
	interfaceForce = -sigmaI * ki * ni; 
}

void Fluid::computeForces(float dt, force_t externalForce)
{
	//TODO - add other forces for now, just add gravity
	#pragma omp parallel for
	for (int i = 0; i < (int)theParticles.size(); i++) 
	{
		glm::vec3 pressureForce = computePressure(dt, i); 
		glm::vec3 viscosityForce = computeViscosity(dt, i); 
		glm::vec3 finalForce;
		glm::vec3 externalForces = externalForce(theParticles.at(i)->getPosition());
		glm::vec3 interfaceForce; 
		glm::vec3 surfaceTension;
			
		computeSurfaceAndInterfaceTension(i, interfaceForce, surfaceTension);
		finalForce = pressureForce + theParticles.at(i)->getDensity()*externalForces + viscosityForce + surfaceTension + interfaceForce; 
		theParticles.at(i)->setForce( finalForce ); 
	}
}


void Fluid::integrate(float dt)
{
	
	//Euler just in case leapfrog is wrong
	#pragma omp parallel for
	for (int i = 0; i < (int)theParticles.size(); i++)
	{
		Particle * p = theParticles.at(i); 
		p->setVelocity(p->getVelocity() + dt * p->getForce() / p->getDensity());
		p->setPostion(p->getPosition() + dt * p->getVelocity());
	}

	return;
	
	//********************************************************************
	//Leapfrog from http://image.diku.dk/projects/media/kelager.06.pdf
	std::vector<glm::vec3> accel1; 
	std::vector<vec3> vel1; 
	for (unsigned int i = 0; i < theParticles.size(); i++)
	{
		glm::vec3 a = theParticles.at(i)->getForce() / theParticles.at(i)->getDensity();
		accel1.push_back(a);
		glm::vec3 velTminusHalf;
		if (i < accel0.size()) {
			velTminusHalf = vel0.at(i) - 0.5f * accel0.at(i) * dt;
		} else {
			velTminusHalf = theParticles.at(i)->getVelocity() - 0.5f * vec3(0.0, -9.8, 0.0) * dt; 
		}
		glm::vec3 velTplusHalf = velTminusHalf + dt * a; 
		theParticles.at(i)->setPostion(theParticles.at(i)->getPosition() + velTplusHalf * dt); 
		theParticles.at(i)->setVelocity((velTminusHalf + velTplusHalf) / 2.0f); 
		vel1.push_back(velTplusHalf); 
	}
	accel0.clear();
	accel0 = accel1; 
	vel0.clear();
	vel0 = vel1;
}

void Fluid::resolveCollisions()
{
	float damping = 1.f;
	// p.velocity = p.velocity - 2.0 * Dot(normal, p.velocity) * normal;
	for (unsigned int i = 0; i < theParticles.size(); i++)
	{
		glm::vec3 pos = theParticles.at(i)->getPosition();
		glm::vec3 vel = theParticles.at(i)->getVelocity();
		bool updated = false;
		//For now, don't handle corners
		if (pos.x < containerMin.x) {
			glm::vec3 normal = glm::vec3(1, 0, 0); 
			updated = true;
			//Normal of wall is 
			pos.x = containerMin.x;
			vel.x *= -damping;
		}
		if (pos.y < containerMin.y) {
			glm::vec3 normal = glm::vec3(0, 1.0, 0); 
			updated = true;
			pos.y = containerMin.y;
			vel.y *= -damping;  
		}
		if (pos.z < containerMin.z) {
			glm::vec3 normal = glm::vec3(0, 0, 1); 
			updated = true;
			pos.z = containerMin.z;
			vel.z *= -damping; 
		}

		if (pos.x > containerMax.x) {
			glm::vec3 normal = glm::vec3(-1, 0, 0); 
			updated = true;
			pos.x = containerMax.x;
			vel.x *= -damping;
		}
		if (pos.y > containerMax.y) {
			glm::vec3 normal = glm::vec3(0, -1.0, 0); 
			updated = true;
			pos.y = containerMax.y;
			vel.y *= -damping;  
		}
		if (pos.z > containerMax.z) {
			glm::vec3 normal = glm::vec3(0, 0, -1.0); 
			updated = true;
			pos.z = containerMax.z;
			vel.z *= -damping; 
		}

		if (updated)
		{
			theParticles.at(i)->setPostion(pos); 
			theParticles.at(i)->setVelocity(vel);
		}
	}
}

inline float Fluid::wPoly6(float r, float h)
{
	if (0 <= r && r <= h) {
		float c = 315.0f / (64.0f * PI * pow(h, 9)); 
		float w = pow(pow(h, 2) - pow(r, 2), 3); 
		return c * w; 
	} else {
		return 0.0;
	}
}

inline glm::vec3 Fluid::wPoly6Grad(glm::vec3 r, float h)
{
	float lr = length( r );
	if (0 <= lr && lr <= h) {
		float lrs = lr*lr-h*h;
		lrs *= lrs;
		float c = 315.0f / (64.0f * PI * pow(h, 9)); 
		float x = -6*r.x*lrs;
		float y = -6*r.y*lrs;
		float z = -6*r.z*lrs;
		vec3 w(x,y,z); 
		return c * w; 
	} else {
		return vec3( 0.0 );
	}
}

inline float Fluid::wPoly6Lap(glm::vec3 r, float h)
{
	float c = 315.0f / (64.0f * PI * pow(h, 9)); 
	float lr = glm::length(r); 
	float lrS = lr * lr; 
	if (0 <= lr && lr <= h) {
		float w = -6 * ( 3 * h * h * h * h - 10 * h * h * lrS + 7 * lrS * lrS);
		return c * w; 
	} else {
		return 0.0f; 
	}
}

//Used for pressure calcs
inline float Fluid::wSpiky(float r, float h)
{
	if (0 <= r && r <= h) {
		float c = 15.0f / (PI * pow(h, 6)); 
		float w = pow (h - r, 3); 
		return c*w;
	} else {
		return 0; 
	}
}
inline vec3 Fluid::wSpikyGrad(vec3 r, float h)
{
	float lr = length( r );
	if (0 <= lr && lr <= h) {
		float c = 15.0f / (PI * pow(h, 6)); 
		float numerator = pow(lr - h, 2);
		float x = -3.f*r.x*numerator/(lr+1e-15f);
		float y = -3.f*r.y*numerator/(lr+1e-15f);
		float z = -3.f*r.z*numerator/(lr+1e-15f);
		vec3 w = vec3( x, y, z );
		return c*w;
	} else {
		return vec3(0); 
	}
}

inline float Fluid::wViscosity(float r, float h)
{
	return 0.0f;
}

inline glm::vec3 Fluid::wViscosityGrad(glm::vec3 r, float h)
{
	return vec3(0.0f);
}

inline float Fluid::wViscosityLap(glm::vec3 r, float h)
{
	float lr = length( r );
	if (0 <= lr && lr <= h) {
		float c = 45.0f / (PI * pow(h, 6)); 
		float w = h - lr;
		return c*w;
	} else {
		return 0; 
	}
}


const std::vector<Particle*>& Fluid::getParticles()
{
	return theParticles; 
}

vec4 Fluid::field( vec3 pos, std::map<int, vec3> colorMap )
{
	std::vector< Particle * >::iterator it;
	float potential = 0;
	vec3 color(0.0f);

	for( int x = -1; x <= 1; x ++ )
	{
		for( int y = -1; y <= 1; y ++ )
		{
			for( int z = -1; z <= 1; z ++ )
			{
				Box * box = container( pos+h*vec3(x,y,z) );
				if( box->frame < frame )
					continue;
				std::vector< Particle * > particles = box->particles;
				for( it = particles.begin(); it != particles.end(); ++it )
				{
					float rho = (*it)->getDensity(); 
					float r = glm::distance( (*it)->getPosition(), pos );
					float weight = (*it)->getMass() / (1e-15f + rho) * wPoly6(r, h);
					color += colorMap[(*it)->getIndex()]*weight;
					potential += weight; 
				}
			}
		}
	}
	return vec4( color, potential );
}