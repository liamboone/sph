#include "fluid.h"
#include <iostream>

//Globals
const float radius = 0.0451f;
const float h = 3*radius;
const float k = 20.0f; //TODO - make this function dependant on the temp
const float PI = 3.14159265f; 
const float mu = 50.f;
const float sigma = 0.6f;
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
const bool useMultiFluidCalcs = !true; 
const bool loadFromMesh = false; 

//Container size 
const vec3 containerMin(-1.0, 0, -1.0);
const vec3 containerMax(1.0, 3, 1.0); 

Fluid::Fluid(void) : container( h, containerMin, containerMax)
{
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
		return -1;

	//Calc the barycentric coord for v
	vec3 q = cross(s, e1);
	v = f * dot(rayDirection, q);

	if (v < 0.0 || u + v > 1.0) {
		return -1;
	}

	//Ray intersection
	//float dotTest = dot(e2, q);
	float t = f * dot(e2, q);
	if (t > 0.00001) {
		return t;
	} else {
		return -1;
	}
}


bool Fluid::rayCubeIntersect(vec3 const& P0) {

	float xMax = theMesh->xmax;
	float xMin = theMesh->xmin;
	float yMax = theMesh->ymax;
	float yMin = theMesh->ymin;
	float zMax = theMesh->zmax;
	float zMin = theMesh->zmin;

	if (P0.x > xMin && P0.x < xMax)
	{
		if (P0.y > yMin && P0.y < yMax)
		{
			if (P0.z > zMin && P0.z < zMax)
			{
				return true; 
			}
		}
	}


	return false;
}

bool Fluid::insideOutside(vec3 p)
{
	//Check if inside the bounding box
	if (rayCubeIntersect(p) == true)
	{
		int numIntersects = 0; 
		vector<vec4> points = *theMesh->getPoints(); 
		for (unsigned int i = 0; i < (*theMesh->getFaces()).size(); i++)
		//Then check if its actually inside the mesh
		{
			vector<int> face = (*theMesh->getFaces()).at(i); 
			glm::vec3 p1(points[face[0]].x, points[face[0]].y, points[face[0]].z);
			glm::vec3 p2(points[face[1]].x, points[face[1]].y, points[face[1]].z);
			glm::vec3 p3(points[face[2]].x, points[face[2]].y, points[face[2]].z);
			if (rayTriangleIntersect(p, vec3(0, 1, 0), p1, p2, p3) != -1) {
				numIntersects++;
			}
		}
		if (numIntersects % 2 == 0)
			return false;
		else 
			return true;
	} 
	
	return false; 
}
//********************************************************************************************
//Particle creation 
//********************************************************************************************
//Creates particles inside of theMesh obj
void Fluid::createParticlesFromMesh()
{
	//Run through the grid & add particles if we are inside the mesh
	for( float x = containerMin.x; x < containerMax.x; x+= 0.1f )
	{
		for( float y = containerMin.y; y < containerMax.y; y+=0.1f )
		{
			for( float z = containerMin.z; z < containerMax.z; z+=0.1f )
			{
				vec3 pos(x, y, z); 
				if (insideOutside(pos) == true) {
					//Make a little random
					float rx = 0.001f*((float) rand() / (RAND_MAX)); 
					float ry = 0.001f*((float) rand() / (RAND_MAX)); 
					float rz = 0.001f*((float) rand() / (RAND_MAX)); 
					pos.x += rx;
					pos.y += ry;
					pos.z += rz;
					Particle * p = new Particle(restDensity, mass, pos, glm::vec3(0.0, 0.0, 0.0));
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
		//Add more dense fluid in cube
		for( float z = -0.4; z < 0.4; z+=0.1f )
		{
			for( float y = 0; y < 1.0; y += 0.1f )
			{
				for( float x = -0.4; x < 0.4; x+=0.1f )
				{
					float rx = 0.01f*((float) rand() / (RAND_MAX)); 
					float ry = 0.01f*((float) rand() / (RAND_MAX)); 
					float rz = 0.01f*((float) rand() / (RAND_MAX)); 

					vec3 pos(x+rx, y+ry, z+rz);
					//Density, mass, position, velocity (particle inputs)
					//**************************************
					//TODO - play with #'s here!! 
					Particle * p = new Particle(600, 4, pos, glm::vec3(0, 0, 0.0));
					p->setColor(vec3(0, 0.2, 1.0)); 
					p->setViscosity(16.f); 
					p->setCi(-0.5); // This is a polar liquid (water)
					//**************************************
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
		for (float z = -0.1; z < 0.1; z+=0.1)
		{
			for (float x = -0.1; x < 0.1; x+=0.1)
			{
				float rx = 0.01*((double) rand() / (RAND_MAX)); 
				float ry = 0.01*((double) rand() / (RAND_MAX)); 
				float rz = 0.01*((double) rand() / (RAND_MAX)); 

				vec3 pos(x+rx, 1.5+ry, z+rz);
				//Density, mass, position, velocity (particle inputs)
				//**************************************
				//TODO - play with #'s here!! 
				Particle * p = new Particle(800, 10, pos, glm::vec3(-0.50, -5.20, 0.0));
				p->setColor(vec3(0.1, 1.0, 0.0)); 
				p->setViscosity(30.f);
				 //This is non-polar liquid (ie oil)
				p->setCi(0.5f);
				//**************************************//
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
		for( float y = 0; y < 0.55; y += 0.1 )
		{
			for( float x = containerMin.x / 1.5f; x < containerMax.x / 1.5f; x += 0.1 )
			{
				for( float z = containerMin.z / 1.5f; z < containerMax.z / 1.5f; z += 0.1 )
				{
					float rx = 0.01f*(float)rand() / RAND_MAX; 
					float ry = 0.01f*(float)rand() / RAND_MAX; 
					float rz = 0.01f*(float)rand() / RAND_MAX; 
					vec3 pos(x+rx, y+ry, z+rz);
					//TODO - mass freaks out when its 0.012 like they give in the paper!!! :(
					Particle * p = new Particle(restDensity+500, mass, pos, glm::vec3(0));
					theParticles.push_back(p);
					p->setColor(vec3(1.0, 0.0, 0.0)); 
					p->setTemp(15.f); 
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
		for( float y = 0.65; y < 1.3; y += 0.1 )
		{
			for( float x = containerMin.x / 1.5f; x < containerMax.x / 1.5f; x += 0.1 )
			{
				for( float z = containerMin.z / 1.5f; z < containerMax.z / 1.5f; z += 0.1 )
				{
					float rx = 0.01f*(float)rand() / RAND_MAX; 
					float ry = 0.01f*(float)rand() / RAND_MAX; 
					float rz = 0.01f*(float)rand() / RAND_MAX; 
					vec3 pos(x+rx, y+ry, z+rz);
					//was 0.006
					//TODO - figure out if we just wanna use the defaults or if we want to change other params
					Particle * p = new Particle(restDensity-500, mass, pos, glm::vec3(0));
					theParticles.push_back(p);
					p->setColor(vec3(0.0, 0.0, 1.0)); 
					p->setTemp(5.f); 
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

	////Heat to the bottom & cool down the top
	//if (frame < 100) {
	//	float max = containerMax.y - 0.3; 
	//	float min = containerMin.y + 0.3; 
	//	for (int i = 0; i < theParticles.size(); i++)
	//	{
	//		if (theParticles.at(i)->getPosition().y < min) {
	//			theParticles.at(i)->setTemp(theParticles.at(i)->getTemp() + 0.1);
	//		} else if (theParticles.at(i)->getPosition().y > max) {
	//			theParticles.at(i)->setTemp(theParticles.at(i)->getTemp() - 0.1);
	//		}
	//	}
	//}
}

//Initialize fluid from some point
void Fluid::addFluid(float dt)
{
	//Density, mass, position, velocity (particle inputs)
	//*
	for (float z = -0.2; z < 0.2; z+=0.1)
    {
		for (float y = -0.2; y < 0.2; y += 0.1 )
		{
			float rx = 0.01*((double) rand() / (RAND_MAX)); 
			float ry = 0.01*((double) rand() / (RAND_MAX)); 
			float rz = 0.01*((double) rand() / (RAND_MAX)); 
			
			vec3 pos(rx+1, y+ry+1, z+rz);
			//Density, mass, position, velocity (particle inputs)
			Particle * p = new Particle(restDensity, mass, pos, glm::vec3(-3, 0, 0.0));
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
	/*/

	for( float y = 0; y < 3; y += 0.1 )
	{
		for( float x = -0.5; x < 0.5; x += 0.1 )
		{
			for( float z = -0.5; z < 0.5; z += 0.1 )
			{
				float rx = 0.01*(float)rand() / RAND_MAX; 
				float ry = 0.01*(float)rand() / RAND_MAX; 
				float rz = 0.01*(float)rand() / RAND_MAX; 
				vec3 pos(x+rx, y+ry, z+rz);
				Particle * p = new Particle(1000, 1, pos, glm::vec3(0));
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
	//*/
}

//**********************************************************************************************
//The main fluid simulation 
//**********************************************************************************************

//Calls all the SPH fns
void Fluid::Update(float dt, glm::vec3& externalForces)
{
	if (loadFromMesh == true && frame == 0) {
		createParticlesFromMesh();
	} else if (loadFromMesh == false && theParticles.size() < 8000 && frame % 4 == 0) {
		addFluid(dt);
		//addMultiFluid();
		//addLavaLamp(); 
		if( frame == -1 )
		{
			for( float y = 0; y < 3; y += 0.1 )
			{
				for( float x = containerMin.x/2; x < containerMax.x/2; x += 0.1 )
				{
					for( float z = containerMin.z/2; z < containerMax.z/2; z += 0.1 )
					{
						float rx = 0.01f*(float)rand() / RAND_MAX; 
						float ry = 0.01f*(float)rand() / RAND_MAX; 
						float rz = 0.01f*(float)rand() / RAND_MAX; 
						vec3 pos(x+rx, y+ry, z+rz);
						Particle * p = new Particle(restDensity-500, mass, pos, glm::vec3(0));
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
	//Additional computations for multiple fluid interaction
	if (useMultiFluidCalcs) {
		//computeDiffusion(dt); 
		//manageAirBubbles();
	}

	computeDensity(dt);
	computeForces(dt, externalForces);
	integrate(dt); 
	resolveCollisions();
	frame++;
	
	container.clear();
	for (int i = 0; i < theParticles.size(); i++)
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
	for (int i = 0; i < theParticles.size(); i++)
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
		/*/
		for (int j = 0 ; j < theParticles.size(); j++)
		{
			//Compute the distance from particle i to j & add it if its within the kernal radius
			if (glm::distance(theParticles.at(i)->getPosition(), theParticles.at(j)->getPosition()) <= h ) {
				theParticles.at(i)->addNeighbor( theParticles.at(j) ); 
			}
		}
		//*/
	}
}

//Finds the density at the current timestep
void Fluid::computeDensity(float dt)
{
	for (int i = 0; i < theParticles.size(); i++)
	{
		float density = 0; //theParticles.at( i ).getRestDensity();
		vec3 pos_i = theParticles.at(i)->getPosition();
		std::vector<Particle*> neighbors = theParticles.at(i)->getNeighbors();
		for (int j = 0; j < neighbors.size(); j++)
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
	float pi;
	if (useMultiFluidCalcs) {
		pi = (theParticles.at(i)->getK() * (theParticles.at(i)->getDensity() - theParticles.at(i)->getRestDensity())); 
	} else {
		pi = (k * (theParticles.at(i)->getDensity() - theParticles.at(i)->getRestDensity())); 
	}
	std::vector<Particle*> neighbors = theParticles.at(i)->getNeighbors();
	for (int j = 0; j < neighbors.size(); j++) {
		float pj;
		if (useMultiFluidCalcs) {
			pj = (neighbors.at(j)->getK() * (neighbors.at(j)->getDensity() - neighbors.at(j)->getRestDensity())); 
		} else {
			pj = (k* (neighbors.at(j)->getDensity() - neighbors.at(j)->getRestDensity())); 
		}
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
	for (int j = 0; j < neighbors.size(); j++)
	{
		glm::vec3 r = theParticles.at(i)->getPosition() - neighbors.at(j)->getPosition();
		glm::vec3 vel = (neighbors.at(j)->getVelocity() - theParticles.at(i)->getVelocity()) / (neighbors.at(j)->getDensity() + 1e-15f);
		if (useMultiFluidCalcs) {
			float muAvg = (neighbors.at(j)->getViscosity() + theParticles.at(i)->getViscosity()) / 2.0f; 
			v += muAvg * neighbors.at(j)->getMass() * vel * wViscosityLap( r, h );
		} else {
			v +=  neighbors.at(j)->getMass() * vel * wViscosityLap( r, h );
		}

	}
	if (useMultiFluidCalcs) 
		return v; 
	else 
		return mu * v; 
}

glm::vec3 Fluid::computeSurfaceTension(float dt, int i)
{
	glm::vec3 n(0.0); 
    float k = 0.0; 
	std::vector<Particle*> neighbors = theParticles.at(i)->getNeighbors();
	for (int j = 0; j < neighbors.size(); j++)
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
	for (int j = 0; j < neighbors.size(); j++)
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

void Fluid::computeForces(float dt, glm::vec3 externalForces)
{
	//TODO - add other forces for now, just add gravity
	for (int i = 0; i < theParticles.size(); i++) 
	{
		glm::vec3 pressureForce = computePressure(dt, i); 
		glm::vec3 viscosityForce = computeViscosity(dt, i); 
		glm::vec3 finalForce;
		if (!useMultiFluidCalcs) {
			glm::vec3 surfaceTension = computeSurfaceTension(dt, i);  //TODO - may want to turn on/off based on multifluid
			finalForce = pressureForce + theParticles.at(i)->getDensity()*externalForces + viscosityForce + surfaceTension;
		} else {
			glm::vec3 interfaceForce; 
			glm::vec3 surfaceTension;
			 
			computeSurfaceAndInterfaceTension(i, interfaceForce, surfaceTension);
			//glm::vec3 buoyancy = computeArtificialBuoyancy(i); 
			finalForce = pressureForce + theParticles.at(i)->getDensity()*externalForces + viscosityForce + surfaceTension + interfaceForce; 
		}
		theParticles.at(i)->setForce( finalForce ); 
	}
}


void Fluid::integrate(float dt)
{
	//Euler just in case leapfrog is wrong
	for (int i = 0; i < theParticles.size(); i++)
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
	for (int i = 0; i < theParticles.size(); i++)
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

	// p.velocity = p.velocity - 2.0 * Dot(normal, p.velocity) * normal;
	for (int i = 0; i < theParticles.size(); i++)
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
			vel.x *= -0.5;
		}
		if (pos.y < containerMin.y) {
			glm::vec3 normal = glm::vec3(0, 1.0, 0); 
			updated = true;
			pos.y = containerMin.y;
			vel.y *= -0.5;  
		}
		if (pos.z < containerMin.z) {
			glm::vec3 normal = glm::vec3(0, 0, 1); 
			updated = true;
			pos.z = containerMin.z;
			vel.z *= -0.5; 
		}

		if (pos.x > containerMax.x) {
			glm::vec3 normal = glm::vec3(-1, 0, 0); 
			updated = true;
			pos.x = containerMax.x;
			vel.x *= -0.5;
		}
		if (pos.y > containerMax.y) {
			glm::vec3 normal = glm::vec3(0, -1.0, 0); 
			updated = true;
			pos.y = containerMax.y;
			vel.y *= -0.5;  
		}
		if (pos.z > containerMax.z) {
			glm::vec3 normal = glm::vec3(0, 0, -1.0); 
			updated = true;
			pos.z = containerMax.z;
			vel.z *= -0.5; 
		}

		if (updated)
		{
			theParticles.at(i)->setPostion(pos); 
			theParticles.at(i)->setVelocity(vel);
		}
	}
}

//********************************************************************************************
//Multiple fluid functions
//********************************************************************************************
void  Fluid::computeDiffusion(float dt)
{
 
	//Temp diffusion constant in the paper -> 0.0001.
	float c = 0.1; 
	for (int i = 0; i < theParticles.size(); i++) 
	{
		float deltaTemp = 0.f;
		for (int j = 0; j < theParticles.at(i)->getNeighbors().size(); j++)
		{
			glm::vec3 r = theParticles.at(i)->getPosition() - theParticles.at(i)->getNeighbors().at(j)->getPosition();
			float mass = theParticles.at(i)->getNeighbors().at(j)->getMass(); 
			float num = theParticles.at(i)->getNeighbors().at(j)->getTemp() - theParticles.at(i)->getTemp(); 
			float avgTemp = num / (theParticles.at(i)->getNeighbors().at(j)->getDensity() + 1e-15f); 
			deltaTemp += mass * avgTemp * wPoly6Lap(r, h); 
		}
		//Euler integration for change
		theParticles.at(i)->setTemp(theParticles.at(i)->getTemp() + c * deltaTemp * dt); 
		
		//Rest Density = alpha / temp
		theParticles.at(i)->setRestDensity(10000.f / theParticles.at(i)->getTemp()); 
	}
}

glm::vec3  Fluid::computeArtificialBuoyancy(int i)
{
	float b = 5.f; //Taken from paper
	glm::vec3 bForce = b * (theParticles.at(i)->getDensity() - theParticles.at(i)->getRestDensity()) * vec3(0, -9.8, 0); 
	//return bForce; 
	return vec3(0.0); 
}

void  Fluid::manageAirBubbles()
{
	//Look for large gradient of cp & positive y component & generate air particles with an offset
	//TODO - not sure what cp is exactly...is it just 1?  Then wouldn't it be the same as cs then??
	for (int i = 0; i < theParticles.size(); i++)
	{
		glm::vec3 csGrad = theParticles.at(i)->getCsGrad(); 
		glm::vec3 cpGrad = theParticles.at(i)->getCpGrad(); 
		//if (glm::length(csGrad) > 10)
		//	cout<<"FOUND GRAD: "<<glm::length(csGrad)<<endl;
		if (glm::length(cpGrad) > 1000 && cpGrad.y > 0)
		{
			vec3 pos =  theParticles.at(i)->getPosition() -0.0003f * csGrad; 
			Particle * p = new Particle(100, mass, pos, theParticles.at(i)->getVelocity() );
			p->setCi(0);
			p->setCs(0); 
			theParticles.push_back(p);
			Box * box = container( pos );
			if( box->frame < frame )
			{
				box->particles.clear();
				box->frame = frame;
			}
			box->particles.push_back( p );
		}

		//Delete any air particles that don't meet min conditions
		//TODO - the second check should be cp!!! 
		if (theParticles.at(i)->getCi() == 0 && theParticles.at(i)->getCs() == 0) {
			if ((length(csGrad) < 0.01 && length(cpGrad.y) > 100) || theParticles.at(i)->getDensity() < 5) {
				int k = 0; 
				theParticles.erase(theParticles.begin() + i); 
				//If we are using leapfrog, also delete the corresponding elements
				//if (accel0.size() >= i) {
				//	accel0.erase(accel0.begin() + i); 
				//	vel0.erase(vel0.begin() + i); 
				//}
			}
		}
	}

}


//********************************************************************************************
//Kernal functions
//********************************************************************************************
float Fluid::wPoly6(float r, float h)
{
	if (0 <= r && r <= h) {
		float c = 315.0f / (64.0f * PI * pow(h, 9)); 
		float w = pow(pow(h, 2) - pow(r, 2), 3); 
		return c * w; 
	} else {
		return 0.0;
	}
}

glm::vec3 Fluid::wPoly6Grad(glm::vec3 r, float h)
{
	float lr = length( r );
	if (0 <= lr && lr <= h) {
		float lrs = lr*lr-h*h;
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

float Fluid::wPoly6Lap(glm::vec3 r, float h)
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
float Fluid::wSpiky(float r, float h)
{
	if (0 <= r && r <= h) {
		float c = 15.0f / (PI * pow(h, 6)); 
		float w = pow (h - r, 3); 
		return c*w;
	} else {
		return 0; 
	}
}
vec3 Fluid::wSpikyGrad(vec3 r, float h)
{
	float lr = length( r );
	if (0 <= lr && lr <= h) {
		float c = 15.0f / (PI * pow(h, 6)); 
		float numerator = pow(lr - h, 2);
		float x = -3.f*r.x*numerator/(lr+1e-15);
		float y = -3.f*r.y*numerator/(lr+1e-15);
		float z = -3.f*r.z*numerator/(lr+1e-15);
		vec3 w = vec3( x, y, z );
		return c*w;
	} else {
		return vec3(0); 
	}
}

float Fluid::wViscosity(float r, float h)
{
	return 0.0f;
}

glm::vec3 Fluid::wViscosityGrad(glm::vec3 r, float h)
{
	return vec3(0.0f);
}

float Fluid::wViscosityLap(glm::vec3 r, float h)
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

float Fluid::field( vec3 pos )
{
	Box * box = container( pos );
	std::vector< Particle * >::iterator it;
	if( box->frame < frame )
		return 0;
	std::vector< Particle * > particles = box->particles.at(0)->getNeighbors();
	float potential = 0;
	for( it = particles.begin(); it != particles.end(); ++it )
	{
		float rho = (*it)->getDensity(); 
		float r = glm::distance( (*it)->getPosition(), pos );
		//fPressure = - sum (mj (tempPi + tempPj) / 2 pj * gradient(W(ri - rj, h))
		potential += (*it)->getMass() / (1e-15f + rho) * wPoly6(r, h); 
	}
	return potential;
}