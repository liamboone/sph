#include "fluid.h"
#include <iostream>

//Globals
const double h = 0.1;
const double k = 300; //TODO - make this function dependant on the temp
const double PI = 3.14159265; 

const bool loadFromMesh = true; 

//Container size 
const vec3 containerMin(-3, 0, -3);
const vec3 containerMax(3, 6, 3); 

Fluid::Fluid(void) : container( 40, 40, 40, vec3( -3, 0, -3 ), vec3( 3, 6, 3 ) )
{
	frame = 0;
	/*srand (time(NULL));*/

	if (loadFromMesh == true)
	{
		theMesh = new obj();
		string file = "..\\stanford_bunny\\bunny_scaled.obj";  
		//string file = "..\\bunny.obj";
		objLoader* loader = new objLoader( file, theMesh );
		theMesh->buildVBOs();
		delete loader;

		cout<<"The mesh size: "<<(*theMesh->getFaces()).size()<<endl; 
	}
}


Fluid::~Fluid(void)
{
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
		for (int i = 0; i < (*theMesh->getFaces()).size(); i++)
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

//Creates particles inside of theMesh obj
void Fluid::createParticlesFromMesh()
{
	//Run through the grid & add particles if we are inside the mesh
	for (float x = containerMin.x; x < containerMax.x; x+= 0.2)
	{
		for (float y = containerMin.y; y < containerMax.y; y+=0.2 )
		{
			for (float z = containerMin.z; z < containerMax.z; z+=0.2)
			{
				vec3 pos(x, y, z); 
				if (insideOutside(pos) == true) {
					//Make a little random
					float rx = 0.001*((double) rand() / (RAND_MAX)); 
					float ry = 0.001*((double) rand() / (RAND_MAX)); 
					float rz = 0.001*((double) rand() / (RAND_MAX)); 
					pos.x += rx;
					pos.y += ry;
					pos.z += rz;
					Particle * p = new Particle(1500, 1, pos, glm::vec3(0.0, 0.0, 0.0));
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


//Initialize fluid from some point
void Fluid::addFluid(float dt)
{
	//Density, mass, position, velocity (particle inputs)
	//*
	for (float z = -0.2; z < 0.2; z+=0.1)
    {
		for (float y = 1.8; y < 2.2; y += 0.1 )
		{
			float rx = 0.01*((double) rand() / (RAND_MAX)); 
			float ry = 0.01*((double) rand() / (RAND_MAX)); 
			float rz = 0.01*((double) rand() / (RAND_MAX)); 
			
			vec3 pos(rx-1.5, y+ry, z+rz);
			//Density, mass, position, velocity (particle inputs)
			Particle * p = new Particle(1500, 1, pos, glm::vec3(-3, 0.0, 0.0));
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
	for( float y = 0; y < 0.5; y += 0.1 )
	{
		for( float x = -0.5; x < 0.5; x += 0.1 )
		{
			for( float z = -0.5; z < 0.5; z += 0.1 )
			{
				float rx = 0.01*(float)rand() / RAND_MAX; 
				float ry = 0.01*(float)rand() / RAND_MAX; 
				float rz = 0.01*(float)rand() / RAND_MAX; 
				vec3 pos(x+rx, y+ry, z+rz);
				Particle * p = new Particle(1500, 1, pos, glm::vec3(0));
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
	} else if (loadFromMesh == false && theParticles.size() < 8000 && frame % 15 == 0) {
		addFluid(dt);
	}
	findNeighbors();
	computeDensity(dt);
	computeForces(dt, externalForces);
	integrate(dt); 
	//computeVelocity(dt);
	//computePosition(dt);
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
	float pi = (k * (theParticles.at(i)->getDensity() - theParticles.at(i)->getRestDensity())); 
	std::vector<Particle*> neighbors = theParticles.at(i)->getNeighbors();
	for (int j = 0; j < neighbors.size(); j++) {
		float pj = (k * (neighbors.at(j)->getDensity() - neighbors.at(j)->getRestDensity())); 
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
		glm::vec3 vel = (neighbors.at(j)->getVelocity() - theParticles.at(i)->getVelocity()) / neighbors.at(j)->getDensity();
		v += neighbors.at(j)->getMass()*vel*wViscosityLap( r, h );
	}
	return 15.f*v; 
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
		float denInv = 1.0 / (neighbors.at(j)->getDensity() + 1e-15f); 
		n += mass * denInv * wPoly6Grad(r, h); 
		k += mass * denInv * wPoly6Lap(r, h); 
	}

	k = -k / (glm::length(n) + 1e-15f); 
	n =  n / (glm::length(n) + 1e-15f);
	return 50.f * k * n; 
}


void Fluid::computeForces(float dt, glm::vec3 externalForces)
{
	//TODO - add other forces for now, just add gravity
	for (int i = 0; i < theParticles.size(); i++) 
	{
		glm::vec3 pressureForce = computePressure(dt, i); 
		glm::vec3 viscosityForce = computeViscosity(dt, i); 
		glm::vec3 surfaceTension = computeSurfaceTension(dt, i);
		glm::vec3 finalForce = pressureForce + theParticles.at(i)->getDensity()*externalForces + viscosityForce + surfaceTension; 
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

	////*********************************************************************
	////Leapfrog 
	//std::vector<glm::vec3> accel1; 
	//for (int i = 0; i < theParticles.size(); i++)
	//{



	//	//Acceleration at v + 1
	//	vec3 a = theParticles.at(i)->getForce() / theParticles.at(i)->getDensity(); 
	//	accel1.push_back(a);

	//	//if (i < 10)
	//	//{
	//	//	cout<<"Particle "<<i<<endl;
	//	//	cout<<"Pos: "<< theParticles.at(i)->getPosition().x<<" "<<theParticles.at(i)->getPosition().y<<" "<<theParticles.at(i)->getPosition().z<<endl;
	//	//	cout<<"Vel: " << theParticles.at(i)->getVelocity().x<<" "<<theParticles.at(i)->getVelocity().y<<" "<<theParticles.at(i)->getVelocity().z<<endl;
	//	//	cout<<"New accel: "<<a.x<<" " <<a.y<<" "<<a.z<<endl;
	//	//	if (accel0.size() >= 10)
	//	//	{
	//	//		cout<<"Old accel: "<<accel0.at(i).x<<" "<<accel0.at(i).y<<" "<<accel0.at(i).z<<endl;
	//	//	}
	//	//	cout<<endl;
	//	//}

	//	//Get v 1/2
	//	vec3 velHalfStep;	
	//	if (i >= accel0.size() || accel0.size() == 0) {
	//		velHalfStep = theParticles.at(i)->getVelocity() + 0.5f * vec3(0, -9.8, 0) * dt; 
	//		theParticles.at(i)->setVelocity(velHalfStep); 
	//	} 
	//	else {
	//		if ( accel0.at(i) != accel0.at(i))
	//			std::cout<<"Error in stored acceleration"<<std::endl;
	//		velHalfStep = theParticles.at(i)->getVelocity() + 0.5f * (accel0.at(i))* dt; 
	//	}

	//	//Get next postion with half step vel
	//	if (i <= accel0.size() && accel0.size() != 0) {
	//		theParticles.at(i)->setVelocity(velHalfStep + 0.5f * a * dt);
	//	}		
	//	theParticles.at(i)->setPostion(theParticles.at(i)->getPosition() + theParticles.at(i)->getVelocity() * dt); 
	//}
	//accel0.clear();
	//accel0 = accel1; 
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
		if (pos.x < -3) {
			glm::vec3 normal = glm::vec3(1, 0, 0); 
			updated = true;
			//Normal of wall is 
			pos.x = -3;
			vel.x *= -0.9;
		}
		if (pos.y < 0) {
			glm::vec3 normal = glm::vec3(0, 1.0, 0); 
			updated = true;
			pos.y = 0;
			vel.y *= -0.9;  
		}
		if (pos.z < -3) {
			glm::vec3 normal = glm::vec3(0, 0, 1); 
			updated = true;
			pos.z = -3;
			vel.z *= -0.9; 
		}

		if (pos.x > 3) {
			glm::vec3 normal = glm::vec3(-1, 0, 0); 
			updated = true;
			pos.x = 3;
			vel.x *= -0.9;
		}
		if (pos.y > 6) {
			glm::vec3 normal = glm::vec3(0, -1.0, 0); 
			updated = true;
			pos.y = 6;
			vel.y *= -0.9;  
		}
		if (pos.z > 3) {
			glm::vec3 normal = glm::vec3(0, 0, -1.0); 
			updated = true;
			pos.z = 3;
			vel.z *= -0.9; 
		}

		if (updated)
		{
			theParticles.at(i)->setPostion(pos); 
			theParticles.at(i)->setVelocity(vel);
		}
	}
}


float Fluid::wPoly6(float r, float h)
{
	if (0 <= r && r <= h) {
		float c = 315.0 / (64.0 * PI * pow(h, 9)); 
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
		float c = 315.0 / (64.0 * PI * pow(h, 9)); 
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
	float c = 315.0 / (64.0 * PI * pow(h, 9)); 
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
		float c = 15.0 / (PI * pow(h, 6)); 
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
		float c = 15.0 / (PI * pow(h, 6)); 
		float numerator = pow(lr - h, 2);
		float x = -3*r.x*numerator/(lr+1e-15);
		float y = -3*r.y*numerator/(lr+1e-15);
		float z = -3*r.z*numerator/(lr+1e-15);
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

