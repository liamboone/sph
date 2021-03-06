#include "world.h"
#include "../ObjCore/objloader.h"

World::World(const std::string& filename) 
{
    loadFromFile(filename);
}

World::~World() 
{
    for (unsigned int i = 0; i < m_shapes.size(); i++)
    {
        delete m_shapes.at(i);
    }
    m_shapes.clear();
}

void World::loadFromFile( const std::string& filename )
{
	//TODO: actually read from a file...
	
	World::Shape * shape = new World::Cube();

	shape->translate( vec3( 0.0f, -0.05f, 0.0f ) );
	shape->scale( vec3( 100.0f, 0.1f, 100.0f ) );
	shape->setColor( vec3( 1.0f ) );
	
	m_shapes.push_back( shape );
	/*
	shape = new World::Mesh("../obj/bunny.obj");

	shape->translate( vec3( 0.0f, -0.05f, 0.0f ) );
	shape->scale( vec3( 100.0f, 0.1f, 100.0f ) );
	shape->setColor( vec3( 1.0f ) );
	
	m_shapes.push_back( shape );
	*/
}

void World::Mesh::initMesh( )
{
	glGenBuffers(1, &vbo);
	glGenBuffers(1, &cbo);
	glGenBuffers(1, &nbo);
	glGenBuffers(1, &ibo);

	obj * meshObj = new obj();
    objLoader* loader = new objLoader( mesh, meshObj );
    meshObj->buildVBOs();
    delete loader;

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, meshObj->getVBOsize() * sizeof(float), meshObj->getVBO(), GL_STATIC_DRAW); 

	glBindBuffer(GL_ARRAY_BUFFER, nbo);
	glBufferData(GL_ARRAY_BUFFER, meshObj->getNBOsize() * sizeof(float), meshObj->getNBO(), GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, meshObj->getIBOsize() * sizeof(unsigned short), meshObj->getIBO(), GL_STATIC_DRAW);
	
	num = meshObj->getIBOsize();

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	delete meshObj;
}

void World::draw( GLuint ploc, GLuint cloc, GLuint nloc, GLint Mloc )
{
	for( unsigned int i = 0; i < m_shapes.size(); i ++ )
	{
		m_shapes.at( i )->draw( ploc, cloc, nloc, Mloc );
	}
}

World::Shape::~Shape()
{
	glDeleteBuffers(1, &vbo);
	glDeleteBuffers(1, &cbo);
	glDeleteBuffers(1, &nbo);
	glDeleteBuffers(1, &ibo);
}

void World::Cube::setColor( vec3 rgb )
{
	color = rgb;
	float* colors = new float[72];

	for( int i = 0; i < 72; i+=3 )
	{
		colors[i  ] = color.r;
		colors[i+1] = color.g;
		colors[i+2] = color.b;
	}

	glBindBuffer(GL_ARRAY_BUFFER, cbo);
	glBufferData(GL_ARRAY_BUFFER, 72 * sizeof(float), colors, GL_STATIC_DRAW);

	delete[] colors;
}

void World::Sphere::setColor( vec3 rgb )
{
	color = rgb;
	float* colors = new float[675];

	for( int i = 0; i < 675; i+=3 )
	{
		colors[i  ] = color.r;
		colors[i+1] = color.g;
		colors[i+2] = color.b;
	}

	glBindBuffer(GL_ARRAY_BUFFER, cbo);
	glBufferData(GL_ARRAY_BUFFER, 675 * sizeof(float), colors, GL_STATIC_DRAW);

	delete[] colors;
}

void World::Shape::draw( GLuint ploc, GLuint cloc, GLuint nloc, GLint Mloc )
{
	//activate our three kinds of information
	glEnableVertexAttribArray(ploc);
	glEnableVertexAttribArray(cloc);
	glEnableVertexAttribArray(nloc);
	
	//bind again
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glVertexAttribPointer(ploc, 4, GL_FLOAT, 0, 0, static_cast<char*>(0));
	
	glBindBuffer(GL_ARRAY_BUFFER, cbo);
	glVertexAttribPointer(cloc, 3, GL_FLOAT, 0, 0, static_cast<char*>(0));
	
	glBindBuffer(GL_ARRAY_BUFFER, nbo);
	glVertexAttribPointer(nloc, 4, GL_FLOAT, 0, 0, static_cast<char*>(0));

	//set the modelview matrix again since it changed
	glUniformMatrix4fv(Mloc, 1, GL_FALSE, &M[0][0]);

	//draw again, even the indices from before are good
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
	glDrawElements(GL_TRIANGLES, num, GL_UNSIGNED_SHORT, 0);
	
	//shut off the information since we're done drawing
	glDisableVertexAttribArray(ploc);
	glDisableVertexAttribArray(cloc);
	glDisableVertexAttribArray(nloc);
}

void World::Cube::initMesh()
{
	//Create the VBOs and IBO we'll be using to render images in OpenGL
	glGenBuffers(1, &vbo);
	glGenBuffers(1, &cbo);
	glGenBuffers(1, &nbo);
	glGenBuffers(1, &ibo);

	//set the number of indicies
	num = 36;

	float* vertices = new float[96];
	float* normals = new float[96];
	unsigned short* indices = new unsigned short[num];
	float* colors = new float[72];
	
	float* verts = vertices;
	float* norms = normals;
	unsigned short* indxs = indices;
	
	//TOP face 
	verts[ 0] = 0.5; verts[ 1] = 0.5; verts[ 2] = 0.5; verts[ 3] = 1;
	verts[ 4] =-0.5; verts[ 5] = 0.5; verts[ 6] = 0.5; verts[ 7] = 1;
	verts[ 8] =-0.5; verts[ 9] = 0.5; verts[10] =-0.5; verts[11] = 1;
	verts[12] = 0.5; verts[13] = 0.5; verts[14] =-0.5; verts[15] = 1;

	norms[ 0] = 0; norms[ 1] = 1; norms[ 2] = 0; norms[ 3] = 0;
	norms[ 4] = 0; norms[ 5] = 1; norms[ 6] = 0; norms[ 7] = 0;
	norms[ 8] = 0; norms[ 9] = 1; norms[10] = 0; norms[11] = 0;
	norms[12] = 0; norms[13] = 1; norms[14] = 0; norms[15] = 0;

	indxs[ 0] = 2; indxs[ 1] = 1; indxs[ 2] = 0; 
	indxs[ 3] = 0; indxs[ 4] = 3; indxs[ 5] = 2; 

	verts += 16;
	norms += 16;
	indxs += 6;
	
	//BOTTOM face verts
	verts[ 0] = 0.5; verts[ 1] =-0.5; verts[ 2] = 0.5; verts[ 3] = 1;
	verts[ 4] = 0.5; verts[ 5] =-0.5; verts[ 6] =-0.5; verts[ 7] = 1;
	verts[ 8] =-0.5; verts[ 9] =-0.5; verts[10] =-0.5; verts[11] = 1;
	verts[12] =-0.5; verts[13] =-0.5; verts[14] = 0.5; verts[15] = 1;

	norms[ 0] = 0; norms[ 1] =-1; norms[ 2] = 0; norms[ 3] = 0;
	norms[ 4] = 0; norms[ 5] =-1; norms[ 6] = 0; norms[ 7] = 0;
	norms[ 8] = 0; norms[ 9] =-1; norms[10] = 0; norms[11] = 0;
	norms[12] = 0; norms[13] =-1; norms[14] = 0; norms[15] = 0;

	indxs[ 0] = 6; indxs[ 1] = 5; indxs[ 2] = 4; 
	indxs[ 3] = 4; indxs[ 4] = 7; indxs[ 5] = 6; 

	verts += 16;
	norms += 16;
	indxs += 6;

	//FRONT face verts
	verts[ 0] = 0.5; verts[ 1] = 0.5; verts[ 2] = 0.5; verts[ 3] = 1;
	verts[ 4] = 0.5; verts[ 5] =-0.5; verts[ 6] = 0.5; verts[ 7] = 1;
	verts[ 8] =-0.5; verts[ 9] =-0.5; verts[10] = 0.5; verts[11] = 1;
	verts[12] =-0.5; verts[13] = 0.5; verts[14] = 0.5; verts[15] = 1;

	norms[ 0] = 0; norms[ 1] = 0; norms[ 2] = 1; norms[ 3] = 0;
	norms[ 4] = 0; norms[ 5] = 0; norms[ 6] = 1; norms[ 7] = 0;
	norms[ 8] = 0; norms[ 9] = 0; norms[10] = 1; norms[11] = 0;
	norms[12] = 0; norms[13] = 0; norms[14] = 1; norms[15] = 0;

	indxs[ 0] =10; indxs[ 1] = 9; indxs[ 2] = 8; 
	indxs[ 3] = 8; indxs[ 4] =11; indxs[ 5] =10; 

	verts += 16;
	norms += 16;
	indxs += 6;

	//BACK face verts
	verts[ 0] = 0.5; verts[ 1] = 0.5; verts[ 2] =-0.5; verts[ 3] = 1;
	verts[ 4] =-0.5; verts[ 5] = 0.5; verts[ 6] =-0.5; verts[ 7] = 1;
	verts[ 8] =-0.5; verts[ 9] =-0.5; verts[10] =-0.5; verts[11] = 1;
	verts[12] = 0.5; verts[13] =-0.5; verts[14] =-0.5; verts[15] = 1;

	norms[ 0] = 0; norms[ 1] = 0; norms[ 2] =-1; norms[ 3] = 0;
	norms[ 4] = 0; norms[ 5] = 0; norms[ 6] =-1; norms[ 7] = 0;
	norms[ 8] = 0; norms[ 9] = 0; norms[10] =-1; norms[11] = 0;
	norms[12] = 0; norms[13] = 0; norms[14] =-1; norms[15] = 0;

	indxs[ 0] =14; indxs[ 1] =13; indxs[ 2] =12; 
	indxs[ 3] =12; indxs[ 4] =15; indxs[ 5] =14; 

	verts += 16;
	norms += 16;
	indxs += 6;

	//RIGHT face verts
	verts[ 0] = 0.5; verts[ 1] = 0.5; verts[ 2] = 0.5; verts[ 3] = 1;
	verts[ 4] = 0.5; verts[ 5] =-0.5; verts[ 6] = 0.5; verts[ 7] = 1;
	verts[ 8] = 0.5; verts[ 9] =-0.5; verts[10] =-0.5; verts[11] = 1;
	verts[12] = 0.5; verts[13] = 0.5; verts[14] =-0.5; verts[15] = 1;

	norms[ 0] = 1; norms[ 1] = 0; norms[ 2] = 0; norms[ 3] = 0;
	norms[ 4] = 1; norms[ 5] = 0; norms[ 6] = 0; norms[ 7] = 0;
	norms[ 8] = 1; norms[ 9] = 0; norms[10] = 0; norms[11] = 0;
	norms[12] = 1; norms[13] = 0; norms[14] = 0; norms[15] = 0;

	indxs[ 0] =16; indxs[ 1] =17; indxs[ 2] =18; 
	indxs[ 3] =18; indxs[ 4] =19; indxs[ 5] =16; 

	verts += 16;
	norms += 16;
	indxs += 6;

	//LEFT face verts
	verts[ 0] =-0.5; verts[ 1] = 0.5; verts[ 2] = 0.5; verts[ 3] = 1;
	verts[ 4] =-0.5; verts[ 5] = 0.5; verts[ 6] =-0.5; verts[ 7] = 1;
	verts[ 8] =-0.5; verts[ 9] =-0.5; verts[10] =-0.5; verts[11] = 1;
	verts[12] =-0.5; verts[13] =-0.5; verts[14] = 0.5; verts[15] = 1;

	norms[ 0] =-1; norms[ 1] = 0; norms[ 2] = 0; norms[ 3] = 0;
	norms[ 4] =-1; norms[ 5] = 0; norms[ 6] = 0; norms[ 7] = 0;
	norms[ 8] =-1; norms[ 9] = 0; norms[10] = 0; norms[11] = 0;
	norms[12] =-1; norms[13] = 0; norms[14] = 0; norms[15] = 0;

	indxs[ 0] =20; indxs[ 1] =21; indxs[ 2] =22; 
	indxs[ 3] =22; indxs[ 4] =23; indxs[ 5] =20;
	
	for( int i = 0; i < 72; i+=3 )
	{
		colors[i  ] = color.r;
		colors[i+1] = color.g;
		colors[i+2] = color.b;
	}

	glBindBuffer(GL_ARRAY_BUFFER, cbo);
	glBufferData(GL_ARRAY_BUFFER, 72 * sizeof(float), colors, GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, 96 * sizeof(float), vertices, GL_STATIC_DRAW); 

	glBindBuffer(GL_ARRAY_BUFFER, nbo);
	glBufferData(GL_ARRAY_BUFFER, 96 * sizeof(float), normals, GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, num * sizeof(unsigned short), indices, GL_STATIC_DRAW);
	
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	delete[] colors;
	delete[] indices;
	delete[] normals;
	delete[] vertices;
}

void World::Sphere::initMesh( )
{
	//Create the VBOs and IBO we'll be using to render images in OpenGL
	glGenBuffers(1, &vbo);
	glGenBuffers(1, &cbo);
	glGenBuffers(1, &nbo);
	glGenBuffers(1, &ibo);

	//set the number of indicies
	num = 1350;

	float* vertices = new float[15*15*4];
	float* normals = new float[15*15*4];
	unsigned short* indices = new unsigned short[num];

	int i = 0;
	int idx = 0;
	for( float theta = 0; theta < 15; theta += 1 )
	{
		for( float azimuth = 0; azimuth < 15; azimuth += 1 )
		{
			float x = sin( (azimuth)/14.0f * 3.14159f )*cos(theta/7.5f * 3.14159f);
			float z = sin( (azimuth)/14.0f * 3.14159f )*sin(theta/7.5f * 3.14159f);
			float y = cos( (azimuth)/14.0f * 3.14159f );
			vertices[i] = x/2; vertices[i+1] = y/2; vertices[i+2] = z/2; vertices[i+3] = 1;
			normals[i] = x; normals[i+1] = y; normals[i+2] = z; normals[i+3] = 0;
			i+=4;
			if( azimuth < 14 )
			{
				indices[idx  ] = (int)theta * 15 + (int)azimuth;
				indices[idx+1] = (int)( theta+1 )%15 * 15 + (int)azimuth + 1;
				indices[idx+2] = (int)theta * 15 + (int)azimuth + 1;
				
				indices[idx+3] = (int)theta * 15 + (int)azimuth;
				indices[idx+4] = (int)( theta+1 )%15 * 15 + (int)azimuth;
				indices[idx+5] = (int)( theta+1 )%15 * 15 + (int)azimuth + 1;
				idx += 6;
			}
		}
	}
	
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, 900 * sizeof(float), vertices, GL_STATIC_DRAW); 

	glBindBuffer(GL_ARRAY_BUFFER, nbo);
	glBufferData(GL_ARRAY_BUFFER, 900 * sizeof(float), normals, GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, 1350 * sizeof(unsigned short), indices, GL_STATIC_DRAW);
	
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	delete[] indices;
	delete[] normals;
	delete[] vertices;
}
