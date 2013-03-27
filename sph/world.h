#ifndef _SPH_WORLD_H_
#define _SPH_WORLD_H_

#include <iostream>
#include <vector>
#include <string>

#include "memleak.h"

// glew
#include "glew.h"

// glm
#include "../glm/glm.hpp"
#include "../glm/gtc/matrix_transform.hpp"

using namespace glm;

class World
{
public:
	World(const std::string& filename);
    virtual ~World();

    void loadFromFile(const std::string& filename);
    virtual void draw( GLuint ploc, GLuint cloc, GLuint nloc, GLint Mloc );

    enum ShapeType {SPHERE, CUBE, CYLINDER};
    class Shape
    {
    public:
        ShapeType getType() const { return type; }
		virtual void draw( GLuint ploc, GLuint cloc, GLuint nloc, GLint Mloc );
		void setColor( float r, float g, float b ) { setColor( vec3( r, g, b ) ); }
		void translate( vec3 t ) { M = glm::translate( M, t ); }
		void scale( vec3 s ) { M = glm::scale( M, s ); }
		void rotate( vec3 axis, float theta ) { M = glm::rotate( M, theta, axis ); }
		
		virtual void setColor( vec3 rgb ) { color = rgb; }
		
    protected:
		unsigned int vbo;
		unsigned int cbo;
		unsigned int nbo;
		unsigned int ibo;
		unsigned int num;

        Shape(ShapeType t) : type(t), M(1.0f), color(1.0f) {}
        ShapeType type;

		mat4 M;
		vec3 color;

		virtual void initMesh() {}
    };

    class Sphere : public Shape
    {
    public:
        Sphere() : Shape(SPHERE), r(1.0) { initMesh(); }
        double r;
    };

    class Cube : public Shape
    {
    public:
        Cube() : Shape(CUBE), hx(0.5), hy(0.5), hz(0.5) { initMesh(); }
        double hx, hy, hz;
		void setColor( vec3 rgb );

	protected:
		void initMesh();
    };

    class Cylinder : public Shape
    {
    public:
        Cylinder() : Shape(CYLINDER), r(1.0) { initMesh(); }
        vec3 start, end;
        double r;
    };

    std::vector<Shape*> m_shapes;
};

#endif