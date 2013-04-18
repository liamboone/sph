#include "display.h"
#include "../glm/gtc/type_ptr.hpp"

#include <algorithm>

#define SHADOW_MAP_SIZE 1024
#define DMAP_SIZE 64

//BEGIN glsl utilities
// from swiftless.com
char* Display::textFileRead(const char* fileName) {
    char* text;
    
    assert(fileName != NULL);

	FILE *file = fopen(fileName, "rt");
	
	assert(file != NULL);
	
	fseek(file, 0, SEEK_END);
	int count = ftell(file);
	rewind(file);
	
	if (count > 0) {
		text = (char*)malloc(sizeof(char) * (count + 1));
		count = fread(text, sizeof(char), count, file);
		text[count] = '\0';	//cap off the string with a terminal symbol, fixed by Cory
	}
	fclose(file);
	
    return text;
}
void Display::printLinkInfoLog(int prog) 
{
	int infoLogLen = 0;
	int charsWritten = 0;
	GLchar *infoLog;

	glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &infoLogLen);

	// should additionally check for OpenGL errors here

	if (infoLogLen > 0)
	{
		infoLog = new GLchar[infoLogLen];
		// error check for fail to allocate memory omitted
		glGetProgramInfoLog(prog,infoLogLen, &charsWritten, infoLog);
		std::cout << "InfoLog:" << std::endl << infoLog << std::endl;
		delete [] infoLog;
	}
}
void Display::printShaderInfoLog(int shader)
{
	int infoLogLen = 0;
	int charsWritten = 0;
	GLchar *infoLog;

	glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLogLen);

	// should additionally check for OpenGL errors here

	if (infoLogLen > 0)
	{
		infoLog = new GLchar[infoLogLen];
		// error check for fail to allocate memory omitted
		glGetShaderInfoLog(shader,infoLogLen, &charsWritten, infoLog);
		std::cout << "InfoLog:" << std::endl << infoLog << std::endl;
		delete [] infoLog;
	}

	// should additionally check for OpenGL errors here
}
//END glsl utilities

Display::Display() 
{
	colorMap[1] = vec3( 1.0f, 0.5f, 0.5f );
	colorMap[2] = vec3( 0.5f, 0.5f, 1.0f );
	colorMap[4] = vec3( 0.5f, 1.0f, 0.5f );
}

Display::~Display()
{
	delete camera;
	delete world;
	
	glDeleteBuffers(1, &fbo);

	//Tear down the shader program in reverse of building it
	/*glDetachShader(shaderProgram, vertexShader);
	glDetachShader(shaderProgram, fragmentShader);
	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);
	glDeleteProgram(shaderProgram);*/
}

void Display::init()
{
	glewInit();
	
	world = new World( "" );
	camera = new Camera();

	// assign locations
	positionLocation = 0;
	shadowPositionLocation = 0;
	normalLocation = 1;
	colorLocation = 2;
	
	//Everybody does this
	glClearColor(0.6f, 0.8f, 1.0f, 1.0f);
	glEnable(GL_DEPTH_TEST);
	glClearDepth(1.0);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_CULL_FACE);
	glFrontFace(GL_CCW);

	initShaders();

	//Get the uniform locations for our shaders, unfortunately they can not be set by us, we have
	//to ask OpenGL for them
	u_modelMatrixLocation		= glGetUniformLocation(shaderProgram, "u_modelMatrix");
	u_projMatrixLocation		= glGetUniformLocation(shaderProgram, "u_projMatrix");
	u_lightPositionLocation		= glGetUniformLocation(shaderProgram, "u_lightPosition");
	u_lightColorLocation		= glGetUniformLocation(shaderProgram, "u_lightColor");
	u_camPositionLocation		= glGetUniformLocation(shaderProgram, "u_camPosition");
	u_shadowMapLocation			= glGetUniformLocation(shaderProgram, "u_shadowMap");
	u_shadowBiasMatrixLocation	= glGetUniformLocation(shaderProgram, "u_shadowProjMatrix");

	u_resolutionLocation		= glGetUniformLocation(raymarchShaderProgram, "u_resolution");
	u_rayCamPositionLocation	= glGetUniformLocation(raymarchShaderProgram, "u_camPosition");
	u_distanceMapLocation		= glGetUniformLocation(raymarchShaderProgram, "u_distanceMap");
	u_cMinLocation				= glGetUniformLocation(raymarchShaderProgram, "u_cMin");
	u_cMaxLocation				= glGetUniformLocation(raymarchShaderProgram, "u_cMax");
	
	u_shadowModelMatrixLocation	= glGetUniformLocation(shadowShaderProgram, "u_modelMatrix");
	u_shadowProjMatrixLocation	= glGetUniformLocation(shadowShaderProgram, "u_projMatrix");
	

	lightCol = vec3( 1.0f, 1.0f, 1.0f );
	lightPos = vec3( -5.0f, 10.0f, 3.0f );

	camera->setViewport( 640, 480 );

	vec3 cpos = camera->getPos();
	mat4 cmat = camera->getMat4();
	mat4 mmat( 1.0f );
	
	glUseProgram(shaderProgram);
	glUniformMatrix4fv(u_modelMatrixLocation, 1, GL_FALSE, &mmat[0][0]);
	glUniformMatrix4fv(u_projMatrixLocation, 1, GL_FALSE, &cmat[0][0]);
	glUniform3f(u_lightPositionLocation, lightPos.x, lightPos.y, lightPos.z );
	glUniform3f(u_lightColorLocation, lightCol.x, lightCol.y, lightCol.z );
	glUniform3f(u_camPositionLocation, cpos.x, cpos.y, cpos.z );

	glUseProgram(raymarchShaderProgram);
	glUniform2i(u_resolutionLocation, camera->getWidth(), camera->getHeight() );
	glUniform3f(u_rayCamPositionLocation, cpos.x, cpos.y, cpos.z );
	
	// The framebuffer
	glGenFramebuffers(1, &fbo);
	glBindFramebuffer(GL_FRAMEBUFFER, fbo);
 
	// Depth texture. Slower than a depth buffer, but you can sample it later in your shader
	glGenTextures(1, &depthTexture);
	glBindTexture(GL_TEXTURE_2D, depthTexture);
	glTexImage2D(GL_TEXTURE_2D, 0,GL_DEPTH_COMPONENT16, SHADOW_MAP_SIZE, SHADOW_MAP_SIZE,
		0,GL_DEPTH_COMPONENT, GL_FLOAT, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
 
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthTexture, 0);

	glDrawBuffer(GL_NONE); // No color buffer is drawn to.
 
	// Always check that our framebuffer is ok
	assert(glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE);
	
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glBindTexture(GL_TEXTURE_2D, 0);
	
	// request 1 texture name from OpenGL
	glGenTextures(1, &fluidTexture);	
	// tell OpenGL we're going to be setting up the texture name it gave us	
	glBindTexture(GL_TEXTURE_3D, fluidTexture);	
	// when this texture needs to be shrunk to fit on small polygons, use linear interpolation of the texels to determine the color
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	// when this texture needs to be magnified to fit on a big polygon, use linear interpolation of the texels to determine the color
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	// we want the texture to repeat over the S axis, so if we specify coordinates out of range we still get textured.
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP);
	// same as above for T axis
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP);
	// same as above for R axis
	glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP);
}

void Display::loadShader( const char* vertFile, const char* fragFile, const char* geomFile, Shader & shader )
{
	shader.vertex = glCreateShader(GL_VERTEX_SHADER);
	shader.fragment = glCreateShader(GL_FRAGMENT_SHADER);
	if( geomFile != NULL ) shader.geometry = glCreateShader(GL_GEOMETRY_SHADER);
	shader.program = glCreateProgram();
	
	//load up the source, compile and link the shader program
	const char* vertSource = textFileRead(vertFile);
	const char* fragSource = textFileRead(fragFile);
	const char* geomSource = NULL;
	if( geomFile != NULL ) geomSource = textFileRead(geomFile);
	glShaderSource(shader.vertex, 1, &vertSource, 0);
	glShaderSource(shader.fragment, 1, &fragSource, 0);
	if( geomFile != NULL ) glShaderSource(shader.geometry, 1, &geomSource, 0);
	glCompileShader(shader.vertex);
	glCompileShader(shader.fragment);
	if( geomFile != NULL ) glCompileShader(shader.geometry);

	//For your convenience, i decided to throw in some compiler/linker output helper functions
	//from CIS 565
	GLint compiled;
	glGetShaderiv(shader.vertex, GL_COMPILE_STATUS, &compiled);
	if (!compiled)
	{
		printShaderInfoLog(shader.vertex);
	} 
	glGetShaderiv(shader.fragment, GL_COMPILE_STATUS, &compiled);
	if (!compiled)
	{
		printShaderInfoLog(shader.fragment);
	} 
	if( geomFile != NULL )
	{
		glGetShaderiv(shader.geometry, GL_COMPILE_STATUS, &compiled);
		if (!compiled)
		{
			printShaderInfoLog(shader.geometry);
		} 
	}
	//finish shader setup
	glAttachShader(shader.program, shader.vertex);
	glAttachShader(shader.program, shader.fragment);
	if( geomFile != NULL ) glAttachShader(shader.program, shader.geometry);
	glLinkProgram(shader.program);
	
	//check for linking success
	GLint linked;
	glGetProgramiv(shader.program,GL_LINK_STATUS, &linked);
	if (!linked) 
	{
		printLinkInfoLog(shader.program);
	}
}

//Thanks to Tiantian Liu @ University of Pennsylvania, 2012
void Display::initShaders()
{
	Shader shader;
	loadShader( "diffuse.vert", "diffuse.frag", NULL, shader );
	shaderProgram = shader.program;
	glBindAttribLocation(shaderProgram, positionLocation, "vs_position");
	glBindAttribLocation(shaderProgram, normalLocation, "vs_normal");
	glBindAttribLocation(shaderProgram, colorLocation, "vs_color");

	loadShader( "shadow.vert", "shadow.frag", NULL, shader );
	shadowShaderProgram = shader.program;
	glBindAttribLocation(shadowShaderProgram, shadowPositionLocation, "vs_position");
	
	loadShader( "fsq.vert", "fsq.frag", "fsq.geom", shader );
	raymarchShaderProgram = shader.program;
}

void Display::draw()
{
	std::vector<Particle*> particles = theFluid->getParticles();
	int w = camera->getWidth();
	int h = camera->getHeight();

	//BEGIN render from light
	glUseProgram(shadowShaderProgram);
	glBindFramebuffer(GL_FRAMEBUFFER, fbo);
	// Clear previous frame values
	glClear( GL_DEPTH_BUFFER_BIT);
	
	//Disable color rendering, we only want to write to the Z-Buffer
	glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE); 
	glCullFace( GL_FRONT );

	camera->setViewport(SHADOW_MAP_SIZE, SHADOW_MAP_SIZE);
	glViewport(0,0,SHADOW_MAP_SIZE,SHADOW_MAP_SIZE);
	vec3 pos = camera->getPos();
	camera->setPos( vec3( lightPos.x, lightPos.y, lightPos.z ) );
	
	mat4 cmat = camera->getMat4();
	
	glUniformMatrix4fv(u_shadowProjMatrixLocation, 1, GL_FALSE, &cmat[0][0]);

	if( pos.y > 0 ) 
		world->draw( shadowPositionLocation, colorLocation, normalLocation, u_shadowModelMatrixLocation );
	//TODO: draw particles
	World::Shape * particle = new World::Cube();
	for (unsigned int i = 0; i < particles.size(); i++)
	{
		int pidx = particles.at( i )->getIndex();
		if( flags & pidx )
		{
			particle->clearMat();
			particle->translate(particles.at(i)->getPosition()); 
			particle->scale( vec3( 0.04 ) );
			particle->draw( shadowPositionLocation, colorLocation, normalLocation, u_shadowModelMatrixLocation );
		}
	}
	//END render from light

	//BEGIN render from camera
	glUseProgram(shaderProgram);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUniform1i( u_shadowMapLocation, 7 );
	glActiveTexture(GL_TEXTURE7);
	glBindTexture(GL_TEXTURE_2D,depthTexture);

	float bias[16] = {	0.5, 0.0, 0.0, 0.0, 
						0.0, 0.5, 0.0, 0.0,
						0.0, 0.0, 0.5, 0.0,
						0.5, 0.5, 0.5, 1.0  };
	
	
	mat4 cameraMatrix = make_mat4( bias ) * cmat;
	glUniformMatrix4fv(u_shadowBiasMatrixLocation, 1, GL_FALSE, &cameraMatrix[0][0]);
	
	camera->setViewport(w, h);
	camera->setPos( pos );
	glViewport(0,0,w,h);
	glCullFace( GL_BACK );
	if( pos.y > 0 ) 
		world->draw( positionLocation, colorLocation, normalLocation, u_modelMatrixLocation );
	//TODO: draw particles
	vec3 red( 1,0,0 );
	vec3 blue( 0,0,1 );

	for (unsigned int i = 0; i < particles.size(); i++)
	{
		int pidx = particles.at( i )->getIndex();
		if( flags & pidx )
		{
			particle->clearMat();
			particle->translate(particles.at(i)->getPosition()); 
			particle->scale( vec3( 0.04 ) );

			if( flags & DFLAG_VEL )
			{
				particle->setColor( 0.1f+glm::clamp( particles.at(i)->getVelocity()*particles.at(i)->getVelocity()/5.0f, vec3(0.0), vec3(1.0) ) );
			}
			else if( flags & DFLAG_TEMP )
			{
				float alpha = ( particles.at( i )->getTemp() - 5 ) / 10;
				particle->setColor(alpha*red + (1-alpha)*blue);
			}
			else
			{
				particle->setColor( colorMap[pidx] );
			}
			particle->draw( positionLocation, colorLocation, normalLocation, u_modelMatrixLocation );
		}
	}	
	delete particle;
	
	//END render from camera
	glBindTexture(GL_TEXTURE_2D,0);
}

void Display::updateCamera()
{
	vec3 cpos = camera->getPos();
	mat4 cmat = camera->getMat4();
	
	glUseProgram(shaderProgram);
	glUniformMatrix4fv(u_projMatrixLocation, 1, GL_FALSE, &cmat[0][0]);
	glUniform3f(u_camPositionLocation, cpos.x, cpos.y, cpos.z );
	
	glUseProgram(raymarchShaderProgram);
	glUniform2i(u_resolutionLocation, camera->getWidth(), camera->getHeight() );
	glUniform3f(u_rayCamPositionLocation, cpos.x, cpos.y, cpos.z );
}

void Display::setFluids(Fluid *fluid)
{
	theFluid = fluid; 
}

void Display::march()
{
	glUseProgram( raymarchShaderProgram );
	glUniform1i( u_distanceMapLocation, 4 );
	glActiveTexture(GL_TEXTURE4);
	glBindTexture(GL_TEXTURE_3D,fluidTexture);
	std::map<int, bool> toCalc;
	std::vector<Particle*> particles = theFluid->getParticles();
	
	vec3 cMin = theFluid->containerMin;
	vec3 cMax = theFluid->containerMax;
	vec3 span = cMax-cMin;
	
	glUniform3f( u_cMinLocation, cMin.x, cMin.y, cMin.z );
	glUniform3f( u_cMaxLocation, cMax.x, cMax.y, cMax.z );

	for (unsigned int pid = 0; pid < particles.size(); pid++)
	{
		vec3 P = (particles.at(pid)->getPosition()-cMin)/span;

		float x = P.x;
		float y = P.y;
		float z = P.z;

		int i,j,k;
	
		i = (int) floor( DMAP_SIZE * x );
		j = (int) floor( DMAP_SIZE * y );
		k = (int) floor( DMAP_SIZE * z );

		int index = i + j*DMAP_SIZE + k*DMAP_SIZE*DMAP_SIZE;

		toCalc[index+1+DMAP_SIZE+DMAP_SIZE*DMAP_SIZE] = true;
		toCalc[index+1+DMAP_SIZE+0					] = true;
		toCalc[index+1+DMAP_SIZE-DMAP_SIZE*DMAP_SIZE] = true;
		toCalc[index+1+0		+DMAP_SIZE*DMAP_SIZE] = true;
		toCalc[index+1+0		+0					] = true;
		toCalc[index+1+0		-DMAP_SIZE*DMAP_SIZE] = true;
		toCalc[index+1-DMAP_SIZE+DMAP_SIZE*DMAP_SIZE] = true;
		toCalc[index+1-DMAP_SIZE+0					] = true;
		toCalc[index+1-DMAP_SIZE-DMAP_SIZE*DMAP_SIZE] = true;

		toCalc[index+0+DMAP_SIZE+DMAP_SIZE*DMAP_SIZE] = true;
		toCalc[index+0+DMAP_SIZE+0					] = true;
		toCalc[index+0+DMAP_SIZE-DMAP_SIZE*DMAP_SIZE] = true;
		toCalc[index+0+0		+DMAP_SIZE*DMAP_SIZE] = true;
		toCalc[index+0+0		+0					] = true;
		toCalc[index+0+0		-DMAP_SIZE*DMAP_SIZE] = true;
		toCalc[index+0-DMAP_SIZE+DMAP_SIZE*DMAP_SIZE] = true;
		toCalc[index+0-DMAP_SIZE+0					] = true;
		toCalc[index+0-DMAP_SIZE-DMAP_SIZE*DMAP_SIZE] = true;

		toCalc[index-1+DMAP_SIZE+DMAP_SIZE*DMAP_SIZE] = true;
		toCalc[index-1+DMAP_SIZE+0					] = true;
		toCalc[index-1+DMAP_SIZE-DMAP_SIZE*DMAP_SIZE] = true;
		toCalc[index-1+0		+DMAP_SIZE*DMAP_SIZE] = true;
		toCalc[index-1+0		+0					] = true;
		toCalc[index-1+0		-DMAP_SIZE*DMAP_SIZE] = true;
		toCalc[index-1-DMAP_SIZE+DMAP_SIZE*DMAP_SIZE] = true;
		toCalc[index-1-DMAP_SIZE+0					] = true;
		toCalc[index-1-DMAP_SIZE-DMAP_SIZE*DMAP_SIZE] = true;
	}

	float * texels = new float[ DMAP_SIZE*DMAP_SIZE*DMAP_SIZE*4 ];
	for( int i = 0; i < DMAP_SIZE*DMAP_SIZE*DMAP_SIZE*4; i +=4 )
	{
		texels[i] = 0.25;
		texels[i+1] = 0.25;
		texels[i+2] = 0.25;
		texels[i+3] = 0;
	}

	std::map<int,bool>::iterator it;
	for( it = toCalc.begin(); it != toCalc.end(); ++ it )
	{
		int i = it->first;
		if( i >= 0 && i < DMAP_SIZE*DMAP_SIZE*DMAP_SIZE )
		{
			float x = (float)( (i)%DMAP_SIZE ) / (float)(DMAP_SIZE-1) * span.x + cMin.x;
			float y = (float)( (i/DMAP_SIZE)%DMAP_SIZE ) / (float)(DMAP_SIZE-1) * span.y + cMin.y;
			float z = (float)( (i/DMAP_SIZE/DMAP_SIZE)%DMAP_SIZE ) / (float)(DMAP_SIZE-1) * span.z + cMin.z;
			vec4 field = theFluid->field( vec3( x, y, z ), colorMap );
			texels[4*i] = field.x;
			texels[4*i+1] = field.y;
			texels[4*i+2] = field.z;
			texels[4*i+3] = field.w;
		}
	}

	glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA, DMAP_SIZE, DMAP_SIZE, DMAP_SIZE, 0, GL_RGBA, GL_FLOAT, texels);
	delete texels;

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//tip the dominoes
	glBegin(GL_POINTS);
	glVertex3f(1.0f, 1.0f, 1.0f);
	glEnd();
}