#include "display.h"
#include "../glm/gtc/type_ptr.hpp"

#define SHADOW_MAP_SIZE 1024

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
	glewInit();
	
	// assign locations
	positionLocation = 0;
	colorLocation = 1;
	normalLocation = 2;
	
	//Everybody does this
	glClearColor(0, 0, 0, 1);
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
	
	u_shadowModelMatrixLocation	= glGetUniformLocation(shadowShaderProgram, "u_modelMatrix");
	u_shadowProjMatrixLocation	= glGetUniformLocation(shadowShaderProgram, "u_projMatrix");

	lightCol = vec3( 1.0f, 1.0f, 1.0f );
	lightPos = vec3( 5.0f, 10.0f, 3.0f );

	world = new World( "" );

	camera = new Camera();
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
}

//Thanks to Tiantian Liu @ University of Pennsylvania, 2012
void Display::initShaders()
{
	//here is stuff for setting up our shaders
	const char* fragFile = "diffuse.frag";
	const char* vertFile = "diffuse.vert";
	vertexShader = glCreateShader(GL_VERTEX_SHADER);
	fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	shaderProgram = glCreateProgram();
	
	const char* fragShadowFile = "shadow.frag";
	const char* vertShadowFile = "shadow.vert";
	vertexShadowShader = glCreateShader(GL_VERTEX_SHADER);
	fragmentShadowShader = glCreateShader(GL_FRAGMENT_SHADER);
	shadowShaderProgram = glCreateProgram();
	
	//load up the source, compile and link the shader program
	const char* vertSource = textFileRead(vertFile);
	const char* fragSource = textFileRead(fragFile);
	glShaderSource(vertexShader, 1, &vertSource, 0);
	glShaderSource(fragmentShader, 1, &fragSource, 0);
	glCompileShader(vertexShader);
	glCompileShader(fragmentShader);
	
	const char* vertShadowSource = textFileRead(vertFile);
	const char* fragShadowSource = textFileRead(fragFile);
	glShaderSource(vertexShadowShader, 1, &vertSource, 0);
	glShaderSource(fragmentShadowShader, 1, &fragSource, 0);
	glCompileShader(vertexShadowShader);
	glCompileShader(fragmentShadowShader);

	//For your convenience, i decided to throw in some compiler/linker output helper functions
	//from CIS 565
	GLint compiled;
	glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &compiled);
	if (!compiled)
	{
		printShaderInfoLog(vertexShader);
	} 
	glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &compiled);
	if (!compiled)
	{
		printShaderInfoLog(fragmentShader);
	} 
	
	//set the attribute locations for our shaders
	glBindAttribLocation(shaderProgram, positionLocation, "vs_position");
	glBindAttribLocation(shaderProgram, normalLocation, "vs_normal");
	glBindAttribLocation(shaderProgram, colorLocation, "vs_color");

	//finish shader setup
	glAttachShader(shaderProgram, vertexShader);
	glAttachShader(shaderProgram, fragmentShader);
	glLinkProgram(shaderProgram);
	
	//check for linking success
	GLint linked;
	glGetProgramiv(shaderProgram,GL_LINK_STATUS, &linked);
	if (!linked) 
	{
		printLinkInfoLog(shaderProgram);
	}

	//compile shadow stuff
	glGetShaderiv(vertexShadowShader, GL_COMPILE_STATUS, &compiled);
	if (!compiled)
	{
		printShaderInfoLog(vertexShadowShader);
	} 
	glGetShaderiv(fragmentShadowShader, GL_COMPILE_STATUS, &compiled);
	if (!compiled)
	{
		printShaderInfoLog(fragmentShadowShader);
	} 
	
	//set the attribute locations for our shaders
	glBindAttribLocation(shaderProgram, shadowPositionLocation, "vs_position");

	//finish shader setup
	glAttachShader(shadowShaderProgram, vertexShadowShader);
	glAttachShader(shadowShaderProgram, fragmentShadowShader);
	glLinkProgram(shadowShaderProgram);
	
	//check for linking success
	glGetProgramiv(shadowShaderProgram,GL_LINK_STATUS, &linked);
	if (!linked) 
	{
		printLinkInfoLog(shadowShaderProgram);
	}
}

void Display::draw()
{
	std::vector<Particle> particles = theFluid->getParticles();
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

	world->draw( positionLocation, colorLocation, normalLocation, u_shadowModelMatrixLocation );
	//TODO: draw particles
	for (int i = 0; i < particles.size(); i++)
	{
		World::Shape * particle = new World::Cube();
		particle->translate(particles.at(i).getPosition()); 
		particle->scale(vec3(0.1));
		particle->draw( positionLocation, colorLocation, normalLocation, u_shadowModelMatrixLocation );
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
	world->draw( positionLocation, colorLocation, normalLocation, u_modelMatrixLocation );
	//TODO: draw particles
	for (int i = 0; i < particles.size(); i++)
	{
		World::Shape * particle = new World::Cube();
		particle->translate(particles.at(i).getPosition()); 
		particle->scale(vec3(0.1));
		particle->draw( positionLocation, colorLocation, normalLocation, u_modelMatrixLocation );
	}
	//END render from camera
}

void Display::updateCamera()
{
	vec3 cpos = camera->getPos();
	mat4 cmat = camera->getMat4();
	
	glUniformMatrix4fv(u_projMatrixLocation, 1, GL_FALSE, &cmat[0][0]);
	glUniform3f(u_camPositionLocation, cpos.x, cpos.y, cpos.z );
}

void Display::setFluids(Fluid *fluid)
{
	theFluid = fluid; 
}
