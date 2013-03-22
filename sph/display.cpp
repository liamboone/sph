#include "display.h"

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

	//Always remember that it doesn't do much good if you don't have OpenGL actually use the shaders
	glUseProgram(shaderProgram);

	lightCol = vec3( 1.0f, 1.0f, 1.0f );
	lightPos = vec3( 5.0f, 10.0f, 3.0f );

	world = new World( "" );

	camera = new Camera();
	camera->setViewport( 640, 480 );

	vec3 cpos = camera->getPos();
	mat4 cmat = camera->getMat4();
	mat4 mmat( 1.0f );

	glUniformMatrix4fv(u_modelMatrixLocation, 1, GL_FALSE, &mmat[0][0]);
	glUniformMatrix4fv(u_projMatrixLocation, 1, GL_FALSE, &cmat[0][0]);
	glUniform3f(u_lightPositionLocation, lightPos.x, lightPos.y, lightPos.z );
	glUniform3f(u_lightColorLocation, lightCol.x, lightCol.y, lightCol.z );
	glUniform3f(u_camPositionLocation, cpos.x, cpos.y, cpos.z );
}

//Thanks to Tiantian Liu @ University of Pennsylvania, 2012
void Display::initShaders()
{
	//here is stuff for setting up our shaders
	const char* fragFile = "pass.frag";
	const char* vertFile = "pass.vert";
	vertexShader = glCreateShader(GL_VERTEX_SHADER);
	fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
	shaderProgram = glCreateProgram();
	
	//load up the source, compile and link the shader program
	const char* vertSource = textFileRead(vertFile);
	const char* fragSource = textFileRead(fragFile);
	glShaderSource(vertexShader, 1, &vertSource, 0);
	glShaderSource(fragmentShader, 1, &fragSource, 0);
	glCompileShader(vertexShader);
	glCompileShader(fragmentShader);

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
}

void Display::draw()
{
	world->draw( positionLocation, colorLocation, normalLocation, u_modelMatrixLocation );
	//TODO: draw particles
}

void Display::updateCamera()
{
	vec3 cpos = camera->getPos();
	mat4 cmat = camera->getMat4();
	
	glUniformMatrix4fv(u_projMatrixLocation, 1, GL_FALSE, &cmat[0][0]);
	glUniform3f(u_camPositionLocation, cpos.x, cpos.y, cpos.z );
}
