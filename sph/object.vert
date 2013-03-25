#version 330

uniform mat4 u_projMatrix;
uniform mat4 u_modelMatrix;

in vec4 vs_position;

void main()
{
	//built-in things to pass down the pipeline
	gl_Position = u_projMatrix * u_modelMatrix * vs_position;
}