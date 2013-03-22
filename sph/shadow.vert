#version 330
 
uniform mat4 u_shadowProjMatrix;

in vec4 vs_position;
 
void main()
{
	gl_Position =  u_shadowProjMatrix * vs_position;
}