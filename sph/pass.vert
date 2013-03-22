#version 330

uniform mat4 u_projMatrix;
uniform mat4 u_modelMatrix;
uniform vec3 u_lightPosition;

in vec4 vs_position;
in vec4 vs_normal;
in vec3 vs_color;

out vec3 fs_position;
out vec3 fs_normal;
out vec3 fs_color;
out vec3 fs_light_vector;

void main()
{
	fs_color = vs_color;
	fs_normal = (transpose(inverse(u_modelMatrix)) * vs_normal).xyz;
	fs_position = (u_modelMatrix * vs_position).xyz;
	
	fs_light_vector = vec3(u_lightPosition - (u_modelMatrix * vs_position).xyz);

	//built-in things to pass down the pipeline
	gl_Position = u_projMatrix * u_modelMatrix * vs_position;
}