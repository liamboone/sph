#version 330

layout(points) in;
layout(triangle_strip, max_vertices = 4) out;

out vec2 fs_coord;

void main() 
{
    gl_Position = vec4( 1.0, 1.0, 0.0, 1.0 );
    fs_coord = vec2( 0.5, 0.5 );
    EmitVertex();

    gl_Position = vec4(-1.0, 1.0, 0.0, 1.0 );
    fs_coord = vec2( -0.5, 0.5 ); 
    EmitVertex();

    gl_Position = vec4( 1.0,-1.0, 0.0, 1.0 );
    fs_coord = vec2( 0.5, -0.5 ); 
    EmitVertex();

    gl_Position = vec4(-1.0,-1.0, 0.0, 1.0 );
    fs_coord = vec2( -0.5, -0.5 ); 
    EmitVertex();

    EndPrimitive(); 
}