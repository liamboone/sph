// codebit originally written by Tiantian Liu @ University of Pennsylvania, 2012
#include "camera.h"

void Camera::init( float _fovy, vec3 _pos, vec3 _target, vec3 _up, float _znear, float _zfar, int _width, int _height )
{
	fovy = _fovy;
	pos = _pos;
	target = _target;
	up = _up;

	zNear = 0.1;
	zFar = 30.0;

	width = 640;
	height = 480;
}

void Camera::setViewport( int w, int h )
{
	width = w;
	height = h;
}

mat4 Camera::getMat4()
{
	mat4 projection = perspective(fovy, static_cast<float>(width) / static_cast<float>(height), zNear, zFar);
	mat4 view = lookAt(pos, target, up);

	return (projection * view);
}

void Camera::zoom( int dz )
{
	vec3 newEye = pos-target;
	newEye = glm::clamp(length(newEye)+dz/10.0f, 1.0f, 50.0f)*normalize(newEye);
	pos = newEye + target;
}

void Camera::orbit( int dx, int dy )
{
	vec3 newEye = pos-target;
	float R = length( newEye );
	newEye = normalize( newEye );
	float zenith = glm::clamp( (float) ( acos( newEye.y ) - dy/180.0f*3.14159f ), 0.1f, 3.04159f );
	float azimuth = atan2( newEye.z, newEye.x ) + dx/180.0f*3.14159f;
	
	newEye.x = R*sin(zenith)*cos(azimuth);
	newEye.y = R*cos(zenith);
	newEye.z = R*sin(zenith)*sin(azimuth);

	pos = newEye+target;
}