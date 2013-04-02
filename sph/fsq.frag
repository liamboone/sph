#version 330 core

uniform vec2 u_resolution;
uniform vec3 u_camPosition;
uniform sampler3D u_distanceMap;

in vec2 fs_coord;

out vec4 outColor;

// 1.047 ~= 60 degrees

void main()
{
	float ratio = float(u_resolution.x) / float(u_resolution.y);
	vec2 coords = vec2( fs_coord.x * ratio, fs_coord.y ) * sin( 1.047 );
	vec3 eye = u_camPosition;
	vec3 look = normalize( -eye );
	vec3 up = vec3( 0,1,0 );
	vec3 right = normalize( cross( look, up ) );
	up = normalize( cross( right, look ) );
	vec3 view = normalize( look + coords.x*right + coords.y*up );

	float tfloor = -eye.y / view.y;

	vec3 t1 = ( vec3( -3, 0, -3 ) - eye ) / view;
	vec3 t2 = ( vec3( 3, 6, 3 ) - eye ) / view;
	
	vec3 tmin = min( t1, t2 );
	vec3 tmax = max( t1, t2 );

	if( tmax.x < tmin.y || tmax.x < tmin.z || tmax.y < tmin.x || 
		tmax.y < tmin.z || tmax.z < tmin.x || tmax.z < tmin.y || eye.y < 0)
	{
		if( tfloor > 0 )
			outColor = vec4( 1 );
		else
			outColor = vec4( view, 1 );
	}
	else
	{
		float t = max( tmin.x, max( tmin.y, tmin.z ) );
		float exit = min( tmax.x, min( tmax.y, tmax.z ) );
		vec3 origin = eye + view * t;
		outColor = vec4( 0.9 );
		for( int i = 0; i < 100; i ++ )
		{
			vec3 dMapCoords = ( origin - vec3( -3, 0, -3 ) ) / 6;
			float dist = texture( u_distanceMap, dMapCoords ).x;
			t += dist * 0.99;
			origin = eye + view * t;
			if( dist < 0.01 )
			{
				outColor = vec4( 1,0,0,1 );
				break;
			}
			if( t > exit )
				break;
		}
	}
}