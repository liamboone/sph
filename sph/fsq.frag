#version 330 core

uniform ivec2 u_resolution;
uniform vec3 u_camPosition;
uniform sampler3D u_distanceMap;

in vec2 fs_coord;

out vec4 outColor;

// 1.047 ~= 60 degrees

#define GRAD_SAMPLE 0.01
#define THRESHOLD 0.4

vec3 castRayToFloor( vec3 pos, vec3 dir, vec3 lpos )
{
	float tfloor = -pos.y / dir.y;
	vec3 origin = pos + dir * tfloor;
		
	vec3 light = normalize(lpos-origin);
	float diffuseTerm = clamp( dot(vec3(0,1,0), light), 0.0, 1.0);

	if( tfloor > 0 )
	{
		vec2 checker = 0.25 + 0.25*floor( mod( origin.xz, 2 ) );
		return vec3( vec3( 0.9*diffuseTerm*(checker.x+checker.y) ) );
	}
	
	return vec3( 0.6, 0.8, 1 );
}

void main()
{
	vec3 lpos = vec3(-5,10,3);
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
		outColor = vec4( castRayToFloor( eye, view, lpos ), 1 );
	}
	else
	{
		float t = max( tmin.x, max( tmin.y, tmin.z ) );
		float exit = min( tmax.x, min( tmax.y, tmax.z ) );

		outColor = vec4( castRayToFloor( eye, view, lpos ), 1 );

		vec3 origin = eye + view * t;
		for( int i = 0; i < 200; i ++ )
		{
			vec3 dMapCoords = ( origin - vec3( -3, 0, -3 ) ) / 6;
			float dist = texture( u_distanceMap, dMapCoords ).x;
			if( dist > THRESHOLD )
			{
				for( int j = 0; j < 10; j ++ )
				{
					t -= 0.01;
					origin = eye + view * t;
					dist = texture( u_distanceMap, dMapCoords ).x;
					if( dist < THRESHOLD ) 
					{
						t += 0.01;
						break;
					}
				}
				origin = eye + view * t;
				vec3 light = normalize(origin-lpos);
				outColor = vec4( 0.5,0.5,1,1 );
				vec3 normal = normalize( vec3(	texture( u_distanceMap, dMapCoords + vec3( GRAD_SAMPLE, 0, 0 ) ).x - texture( u_distanceMap, dMapCoords - vec3( GRAD_SAMPLE, 0, 0 ) ).x,
												texture( u_distanceMap, dMapCoords + vec3( 0, GRAD_SAMPLE, 0 ) ).x - texture( u_distanceMap, dMapCoords - vec3( 0, GRAD_SAMPLE, 0 ) ).x,
												texture( u_distanceMap, dMapCoords + vec3( 0, 0, GRAD_SAMPLE ) ).x - texture( u_distanceMap, dMapCoords - vec3( 0, 0, GRAD_SAMPLE ) ).x));
				float diffuseTerm = clamp( dot(normal, light), 0.0, 1.0);
				float specularTerm = clamp( pow( max( 0.0, dot( reflect(-light, normal), view)), 50.0), 0.0, 1.0);
				vec4 refcolor = vec4( castRayToFloor( origin, normalize(refract( normalize(view), -normal, 1.15 )), lpos ), 1 );
				outColor = refcolor*0.9 + 0.1*vec4( 0.5,0.5,1,1 )*diffuseTerm + specularTerm;
				break;
			}
			t += 0.05;//min( dist*0.8, 0.1 );
			origin = eye + view * t;
			if( t > exit )
				break;
		}
	}
}