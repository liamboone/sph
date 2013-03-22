#version 330 

uniform vec3 u_lightColor;
uniform vec3 u_camPosition;
uniform sampler2D u_shadowMap;

//these are the interpolated values out of the rasterizer, so you can't know
//their specific values without knowing the vertices that contributed to them
in vec3 fs_normal;
in vec3 fs_light_vector;
in vec3 fs_color;
in vec3 fs_position;
in float fs_isLine;
in vec4 fs_shadow;

out vec4 out_Color;

void main() {
    //base colors for materials
    vec4 diffuseColor = vec4(fs_color, 1.0);
	vec4 lightColor = vec4(u_lightColor, 1.0);
    
	vec3 view = normalize( u_camPosition - fs_position );
	vec3 light = normalize(fs_light_vector);
	vec3 normal = normalize(fs_normal);

	vec4 shadow = fs_shadow / fs_shadow.w;
	shadow.z += 0.000005;

	float ldist = texture2D( u_shadowMap, shadow.st ).z;

	float shadowTerm = 1.0;
	if( fs_shadow.w > 0 ) shadowTerm = ( ldist < shadow.z ? 0.5 : 1.0 );

    //calculate diffuse term and clamp to the range [0, 1]
    float diffuseTerm = clamp(dot(normal, light), 0.0, 1.0);
    
    //calculate ambient term
    float ambientTerm = 0.2;
    
    out_Color = lightColor * diffuseColor * clamp( diffuseTerm + ambientTerm, 0.0, 1.0 );

	if ( diffuseTerm > 0.1 && ldist > shadow.z ) 
		out_Color += lightColor * pow( max( 0.0, dot( reflect(-light, normal), view)), 50.0);

	out_Color *= shadowTerm;

	if ( fs_isLine > 0 )
		out_Color = diffuseColor; 
}