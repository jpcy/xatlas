$input v_normal, v_texcoord0

#include <bgfx_shader.sh>
#include "shared.h"

SAMPLER2D(s_diffuse, 0);
SAMPLER2D(s_emission, 1);
SAMPLER2D(s_lightmap, 2);

uniform vec4 u_diffuse;
uniform vec4 u_emission;
uniform vec4 u_lightDir;
uniform vec4 u_shade_diffuse_emission;
#define u_shadeType uint(u_shade_diffuse_emission.x)
#define u_diffuseType uint(u_shade_diffuse_emission.y)
#define u_emissionType uint(u_shade_diffuse_emission.z)

void main()
{
	vec4 emission = u_emission;
	if (u_emissionType == EMISSION_TEXTURE)
		emission.rgb = texture2D(s_emission, v_texcoord0.xy).rgb;
	vec4 color = vec4(0.0, 0.0, 0.0, 1.0);
	if (emission.r > 0.0 || emission.g > 0.0 || emission.b > 0.0)
		color = emission;
	else
	{
		vec4 diffuse = u_diffuse;
		if (u_diffuseType == DIFFUSE_TEXTURE)
			diffuse *= texture2D(s_diffuse, v_texcoord0.xy);
		if (u_shadeType == SHADE_FLAT)
			color.rgb = diffuse.rgb * (dot(v_normal, u_lightDir.xyz) * 0.5 + 0.5); // half lambert
		else if (u_shadeType == SHADE_LIGHTMAP)
			color.rgb = diffuse.rgb * texture2D(s_lightmap, v_texcoord0.zw).rgb;
		else if (u_shadeType == SHADE_LIGHTMAP_ONLY)
			color.rgb = texture2D(s_lightmap, v_texcoord0.zw).rgb;
		color.a = diffuse.a;
	}
	gl_FragColor = color;
}
