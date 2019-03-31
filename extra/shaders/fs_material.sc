$input v_normal, v_texcoord0

#include <bgfx_shader.sh>
#include "shared.h"

SAMPLER2D(u_diffuseSampler, 0);
SAMPLER2D(u_emissionSampler, 1);
SAMPLER2D(u_lightmapSampler, 2);

uniform vec4 u_diffuse;
uniform vec4 u_emission;
uniform vec4 u_lightDir;
uniform vec4 u_shade_diffuse_emission;
#define u_shadeType uint(u_shade_diffuse_emission.x)
#define u_diffuseType uint(u_shade_diffuse_emission.y)
#define u_emissionType uint(u_shade_diffuse_emission.z)

void main()
{
	vec3 diffuse = u_diffuse.rgb;
	if (u_diffuseType == DIFFUSE_TEXTURE)
		diffuse *= texture2D(u_diffuseSampler, v_texcoord0.xy).rgb;
	vec3 color = vec3_splat(0.0);
	if (u_shadeType == SHADE_FLAT)
		color = diffuse * (dot(v_normal, u_lightDir.xyz) * 0.5 + 0.5); // half lambert
	else if (u_shadeType == SHADE_EMISSIVE)
	{
		vec3 emission = u_emission.rgb;
		if (u_emissionType == EMISSION_TEXTURE)
			emission *= texture2D(u_emissionSampler, v_texcoord0.xy).rgb;
		color = emission;
	}
	else if (u_shadeType == SHADE_LIGHTMAP)
		color = diffuse * texture2D(u_lightmapSampler, v_texcoord0.zw).rgb;
	gl_FragColor = vec4(color, 1.0);
}
