$input v_normal, v_texcoord0

#include <bgfx_shader.sh>
#include "shared.h"

SAMPLER2D(u_lightmapSampler, 0);

uniform vec4 u_diffuse;
uniform vec4 u_emission;
uniform vec4 u_lightDir_shadeType;
#define u_lightDir u_lightDir_shadeType.xyz
#define u_shadeType uint(u_lightDir_shadeType.w)

void main()
{
	if (u_shadeType == SHADE_FLAT) {
		vec3 color = u_diffuse.rgb * (dot(v_normal, u_lightDir) * 0.5 + 0.5); // half lambert
		gl_FragColor = vec4(color, u_diffuse.a);
	} else if (u_shadeType == SHADE_EMISSIVE) {
		gl_FragColor = vec4(u_emission.rgb, 1.0);
	} else {
		vec3 color = u_diffuse.rgb * texture2D(u_lightmapSampler, v_texcoord0.zw).rgb;
		gl_FragColor = vec4(color, u_diffuse.a);
	}
}
