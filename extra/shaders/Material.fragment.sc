$input v_normal, v_texcoord0

#include <bgfx_shader.sh>

SAMPLER2D(u_lightmap, 0);

uniform vec4 u_diffuse;
uniform vec4 u_emission;
uniform vec4 u_lightDir_shadeType;
#define u_lightDir u_lightDir_shadeType.xyz
#define u_shadeType uint(u_lightDir_shadeType.w)

void main()
{
	if (u_shadeType == 0u) {
		vec3 color = u_diffuse.rgb * (dot(v_normal, u_lightDir) * 0.5 + 0.5) + u_emission.rgb; // half lambert
		gl_FragColor = vec4(color, u_diffuse.a);
	} else if (u_emission.r > 0.0 || u_emission.g > 0.0 || u_emission.b > 0.0) {
		gl_FragColor = vec4(u_emission.rgb, 1.0);
	} else {
		vec3 color = u_diffuse.rgb * texture2D(u_lightmap, v_texcoord0.zw).rgb;
		gl_FragColor = vec4(color, u_diffuse.a);
	}
}
