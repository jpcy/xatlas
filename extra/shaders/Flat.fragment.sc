$input v_normal

#include <bgfx_shader.sh>

uniform vec4 u_diffuse;
uniform vec4 u_emission;
uniform vec4 u_lightDir;

void main()
{
	vec3 color = u_diffuse.rgb * (dot(v_normal, u_lightDir.xyz) * 0.5 + 0.5) + u_emission.rgb; // half lambert
	gl_FragColor = vec4(color, u_diffuse.a);
}
