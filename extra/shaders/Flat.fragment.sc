$input v_normal

#include <bgfx_shader.sh>

uniform vec4 u_color;
uniform vec4 u_lightDir;

void main()
{
	vec3 diffuse = u_color.rgb * (dot(v_normal, u_lightDir.xyz) * 0.5 + 0.5); // half lambert
	gl_FragColor = vec4(diffuse, u_color.a);
}
