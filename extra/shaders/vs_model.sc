$input a_color0, a_normal, a_position, a_texcoord0
$output v_color0, v_normal, v_texcoord0

#include <bgfx_shader.sh>

void main()
{
	v_color0 = a_color0;
	v_normal = normalize(mul(u_model[0], vec4(a_normal, 0.0)).xyz);
	v_texcoord0 = a_texcoord0;
	gl_Position = mul(u_modelViewProj, vec4(a_position, 1.0));
}
