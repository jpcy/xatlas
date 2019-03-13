$input a_normal, a_position
$output v_normal

#include <bgfx_shader.sh>

void main()
{
	v_normal = normalize(mul(u_model[0], vec4(a_normal, 0.0)).xyz);
	gl_Position = mul(u_modelViewProj, vec4(a_position, 1.0));
}
