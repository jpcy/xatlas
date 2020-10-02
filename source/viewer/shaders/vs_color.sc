$input a_color0, a_position
$output v_color0

#include <bgfx_shader.sh>

void main()
{
	v_color0 = a_color0;
	gl_Position = mul(u_modelViewProj, vec4(a_position, 1.0));
}
