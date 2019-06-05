$input a_position, a_texcoord0
$output v_barycentric

#include <bgfx_shader.sh>

void main()
{
	v_barycentric = a_texcoord0.xyz;
	gl_Position = mul(u_modelViewProj, vec4(a_position, 1.0));
}
