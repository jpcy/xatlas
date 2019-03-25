$input a_texcoord0
$output v_texcoord0

#include <bgfx_shader.sh>

uniform vec4 u_textureSize_cellSize;

void main()
{
	v_texcoord0 = a_texcoord0;
	gl_Position = mul(u_modelViewProj, vec4(a_texcoord0.zw * u_textureSize_cellSize.xy, 0.0, 1.0));
}
