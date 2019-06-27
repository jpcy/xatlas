$input v_color0, v_texcoord0

#include <bgfx_shader.sh>

uniform vec4 u_textureSize_cellSize;

void main()
{
	uint x = uint(v_texcoord0.z * u_textureSize_cellSize.x);
	uint y = uint(v_texcoord0.w * u_textureSize_cellSize.y);
	uint cellSize = uint(u_textureSize_cellSize.z);
	float scale = 1.0;
	if (cellSize > 0u)
		scale = (x / cellSize % 2u) != (y / cellSize % 2u) ? 0.75 : 1.0;
	gl_FragColor = vec4(v_color0.rgb * scale, v_color0.a);
}
