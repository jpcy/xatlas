$input v_color0, v_texcoord0

#include <bgfx_shader.sh>
#include "shared.h"

uniform vec4 u_textureSize_cellSize;
uniform vec4 u_colorChartType;

void main()
{
	uint x = uint(v_texcoord0.z * u_textureSize_cellSize.x);
	uint y = uint(v_texcoord0.w * u_textureSize_cellSize.y);
	uint cellSize = uint(u_textureSize_cellSize.z);
	float scale = 1.0;
	if (cellSize > 0u)
		scale = (x / cellSize % 2u) != (y / cellSize % 2u) ? 0.75 : 1.0;
	uint colorChartType = uint(u_colorChartType.x);
	uint chartType = uint(v_color0.a + 0.5);
	vec3 rgb;
	if (colorChartType == CHART_TYPE_ANY || colorChartType == chartType)
		rgb = v_color0.rgb * scale;
	else
		rgb = vec3_splat(0.75);
	gl_FragColor = vec4(rgb, 1.0);
}
