$input v_normal, v_texcoord0

#include <bgfx_shader.sh>
#include "shared.h"

SAMPLER2D(s_diffuse, 0);
SAMPLER2D(s_emission, 1);
SAMPLER2D(s_lightmap, 2);
SAMPLER2D(s_faceData, 3);

uniform vec4 u_diffuse;
uniform vec4 u_emission;
uniform vec4 u_lightDir;
uniform vec4 u_shade_overlay_diffuse_emission;
#define u_shadeType uint(u_shade_overlay_diffuse_emission.x)
#define u_overlayType uint(u_shade_overlay_diffuse_emission.y)
#define u_diffuseType uint(u_shade_overlay_diffuse_emission.z)
#define u_emissionType uint(u_shade_overlay_diffuse_emission.w)
uniform vec4 u_textureSize_cellSize2;
uniform vec4 u_overlayOpacity_colorChartType;
#define u_overlayOpacity u_overlayOpacity_colorChartType.x
#define u_colorChartType uint(u_overlayOpacity_colorChartType.y)
uniform vec4 u_meshColor_primitiveIdStart;
#define u_meshColor u_meshColor_primitiveIdStart.rgb
#define u_primitiveIdStart uint(u_meshColor_primitiveIdStart.w)

void main()
{
	vec4 emission = u_emission;
	if (u_emissionType == EMISSION_TEXTURE)
		emission.rgb = texture2D(s_emission, v_texcoord0.xy).rgb;
	vec4 color = vec4(0.0, 0.0, 0.0, 1.0);
	if (emission.r > 0.0 || emission.g > 0.0 || emission.b > 0.0)
		color = emission;
	else
	{
		vec4 diffuse = u_diffuse;
		if (u_diffuseType == DIFFUSE_TEXTURE)
			diffuse *= texture2D(s_diffuse, v_texcoord0.xy);
		if (u_shadeType == SHADE_FLAT)
			color.rgb = diffuse.rgb * (dot(v_normal, u_lightDir.xyz) * 0.5 + 0.5); // half lambert
		else if (u_shadeType == SHADE_LIGHTMAP)
			color.rgb = diffuse.rgb * texture2D(s_lightmap, v_texcoord0.zw).rgb;
		else if (u_shadeType == SHADE_LIGHTMAP_ONLY)
			color.rgb = texture2D(s_lightmap, v_texcoord0.zw).rgb;
		color.a = diffuse.a;
	}
	if (u_overlayType != OVERLAY_NONE)
	{
		vec3 overlayColor = vec3_splat(1.0);
		if (u_overlayType == OVERLAY_CHART)
		{
			int u = int((u_primitiveIdStart + uint(gl_PrimitiveID)) % FACE_DATA_TEXTURE_WIDTH);
			int v = int((u_primitiveIdStart + uint(gl_PrimitiveID)) / FACE_DATA_TEXTURE_WIDTH);
			vec4 faceData = texelFetch(s_faceData, ivec2(u, v), 0);
			uint x = uint(v_texcoord0.z * u_textureSize_cellSize2.x);
			uint y = uint(v_texcoord0.w * u_textureSize_cellSize2.y);
			uint cellSize = uint(u_textureSize_cellSize2.z);
			float scale = 1.0;
			if (cellSize > 0u)
				scale = (x / cellSize % 2u) != (y / cellSize % 2u) ? 0.75 : 1.0;
			uint chartType = uint(faceData.a + 0.5);
			if (u_colorChartType == CHART_TYPE_ANY || u_colorChartType == chartType)
				overlayColor = faceData.rgb * scale;
			else
				overlayColor = vec3_splat(0.75) * scale;
		}
		else if (u_overlayType == OVERLAY_MESH)
			overlayColor = u_meshColor;
		else if (u_overlayType == OVERLAY_STRETCH)
		{
			int u = int((u_primitiveIdStart + uint(gl_PrimitiveID)) % FACE_DATA_TEXTURE_WIDTH);
			int v = int((u_primitiveIdStart + uint(gl_PrimitiveID)) / FACE_DATA_TEXTURE_WIDTH);
			vec4 faceData = texelFetch(s_faceData, ivec2(u, v), 0);
			overlayColor = mix(vec3(0.1, 0.9, 0.1), vec3(0.9, 0.1, 0.01), faceData.x);
		}
		color.rgb = color.rgb * (1.0 - u_overlayOpacity) + overlayColor * u_overlayOpacity;
	}
	gl_FragColor = color;
}
