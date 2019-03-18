#include <bgfx_compute.sh>

FRAMEBUFFER_UIMAGE2D_RW(u_rayBundleHeaderSampler, r32ui, 0);
FRAMEBUFFER_UIMAGE2D_RW(u_rayBundleDataSampler, rgba32ui, 1);
FRAMEBUFFER_IMAGE2D_RW(u_lightmapSampler, rgba32f, 2);

uniform vec4 u_rayBundleDataResolution;
#define u_dataResolution uint(u_rayBundleDataResolution.x)

ivec2 dataUv(uint offset, uint pixel)
{
	return ivec2((offset * 3u + pixel) % u_dataResolution, (offset * 3u + pixel) / u_dataResolution);
}

struct Node
{
	vec3 color;
	float depth;
	vec2 texcoord;
};

#define MAX_NODES 16

void main()
{
#if BGFX_SHADER_LANGUAGE_GLSL
	ivec2 uv = ivec2(gl_FragCoord.xy);
	uint offset = imageLoad(u_rayBundleHeaderSampler, uv).x;
	if (offset != 0xffffffff) {
		Node nodes[MAX_NODES];
		uint numNodes = 0u;
		while (offset != 0xffffffff && numNodes < MAX_NODES) {
			uvec4 color_offset = imageLoad(u_rayBundleDataSampler, dataUv(offset, 0u));
			uvec4 normal_depth = imageLoad(u_rayBundleDataSampler, dataUv(offset, 1u));
			uvec4 texcoord = imageLoad(u_rayBundleDataSampler, dataUv(offset, 2u));
			nodes[numNodes].color.r = uintBitsToFloat(color_offset.r);
			nodes[numNodes].color.g = uintBitsToFloat(color_offset.g);
			nodes[numNodes].color.b = uintBitsToFloat(color_offset.b);
			nodes[numNodes].depth = uintBitsToFloat(normal_depth.w);
			nodes[numNodes].texcoord = vec2(uintBitsToFloat(texcoord.x), uintBitsToFloat(texcoord.y));
			offset = color_offset.w;
			numNodes++;
		}
		for (uint i = 0u; i < numNodes; i++) {
			for (uint j = i + 1u; j < numNodes; j++) {
				if (nodes[i].depth > nodes[j].depth) {
					Node temp = nodes[i];
					nodes[i] = nodes[j];
					nodes[j] = temp;
				}
			}
		}
		vec3 color = vec3_splat(0.0);
		for (uint j = 0u; j < numNodes; j++) {
			color = mix(color, nodes[j].color, 0.5);
		}
		imageStore(u_lightmapSampler, uv, vec4(color, 1.0));
	} else {
		imageStore(u_lightmapSampler, uv, vec4_splat(0.0));
	}
#endif
}
