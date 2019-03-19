#include <bgfx_compute.sh>

FRAMEBUFFER_UIMAGE2D_RW(u_rayBundleHeaderSampler, r32ui, 0);
FRAMEBUFFER_UIMAGE2D_RW(u_rayBundleDataSampler, rgba32ui, 1);
FRAMEBUFFER_IMAGE2D_RW(u_lightmapSampler, rgba32f, 2);

uniform vec4 u_lightmapSize_dataSize;
#define u_lightmapSize u_lightmapSize_dataSize.xy
#define u_dataSize uint(u_lightmapSize_dataSize.z)
uniform vec4 u_rayNormal;

ivec2 dataUv(uint offset, uint pixel)
{
	return ivec2((offset * 3u + pixel) % u_dataSize, (offset * 3u + pixel) / u_dataSize);
}

#if BGFX_SHADER_LANGUAGE_GLSL
void setLuxel(vec2 texCoord, vec3 color) {
	ivec2 uv = ivec2(texCoord * u_lightmapSize);
	if (uv.x > 0 && uv.y > 0) {
		// Should be imageAtomicAdd, but that only works with ints.
		vec4 current = imageLoad(u_lightmapSampler, uv);
		// https://blog.demofox.org/2016/08/23/incremental-averaging/
		imageStore(u_lightmapSampler, uv, vec4(mix(current.rgb, color.rgb, 1.0 / (current.a + 1.0)), current.a + 1.0));
	}
}
#endif

struct Node
{
	vec3 color;
	float depth;
	vec3 normal;
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
			nodes[numNodes].normal.x = uintBitsToFloat(normal_depth.x);
			nodes[numNodes].normal.y = uintBitsToFloat(normal_depth.y);
			nodes[numNodes].normal.z = uintBitsToFloat(normal_depth.z);
			nodes[numNodes].depth = uintBitsToFloat(normal_depth.w);
			nodes[numNodes].texcoord = vec2(uintBitsToFloat(texcoord.x), uintBitsToFloat(texcoord.y));
			offset = color_offset.w;
			numNodes++;
		}
		// need at least 2 nodes to transfer radiance
		if (numNodes <= 1)
			return;
		for (uint i = 0u; i < numNodes; i++) {
			for (uint j = i + 1u; j < numNodes; j++) {
				if (nodes[i].depth > nodes[j].depth || (nodes[i].depth == nodes[j].depth && dot(nodes[i].normal, u_rayNormal.xyz) < 0.0)) {
					Node temp = nodes[i];
					nodes[i] = nodes[j];
					nodes[j] = temp;
				}
			}
		}
		float brdf = 1.0;
		for (uint j = 0; j < numNodes - 1; j += 2) {
			// n1 to n2
			bool n2Forward = dot(nodes[j + 1].normal, u_rayNormal.xyz) > 0;
			float d = dot(nodes[j + 1].normal, n2Forward ? u_rayNormal.xyz : -u_rayNormal.xyz);
			if (d > 0)
				setLuxel(nodes[j + 1].texcoord, brdf * nodes[j + 0].color * d);
			// n2 to n1
			bool n1Forward = dot(nodes[j + 0].normal, u_rayNormal.xyz) > 0;
			d = dot(nodes[j + 0].normal, n1Forward ? u_rayNormal.xyz : -u_rayNormal.xyz);
			if (d > 0)
				setLuxel(nodes[j + 0].texcoord, brdf * nodes[j + 1].color * d);
		}
		/*for (uint j = 0u; j < numNodes; j++) {
			imageStore(u_lightmapSampler, ivec2(nodes[j].texcoord * u_lightmapSize), vec4(nodes[j].color, 1.0));
		}*/
	}
#endif
}
