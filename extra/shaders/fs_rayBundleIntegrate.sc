#include <bgfx_compute.sh>

FRAMEBUFFER_UIMAGE2D_RW(u_rayBundleHeaderSampler, r32ui, 1);

#if BGFX_SHADER_LANGUAGE_GLSL
FRAMEBUFFER_UIMAGE2D_RW(u_rayBundleDataSampler, rgba32ui, 2);
#else
USAMPLER2D(u_rayBundleDataSampler, 2);
#endif

FRAMEBUFFER_UIMAGE2D_RW(u_rayBundleLightmapSampler, r32ui, 3);

uniform vec4 u_lightmapSize_dataSize;
#define u_lightmapSize u_lightmapSize_dataSize.xy
#define u_dataSize uint(u_lightmapSize_dataSize.z)
uniform vec4 u_rayNormal;
uniform vec4 u_skyColor_enabled;
#define u_skyColor u_skyColor_enabled.rgb
#define u_skyEnabled uint(u_skyColor_enabled.w)

ivec2 rayBundleDataUv(uint offset, uint pixel)
{
	return ivec2((offset * 3u + pixel) % u_dataSize, (offset * 3u + pixel) / u_dataSize);
}

ivec2 rayBundleLightmapDataUv(vec2 uv, uint pixel)
{
	return ivec2(uint(uv.x * u_lightmapSize.x) * 4u + pixel, uint(uv.y * u_lightmapSize.y));
}

struct Node
{
	vec3 color;
	float depth;
	vec3 normal;
	vec2 texcoord;
};

#define MAX_NODES 64

void main()
{
	uint offset = imageLoad(u_rayBundleHeaderSampler, ivec2(gl_FragCoord.xy)).x;
	if (offset != 0xffffffff) {
		Node nodes[MAX_NODES];
		uint numNodes = 0u;
		while (offset != 0xffffffff && numNodes < MAX_NODES) {
#if BGFX_SHADER_LANGUAGE_GLSL
			uvec4 color_offset = imageLoad(u_rayBundleDataSampler, rayBundleDataUv(offset, 0u));
			uvec4 normal_depth = imageLoad(u_rayBundleDataSampler, rayBundleDataUv(offset, 1u));
			uvec4 texcoord = imageLoad(u_rayBundleDataSampler, rayBundleDataUv(offset, 2u));
#else
			uvec4 color_offset = texelFetch(u_rayBundleDataSampler, rayBundleDataUv(offset, 0u), 0);
			uvec4 normal_depth = texelFetch(u_rayBundleDataSampler, rayBundleDataUv(offset, 1u), 0);
			uvec4 texcoord = texelFetch(u_rayBundleDataSampler, rayBundleDataUv(offset, 2u), 0);
#endif
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
		vec3 nodeRadiance[MAX_NODES];
		for (uint i = 0u; i < numNodes; i++) {
			for (uint j = i + 1u; j < numNodes; j++) {
				if (nodes[i].depth > nodes[j].depth || (nodes[i].depth == nodes[j].depth && dot(nodes[i].normal, u_rayNormal.xyz) > 0.0)) {
					Node temp = nodes[i];
					nodes[i] = nodes[j];
					nodes[j] = temp;
				}
			}
			nodeRadiance[i] = vec3_splat(0.0);
		}
		// need at least 2 nodes to transfer radiance
		if (numNodes >= 2u) {
			float brdf = 1.0;
			for (uint i = 0u; i < numNodes - 1u; i++) {
				float n1cosTheta = dot(nodes[i + 0u].normal, u_rayNormal.xyz);
				float n2cosTheta = dot(nodes[i + 1u].normal, -u_rayNormal.xyz);
				if (n1cosTheta > 0.0 && n2cosTheta > 0.0) {
					float cosTheta = n1cosTheta * n2cosTheta;
					nodeRadiance[i + 0u] = brdf * nodes[i + 1u].color * cosTheta;
					nodeRadiance[i + 1u] = brdf * nodes[i + 0u].color * cosTheta;
				}
			}
		}
		if (u_skyEnabled != 0u && numNodes > 0u) {
			if (dot(nodes[0u].normal, -u_rayNormal.xyz) > 0.0)
				nodeRadiance[0u] = u_skyColor;
			if (dot(nodes[numNodes - 1u].normal, u_rayNormal.xyz) > 0.0)
				nodeRadiance[numNodes - 1u] = u_skyColor;
		}
		for (uint j = 0u; j < numNodes; j++) {
			vec3 color = nodeRadiance[j];
			vec2 uv = nodes[j].texcoord;
			if (uv.x > 0.0 && uv.y > 0.0) {
				imageAtomicAdd(u_rayBundleLightmapSampler, rayBundleLightmapDataUv(uv, 0u), uint(color.r * 255.0));
				imageAtomicAdd(u_rayBundleLightmapSampler, rayBundleLightmapDataUv(uv, 1u), uint(color.g * 255.0));
				imageAtomicAdd(u_rayBundleLightmapSampler, rayBundleLightmapDataUv(uv, 2u), uint(color.b * 255.0));
				imageAtomicAdd(u_rayBundleLightmapSampler, rayBundleLightmapDataUv(uv, 3u), 1u);
			}
		}
	}
}
