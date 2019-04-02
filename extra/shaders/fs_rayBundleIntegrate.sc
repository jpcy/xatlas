#include <bgfx_compute.sh>
#include "rayBundleFragment.sh"

#define DEBUG_RAY_BUNDLE 0

FRAMEBUFFER_UIMAGE2D_RO(u_rayBundleHeaderSampler, r32ui, 0);

#if BGFX_SHADER_LANGUAGE_GLSL
FRAMEBUFFER_UIMAGE2D_RO(u_rayBundleDataSampler, rgba32ui, 1);
#else
USAMPLER2D(u_rayBundleDataSampler, 1);
#endif

FRAMEBUFFER_UIMAGE2D_RW(u_rayBundleLightmapSampler, r32ui, 2);

#if DEBUG_RAY_BUNDLE && BGFX_SHADER_LANGUAGE_GLSL
FRAMEBUFFER_IMAGE2D_RW(u_rayBundleDebugIntegrateSampler, rgba8, 3);
#endif

uniform vec4 u_lightmapSize_dataSize;
#define u_lightmapSize u_lightmapSize_dataSize.xy
#define u_dataSize uint(u_lightmapSize_dataSize.z)
uniform vec4 u_rayNormal;
uniform vec4 u_skyColor_enabled;
#define u_skyColor u_skyColor_enabled.rgb
#define u_skyEnabled uint(u_skyColor_enabled.w)

ivec2 rayBundleDataUv(uint offset, uint pixel)
{
	return ivec2((offset * 2u + pixel) % u_dataSize, (offset * 2u + pixel) / u_dataSize);
}

ivec2 rayBundleLightmapDataUv(vec2 uv, uint pixel)
{
	return ivec2(uint(uv.x) * 4u + pixel, uint(uv.y));
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
			RayBundleFragmentData fragData;
#if BGFX_SHADER_LANGUAGE_GLSL
			fragData.data0 = imageLoad(u_rayBundleDataSampler, rayBundleDataUv(offset, 0u));
			fragData.data1 = imageLoad(u_rayBundleDataSampler, rayBundleDataUv(offset, 1u));
#else
			fragData.data0 = texelFetch(u_rayBundleDataSampler, rayBundleDataUv(offset, 0u), 0);
			fragData.data1 = texelFetch(u_rayBundleDataSampler, rayBundleDataUv(offset, 1u), 0);
#endif
			RayBundleFragment frag = decodeRayBundleFragment(fragData);
			nodes[numNodes].color = frag.color;
			nodes[numNodes].normal = frag.normal;
			nodes[numNodes].depth = frag.depth;
			nodes[numNodes].texcoord = frag.texcoord;
#if BGFX_SHADER_LANGUAGE_GLSL
#if DEBUG_RAY_BUNDLE
			//imageStore(u_rayBundleDebugIntegrateSampler, ivec2(nodes[numNodes].texcoord.x, nodes[numNodes].texcoord.y), vec4(1.0, 0.0, 1.0, 1.0));
			//imageStore(u_rayBundleDebugIntegrateSampler, ivec2(gl_FragCoord.xy), vec4(1.0, 0.0, 1.0, 1.0));
#endif
#endif
			offset = frag.offset;
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
#if BGFX_SHADER_LANGUAGE_GLSL
#if DEBUG_RAY_BUNDLE
				imageStore(u_rayBundleDebugIntegrateSampler, ivec2(uv.x, uv.y), vec4(1.0, 0.0, 1.0, 1.0));		
#endif
#endif
			}
		}
	}
}
