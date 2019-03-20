#include <bgfx_compute.sh>

FRAMEBUFFER_UIMAGE2D_RW(u_rayBundleHeaderSampler, r32ui, 0);
FRAMEBUFFER_UIMAGE2D_RW(u_rayBundleDataSampler, rgba32ui, 1);
FRAMEBUFFER_IMAGE2D_RW(u_lightmap0Sampler, rgba32f, 2);

uniform vec4 u_lightmapSize_dataSize;
#define u_lightmapSize u_lightmapSize_dataSize.xy
#define u_dataSize uint(u_lightmapSize_dataSize.z)
uniform vec4 u_rayNormal;
uniform vec4 u_skyColor_enabled;
#define u_skyColor u_skyColor_enabled.rgb
#define u_skyEnabled uint(u_skyColor_enabled.w)

ivec2 dataUv(uint offset, uint pixel)
{
	return ivec2((offset * 3u + pixel) % u_dataSize, (offset * 3u + pixel) / u_dataSize);
}

#if BGFX_SHADER_LANGUAGE_GLSL
void setLuxel(vec2 texCoord, vec3 color) {
	ivec2 uv = ivec2(texCoord * u_lightmapSize);
	if (uv.x > 0 && uv.y > 0)
		imageStore(u_lightmap0Sampler, uv, vec4(color.rgb, 1.0));
}
#endif

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
		for (uint i = 0u; i < numNodes; i++) {
			for (uint j = i + 1u; j < numNodes; j++) {
				if (nodes[i].depth > nodes[j].depth || (nodes[i].depth == nodes[j].depth && dot(nodes[i].normal, u_rayNormal.xyz) < 0.0)) {
					Node temp = nodes[i];
					nodes[i] = nodes[j];
					nodes[j] = temp;
				}
			}
		}
		if (u_skyEnabled != 0u && numNodes > 0u) {
			float d = dot(nodes[0u].normal, u_rayNormal.xyz);
			if (d > 0.0)
				setLuxel(nodes[0u].texcoord, u_skyColor * d);
			if (numNodes > 1u) {
				float d = dot(nodes[numNodes - 1u].normal, -u_rayNormal.xyz);
				if (d > 0.0)
					setLuxel(nodes[numNodes - 1u].texcoord, u_skyColor * d);
			}
		}
		// need at least 2 nodes to transfer radiance
		else if (numNodes >= 2u) {
			float brdf = 1.0;
			for (uint j = 0u; j < numNodes - 1u; j += 2u) {
				// n1 to n2
				bool n2Forward = dot(nodes[j + 1u].normal, u_rayNormal.xyz) > 0.0;
				float d = dot(nodes[j + 1u].normal, n2Forward ? u_rayNormal.xyz : -u_rayNormal.xyz);
				if (d > 0.0)
					setLuxel(nodes[j + 1u].texcoord, brdf * nodes[j + 0u].color * d);
				// n2 to n1
				bool n1Forward = dot(nodes[j + 0u].normal, u_rayNormal.xyz) > 0.0;
				d = dot(nodes[j + 0u].normal, n1Forward ? u_rayNormal.xyz : -u_rayNormal.xyz);
				if (d > 0.0)
					setLuxel(nodes[j + 0u].texcoord, brdf * nodes[j + 1u].color * d);
			}
		}
	}
#endif
}
