$input v_normal, v_texcoord0

#include <bgfx_compute.sh>

#define DEBUG_RAY_BUNDLE 0

FRAMEBUFFER_UIMAGE2D_RW(u_atomicCounterSampler, r32ui, 0);
FRAMEBUFFER_UIMAGE2D_RW(u_rayBundleHeaderSampler, r32ui, 1);
FRAMEBUFFER_UIMAGE2D_RW(u_rayBundleDataSampler, rgba32ui, 2);

#if DEBUG_RAY_BUNDLE
FRAMEBUFFER_IMAGE2D_RW(u_rayBundleDebugWriteSampler, rgba8, 3);
#endif

uniform vec4 u_diffuse;
uniform vec4 u_emission;
uniform vec4 u_lightmapSize_dataSize;
#define u_dataSize uint(u_lightmapSize_dataSize.z)

ivec2 rayBundleDataUv(uint offset, uint pixel)
{
	return ivec2((offset * 3u + pixel) % u_dataSize, (offset * 3u + pixel) / u_dataSize);
}

void main()
{
	vec3 color = u_emission.rgb;
	uint newOffset = imageAtomicAdd(u_atomicCounterSampler, ivec2(0, 0), 1u);
	if (newOffset >= u_dataSize * u_dataSize * 3u) {
		discard;
		return;
	}
	uint oldOffset = imageAtomicExchange(u_rayBundleHeaderSampler, ivec2(gl_FragCoord.xy), newOffset);
	uvec4 color_offset;
	color_offset.x = floatBitsToUint(color.r);
	color_offset.y = floatBitsToUint(color.g);
	color_offset.z = floatBitsToUint(color.b);
	color_offset.w = oldOffset;
	uvec4 normal_depth;
	normal_depth.x = floatBitsToUint(v_normal.r);
	normal_depth.y = floatBitsToUint(v_normal.g);
	normal_depth.z = floatBitsToUint(v_normal.b);
	normal_depth.w = floatBitsToUint(gl_FragCoord.z);
	uvec4 texcoord;
	texcoord.x = floatBitsToUint(v_texcoord0.z);
	texcoord.y = floatBitsToUint(v_texcoord0.w);
	texcoord.z = 0u;
	texcoord.w = 0u;
	imageStore(u_rayBundleDataSampler, rayBundleDataUv(newOffset, 0u), color_offset);
	imageStore(u_rayBundleDataSampler, rayBundleDataUv(newOffset, 1u), normal_depth);
	imageStore(u_rayBundleDataSampler, rayBundleDataUv(newOffset, 2u), texcoord);
#if DEBUG_RAY_BUNDLE
	//imageStore(u_rayBundleDebugWriteSampler, ivec2(gl_FragCoord.xy), vec4(vec3_splat(float(newOffset) / (1024.0 * 1024.0)), 1.0));
	imageStore(u_rayBundleDebugWriteSampler, ivec2(gl_FragCoord.xy), vec4(1.0, 0.0, 1.0, 1.0));
	//imageStore(u_rayBundleDebugWriteSampler, ivec2(v_texcoord0.z * u_lightmapSize_dataSize.x, v_texcoord0.w * u_lightmapSize_dataSize.y), vec4(1.0, 0.0, 1.0, 1.0));
#endif
}
