$input v_normal

#include <bgfx_compute.sh>

FRAMEBUFFER_UIMAGE2D_RW(u_atomicCounterSampler, r32ui, 0);
FRAMEBUFFER_UIMAGE2D_RW(u_rayBundleHeaderSampler, r32ui, 1);
FRAMEBUFFER_UIMAGE2D_RW(u_rayBundleDataSampler, rgba32ui, 2);

uniform vec4 u_diffuse;
uniform vec4 u_emission;
uniform vec4 u_lightDir;
uniform vec4 u_rayBundleDataResolution;
#define u_dataResolution uint(u_rayBundleDataResolution.x)

ivec2 dataUv(uint offset, uint pixel)
{
	return ivec2((offset * 2u + pixel) % u_dataResolution, (offset * 2u + pixel) / u_dataResolution);
}

void main()
{
	ivec2 uv = ivec2(gl_FragCoord.xy);
	vec3 color = u_diffuse.rgb * (dot(v_normal, u_lightDir.xyz) * 0.5 + 0.5) + u_emission.rgb; // half lambert
#if BGFX_SHADER_LANGUAGE_GLSL
	uint newOffset = imageAtomicAdd(u_atomicCounterSampler, ivec2(0, 0), 1u);
	if (newOffset >= u_dataResolution * u_dataResolution * 2u) {
		discard;
		return;
	}
	uint oldOffset = imageAtomicExchange(u_rayBundleHeaderSampler, uv, newOffset);
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
	imageStore(u_rayBundleDataSampler, dataUv(newOffset, 0u), color_offset);
	imageStore(u_rayBundleDataSampler, dataUv(newOffset, 1u), normal_depth);
#endif
}
