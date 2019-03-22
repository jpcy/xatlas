#include <bgfx_compute.sh>

FRAMEBUFFER_UIMAGE2D_RW(u_rayBundleLightmapSampler, r32ui, 1);
FRAMEBUFFER_IMAGE2D_RW(u_lightmapSampler, rgba32f, 2);

void main()
{
#if BGFX_SHADER_LANGUAGE_GLSL
	ivec2 uv = ivec2(gl_FragCoord.xy);
	uint r = imageLoad(u_rayBundleLightmapSampler, ivec2(uv.x * 4u + 0u, uv.y)).r;
	uint g = imageLoad(u_rayBundleLightmapSampler, ivec2(uv.x * 4u + 1u, uv.y)).r;
	uint b = imageLoad(u_rayBundleLightmapSampler, ivec2(uv.x * 4u + 2u, uv.y)).r;
	uint a = imageLoad(u_rayBundleLightmapSampler, ivec2(uv.x * 4u + 3u, uv.y)).r;
	if (a > 0u) {
		vec4 color = vec4(float(r) / 255.0 / float(a), float(g) / 255.0 / float(a), float(b) / 255.0 / float(a), 1.0);
		imageStore(u_lightmapSampler, uv, color);
	}
#endif
}
