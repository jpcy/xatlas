#include <bgfx_compute.sh>

FRAMEBUFFER_UIMAGE2D_RW(u_rayBundleLightmapSampler, r32ui, 0);
FRAMEBUFFER_IMAGE2D_RW(u_lightmapSampler, rgba32f, 1);
uniform vec4 u_clearLightmaps;
#define u_clearLightmap0 uint(u_clearLightmaps.x)
#define u_clearLightmap1 uint(u_clearLightmaps.y)

void main()
{
	ivec2 uv = ivec2(gl_FragCoord.xy);
	if (u_clearLightmap0 != 0u) {
		imageStore(u_rayBundleLightmapSampler, ivec2(uv.x * 4u + 0u, uv.y), uvec4_splat(0u));
		imageStore(u_rayBundleLightmapSampler, ivec2(uv.x * 4u + 1u, uv.y), uvec4_splat(0u));
		imageStore(u_rayBundleLightmapSampler, ivec2(uv.x * 4u + 2u, uv.y), uvec4_splat(0u));
		imageStore(u_rayBundleLightmapSampler, ivec2(uv.x * 4u + 3u, uv.y), uvec4_splat(0u));
	}
	if (u_clearLightmap1 != 0u)
		imageStore(u_lightmapSampler, uv, vec4_splat(0.0));
}
