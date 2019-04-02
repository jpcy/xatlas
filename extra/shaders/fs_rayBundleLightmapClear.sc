#include <bgfx_compute.sh>

FRAMEBUFFER_UIMAGE2D_RW(s_rayBundleLightmap, r32ui, 0);

void main()
{
	ivec2 uv = ivec2(gl_FragCoord.xy);
	imageStore(s_rayBundleLightmap, ivec2(uv.x * 4u + 0u, uv.y), uvec4_splat(0u));
	imageStore(s_rayBundleLightmap, ivec2(uv.x * 4u + 1u, uv.y), uvec4_splat(0u));
	imageStore(s_rayBundleLightmap, ivec2(uv.x * 4u + 2u, uv.y), uvec4_splat(0u));
	imageStore(s_rayBundleLightmap, ivec2(uv.x * 4u + 3u, uv.y), uvec4_splat(0u));
}
