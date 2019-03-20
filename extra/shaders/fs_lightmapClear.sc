#include <bgfx_compute.sh>

FRAMEBUFFER_IMAGE2D_RW(u_lightmap0Sampler, rgba32f, 1);
FRAMEBUFFER_IMAGE2D_RW(u_lightmap1Sampler, rgba32f, 2);
FRAMEBUFFER_IMAGE2D_RW(u_lightmap2Sampler, rgba32f, 3);
uniform vec4 u_clearLightmaps;
#define u_clearLightmap0 uint(u_clearLightmaps.x)
#define u_clearLightmap1 uint(u_clearLightmaps.y)
#define u_clearLightmap2 uint(u_clearLightmaps.z)

void main()
{
	ivec2 uv = ivec2(gl_FragCoord.xy);
	if (u_clearLightmap0 != 0u)
		imageStore(u_lightmap0Sampler, uv, vec4_splat(0.0));
	if (u_clearLightmap1 != 0u)
		imageStore(u_lightmap1Sampler, uv, vec4_splat(0.0));
	if (u_clearLightmap2 != 0u)
		imageStore(u_lightmap2Sampler, uv, vec4(0.0, 0.0, 0.0, 1.0));
}
