#include <bgfx_compute.sh>

FRAMEBUFFER_IMAGE2D_RW(u_lightmapSampler, rgba32f, 1);

void main()
{
	imageStore(u_lightmapSampler, ivec2(gl_FragCoord.xy), vec4_splat(0.0));
}
