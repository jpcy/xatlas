#include <bgfx_compute.sh>

FRAMEBUFFER_IMAGE2D_RW(u_rayBundleColorSampler, rgba8, 1);

void main()
{
	imageStore(u_rayBundleColorSampler, ivec2(gl_FragCoord.xy), vec4(0.0, 0.0, 0.0, 0.0));
}
