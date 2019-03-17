#include <bgfx_compute.sh>

FRAMEBUFFER_UIMAGE2D_RW(u_rayBundleHeaderSampler, r32ui, 1);

void main()
{
	imageStore(u_rayBundleHeaderSampler, ivec2(gl_FragCoord.xy), ivec4(0xffffffff, 0, 0, 0));
}
