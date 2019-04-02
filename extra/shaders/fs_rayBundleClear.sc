#include <bgfx_compute.sh>

FRAMEBUFFER_UIMAGE2D_RW(s_rayBundleHeader, r32ui, 0);

void main()
{
	imageStore(s_rayBundleHeader, ivec2(gl_FragCoord.xy), ivec4(0xffffffff, 0, 0, 0));
}
