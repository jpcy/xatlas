#include <bgfx_compute.sh>

FRAMEBUFFER_UIMAGE2D_RW(s_atomicCounter, r32ui, 0);

void main()
{
    imageStore(s_atomicCounter, ivec2(gl_FragCoord.xy), ivec4(0, 0, 0, 0));
}
