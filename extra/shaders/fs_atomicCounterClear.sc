#include <bgfx_compute.sh>

FRAMEBUFFER_UIMAGE2D_RW(u_atomicCounterSampler, r32ui, 0);

void main()
{
    imageStore(u_atomicCounterSampler, ivec2(gl_FragCoord.xy), ivec4(0, 0, 0, 0));
}
