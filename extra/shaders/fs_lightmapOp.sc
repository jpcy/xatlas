#include <bgfx_compute.sh>
#include "shared.h"

FRAMEBUFFER_IMAGE2D_RW(u_lightmapSampler, rgba32f, 0);
uniform vec4 u_lightmapOp;

void main()
{
	ivec2 uv = ivec2(gl_FragCoord.xy);
	uint op = uint(u_lightmapOp.x);
	if (op == LIGHTMAP_OP_CLEAR)
		imageStore(u_lightmapSampler, uv, vec4_splat(0.0));
}
