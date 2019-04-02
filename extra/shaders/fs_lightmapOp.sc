#include <bgfx_compute.sh>
#include "shared.h"

FRAMEBUFFER_IMAGE2D_RW(s_lightmapCurrPass, rgba32f, 0);
FRAMEBUFFER_IMAGE2D_RW(s_lightmapPrevPass, rgba32f, 1);
FRAMEBUFFER_IMAGE2D_RW(s_lightmapSum, rgba32f, 2);
FRAMEBUFFER_IMAGE2D_RW(s_lightmap, rgba32f, 3);
uniform vec4 u_lightmapOp_bounceWeight;
#define u_op uint(u_lightmapOp_bounceWeight.x)
#define u_bounceWeight u_lightmapOp_bounceWeight.y

void main()
{
	ivec2 uv = ivec2(gl_FragCoord.xy);
	if (u_op == LIGHTMAP_OP_CLEAR_CURR) {
		imageStore(s_lightmapCurrPass, uv, vec4_splat(0.0));
	}
	else if (u_op == LIGHTMAP_OP_CLEAR_SUM) {
		imageStore(s_lightmapSum, uv, vec4_splat(0.0));
	}
	else if (u_op == LIGHTMAP_OP_WRITE_FINAL) {
		// final = sum + current
		vec4 sum = imageLoad(s_lightmapSum, uv);
		vec4 curr = imageLoad(s_lightmapCurrPass, uv) * u_bounceWeight;
		imageStore(s_lightmap, uv, sum + curr);
	}
	else if (u_op == LIGHTMAP_OP_FINISH_PASS) {
		// previous = current
		// sum += current
		vec4 curr = imageLoad(s_lightmapCurrPass, uv) * u_bounceWeight;
		imageStore(s_lightmapPrevPass, uv, curr);
		vec4 sum = imageLoad(s_lightmapSum, uv);
		imageStore(s_lightmapSum, uv, sum + curr);
	}
}
