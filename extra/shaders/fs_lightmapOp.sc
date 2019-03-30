#include <bgfx_compute.sh>
#include "shared.h"

FRAMEBUFFER_IMAGE2D_RW(u_lightmapCurrPassSampler, rgba32f, 0);
FRAMEBUFFER_IMAGE2D_RW(u_lightmapPrevPassSampler, rgba32f, 1);
FRAMEBUFFER_IMAGE2D_RW(u_lightmapSumSampler, rgba32f, 2);
FRAMEBUFFER_IMAGE2D_RW(u_lightmapSampler, rgba32f, 3);
uniform vec4 u_lightmapOp_bounceWeight;
#define u_op uint(u_lightmapOp_bounceWeight.x)
#define u_bounceWeight u_lightmapOp_bounceWeight.y

void main()
{
	ivec2 uv = ivec2(gl_FragCoord.xy);
	if (u_op == LIGHTMAP_OP_CLEAR_CURR) {
		imageStore(u_lightmapCurrPassSampler, uv, vec4_splat(0.0));
	}
	else if (u_op == LIGHTMAP_OP_CLEAR_SUM) {
		imageStore(u_lightmapSumSampler, uv, vec4_splat(0.0));
	}
	else if (u_op == LIGHTMAP_OP_WRITE_FINAL) {
		// final = sum + current
		vec4 sum = imageLoad(u_lightmapSumSampler, uv);
		vec4 curr = imageLoad(u_lightmapCurrPassSampler, uv) * u_bounceWeight;
		imageStore(u_lightmapSampler, uv, sum + curr);
	}
	else if (u_op == LIGHTMAP_OP_FINISH_PASS) {
		// previous = current
		// sum += current
		vec4 curr = imageLoad(u_lightmapCurrPassSampler, uv) * u_bounceWeight;
		imageStore(u_lightmapPrevPassSampler, uv, curr);
		vec4 sum = imageLoad(u_lightmapSumSampler, uv);
		imageStore(u_lightmapSumSampler, uv, sum + curr);
	}
}
