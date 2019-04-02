#include <bgfx_compute.sh>

FRAMEBUFFER_UIMAGE2D_RO(s_rayBundleLightmap, r32ui, 0);
FRAMEBUFFER_IMAGE2D_RW(s_lightmap, rgba32f, 1);

vec4 toGamma(vec4 v)
{
	return vec4(pow(abs(v.rgb), vec3_splat(1.0 / 2.2)), v.a);
}

void main()
{
	ivec2 uv = ivec2(gl_FragCoord.xy);
	uint r = imageLoad(s_rayBundleLightmap, ivec2(uv.x * 4u + 0u, uv.y)).r;
	uint g = imageLoad(s_rayBundleLightmap, ivec2(uv.x * 4u + 1u, uv.y)).r;
	uint b = imageLoad(s_rayBundleLightmap, ivec2(uv.x * 4u + 2u, uv.y)).r;
	uint a = imageLoad(s_rayBundleLightmap, ivec2(uv.x * 4u + 3u, uv.y)).r;
	if (a > 0u) {
		float sum = float(a);
		vec4 color = vec4(float(r) / 255.0 / sum, float(g) / 255.0 / sum, float(b) / 255.0 / sum, 1.0);
		imageStore(s_lightmap, uv, toGamma(color));
	}
}
