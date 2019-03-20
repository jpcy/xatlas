#include <bgfx_compute.sh>

FRAMEBUFFER_IMAGE2D_RW(u_lightmap0Sampler, rgba32f, 1);
FRAMEBUFFER_IMAGE2D_RW(u_lightmap1Sampler, rgba32f, 2);

void main()
{
#if BGFX_SHADER_LANGUAGE_GLSL
	ivec2 uv = ivec2(gl_FragCoord.xy);
	vec4 color0 = imageLoad(u_lightmap0Sampler, uv);
	if (color0.a > 0.0) {
		vec4 color1 = imageLoad(u_lightmap1Sampler, uv);
		imageStore(u_lightmap1Sampler, uv, color0 + color1);
	}
#endif
}
