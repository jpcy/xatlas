#include <bgfx_compute.sh>

FRAMEBUFFER_IMAGE2D_RW(u_lightmap0Sampler, rgba32f, 1);
FRAMEBUFFER_IMAGE2D_RW(u_lightmap1Sampler, rgba32f, 2);

void main()
{
#if BGFX_SHADER_LANGUAGE_GLSL
	ivec2 uv = ivec2(gl_FragCoord.xy);
	vec4 color = imageLoad(u_lightmap0Sampler, uv);
	if (color.a > 0.0)
		imageStore(u_lightmap1Sampler, uv, vec4(color.rgb / color.a, 1.0));
#endif
}
