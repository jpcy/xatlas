#include <bgfx_compute.sh>

FRAMEBUFFER_IMAGE2D_RW(u_lightmapSampler, rgba32f, 1);

void main()
{
#if BGFX_SHADER_LANGUAGE_GLSL
	ivec2 uv = ivec2(gl_FragCoord.xy);
	vec4 color = imageLoad(u_lightmapSampler, uv);
	if (color.a <= 0.0)
		return;
	float ia = 1.0 / color.a;
	color.r *= ia;
	color.g *= ia;
	color.b *= ia;
	color.a = 1.0;
	imageStore(u_lightmapSampler, uv, color);
#endif
}
