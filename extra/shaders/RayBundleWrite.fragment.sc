$input v_normal

#include <bgfx_compute.sh>

FRAMEBUFFER_UIMAGE2D_RW(u_atomicCounterSampler, r32ui, 0);
FRAMEBUFFER_IMAGE2D_RW(u_rayBundleColorSampler, rgba8, 1);

uniform vec4 u_diffuse;
uniform vec4 u_emission;
uniform vec4 u_lightDir;

void main()
{
	ivec2 uv = ivec2(gl_FragCoord.xy);
	vec3 color = u_diffuse.rgb * (dot(v_normal, u_lightDir.xyz) * 0.5 + 0.5) + u_emission.rgb; // half lambert
/*#if BGFX_SHADER_LANGUAGE_GLSL
	uint value = imageAtomicAdd(u_atomicCounterSampler, ivec2(0, 0), 1u);
	imageStore(u_rayBundleColorSampler, uv, vec4(float(value % 1000u) / 1000.0, 0.0, 0.0, 1.0));
#else*/
	imageStore(u_rayBundleColorSampler, uv, vec4(color, 1.0));
//#endif
}
