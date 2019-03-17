#include <bgfx_compute.sh>

USAMPLER2D(u_rayBundleHeaderSampler, 0);
USAMPLER2D(u_rayBundleDataSampler, 1);

uniform vec4 u_rayBundleDataResolution;
#define u_dataResolution uint(u_rayBundleDataResolution.x)

ivec2 dataUv(uint offset, uint pixel)
{
	return ivec2((offset * 2u + pixel) % u_dataResolution, (offset * 2u + pixel) / u_dataResolution);
}

void main()
{
	ivec2 uv = ivec2(gl_FragCoord.xy);
	uint offset = texelFetch(u_rayBundleHeaderSampler, uv, 0).x;
	if (offset != 0xffffffff) {
		uvec4 color_offset = texelFetch(u_rayBundleDataSampler, dataUv(offset, 0u), 0);
		uvec4 normal_depth = texelFetch(u_rayBundleDataSampler, dataUv(offset, 1u), 0);
		gl_FragColor = vec4(uintBitsToFloat(color_offset.r), uintBitsToFloat(color_offset.g), uintBitsToFloat(color_offset.b), 1.0);
	} else {
		gl_FragColor = vec4_splat(0.0);
	}
}
