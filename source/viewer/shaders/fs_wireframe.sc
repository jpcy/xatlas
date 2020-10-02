$input v_barycentric

/*
 * Copyright 2016 Dario Manesku. All rights reserved.
 * License: https://github.com/bkaradzic/bgfx#license-bsd-2-clause
 */

#include <bgfx_shader.sh>

uniform vec4 u_color;
uniform vec4 u_thickness;

void main()
{
	vec3  color   = u_color.rgb;
	float opacity = u_color.a;
	float thickness = u_thickness.x;

	vec3 fw = abs(dFdx(v_barycentric)) + abs(dFdy(v_barycentric));
	vec3 val = smoothstep(vec3_splat(0.0), fw*thickness, v_barycentric);
	float edge = min(min(val.x, val.y), val.z); // Gets to 0.0 around the edges.

	vec4 rgba = vec4(color, (1.0-edge)*opacity);
	gl_FragColor = rgba;
}
