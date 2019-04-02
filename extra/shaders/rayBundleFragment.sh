struct RayBundleFragment
{
	vec3 color;
	uint offset;
	vec3 normal;
	float depth;
	vec2 texcoord;
};

struct RayBundleFragmentData
{
	uvec4 data0;
	uvec4 data1;
};

RayBundleFragmentData encodeRayBundleFragment(vec3 color, uint offset, vec3 normal, float depth, vec2 texcoord)
{
	RayBundleFragmentData data;
	data.data0.x = packHalf2x16(color.rg);
	data.data0.y = packHalf2x16(vec2(color.b, depth));
	data.data0.z = offset;
	data.data0.w = packHalf2x16(normal.xy);
	data.data1.x = packHalf2x16(vec2(normal.z, 0.0));
	data.data1.y = packHalf2x16(texcoord);
	data.data1.z = 0u;
	data.data1.w = 0u;
	return data;
}

RayBundleFragment decodeRayBundleFragment(RayBundleFragmentData data)
{
	RayBundleFragment frag;
	frag.color.rg = unpackHalf2x16(data.data0.x);
	vec2 temp = unpackHalf2x16(data.data0.y);
	frag.color.b = temp.x;
	frag.depth = temp.y;
	frag.offset = data.data0.z;
	frag.normal.xy = unpackHalf2x16(data.data0.w);
	frag.normal.z = unpackHalf2x16(data.data1.x).x;
	frag.texcoord = unpackHalf2x16(data.data1.y);
	return frag;
}
