static const uint8_t fs_material_gl[1453] =
{
	0x46, 0x53, 0x48, 0x06, 0x91, 0xb0, 0x36, 0xc4, 0x00, 0x00, 0x00, 0x00, 0x07, 0x00, 0x09, 0x73, // FSH...6........s
	0x5f, 0x64, 0x69, 0x66, 0x66, 0x75, 0x73, 0x65, 0x00, 0x01, 0x00, 0x00, 0x01, 0x00, 0x0a, 0x73, // _diffuse.......s
	0x5f, 0x65, 0x6d, 0x69, 0x73, 0x73, 0x69, 0x6f, 0x6e, 0x00, 0x01, 0x00, 0x00, 0x01, 0x00, 0x0a, // _emission.......
	0x73, 0x5f, 0x6c, 0x69, 0x67, 0x68, 0x74, 0x6d, 0x61, 0x70, 0x00, 0x01, 0x00, 0x00, 0x01, 0x00, // s_lightmap......
	0x09, 0x75, 0x5f, 0x64, 0x69, 0x66, 0x66, 0x75, 0x73, 0x65, 0x02, 0x01, 0x00, 0x00, 0x01, 0x00, // .u_diffuse......
	0x0a, 0x75, 0x5f, 0x65, 0x6d, 0x69, 0x73, 0x73, 0x69, 0x6f, 0x6e, 0x02, 0x01, 0x00, 0x00, 0x01, // .u_emission.....
	0x00, 0x0a, 0x75, 0x5f, 0x6c, 0x69, 0x67, 0x68, 0x74, 0x44, 0x69, 0x72, 0x02, 0x01, 0x00, 0x00, // ..u_lightDir....
	0x01, 0x00, 0x18, 0x75, 0x5f, 0x73, 0x68, 0x61, 0x64, 0x65, 0x5f, 0x64, 0x69, 0x66, 0x66, 0x75, // ...u_shade_diffu
	0x73, 0x65, 0x5f, 0x65, 0x6d, 0x69, 0x73, 0x73, 0x69, 0x6f, 0x6e, 0x02, 0x01, 0x00, 0x00, 0x01, // se_emission.....
	0x00, 0x17, 0x05, 0x00, 0x00, 0x69, 0x6e, 0x20, 0x76, 0x65, 0x63, 0x33, 0x20, 0x76, 0x5f, 0x6e, // .....in vec3 v_n
	0x6f, 0x72, 0x6d, 0x61, 0x6c, 0x3b, 0x0a, 0x69, 0x6e, 0x20, 0x76, 0x65, 0x63, 0x34, 0x20, 0x76, // ormal;.in vec4 v
	0x5f, 0x74, 0x65, 0x78, 0x63, 0x6f, 0x6f, 0x72, 0x64, 0x30, 0x3b, 0x0a, 0x75, 0x6e, 0x69, 0x66, // _texcoord0;.unif
	0x6f, 0x72, 0x6d, 0x20, 0x73, 0x61, 0x6d, 0x70, 0x6c, 0x65, 0x72, 0x32, 0x44, 0x20, 0x73, 0x5f, // orm sampler2D s_
	0x64, 0x69, 0x66, 0x66, 0x75, 0x73, 0x65, 0x3b, 0x0a, 0x75, 0x6e, 0x69, 0x66, 0x6f, 0x72, 0x6d, // diffuse;.uniform
	0x20, 0x73, 0x61, 0x6d, 0x70, 0x6c, 0x65, 0x72, 0x32, 0x44, 0x20, 0x73, 0x5f, 0x65, 0x6d, 0x69, //  sampler2D s_emi
	0x73, 0x73, 0x69, 0x6f, 0x6e, 0x3b, 0x0a, 0x75, 0x6e, 0x69, 0x66, 0x6f, 0x72, 0x6d, 0x20, 0x73, // ssion;.uniform s
	0x61, 0x6d, 0x70, 0x6c, 0x65, 0x72, 0x32, 0x44, 0x20, 0x73, 0x5f, 0x6c, 0x69, 0x67, 0x68, 0x74, // ampler2D s_light
	0x6d, 0x61, 0x70, 0x3b, 0x0a, 0x75, 0x6e, 0x69, 0x66, 0x6f, 0x72, 0x6d, 0x20, 0x76, 0x65, 0x63, // map;.uniform vec
	0x34, 0x20, 0x75, 0x5f, 0x64, 0x69, 0x66, 0x66, 0x75, 0x73, 0x65, 0x3b, 0x0a, 0x75, 0x6e, 0x69, // 4 u_diffuse;.uni
	0x66, 0x6f, 0x72, 0x6d, 0x20, 0x76, 0x65, 0x63, 0x34, 0x20, 0x75, 0x5f, 0x65, 0x6d, 0x69, 0x73, // form vec4 u_emis
	0x73, 0x69, 0x6f, 0x6e, 0x3b, 0x0a, 0x75, 0x6e, 0x69, 0x66, 0x6f, 0x72, 0x6d, 0x20, 0x76, 0x65, // sion;.uniform ve
	0x63, 0x34, 0x20, 0x75, 0x5f, 0x6c, 0x69, 0x67, 0x68, 0x74, 0x44, 0x69, 0x72, 0x3b, 0x0a, 0x75, // c4 u_lightDir;.u
	0x6e, 0x69, 0x66, 0x6f, 0x72, 0x6d, 0x20, 0x76, 0x65, 0x63, 0x34, 0x20, 0x75, 0x5f, 0x73, 0x68, // niform vec4 u_sh
	0x61, 0x64, 0x65, 0x5f, 0x64, 0x69, 0x66, 0x66, 0x75, 0x73, 0x65, 0x5f, 0x65, 0x6d, 0x69, 0x73, // ade_diffuse_emis
	0x73, 0x69, 0x6f, 0x6e, 0x3b, 0x0a, 0x76, 0x6f, 0x69, 0x64, 0x20, 0x6d, 0x61, 0x69, 0x6e, 0x20, // sion;.void main 
	0x28, 0x29, 0x0a, 0x7b, 0x0a, 0x20, 0x20, 0x76, 0x65, 0x63, 0x34, 0x20, 0x63, 0x6f, 0x6c, 0x6f, // ().{.  vec4 colo
	0x72, 0x5f, 0x31, 0x3b, 0x0a, 0x20, 0x20, 0x76, 0x65, 0x63, 0x34, 0x20, 0x65, 0x6d, 0x69, 0x73, // r_1;.  vec4 emis
	0x73, 0x69, 0x6f, 0x6e, 0x5f, 0x32, 0x3b, 0x0a, 0x20, 0x20, 0x65, 0x6d, 0x69, 0x73, 0x73, 0x69, // sion_2;.  emissi
	0x6f, 0x6e, 0x5f, 0x32, 0x20, 0x3d, 0x20, 0x75, 0x5f, 0x65, 0x6d, 0x69, 0x73, 0x73, 0x69, 0x6f, // on_2 = u_emissio
	0x6e, 0x3b, 0x0a, 0x20, 0x20, 0x69, 0x66, 0x20, 0x28, 0x28, 0x75, 0x69, 0x6e, 0x74, 0x28, 0x75, // n;.  if ((uint(u
	0x5f, 0x73, 0x68, 0x61, 0x64, 0x65, 0x5f, 0x64, 0x69, 0x66, 0x66, 0x75, 0x73, 0x65, 0x5f, 0x65, // _shade_diffuse_e
	0x6d, 0x69, 0x73, 0x73, 0x69, 0x6f, 0x6e, 0x2e, 0x7a, 0x29, 0x20, 0x3d, 0x3d, 0x20, 0x31, 0x75, // mission.z) == 1u
	0x29, 0x29, 0x20, 0x7b, 0x0a, 0x20, 0x20, 0x20, 0x20, 0x65, 0x6d, 0x69, 0x73, 0x73, 0x69, 0x6f, // )) {.    emissio
	0x6e, 0x5f, 0x32, 0x2e, 0x78, 0x79, 0x7a, 0x20, 0x3d, 0x20, 0x74, 0x65, 0x78, 0x74, 0x75, 0x72, // n_2.xyz = textur
	0x65, 0x20, 0x28, 0x73, 0x5f, 0x65, 0x6d, 0x69, 0x73, 0x73, 0x69, 0x6f, 0x6e, 0x2c, 0x20, 0x76, // e (s_emission, v
	0x5f, 0x74, 0x65, 0x78, 0x63, 0x6f, 0x6f, 0x72, 0x64, 0x30, 0x2e, 0x78, 0x79, 0x29, 0x2e, 0x78, // _texcoord0.xy).x
	0x79, 0x7a, 0x3b, 0x0a, 0x20, 0x20, 0x7d, 0x3b, 0x0a, 0x20, 0x20, 0x63, 0x6f, 0x6c, 0x6f, 0x72, // yz;.  };.  color
	0x5f, 0x31, 0x20, 0x3d, 0x20, 0x76, 0x65, 0x63, 0x34, 0x28, 0x30, 0x2e, 0x30, 0x2c, 0x20, 0x30, // _1 = vec4(0.0, 0
	0x2e, 0x30, 0x2c, 0x20, 0x30, 0x2e, 0x30, 0x2c, 0x20, 0x31, 0x2e, 0x30, 0x29, 0x3b, 0x0a, 0x20, // .0, 0.0, 1.0);. 
	0x20, 0x69, 0x66, 0x20, 0x28, 0x28, 0x28, 0x28, 0x65, 0x6d, 0x69, 0x73, 0x73, 0x69, 0x6f, 0x6e, //  if ((((emission
	0x5f, 0x32, 0x2e, 0x78, 0x20, 0x3e, 0x20, 0x30, 0x2e, 0x30, 0x29, 0x20, 0x7c, 0x7c, 0x20, 0x28, // _2.x > 0.0) || (
	0x65, 0x6d, 0x69, 0x73, 0x73, 0x69, 0x6f, 0x6e, 0x5f, 0x32, 0x2e, 0x79, 0x20, 0x3e, 0x20, 0x30, // emission_2.y > 0
	0x2e, 0x30, 0x29, 0x29, 0x20, 0x7c, 0x7c, 0x20, 0x28, 0x65, 0x6d, 0x69, 0x73, 0x73, 0x69, 0x6f, // .0)) || (emissio
	0x6e, 0x5f, 0x32, 0x2e, 0x7a, 0x20, 0x3e, 0x20, 0x30, 0x2e, 0x30, 0x29, 0x29, 0x29, 0x20, 0x7b, // n_2.z > 0.0))) {
	0x0a, 0x20, 0x20, 0x20, 0x20, 0x63, 0x6f, 0x6c, 0x6f, 0x72, 0x5f, 0x31, 0x20, 0x3d, 0x20, 0x65, // .    color_1 = e
	0x6d, 0x69, 0x73, 0x73, 0x69, 0x6f, 0x6e, 0x5f, 0x32, 0x3b, 0x0a, 0x20, 0x20, 0x7d, 0x20, 0x65, // mission_2;.  } e
	0x6c, 0x73, 0x65, 0x20, 0x7b, 0x0a, 0x20, 0x20, 0x20, 0x20, 0x76, 0x65, 0x63, 0x34, 0x20, 0x64, // lse {.    vec4 d
	0x69, 0x66, 0x66, 0x75, 0x73, 0x65, 0x5f, 0x33, 0x3b, 0x0a, 0x20, 0x20, 0x20, 0x20, 0x64, 0x69, // iffuse_3;.    di
	0x66, 0x66, 0x75, 0x73, 0x65, 0x5f, 0x33, 0x20, 0x3d, 0x20, 0x75, 0x5f, 0x64, 0x69, 0x66, 0x66, // ffuse_3 = u_diff
	0x75, 0x73, 0x65, 0x3b, 0x0a, 0x20, 0x20, 0x20, 0x20, 0x69, 0x66, 0x20, 0x28, 0x28, 0x75, 0x69, // use;.    if ((ui
	0x6e, 0x74, 0x28, 0x75, 0x5f, 0x73, 0x68, 0x61, 0x64, 0x65, 0x5f, 0x64, 0x69, 0x66, 0x66, 0x75, // nt(u_shade_diffu
	0x73, 0x65, 0x5f, 0x65, 0x6d, 0x69, 0x73, 0x73, 0x69, 0x6f, 0x6e, 0x2e, 0x79, 0x29, 0x20, 0x3d, // se_emission.y) =
	0x3d, 0x20, 0x31, 0x75, 0x29, 0x29, 0x20, 0x7b, 0x0a, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x64, // = 1u)) {.      d
	0x69, 0x66, 0x66, 0x75, 0x73, 0x65, 0x5f, 0x33, 0x20, 0x3d, 0x20, 0x28, 0x75, 0x5f, 0x64, 0x69, // iffuse_3 = (u_di
	0x66, 0x66, 0x75, 0x73, 0x65, 0x20, 0x2a, 0x20, 0x74, 0x65, 0x78, 0x74, 0x75, 0x72, 0x65, 0x20, // ffuse * texture 
	0x28, 0x73, 0x5f, 0x64, 0x69, 0x66, 0x66, 0x75, 0x73, 0x65, 0x2c, 0x20, 0x76, 0x5f, 0x74, 0x65, // (s_diffuse, v_te
	0x78, 0x63, 0x6f, 0x6f, 0x72, 0x64, 0x30, 0x2e, 0x78, 0x79, 0x29, 0x29, 0x3b, 0x0a, 0x20, 0x20, // xcoord0.xy));.  
	0x20, 0x20, 0x7d, 0x3b, 0x0a, 0x20, 0x20, 0x20, 0x20, 0x69, 0x66, 0x20, 0x28, 0x28, 0x75, 0x69, //   };.    if ((ui
	0x6e, 0x74, 0x28, 0x75, 0x5f, 0x73, 0x68, 0x61, 0x64, 0x65, 0x5f, 0x64, 0x69, 0x66, 0x66, 0x75, // nt(u_shade_diffu
	0x73, 0x65, 0x5f, 0x65, 0x6d, 0x69, 0x73, 0x73, 0x69, 0x6f, 0x6e, 0x2e, 0x78, 0x29, 0x20, 0x3d, // se_emission.x) =
	0x3d, 0x20, 0x75, 0x69, 0x6e, 0x74, 0x28, 0x30, 0x29, 0x29, 0x29, 0x20, 0x7b, 0x0a, 0x20, 0x20, // = uint(0))) {.  
	0x20, 0x20, 0x20, 0x20, 0x63, 0x6f, 0x6c, 0x6f, 0x72, 0x5f, 0x31, 0x2e, 0x78, 0x79, 0x7a, 0x20, //     color_1.xyz 
	0x3d, 0x20, 0x28, 0x64, 0x69, 0x66, 0x66, 0x75, 0x73, 0x65, 0x5f, 0x33, 0x2e, 0x78, 0x79, 0x7a, // = (diffuse_3.xyz
	0x20, 0x2a, 0x20, 0x28, 0x28, 0x0a, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x64, 0x6f, //  * ((.        do
	0x74, 0x20, 0x28, 0x76, 0x5f, 0x6e, 0x6f, 0x72, 0x6d, 0x61, 0x6c, 0x2c, 0x20, 0x75, 0x5f, 0x6c, // t (v_normal, u_l
	0x69, 0x67, 0x68, 0x74, 0x44, 0x69, 0x72, 0x2e, 0x78, 0x79, 0x7a, 0x29, 0x0a, 0x20, 0x20, 0x20, // ightDir.xyz).   
	0x20, 0x20, 0x20, 0x20, 0x2a, 0x20, 0x30, 0x2e, 0x35, 0x29, 0x20, 0x2b, 0x20, 0x30, 0x2e, 0x35, //     * 0.5) + 0.5
	0x29, 0x29, 0x3b, 0x0a, 0x20, 0x20, 0x20, 0x20, 0x7d, 0x20, 0x65, 0x6c, 0x73, 0x65, 0x20, 0x7b, // ));.    } else {
	0x0a, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x69, 0x66, 0x20, 0x28, 0x28, 0x75, 0x69, 0x6e, 0x74, // .      if ((uint
	0x28, 0x75, 0x5f, 0x73, 0x68, 0x61, 0x64, 0x65, 0x5f, 0x64, 0x69, 0x66, 0x66, 0x75, 0x73, 0x65, // (u_shade_diffuse
	0x5f, 0x65, 0x6d, 0x69, 0x73, 0x73, 0x69, 0x6f, 0x6e, 0x2e, 0x78, 0x29, 0x20, 0x3d, 0x3d, 0x20, // _emission.x) == 
	0x31, 0x75, 0x29, 0x29, 0x20, 0x7b, 0x0a, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x63, // 1u)) {.        c
	0x6f, 0x6c, 0x6f, 0x72, 0x5f, 0x31, 0x2e, 0x78, 0x79, 0x7a, 0x20, 0x3d, 0x20, 0x28, 0x64, 0x69, // olor_1.xyz = (di
	0x66, 0x66, 0x75, 0x73, 0x65, 0x5f, 0x33, 0x2e, 0x78, 0x79, 0x7a, 0x20, 0x2a, 0x20, 0x74, 0x65, // ffuse_3.xyz * te
	0x78, 0x74, 0x75, 0x72, 0x65, 0x20, 0x28, 0x73, 0x5f, 0x6c, 0x69, 0x67, 0x68, 0x74, 0x6d, 0x61, // xture (s_lightma
	0x70, 0x2c, 0x20, 0x76, 0x5f, 0x74, 0x65, 0x78, 0x63, 0x6f, 0x6f, 0x72, 0x64, 0x30, 0x2e, 0x7a, // p, v_texcoord0.z
	0x77, 0x29, 0x2e, 0x78, 0x79, 0x7a, 0x29, 0x3b, 0x0a, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x7d, // w).xyz);.      }
	0x20, 0x65, 0x6c, 0x73, 0x65, 0x20, 0x7b, 0x0a, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, //  else {.        
	0x69, 0x66, 0x20, 0x28, 0x28, 0x75, 0x69, 0x6e, 0x74, 0x28, 0x75, 0x5f, 0x73, 0x68, 0x61, 0x64, // if ((uint(u_shad
	0x65, 0x5f, 0x64, 0x69, 0x66, 0x66, 0x75, 0x73, 0x65, 0x5f, 0x65, 0x6d, 0x69, 0x73, 0x73, 0x69, // e_diffuse_emissi
	0x6f, 0x6e, 0x2e, 0x78, 0x29, 0x20, 0x3d, 0x3d, 0x20, 0x32, 0x75, 0x29, 0x29, 0x20, 0x7b, 0x0a, // on.x) == 2u)) {.
	0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x63, 0x6f, 0x6c, 0x6f, 0x72, 0x5f, //           color_
	0x31, 0x2e, 0x78, 0x79, 0x7a, 0x20, 0x3d, 0x20, 0x74, 0x65, 0x78, 0x74, 0x75, 0x72, 0x65, 0x20, // 1.xyz = texture 
	0x28, 0x73, 0x5f, 0x6c, 0x69, 0x67, 0x68, 0x74, 0x6d, 0x61, 0x70, 0x2c, 0x20, 0x76, 0x5f, 0x74, // (s_lightmap, v_t
	0x65, 0x78, 0x63, 0x6f, 0x6f, 0x72, 0x64, 0x30, 0x2e, 0x7a, 0x77, 0x29, 0x2e, 0x78, 0x79, 0x7a, // excoord0.zw).xyz
	0x3b, 0x0a, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x7d, 0x3b, 0x0a, 0x20, 0x20, 0x20, // ;.        };.   
	0x20, 0x20, 0x20, 0x7d, 0x3b, 0x0a, 0x20, 0x20, 0x20, 0x20, 0x7d, 0x3b, 0x0a, 0x20, 0x20, 0x20, //    };.    };.   
	0x20, 0x63, 0x6f, 0x6c, 0x6f, 0x72, 0x5f, 0x31, 0x2e, 0x77, 0x20, 0x3d, 0x20, 0x64, 0x69, 0x66, //  color_1.w = dif
	0x66, 0x75, 0x73, 0x65, 0x5f, 0x33, 0x2e, 0x77, 0x3b, 0x0a, 0x20, 0x20, 0x7d, 0x3b, 0x0a, 0x20, // fuse_3.w;.  };. 
	0x20, 0x67, 0x6c, 0x5f, 0x46, 0x72, 0x61, 0x67, 0x43, 0x6f, 0x6c, 0x6f, 0x72, 0x20, 0x3d, 0x20, //  gl_FragColor = 
	0x63, 0x6f, 0x6c, 0x6f, 0x72, 0x5f, 0x31, 0x3b, 0x0a, 0x7d, 0x0a, 0x0a, 0x00,                   // color_1;.}...
};
