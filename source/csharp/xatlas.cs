/*
MIT License

Copyright (c) 2018-2020 Jonathan Young

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
/*
thekla_atlas
MIT License
https://github.com/Thekla/thekla_atlas
Copyright (c) 2013 Thekla, Inc
Copyright NVIDIA Corporation 2006 -- Ignacio Castano <icastano@nvidia.com>
*/
using System;
using System.Runtime.InteropServices;

namespace Xatlas
{
	public static partial class xatlas
	{
		const string DllName = "xatlas.dll";
		
		public enum ChartType : uint
		{
			Planar,
			Ortho,
			LSCM,
			Piecewise,
			Invalid
		}
		
		public unsafe struct Chart
		{
			public uint *faceArray;
			public uint atlasIndex;
			public uint faceCount;
			public ChartType type;
			public uint material;
		}
		
		public unsafe struct Vertex
		{
			public int atlasIndex;
			public int chartIndex;
			public fixed float uv[2];
			public uint xref;
		}

		public unsafe struct Mesh
		{
			public Chart *chartArray;
			public uint *indexArray;
			public Vertex *vertexArray;
			public uint chartCount;
			public uint indexCount;
			public uint vertexCount;
		}
		
		const uint kImageChartIndexMask = 0x1FFFFFFF;
		const uint kImageHasChartIndexBit = 0x80000000;
		const uint kImageIsBilinearBit = 0x40000000;
		const uint kImageIsPaddingBit = 0x20000000;
		
		public unsafe struct Atlas
		{
			public uint *image;
			public Mesh *meshes;
			public float *utilization;
			public uint width;
			public uint height;
			public uint atlasCount;
			public uint chartCount;
			public uint meshCount;
			public float texelsPerUnit;
		}
		
		public enum IndexFormat : uint
		{
			UInt16,
			UInt32
		}
		
		public unsafe struct MeshDecl
		{
			public void *vertexPositionData;
			public void *vertexNormalData;
			public void *vertexUvData;
			public void *indexData;
			public bool *faceIgnoreData;
			public uint *faceMaterialData;
			public byte *faceVertexCount;
			public uint vertexCount;
			public uint vertexPositionStride;
			public uint vertexNormalStride;
			public uint vertexUvStride;
			public uint indexCount;
			public int indexOffset;
			public uint faceCount;
			public IndexFormat indexFormat;
			public float epsilon;
		}
		
		public enum AddMeshError : uint
		{
			Success,
			Error,
			IndexOutOfRange,
			InvalidFaceVertexCount,
			InvalidIndexCount
		}
		
		public unsafe struct UvMeshDecl
		{
			public void *vertexUvData;
			public void *indexData;
			public uint *faceMaterialData;
			public uint vertexCount;
			public uint vertexStride;
			public uint indexCount;
			public int indexOffset;
			public IndexFormat indexFormat;
		}

		public unsafe delegate void ParameterizeFunc(float *positions, float *texcoords, uint vertexCount, uint *indices, uint indexCount);

		public unsafe struct ChartOptions
		{
			public ParameterizeFunc paramFunc;
			public float maxChartArea;
			public float maxBoundaryLength;
			public float normalDeviationWeight;
			public float roundnessWeight;
			public float straightnessWeight;
			public float normalSeamWeight;
			public float textureSeamWeight;
			public float maxCost;
			public uint maxIterations;
			public bool useInputMeshUvs;
			public bool fixWinding;
		}

		public unsafe struct PackOptions
		{
			public uint maxChartSize;
			public uint padding;
			public float texelsPerUnit;
			public uint resolution;
			public bool bilinear;
			public bool blockAlign;
			public bool bruteForce;
			public bool createImage;
			public bool rotateChartsToAxis;
			public bool rotateCharts;
		}
		
		public enum ProgressCategory : uint
		{
			AddMesh,
			ComputeCharts,
			PackCharts,
			BuildOutputMeshes
		}

        public unsafe delegate bool ProgressFunc(ProgressCategory category, int progress, void *userData);
        //public unsafe delegate int PrintFunc(const char *, ...);
		
		[DllImport(DllName, EntryPoint="xatlasCreate", CallingConvention=CallingConvention.Cdecl)]
		public static extern unsafe Atlas *Create();
		
		[DllImport(DllName, EntryPoint="xatlasDestroy", CallingConvention=CallingConvention.Cdecl)]
		public static extern unsafe void Destroy(Atlas *atlas);
		
		[DllImport(DllName, EntryPoint="xatlasAddMesh", CallingConvention=CallingConvention.Cdecl)]
		public static extern unsafe AddMeshError AddMesh(Atlas *atlas, MeshDecl *meshDecl, uint meshCountHint);
		
		[DllImport(DllName, EntryPoint="xatlasAddMeshJoin", CallingConvention=CallingConvention.Cdecl)]
		public static extern unsafe void AddMeshJoin(Atlas *atlas);
		
		[DllImport(DllName, EntryPoint="xatlasAddUvMesh", CallingConvention=CallingConvention.Cdecl)]
		public static extern unsafe AddMeshError AddUvMesh(Atlas *atlas, UvMeshDecl *decl);
		
		[DllImport(DllName, EntryPoint="xatlasComputeCharts", CallingConvention=CallingConvention.Cdecl)]
		public static extern unsafe void ComputeCharts(Atlas *atlas, IntPtr chartOptions);
		
		[DllImport(DllName, EntryPoint="xatlasPackCharts", CallingConvention=CallingConvention.Cdecl)]
		public static extern unsafe void PackCharts(Atlas *atlas, PackOptions *packOptions);
		
		[DllImport(DllName, EntryPoint="xatlasGenerate", CallingConvention=CallingConvention.Cdecl)]
		public static extern unsafe void Generate(Atlas *atlas, IntPtr chartOptions, PackOptions *packOptions);

		[DllImport(DllName, EntryPoint="xatlasSetProgressCallback", CallingConvention=CallingConvention.Cdecl)]
		public static extern unsafe void SetProgressCallback(Atlas *atlas, ProgressFunc progressFunc, void *progressUserData);

		//[DllImport(DllName, EntryPoint="xatlasSetPrint", CallingConvention=CallingConvention.Cdecl)]
		//public static extern unsafe void SetPrint(PrintFunc print, bool verbose);

		[DllImport(DllName, EntryPoint="xatlasAddMeshErrorString", CallingConvention=CallingConvention.Cdecl)]
		public static extern unsafe IntPtr AddMeshErrorString(AddMeshError error);

		[DllImport(DllName, EntryPoint="xatlasProgressCategoryString", CallingConvention=CallingConvention.Cdecl)]
		public static extern unsafe IntPtr ProgressCategoryString(ProgressCategory category);

		[DllImport(DllName, EntryPoint="xatlasMeshDeclInit", CallingConvention=CallingConvention.Cdecl)]
		public static extern unsafe void MeshDeclInit(MeshDecl *meshDecl);

		[DllImport(DllName, EntryPoint="xatlasUvMeshDeclInit", CallingConvention=CallingConvention.Cdecl)]
		public static extern unsafe void UvMeshDeclInit(UvMeshDecl *uvMeshDecl);

		[DllImport(DllName, EntryPoint="xatlasChartOptionsInit", CallingConvention=CallingConvention.Cdecl)]
		public static extern unsafe void ChartOptionsInit(IntPtr chartOptions);

		[DllImport(DllName, EntryPoint="xatlasPackOptionsInit", CallingConvention=CallingConvention.Cdecl)]
		public static extern unsafe void PackOptionsInit(PackOptions *packOptions);
	}
}
