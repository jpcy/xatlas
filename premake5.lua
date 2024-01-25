newoption
{
	trigger = "asan",
	description = "Enable Clang AddressSanitizer"
}

newoption
{
	trigger = "msan",
	description = "Enable Clang MemorySanitizer"
}

newoption
{
	trigger = "ubsan",
	description = "Enable Clang UndefinedBehaviorSanitizer"
}

dofile("source/viewer/shaders.lua")

if _ACTION == nil then
	return
end

local asanEnabled = false
local msanEnabled = false
local ubsanEnabled = false
if _ACTION == "gmake" and _OPTIONS["cc"] == "clang" then
    if _OPTIONS["asan"] then
        asanEnabled = true
    	sanitizerEnabled = true
    end
	if _OPTIONS["msan"] then
        msanEnabled = true
    	sanitizerEnabled = true
    end
    if _OPTIONS["ubsan"] then
        ubsanEnabled = true
    	sanitizerEnabled = true
    end
end
local sanitizerEnabled = asanEnabled or msanEnabled or ubsanEnabled

function sanitizer()
	if asanEnabled then
		buildoptions { "-fsanitize=address" }
		linkoptions { "-fsanitize=address" }
	end
	if msanEnabled then
		buildoptions { "-fsanitize=memory" }
		linkoptions { "-fsanitize=memory" }
	end
    if ubsanEnabled then
		buildoptions { "-fsanitize=undefined" }
		linkoptions { "-fsanitize=undefined" }
	end
    if sanitizerEnabled then
		buildoptions { "-fno-omit-frame-pointer" }
	end
end

solution "xatlas"
	configurations { "Release", "Debug" }
	if _OPTIONS["cc"] ~= nil then
		location(path.join("build", _ACTION) .. "_" .. _OPTIONS["cc"])
	else
		location(path.join("build", _ACTION))
	end
	platforms { "x86_64", "x86" }
	startproject "viewer"
	filter "platforms:x86"
		architecture "x86"
	filter "platforms:x86_64"
		architecture "x86_64"
	filter "configurations:Debug*"
		defines { "_DEBUG" }
		optimize "Debug"
		symbols "On"
	filter "configurations:Release"
		defines "NDEBUG"
		optimize "Full"
	filter {}
	if sanitizerEnabled then
		optimize "Off"
		symbols "On"
	end
    linkoptions { "-fopenmp" }
	sanitizer()
	
local XATLAS_DIR = "source/xatlas"
		
project "xatlas_static"
	kind "StaticLib"
	language "C++"
	cppdialect "C++11"
	exceptionhandling "Off"
	rtti "Off"
	openmp "On"
	warnings "Extra"
	enablewarnings { "shadow" }
	files { path.join(XATLAS_DIR, "xatlas.*") }
	sanitizer()
	
project "xatlas"
	kind "SharedLib"
	language "C++"
	cppdialect "C++11"
	exceptionhandling "Off"
	rtti "Off"
	openmp "On"
	warnings "Extra"
	enablewarnings { "shadow" }
	defines { "XATLAS_C_API=1", "XATLAS_EXPORT_API=1" }
	files { path.join(XATLAS_DIR, "*.*") }
	sanitizer()
	
local THIRDPARTY_DIR = "source/thirdparty"
local BGFX_DIR = path.join(THIRDPARTY_DIR, "bgfx")
local BIMG_DIR = path.join(THIRDPARTY_DIR, "bimg")
local BX_DIR = path.join(THIRDPARTY_DIR, "bx")
local EIGEN_DIR = path.join(THIRDPARTY_DIR, "eigen")
local EMBREE_DIR = path.join(THIRDPARTY_DIR, "embree3")
local ENKITS_DIR = path.join(THIRDPARTY_DIR, "enkiTS")
local GLFW_DIR = path.join(THIRDPARTY_DIR, "glfw")
local IGL_DIR = path.join(THIRDPARTY_DIR, "libigl")
local MIMALLOC_DIR = path.join(THIRDPARTY_DIR, "mimalloc")
local NATIVEFILEDIALOG_DIR = path.join(THIRDPARTY_DIR, "nativefiledialog")
local OIDN_DIR = path.join(THIRDPARTY_DIR, "oidn")
local OPENFBX_DIR = path.join(THIRDPARTY_DIR, "OpenFBX")

project "test"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++11"
	exceptionhandling "On"
	rtti "Off"
	openmp "On"
	warnings "Extra"
	sanitizer()
	includedirs { XATLAS_DIR, THIRDPARTY_DIR }
	files "source/test/test.cpp"
	links { "tiny_obj_loader", "xatlas_static" }
	filter "action:vs*"
		files "source/xatlas.natvis"
	filter "system:linux"
		links { "pthread" }

project "viewer"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++14"
	exceptionhandling "Off"
	rtti "Off"
	warnings "Extra"
	sanitizer()
	files { "source/viewer/*.*", "source/viewer/shaders/*.*" }
	includedirs
	{
		XATLAS_DIR,
		THIRDPARTY_DIR,
		path.join(BGFX_DIR, "include"),
		path.join(BX_DIR, "include"),
		EIGEN_DIR,
		EMBREE_DIR,
		ENKITS_DIR,
		path.join(GLFW_DIR, "include"),
		path.join(IGL_DIR, "include"),
		path.join(MIMALLOC_DIR, "include"),
		path.join(NATIVEFILEDIALOG_DIR, "include"),
		path.join(OIDN_DIR, "include"),
		OPENFBX_DIR
	}
	links { "bgfx", "bimg", "bx", "cgltf", "enkiTS", "glfw", "imgui", "mimalloc", "nativefiledialog", "objzero", "OpenFBX", "stb_image", "stb_image_resize", "xatlas_static" }
	filter "system:windows"
		links { "bcrypt", "gdi32", "ole32", "psapi", "uuid"}
	filter "system:linux"
		links { "dl", "GL", "gtk-3", "gobject-2.0", "glib-2.0", "pthread", "X11", "Xcursor", "Xinerama", "Xrandr" }
	filter "action:vs*"
		files "source/xatlas.natvis"
		includedirs { path.join(BX_DIR, "include/compat/msvc") }
	filter { "system:windows", "action:gmake" }
		includedirs { path.join(BX_DIR, "include/compat/mingw") }
		
group "examples"
local EXAMPLES_DIR = "source/examples"

project "example"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++11"
	exceptionhandling "Off"
	rtti "Off"
	warnings "Extra"
	sanitizer()
	files { path.join(EXAMPLES_DIR, "example.cpp") }
	includedirs { XATLAS_DIR, THIRDPARTY_DIR }
	links { "stb_image_write", "tiny_obj_loader", "xatlas_static" }
	filter "action:vs*"
		files "source/xatlas.natvis"
	filter "system:linux"
		links { "pthread" }
		
project "example_c99"
	kind "ConsoleApp"
	language "C"
	cdialect "C99"
	warnings "Extra"
	sanitizer()
	files { path.join(EXAMPLES_DIR, "example_c99.c") }
	includedirs { XATLAS_DIR, THIRDPARTY_DIR }
	links { "objzero", "xatlas" }
	filter "system:linux"
		links { "m", "pthread" }
		
project "example_uvmesh"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++11"
	exceptionhandling "Off"
	rtti "Off"
	warnings "Extra"
	sanitizer()
	files { path.join(EXAMPLES_DIR, "example_uvmesh.cpp") }
	includedirs { XATLAS_DIR, THIRDPARTY_DIR }
	links { "stb_image_write", "tiny_obj_loader", "xatlas_static" }
	filter "action:vs*"
		files "source/xatlas.natvis"
	filter "system:linux"
		links { "pthread" }

group "thirdparty"

project "bgfx"
	kind "StaticLib"
	language "C++"
	cppdialect "C++14"
	exceptionhandling "Off"
	rtti "Off"
	sanitizer()
	defines	{ "__STDC_FORMAT_MACROS" }
	files
	{
		path.join(BGFX_DIR, "include/bgfx/**.h"),
		path.join(BGFX_DIR, "src/*.cpp"),
		path.join(BGFX_DIR, "src/*.h")
	}
	excludes
	{
		path.join(BGFX_DIR, "src/amalgamated.cpp")
	}
	includedirs
	{
		path.join(BX_DIR, "include"),
		path.join(BIMG_DIR, "include"),
		path.join(BIMG_DIR, "3rdparty"),
		path.join(BIMG_DIR, "3rdparty/astc-codec/include"),
		path.join(BIMG_DIR, "3rdparty/iqa/include"),
		path.join(BGFX_DIR, "include"),
		path.join(BGFX_DIR, "3rdparty"),
		path.join(BGFX_DIR, "3rdparty/dxsdk/include"),
		path.join(BGFX_DIR, "3rdparty/khronos")
	}
	filter "configurations:Debug"
		defines "BGFX_CONFIG_DEBUG=1"
	filter "action:vs*"
		defines { "_CRT_SECURE_NO_WARNINGS" }
		includedirs { path.join(BX_DIR, "include/compat/msvc") }
		excludes
		{
			path.join(BGFX_DIR, "src/glcontext_glx.cpp"),
			path.join(BGFX_DIR, "src/glcontext_egl.cpp")
		}
	filter { "system:windows", "action:gmake" }
		includedirs { path.join(BX_DIR, "include/compat/mingw") }
		
project "bimg"
	kind "StaticLib"
	language "C++"
	cppdialect "C++14"
	exceptionhandling "Off"
	rtti "Off"
	sanitizer()
	defines
	{
		"BIMG_DECODE_ENABLE=0"
	}
	files
	{
		path.join(BIMG_DIR, "include/bimg/*.h"),
		path.join(BIMG_DIR, "src/*.cpp"),
		path.join(BIMG_DIR, "src/*.h")
	}
	includedirs
	{
		path.join(BX_DIR, "include"),
		path.join(BIMG_DIR, "include")
	}
	filter "action:vs*"
		defines { "_CRT_SECURE_NO_WARNINGS" }
		includedirs { path.join(BX_DIR, "include/compat/msvc") }
	filter { "system:windows", "action:gmake" }
		includedirs { path.join(BX_DIR, "include/compat/mingw") }

project "bx"
	kind "StaticLib"
	language "C++"
	cppdialect "C++14"
	exceptionhandling "Off"
	rtti "Off"
	sanitizer()
	defines	{ "__STDC_FORMAT_MACROS" }
	files
	{
		path.join(BX_DIR, "include/bx/*.h"),
		path.join(BX_DIR, "include/bx/inline/*.inl"),
		path.join(BX_DIR, "include/tinystl/*.h"),
		path.join(BX_DIR, "src/*.cpp")
	}
	excludes
	{
		path.join(BX_DIR, "src/amalgamated.cpp"),
		path.join(BX_DIR, "src/crtnone.cpp")
	}
	includedirs
	{
		path.join(BX_DIR, "3rdparty"),
		path.join(BX_DIR, "include")
	}
	filter "action:vs*"
		defines { "_CRT_SECURE_NO_WARNINGS" }
		includedirs { path.join(BX_DIR, "include/compat/msvc") }
	filter { "system:windows", "action:gmake" }
		includedirs { path.join(BX_DIR, "include/compat/mingw") }

project "cgltf"
	kind "StaticLib"
	language "C"
	sanitizer()
	files(path.join(THIRDPARTY_DIR, "cgltf.*"))
	filter "action:vs*"
		defines { "_CRT_SECURE_NO_WARNINGS" }
		
project "enkiTS"
	kind "StaticLib"
	language "C++"
	cppdialect "C++11"
	exceptionhandling "Off"
	rtti "Off"
	sanitizer()
	files(path.join(ENKITS_DIR, "*.*"))

project "glfw"
	kind "StaticLib"
	language "C"
	sanitizer()
	files
	{
		path.join(GLFW_DIR, "include/*.h"),
		path.join(GLFW_DIR, "src/context.c"),
		path.join(GLFW_DIR, "src/egl_context.c"),
		path.join(GLFW_DIR, "src/init.c"),
		path.join(GLFW_DIR, "src/input.c"),
		path.join(GLFW_DIR, "src/monitor.c"),
		path.join(GLFW_DIR, "src/osmesa_context.c"),
		path.join(GLFW_DIR, "src/vulkan.c"),
		path.join(GLFW_DIR, "src/window.c"),
	}
	includedirs { path.join(GLFW_DIR, "include") }
	filter "system:windows"
		defines "_GLFW_WIN32"
		files
		{
			path.join(GLFW_DIR, "src/win32_init.c"),
			path.join(GLFW_DIR, "src/win32_joystick.c"),
			path.join(GLFW_DIR, "src/win32_monitor.c"),
			path.join(GLFW_DIR, "src/win32_thread.c"),
			path.join(GLFW_DIR, "src/win32_time.c"),
			path.join(GLFW_DIR, "src/win32_window.c"),
			path.join(GLFW_DIR, "src/wgl_context.c")
		}
	filter "system:linux"
		defines "_GLFW_X11"
		files
		{
			path.join(GLFW_DIR, "src/glx_context.c"),
			path.join(GLFW_DIR, "src/linux*.c"),
			path.join(GLFW_DIR, "src/posix*.c"),
			path.join(GLFW_DIR, "src/x11*.c"),
			path.join(GLFW_DIR, "src/xkb*.c")
		}
	filter "action:vs*"
		defines { "_CRT_SECURE_NO_WARNINGS" }
	filter {}
	
project "imgui"
	kind "StaticLib"
	language "C++"
	exceptionhandling "Off"
	rtti "Off"
	sanitizer()
	files(path.join(THIRDPARTY_DIR, "imgui/*.*"))
	
project "mimalloc"
	kind "StaticLib"
	language "C"
	sanitizer()
	includedirs(path.join(MIMALLOC_DIR, "include"))
	files(path.join(MIMALLOC_DIR, "src/*.*"))
	excludes
	{
		path.join(MIMALLOC_DIR, "src/alloc-override*"),
		path.join(MIMALLOC_DIR, "src/page-queue.c"),
		path.join(MIMALLOC_DIR, "src/static.c")
	}
	
project "nativefiledialog"
	kind "StaticLib"
	language "C++"
	exceptionhandling "Off"
	rtti "Off"
	sanitizer()
	includedirs(path.join(NATIVEFILEDIALOG_DIR, "include"))
	files(path.join(NATIVEFILEDIALOG_DIR, "nfd_common.*"))
	filter "system:windows"
		files(path.join(NATIVEFILEDIALOG_DIR, "nfd_win.cpp"))
	filter "system:linux"
		files(path.join(NATIVEFILEDIALOG_DIR, "nfd_gtk.c"))
		buildoptions(os.outputof("pkg-config --cflags gtk+-3.0"))
	filter "action:vs*"
		defines { "_CRT_SECURE_NO_WARNINGS" }
	
project "objzero"
	kind "StaticLib"
	language "C"
	cdialect "C99"
	sanitizer()
	files(path.join(THIRDPARTY_DIR, "objzero/objzero.*"))
	
project "OpenFBX"
	kind "StaticLib"
	language "C++"
	cppdialect "C++14"
	exceptionhandling "Off"
	rtti "Off"
	sanitizer()
	files(path.join(OPENFBX_DIR, "*.*"))
	
project "stb_image"
	kind "StaticLib"
	language "C"
	sanitizer()
	files(path.join(THIRDPARTY_DIR, "stb_image.*"))
	
project "stb_image_resize"
	kind "StaticLib"
	language "C"
	sanitizer()
	files(path.join(THIRDPARTY_DIR, "stb_image_resize.*"))
	
project "stb_image_write"
	kind "StaticLib"
	language "C"
	sanitizer()
	files(path.join(THIRDPARTY_DIR, "stb_image_write.*"))
	
project "tiny_obj_loader"
	kind "StaticLib"
	language "C++"
	exceptionhandling "Off"
	rtti "Off"
	sanitizer()
	files(path.join(THIRDPARTY_DIR, "tiny_obj_loader.*"))

