solution "xatlas"
	configurations { "Release", "Debug" }
	if _OPTIONS["cc"] ~= nil then
		location(path.join("build", _ACTION) .. "_" .. _OPTIONS["cc"])
	else
		location(path.join("build", _ACTION))
	end
	if os.is64bit() and not os.istarget("windows") then
		platforms { "x86_64", "x86" }
	else
		platforms { "x86", "x86_64" }
	end
	startproject "example"
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
		
project "xatlas"
	kind "StaticLib"
	language "C++"
	cppdialect "C++11"
	exceptionhandling "Off"
	rtti "Off"
	warnings "Extra"
	files { "xatlas.cpp", "xatlas.h" }

project "example"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++11"
	exceptionhandling "Off"
	rtti "Off"
	warnings "Extra"
	files "extra/example.cpp"
	links { "stb_image_write", "tiny_obj_loader", "xatlas" }

project "test"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++11"
	exceptionhandling "Off"
	rtti "Off"
	warnings "Extra"
	files "extra/test.cpp"
	links { "tiny_obj_loader", "xatlas" }

local THIRDPARTY_DIR = "extra/thirdparty"
local GLFW_DIR = path.join(THIRDPARTY_DIR, "glfw")

project "viewer"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++11"
	exceptionhandling "Off"
	rtti "Off"
	warnings "Extra"
	files "extra/viewer.cpp"
	includedirs
	{
		THIRDPARTY_DIR,
		path.join(GLFW_DIR, "include")
	}
	links { "flextGL", "glfw", "HandmadeMath", "imgui", "nativefiledialog", "objzero", "xatlas" }
	filter "system:windows"
		links { "gdi32", "ole32", "opengl32", "uuid" }
	filter "system:linux"
		links { "dl", "GL", "gtk-3", "gobject-2.0", "glib-2.0", "pthread", "X11", "Xcursor", "Xinerama", "Xrandr" }

group "thirdparty"

project "flextGL"
	kind "StaticLib"
	language "C"
	files(path.join(THIRDPARTY_DIR, "flextGL.*"))
	filter "toolset:clang or gcc"
		buildoptions "-Wno-incompatible-pointer-types"

project "glfw"
	kind "StaticLib"
	language "C"
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
	
project "HandmadeMath"
	kind "StaticLib"
	language "C"
	files(path.join(THIRDPARTY_DIR, "HandmadeMath.*"))
	
project "imgui"
	kind "StaticLib"
	language "C++"
	exceptionhandling "Off"
	rtti "Off"
	files(path.join(THIRDPARTY_DIR, "imgui/*.*"))
	
project "nativefiledialog"
	kind "StaticLib"
	language "C++"
	exceptionhandling "Off"
	rtti "Off"
	files(path.join(THIRDPARTY_DIR, "nativefiledialog/nfd_common.*"))
	filter "system:windows"
		files(path.join(THIRDPARTY_DIR, "nativefiledialog/nfd_win.cpp"))
	filter "system:linux"
		files(path.join(THIRDPARTY_DIR, "nativefiledialog/nfd_gtk.c"))
		buildoptions(os.outputof("pkg-config --cflags gtk+-3.0"))
	filter "action:vs*"
		defines { "_CRT_SECURE_NO_WARNINGS" }
	
project "objzero"
	kind "StaticLib"
	language "C"
	cdialect "C99"
	files(path.join(THIRDPARTY_DIR, "objzero/objzero.*"))
	
project "stb_image_write"
	kind "StaticLib"
	language "C"
	files(path.join(THIRDPARTY_DIR, "stb_image_write.*"))
	
project "tiny_obj_loader"
	kind "StaticLib"
	language "C++"
	exceptionhandling "Off"
	rtti "Off"
	files(path.join(THIRDPARTY_DIR, "tiny_obj_loader.*"))
