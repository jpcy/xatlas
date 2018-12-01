local BUILD_DIR = path.join("build", _ACTION)

solution "xatlas"
	configurations { "Release", "Debug" }
	location(BUILD_DIR)
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
	filter { "configurations:Release" }
		defines "NDEBUG"
		optimize "Full"
		
project "xatlas"
	kind "StaticLib"
	language "C++"
	cppdialect "C++98"
	rtti "Off"
	warnings "Extra"
	files { "xatlas.cpp", "xatlas.h" }

project "example"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++11"
	rtti "Off"
	warnings "Extra"
	targetprefix ""
	files { "example/*.cpp" }
	links { "xatlas" }
	filter { "system:windows", "action:gmake", "platforms:x86" }
		gccprefix "i686-w64-mingw32-"
	filter { "system:windows", "action:gmake", "platforms:x86_64" }
		gccprefix "x86_64-w64-mingw32-"
