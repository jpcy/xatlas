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
	filter { "configurations:Release" }
		defines "NDEBUG"
		optimize "Full"
		
project "xatlas"
	kind "StaticLib"
	language "C++"
	cppdialect "C++98"
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
	files "example.cpp"
	links "xatlas"

project "test"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++11"
	exceptionhandling "Off"
	rtti "Off"
	warnings "Extra"
	files "test.cpp"
	links "xatlas"
