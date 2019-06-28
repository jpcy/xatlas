dofile("extra/shaders.lua")

if _ACTION == nil then
	return
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
		
project "xatlas"
	kind "StaticLib"
	language "C++"
	cppdialect "C++11"
	exceptionhandling "Off"
	rtti "Off"
	warnings "Extra"
	files { "xatlas.cpp", "xatlas.h" }
	
dofile("extra/projects.lua")
