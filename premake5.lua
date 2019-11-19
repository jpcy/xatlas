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

dofile("extra/shaders.lua")

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
	sanitizer()
		
project "xatlas"
	kind "StaticLib"
	language "C++"
	cppdialect "C++11"
	exceptionhandling "Off"
	rtti "Off"
	warnings "Extra"
	files { "xatlas.cpp", "xatlas.h" }
	sanitizer()
	
dofile("extra/projects.lua")
