solution "atlas"
	configurations { "Release", "Debug" }
	location "build"
	
	if os.is64bit() and not os.is("windows") then
		platforms { "x86_64", "x86" }
	else
		platforms { "x86", "x86_64" }
	end
		
	startproject "atlas"
	
	configuration "platforms:x86"
		architecture "x86"
		
	configuration "platforms:x86_64"
		architecture "x86_64"
	
	configuration "Debug"
		optimize "Debug"
		defines { "_DEBUG" }
		flags "Symbols"
		
	configuration "Release"
		optimize "Full"
		defines "NDEBUG"
		
	configuration "Debug"
		targetdir "build/bin_debug"
				
	configuration "Release"
		targetdir "build/bin_release"
		
	configuration "vs*"
		defines { "_CRT_SECURE_NO_DEPRECATE" }
		
	configuration { "vs*", "x86_64" }
		defines { "_WIN64", "__WIN64__" }

project "atlas"
	kind "ConsoleApp"
	language "C++"
	rtti "Off"
	targetprefix ""
	
	files
	{
		"xatlas.h",
		"thekla/*.cpp"
	}
	
	includedirs
	{
		"extern/tinyobj",
	}
		
	configuration { "windows", "gmake", "x86" }
		gccprefix "i686-w64-mingw32-"
		
	configuration { "windows", "gmake", "x86_64" }
		gccprefix "x86_64-w64-mingw32-"
