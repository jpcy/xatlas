@echo off
premake5.exe --file=../premake5.lua gmake
premake5.exe --file=../premake5.lua --cc=clang gmake
premake5.exe --file=../premake5.lua vs2015
premake5.exe --file=../premake5.lua vs2017
premake5.exe --file=../premake5.lua vs2019
pause
