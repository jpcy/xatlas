@echo off
premake5.exe gmake
premake5.exe --cc=clang gmake
premake5.exe vs2017
pause
