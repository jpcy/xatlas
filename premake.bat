@echo off
"bin/premake5.exe" gmake
"bin/premake5.exe" --cc=clang gmake
"bin/premake5.exe" vs2017
pause
