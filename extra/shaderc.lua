--[[
type
renderer
inputFilename
includeDirs
varyingFilename
outputFilename
bin2c
variableName
--]]
function compileShader(args)
	local command = nil
	if os.ishost("windows") then
		command = "bin\\shaderc.exe"
	elseif os.ishost("linux") then
		command = "`./bin/shaderc64"
	end
	command = command .. string.format(" -i \"%s\" -f \"%s\" -o \"%s\" --varyingdef \"%s\" --type %s", args.includeDirs, args.inputFilename, args.outputFilename, args.varyingFilename, args.type)
	if args.bin2c then
		command = command .. string.format(" --bin2c \"%s\"", args.variableName)
	end
	if args.renderer == "gl" then
		command = command .. " --platform linux -p 130"
	elseif args.renderer == "gles" then
		command = command .. " --platform asm.js"
	elseif args.renderer == "d3d9" or args.renderer == "d3d11" then
		command = command .. " --platform windows"
		if args.type == "fragment" then
			command = command .. " -p ps_"
		else
			command = command .. " -p vs_"
		end
		if args.renderer == "d3d9" then
			command = command .. "3_0"
		elseif args.renderer == "d3d11" then
			command = command .. "4_0"
		end
		command = command .. " -O 3 --Werror"
	end
	if os.ishost("linux") then
		command = command .. "`"
	end
	os.execute(command)
end
