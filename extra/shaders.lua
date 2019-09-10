local BASE_PATH = path.getabsolute("..")
local SHADERS_DIR = path.join(BASE_PATH, "extra", "shaders")
local SHADERS_BIN_DIR = path.join(BASE_PATH, "extra", "shaders_bin")

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
		command = path.join(BASE_PATH, "bin", "shaderc.exe")
	elseif os.ishost("linux") then
		command = path.join(BASE_PATH, "bin", "shaderc")
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
			command = command .. "5_0"
		end
		command = command .. " -O 3 --Werror"
	elseif args.renderer == "vk" then
		command = command .. " --platform linux -p spirv"
	end
	os.execute(command)
end

newaction
{
	trigger = "shaders",
	description = "Compile shaders",
	onStart = function()
		local shaders =
		{
			"fs_blit",
			"fs_chart",
			"fs_color",
			"fs_gui",
			"fs_material",
			"fs_wireframe",
			"vs_blit",
			"vs_chart",
			"vs_color",
			"vs_gui",
			"vs_model",
			"vs_wireframe"
		}
		local renderers = nil
		if os.ishost("windows") then
			renderers = { "d3d11", "gl", "vk" }
		else
			renderers = { "gl", "vk" }
		end
		for _,renderer in pairs(renderers) do
			os.mkdir(path.join(SHADERS_BIN_DIR, renderer))
		end
		pcall(function()
			for _,shader in pairs(shaders) do
				for _,renderer in pairs(renderers) do
					io.write("Compiling " .. shader .. " " .. renderer .. "\n")
					io.flush()
					local shaderType = "vertex"
					if shader:sub(0, 2) == "fs" then shaderType = "fragment" end
					compileShader(
					{
						type = shaderType,
						renderer = renderer,
						inputFilename = path.join(SHADERS_DIR, shader) .. ".sc",
						includeDirs = path.join(BASE_PATH, "extra/thirdparty/bgfx/src"),
						varyingFilename = path.join(SHADERS_DIR, "varying.def.sc"),
						outputFilename = path.join(SHADERS_BIN_DIR, renderer, shader) .. ".h",
						bin2c = true,
						variableName = shader .. "_" .. renderer
					})
				end
			end
			-- Write a header file that includes all the shader headers.
			local filename = path.join(SHADERS_BIN_DIR, "shaders.h")
			io.write("Writing " .. filename .. "\n")
			io.flush()
			local file = assert(io.open(filename, "w"))
			for _,renderer in pairs(renderers) do
				for _,shader in pairs(shaders) do
					file:write(string.format("#include \"%s.h\"\n", path.join(renderer, shader)))
				end
			end
			file:close()
		end)
	end
}
