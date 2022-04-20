workspace "NucKage"
	configurations {
		"Release",
		"Debug"
	}

project "NucKage"
	kind "ConsoleApp"
	language "C++"
	targetdir "bin"
	objdir "objs"
	cppdialect "C++17"
	location "./"

	files {
		"src/**.cpp",
		"src/**.h"
	}

	--User specified path to ROOT CERN libraries--
	ROOTIncludepath = "/usr/include/root"
	ROOTLibpath = "/usr/lib64/root"

	includedirs {
		"src"
	}

	sysincludedirs {
		ROOTIncludepath
	}

	libdirs {
		ROOTLibpath
	}

	links {
		"Gui", "Core", "Imt", "RIO", "Net", "Hist", 
		"Graf", "Graf3d", "Gpad", "ROOTDataFrame", "ROOTVecOps",
		"Tree", "TreePlayer", "Rint", "Postscript", "Matrix",
		"Physics", "MathCore", "Thread", "MultiProc", "m", "dl"
	}

	filter "system:macosx or linux"
		linkoptions {
			"-pthread"
		}

	filter "configurations:Debug"
		symbols "On"

	filter "configurations:Release"
		optimize "On"