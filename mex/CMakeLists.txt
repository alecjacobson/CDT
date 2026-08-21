# Matlab mex bindings, built the way gptoolbox builds its own: find matlab,
# build each source as a shared library with the platform's mex suffix and no
# "lib" prefix, and drop the result next to the sources so that adding this
# directory to the matlab path is all it takes.

find_package(Matlab REQUIRED COMPONENTS MEX_COMPILER MX_LIBRARY)

# These are directory scoped, so they do not leak back into the parent build.
set(CMAKE_SHARED_LIBRARY_SUFFIX ".${Matlab_MEX_EXTENSION}")
set(CMAKE_SHARED_LIBRARY_PREFIX)
set(CMAKE_SHARED_LIBRARY_CXX_FLAGS) # mex does not accept -fPIC
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INSTALL_RPATH "\$ORIGIN")

set(CDT_MEX_SOURCES cdt.cpp)

foreach(MEX_SOURCE ${CDT_MEX_SOURCES})
	get_filename_component(MEX_NAME ${MEX_SOURCE} NAME_WE)
	# The target cannot be called `cdt`: that name belongs to the command line
	# tool. Only the file on disk has to match the matlab function name.
	add_library(mex_${MEX_NAME} SHARED ${MEX_SOURCE})
	set_target_properties(mex_${MEX_NAME} PROPERTIES
		OUTPUT_NAME ${MEX_NAME}
		LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
		RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
	if(WIN32)
		set_target_properties(mex_${MEX_NAME} PROPERTIES LINK_FLAGS "/export:mexFunction")
	endif()
	target_include_directories(mex_${MEX_NAME} PUBLIC ${Matlab_INCLUDE_DIRS})
	# Deliberately not ${Matlab_LIBRARIES}: that also pulls in libMatlabEngine
	# and libMatlabDataArray, which live outside the directories matlab
	# searches when it loads a mex file, so the result fails to load.
	target_link_libraries(mex_${MEX_NAME} PRIVATE
		cdt_core ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY})
endforeach()
