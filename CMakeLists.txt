cmake_minimum_required(VERSION 3.16)
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()

# Set the target architecture.
# All modern x86/x64 processors support AVX2.
# Older x86/x64 processors may support SSE2 but not AVX2.
# Very old x86/x64 processors, or non x86/x64 
# processors, do not support any of the two.
# Modern ARM processors support NEON.

# --- Uncomment max one of the following lines ---

#							#Generic build
#set(ENABLE_SSE2 True)		#Build for x86-64 with SSE2
set(ENABLE_AVX2 True)		#Build for x86-64 with AVX2
#set(ENABLE_ARMNEON True)	#Build for ARM with NEON

# --------- DON NOT EDIT BELOW THIS LINE! --------

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Use CMake -DLGPL=ON if you opt for the LGPL license
option(LGPL "LGPL" ON)

# set the project name
project(cdt)

option(CDT_BUILD_EXECUTABLE "Build the cdt command line tool" ON)
option(CDT_BUILD_PYTHON "Build the nanobind python bindings" OFF)
option(CDT_BUILD_MATLAB "Build the matlab mex bindings" OFF)

# Resolve the target architecture. This must happen after project(), which is
# what populates APPLE and CMAKE_SYSTEM_PROCESSOR.
if(CMAKE_SYSTEM_PROCESSOR MATCHES "^(arm64|aarch64)$" OR CMAKE_OSX_ARCHITECTURES MATCHES "arm64")
	# On ARM the x86 flavours are meaningless. nfg's numerics.h keys off the
	# compiler-defined __ARM_NEON and includes simde's <x86/avx2.h>, so simde
	# must be fetched whether or not the user asked for NEON above.
	set(ENABLE_SSE2 False)
	set(ENABLE_AVX2 False)
	set(ENABLE_ARMNEON True)
elseif(APPLE)
	# disable SSE and AVX on Intel MacOS
	set(ENABLE_SSE2 False)
	set(ENABLE_AVX2 False)
endif()

# AVX includes SSE -> make that explicit
if(ENABLE_AVX2)
set(ENABLE_SSE2 True)
endif()

# Check that no incompatible architectures are slelected
if(ENABLE_ARMNEON AND ENABLE_SSE2)
message(FATAL_ERROR "ARMNEON and SSE2/AVX2 cannot be both set!")
endif()

# Dependencies
include(FetchContent)

if(EXISTS ../.local.cmake)
	set(nfg_SOURCE_DIR "../NFG/")
	set(indirect_predicates_SOURCE_DIR "../Indirect_Predicates/")
else()
	FetchContent_Declare(
		nfg
		URL https://github.com/MarcoAttene/nfg/archive/refs/heads/master.zip)
	FetchContent_Populate(nfg)
	message(STATUS "NFG is available in " ${nfg_SOURCE_DIR})

	FetchContent_Declare(
		Indirect_Predicates
		URL https://github.com/MarcoAttene/Indirect_Predicates/archive/refs/heads/master.zip)
	FetchContent_Populate(Indirect_Predicates)
	message(STATUS "Indirect_Predicates is available in " ${indirect_predicates_SOURCE_DIR})

endif()

if(ENABLE_ARMNEON)
	FetchContent_Declare(
		simde
		URL https://github.com/simd-everywhere/simde/archive/refs/heads/master.zip)
	FetchContent_Populate(simde)
	message(STATUS "simde is available in " ${simde_SOURCE_DIR})
endif()

# nfg's ip_error() prints to stderr and calls exit(0). That is tolerable for a
# command line tool but would take down the whole interpreter when CDT runs
# inside python or matlab. Serve a copy of the header that throws instead,
# placed ahead of the original on the include path, rather than editing the
# fetched (or locally checked out) sources in place.
set(CDT_PATCHED_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/patched_include")
file(READ "${nfg_SOURCE_DIR}/include/numerics.h" CDT_NUMERICS_H)
set(CDT_IP_ERROR_EXITING "inline void ip_error(const char* msg)\n{\n\tstd::cerr << msg;\n\texit(0);\n}")
set(CDT_IP_ERROR_THROWING "#include <stdexcept>\ninline void ip_error(const char* msg)\n{\n\tthrow std::runtime_error(msg);\n}")
string(REPLACE "${CDT_IP_ERROR_EXITING}" "${CDT_IP_ERROR_THROWING}" CDT_NUMERICS_H_PATCHED "${CDT_NUMERICS_H}")
if(CDT_NUMERICS_H_PATCHED STREQUAL CDT_NUMERICS_H)
	message(WARNING "Could not make nfg's ip_error() throw: a fatal error inside CDT will "
		"terminate the host process instead of raising an exception.")
endif()
file(WRITE "${CDT_PATCHED_INCLUDE_DIR}/numerics.h" "${CDT_NUMERICS_H_PATCHED}")

# The algorithm, with no dependency on the command line front end. Both the
# executable and the language bindings build on top of this.
add_library(cdt_core STATIC
	src/giftWrap.cpp
	src/PLC.cpp
	src/delaunay.cpp
	src/cdt.cpp
)
# The bindings link this into a shared object.
set_property(TARGET cdt_core PROPERTY POSITION_INDEPENDENT_CODE ON)

target_include_directories(cdt_core PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}/src
	${CDT_PATCHED_INCLUDE_DIR}
	${nfg_SOURCE_DIR}/include
	${indirect_predicates_SOURCE_DIR}/include
)

if(ENABLE_ARMNEON)
	target_include_directories(cdt_core PUBLIC ${simde_SOURCE_DIR}/simde)
	target_compile_options(cdt_core PUBLIC -DSIMDE_ENABLE_NATIVE_ALIASES)
endif()

# Compiler-specific options. These are PUBLIC because the exact predicates are
# header only: anything that includes them needs the same numerics.
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
	# grant IEEE 754 compliance; remove security controls; turn off annoying warnings
	target_compile_options(cdt_core PUBLIC "/fp:strict;/GS-;/D _CRT_SECURE_NO_WARNINGS")
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(cdt_core PUBLIC "/arch:AVX2")
	elseif(ENABLE_SSE2)
		target_compile_options(cdt_core PUBLIC "/arch:SSE2")
	elseif(ENABLE_ARMNEON)
		target_compile_options(cdt_core PUBLIC "/arch:ARMV8.0")
	endif()
else()
	# set max optimization
	target_compile_options(cdt_core PUBLIC -O3)
	# grant IEEE 754 compliance
	target_compile_options(cdt_core PUBLIC -frounding-math)
	# the exact predicates rely on error-free transformations (e.g. Dekker
	# two-product): forbid the compiler from contracting a*b+c into an FMA, which
	# rounds once instead of twice and silently breaks those identities. This
	# mirrors the /fp:strict used on the MSVC path.
	target_compile_options(cdt_core PUBLIC -ffp-contract=off)
	# set target architecture
	if(ENABLE_AVX2)
		target_compile_options(cdt_core PUBLIC "-mavx2")
		target_compile_options(cdt_core PUBLIC "-mfma")
	elseif(ENABLE_SSE2)
		target_compile_options(cdt_core PUBLIC "-msse2")
	elseif(ENABLE_ARMNEON)
		target_compile_options(cdt_core PUBLIC "-march=armv8-a+simd")
	endif()
endif()

if(CDT_BUILD_EXECUTABLE)
	add_executable(${PROJECT_NAME} src/main.cpp)
	target_link_libraries(${PROJECT_NAME} PRIVATE cdt_core)
endif()

if(CDT_BUILD_PYTHON)
	add_subdirectory(python)
endif()

if(CDT_BUILD_MATLAB)
	add_subdirectory(mex)
endif()
