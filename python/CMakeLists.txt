# nanobind bindings, laid out the same way as the libigl python bindings: one
# extension module per package, imported by the package's __init__.py.

find_package(Python 3.8
	REQUIRED COMPONENTS Interpreter Development.Module
	OPTIONAL_COMPONENTS Development.SABIModule)

include(FetchContent)
FetchContent_Declare(
	nanobind
	GIT_REPOSITORY https://github.com/wjakob/nanobind.git
	GIT_TAG        v2.7.0
)
FetchContent_MakeAvailable(nanobind)

nanobind_add_module(_cdt STABLE_ABI NB_STATIC cdt_module.cpp)
target_link_libraries(_cdt PRIVATE cdt_core)

# So that the package can be imported straight out of the source tree, without
# installing a wheel first.
set_target_properties(_cdt PROPERTIES
	LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/cdt")

# Where scikit-build-core picks it up when building a wheel.
install(TARGETS _cdt LIBRARY DESTINATION cdt)
