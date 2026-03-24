# Set the minimum CMake version required to build the project.
cmake_minimum_required( VERSION 3.6 )

# Silence some warnings on macOS with new CMake versions.
if( NOT ${CMAKE_VERSION} VERSION_LESS 3.9 )
   cmake_policy( SET CMP0068 NEW )
endif()

set(CMAKE_INSTALL_PREFIX "" CACHE PATH "Path to the install directory")
# Use a relative RPATH at install time to avoid embedding absolute, system-specific
# library paths (which can create cycles when multiple external installs are involved).
# Use the literal string "$ORIGIN" so the installed shared object resolves libraries
# relative to its own location at runtime.
set(CMAKE_INSTALL_RPATH "\$ORIGIN")
# Avoid automatically adding link-path-derived directories into the install RPATH
# (this is what causes CMake to try to order /eos/... and /cvmfs/... entries and
# report a cycle). We'll control the install RPATH explicitly above.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH OFF)
set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)

# Set the project's name and version.
project( PartonMatcher )

# Set up the "C++ version" to use.
set( CMAKE_CXX_STANDARD_REQUIRED 17 CACHE STRING
   "Minimum C++ standard required for the build" )
set( CMAKE_CXX_STANDARD 17 CACHE STRING
   "C++ standard to use for the build" )
set( CMAKE_CXX_EXTENSIONS FALSE CACHE BOOL
   "(Dis)allow using compiler extensions" )

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Option to skip ROOT dictionary generation (useful if rootcling is incompatible)
option(SKIP_ROOT_DICTIONARY "Skip ROOT_GENERATE_DICTIONARY step" OFF)

# Specify the install locations for libraries and binaries.
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin )
set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set( CMAKE_INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib )  # Needed by ROOT_GENERATE_DICTIONARY()

# Set the warning flag(s) to use.
set( CMAKE_CXX_FLAGS "-Wall -Wextra -Wshadow -pedantic -O3 -g" )

set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

# Silence boost warnings
add_definitions(-DBOOST_BIND_GLOBAL_PLACEHOLDERS)
add_definitions(-DBOOST_ALLOW_DEPRECATED_HEADERS)

# Add ROOT system directory and require ROOT.
find_package( ROOT 6.28.00 REQUIRED COMPONENTS Core Hist RIO Tree ROOTDataFrame ROOTVecOps)

find_package( FastFrames REQUIRED )

file(GLOB SOURCES "Root/*.cc")
file(GLOB LIB_HEADERS "PartonMatcher/*.h")
if ( NOT ONNXRUNTIME_ROOTDIR STREQUAL "")
   add_definitions(-DONNXRUNTIME_AVAILABLE)
endif()

# Build the shared library.
add_library( PartonMatcher SHARED ${LIB_HEADERS} ${SOURCES} )
target_include_directories( PartonMatcher PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}> $<INSTALL_INTERFACE:> ${FastFrames_INCLUDE_DIRS} ${ROOT_INCLUDE_DIRS} ${ONNXRUNTIME_INCLUDE_DIRS})
find_library(FastFrames_LIBRARY FastFrames ${CMAKE_SOURCE_DIR}/install)
target_link_libraries( PartonMatcher PUBLIC ${FastFrames_LIBRARY} ${ROOT_LIBRARIES} ${ONNXRUNTIME_LIBRARIES})
set_property( TARGET PartonMatcher
   PROPERTY PUBLIC_HEADER ${LIB_HEADERS} )
target_include_directories(PartonMatcher PUBLIC ${CMAKE_CURRENT_LIST_DIR} )
target_compile_options(PartonMatcher PUBLIC -Wno-shadow)

set(SETUP ${CMAKE_CURRENT_BINARY_DIR}/setup.sh)
file(WRITE ${SETUP} "#!/bin/bash\n")
file(APPEND ${SETUP} "# this is an auto-generated setup script\n" )
file(APPEND ${SETUP} "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${CMAKE_INSTALL_PREFIX}/lib\n")

# Install the libraries.
install( TARGETS PartonMatcher
   EXPORT PartonMatcherTargets
   ARCHIVE DESTINATION lib
   LIBRARY DESTINATION lib
   PUBLIC_HEADER DESTINATION include/PartonMatcher )

if(NOT SKIP_ROOT_DICTIONARY)
   ROOT_GENERATE_DICTIONARY(PartonMatcher_dict PartonMatcher/PartonMatcher.h MODULE PartonMatcher LINKDEF Root/LinkDef.h)
   install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lib/lib${PROJECT_NAME}_rdict.pcm ${CMAKE_CURRENT_BINARY_DIR}/lib/lib${PROJECT_NAME}.rootmap DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
else()
   message(STATUS "SKIP_ROOT_DICTIONARY=ON: not generating/installing ROOT dictionaries")
endif()
