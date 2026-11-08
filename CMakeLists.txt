cmake_minimum_required(VERSION 3.24)
# The version comes from the release tags; 0.1.0 is for a checkout without any.
# Builds without the git history (a Docker context) pass it: -DPLANELOC_VERSION=1.2.0
include(cmake/version_from_tags.cmake)
set(PLANELOC_VERSION "" CACHE STRING "Version to report instead of the one from git tags")
if(PLANELOC_VERSION)
  string(REGEX MATCH "^v?([0-9]+\\.[0-9]+\\.[0-9]+)" _ "${PLANELOC_VERSION}")
  if(NOT CMAKE_MATCH_1)
    message(FATAL_ERROR "PLANELOC_VERSION=${PLANELOC_VERSION} is not MAJOR.MINOR.PATCH[-PRERELEASE]")
  endif()
  set(_pl_release "${CMAKE_MATCH_1}")
  string(REGEX REPLACE "^v" "" _pl_build "${PLANELOC_VERSION}")
else()
  planeloc_version_from_tags("${CMAKE_CURRENT_SOURCE_DIR}" 0.1.0 _pl_build _pl_release)
endif()
project(planeloc
  VERSION ${_pl_release}
  DESCRIPTION "Geolocate a recording from aircraft Doppler curves and ADS-B history"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(PLANELOC_BUILD_TESTS "Build unit and end-to-end tests (fetches Catch2)" ON)
option(PLANELOC_NATIVE "Tune for the build machine (-march/-mcpu=native)" OFF)
option(PLANELOC_WERROR "Treat compiler warnings in planeloc's own targets as errors" OFF)
option(PLANELOC_WITH_QGIS "Link libqgis_core and enable native .qgz export" OFF)
option(PLANELOC_CUDA "Run spectrograms, Doppler fits, candidate search and scoring on an NVIDIA GPU" OFF)
option(PLANELOC_PROF_PHASES "Count GPU kernel phase cycles (profiling only)" OFF)
set(PLANELOC_SANITIZE "" CACHE STRING
  "Semicolon-separated sanitizers, e.g. address;undefined (use with Debug)")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# Warnings apply to planeloc's own targets only, not to fetched dependencies.
function(planeloc_warnings target)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 $<$<BOOL:${PLANELOC_WERROR}>:/WX>)
  else()
    target_compile_options(${target} PRIVATE
      "$<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic;-Wshadow>"
      "$<$<AND:$<COMPILE_LANGUAGE:CXX>,$<BOOL:${PLANELOC_WERROR}>>:-Werror>")
  endif()
endfunction()

if(PLANELOC_NATIVE)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-mcpu=native PLANELOC_HAS_MCPU_NATIVE)
  if(PLANELOC_HAS_MCPU_NATIVE)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-mcpu=native>)
  else()
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-march=native>)
  endif()
endif()

if(PLANELOC_SANITIZE)
  string(REPLACE ";" "," _san "${PLANELOC_SANITIZE}")
  add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-fsanitize=${_san};-fno-omit-frame-pointer>")
  add_link_options(-fsanitize=${_san})
endif()

find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

# ---------------------------------------------------------------------------
# The version of this build, regenerated from git at every build.
# ---------------------------------------------------------------------------
set(_version_cpp ${PROJECT_BINARY_DIR}/generated/version.cpp)
add_custom_target(planeloc_git_version
  COMMAND ${CMAKE_COMMAND} -DSOURCE_DIR=${PROJECT_SOURCE_DIR}
          -DBASE=${PROJECT_VERSION} -DFIXED=${_pl_build}
          -DPINNED=${PLANELOC_VERSION} -DOUT=${_version_cpp}
          -P ${PROJECT_SOURCE_DIR}/cmake/git_version.cmake
  BYPRODUCTS ${_version_cpp}
  COMMENT "planeloc version from git"
  VERBATIM)

# ---------------------------------------------------------------------------
# Everything but main() is a library, so the tests link the same code.
# ---------------------------------------------------------------------------
add_library(planeloc_core STATIC
  src/adsb.cpp src/aircraft_types.cpp src/dem.cpp src/doppler.cpp src/dsp.cpp
  src/export.cpp src/geo.cpp src/io.cpp src/mcp.cpp src/qgis_native.cpp
  src/search.cpp src/simulate.cpp src/update.cpp src/util.cpp
  ${_version_cpp})
add_dependencies(planeloc_core planeloc_git_version)
target_include_directories(planeloc_core PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(planeloc_core PUBLIC Threads::Threads ZLIB::ZLIB)
planeloc_warnings(planeloc_core)
# -ffast-math for the spectrogram and scoring loops, but NaN and infinity
# stay meaningful: they mark missing samples and rejected fits.
target_compile_options(planeloc_core PRIVATE
  "$<$<AND:$<CONFIG:Release,RelWithDebInfo>,$<COMPILE_LANGUAGE:CXX>>:-ffast-math;-fno-finite-math-only>"
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-O3>)

if(PLANELOC_CUDA)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    # For portable builds: -DCMAKE_CUDA_ARCHITECTURES="75-real;86-real;89-real;90"
    set(CMAKE_CUDA_ARCHITECTURES native)
  endif()
  if(NOT CMAKE_CUDA_COMPILER AND NOT DEFINED ENV{CUDACXX} AND EXISTS /usr/local/cuda/bin/nvcc)
    set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)  # nvcc is often not on PATH
  endif()
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  set(CMAKE_CUDA_STANDARD 20)
  target_sources(planeloc_core PRIVATE src/gpu.cu)
  set_target_properties(planeloc_core PROPERTIES CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
  set_source_files_properties(src/gpu.cu PROPERTIES
    COMPILE_OPTIONS "-lineinfo;--expt-relaxed-constexpr")
  target_compile_definitions(planeloc_core PUBLIC PLANELOC_CUDA
    $<$<BOOL:${PLANELOC_PROF_PHASES}>:PL_PROF_PHASES>)
  target_link_libraries(planeloc_core PUBLIC CUDA::cudart_static CUDA::cufft)
  message(STATUS "planeloc: CUDA enabled (archs: ${CMAKE_CUDA_ARCHITECTURES})")
else()
  target_sources(planeloc_core PRIVATE src/gpu_stub.cpp)
endif()

if(PLANELOC_WITH_QGIS)
  # QGIS 3.x (Qt5) or QGIS 3.40+/4 (Qt6). Override QGIS_PREFIX if installed elsewhere.
  set(QGIS_PREFIX "/usr" CACHE PATH "QGIS install prefix")
  find_path(QGIS_INCLUDE_DIR qgsapplication.h
    PATHS ${QGIS_PREFIX}/include/qgis ${QGIS_PREFIX}/include
          /Applications/QGIS.app/Contents/Frameworks/qgis_core.framework/Headers REQUIRED)
  find_library(QGIS_CORE_LIBRARY NAMES qgis_core
    PATHS ${QGIS_PREFIX}/lib ${QGIS_PREFIX}/lib64 /Applications/QGIS.app/Contents/Frameworks REQUIRED)
  find_package(Qt6 COMPONENTS Core Gui Xml Widgets QUIET)
  if(Qt6_FOUND)
    set(QT_LIBS Qt6::Core Qt6::Gui Qt6::Xml Qt6::Widgets)
  else()
    find_package(Qt5 COMPONENTS Core Gui Xml Widgets REQUIRED)
    set(QT_LIBS Qt5::Core Qt5::Gui Qt5::Xml Qt5::Widgets)
  endif()
  target_compile_definitions(planeloc_core PUBLIC PLANELOC_WITH_QGIS)
  target_include_directories(planeloc_core SYSTEM PRIVATE ${QGIS_INCLUDE_DIR})
  target_link_libraries(planeloc_core PUBLIC ${QGIS_CORE_LIBRARY} ${QT_LIBS})
  message(STATUS "planeloc: native QGIS export enabled (${QGIS_CORE_LIBRARY})")
endif()

add_executable(planeloc src/main.cpp)
target_link_libraries(planeloc PRIVATE planeloc_core)
planeloc_warnings(planeloc)

include(GNUInstallDirs)
install(TARGETS planeloc RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

if(PLANELOC_BUILD_TESTS)
  include(FetchContent)
  set(FETCHCONTENT_QUIET ON)
  # Offline builds can set FETCHCONTENT_SOURCE_DIR_CATCH2 to a checkout.
  FetchContent_Declare(catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.5.4
    GIT_SHALLOW TRUE)
  FetchContent_MakeAvailable(catch2)
  enable_testing()
  add_subdirectory(tests)
endif()
