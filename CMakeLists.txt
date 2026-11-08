cmake_minimum_required(VERSION 3.16)
project(planeloc LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(PLANELOC_NATIVE "Optimise for the build machine (-march=native)" ON)
option(PLANELOC_WITH_QGIS "Link libqgis_core and enable native .qgz export" OFF)
option(PLANELOC_CUDA "Build CUDA kernels for spectrograms, Doppler fits, candidate search and scoring" OFF)

find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)

add_executable(planeloc
  src/main.cpp src/util.cpp src/geo.cpp src/dsp.cpp src/io.cpp src/dem.cpp
  src/adsb.cpp src/aircraft_types.cpp src/doppler.cpp src/search.cpp
  src/export.cpp src/qgis_native.cpp src/simulate.cpp src/update.cpp src/mcp.cpp)

if(PLANELOC_CUDA)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES native)  # e.g. -DCMAKE_CUDA_ARCHITECTURES="75;86;89" for portable builds
  endif()
  if(NOT CMAKE_CUDA_COMPILER AND NOT DEFINED ENV{CUDACXX} AND EXISTS /usr/local/cuda/bin/nvcc)
    set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)  # nvcc is often not on PATH
  endif()
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)
  set(CMAKE_CUDA_STANDARD 20)
  target_sources(planeloc PRIVATE src/gpu.cu)
  set_target_properties(planeloc PROPERTIES CUDA_ARCHITECTURES "${CMAKE_CUDA_ARCHITECTURES}")
  set_source_files_properties(src/gpu.cu PROPERTIES COMPILE_OPTIONS "-lineinfo;--expt-relaxed-constexpr")
  target_compile_definitions(planeloc PRIVATE PLANELOC_CUDA)
  target_link_libraries(planeloc PRIVATE CUDA::cudart_static CUDA::cufft)
  message(STATUS "planeloc: CUDA enabled (archs: ${CMAKE_CUDA_ARCHITECTURES})")
else()
  target_sources(planeloc PRIVATE src/gpu_stub.cpp)
endif()

target_compile_options(planeloc PRIVATE $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -Wno-unused-parameter>
  $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-O3 -ffast-math -fno-finite-math-only>)
if(PLANELOC_NATIVE)
  target_compile_options(planeloc PRIVATE $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:-march=native>)
endif()
target_link_libraries(planeloc PRIVATE Threads::Threads ZLIB::ZLIB)

if(PLANELOC_WITH_QGIS)
  # QGIS 3.x (Qt5) or QGIS 3.40+/4 (Qt6). Override QGIS_PREFIX if installed elsewhere.
  set(QGIS_PREFIX "/usr" CACHE PATH "QGIS install prefix")
  find_path(QGIS_INCLUDE_DIR qgsapplication.h PATHS ${QGIS_PREFIX}/include/qgis ${QGIS_PREFIX}/include
            /Applications/QGIS.app/Contents/Frameworks/qgis_core.framework/Headers REQUIRED)
  find_library(QGIS_CORE_LIBRARY NAMES qgis_core PATHS ${QGIS_PREFIX}/lib ${QGIS_PREFIX}/lib64
               /Applications/QGIS.app/Contents/Frameworks REQUIRED)
  find_package(Qt6 COMPONENTS Core Gui Xml Widgets QUIET)
  if(Qt6_FOUND)
    set(QT_LIBS Qt6::Core Qt6::Gui Qt6::Xml Qt6::Widgets)
  else()
    find_package(Qt5 COMPONENTS Core Gui Xml Widgets REQUIRED)
    set(QT_LIBS Qt5::Core Qt5::Gui Qt5::Xml Qt5::Widgets)
  endif()
  target_compile_definitions(planeloc PRIVATE PLANELOC_WITH_QGIS)
  target_include_directories(planeloc PRIVATE ${QGIS_INCLUDE_DIR})
  target_link_libraries(planeloc PRIVATE ${QGIS_CORE_LIBRARY} ${QT_LIBS})
  message(STATUS "planeloc: native QGIS export enabled (${QGIS_CORE_LIBRARY})")
endif()

install(TARGETS planeloc RUNTIME DESTINATION bin)
