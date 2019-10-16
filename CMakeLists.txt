cmake_minimum_required (VERSION 3.9)

# Use new policies on versions of CMake that have been tested
cmake_policy(VERSION 3.9.3...3.15.4)

# Let us use, e.g., PNG_ROOT=/path/to/libpng, in the environment to
# tell find_package() where to look. Usefull on clusters.
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()


project (rescal-snow
  DESCRIPTION "Simulating snow self-organization with cellular automata"
  LANGUAGES C
)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
find_package(Sanitizers)

set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(PNG 1.6 REQUIRED)
find_package(Threads REQUIRED)

set(COMMON_SOURCES
  src/cells.c
  src/lgca.c
  src/doublets.c
  src/format.c
  src/models.c
  src/param.c
  src/simul.c
  src/surface.c
  src/space.c
  src/trace.c
  src/transitions.c
  src/view.c
)

add_executable(rescal
  ${COMMON_SOURCES}
  src/callbacks.c
  src/drand48.c
  src/entry.c
  src/rescal.c
)

add_executable(genesis
  src/format.c
  src/param.c
  src/genesis.c
)

add_executable(regenesis
  src/format.c
  src/regenesis.c
)

add_executable(bin2png
  ${COMMON_SOURCES}
  src/bin2png.c
)

add_executable(csp2png
  ${COMMON_SOURCES}
  src/csp2png.c
)

add_executable(cspinfo
  src/format.c
  src/cspinfo.c
)

add_executable(csp2bin
  src/format.c
  src/csp2bin.c
)

function(setup_target target)
  set_property(TARGET ${target} PROPERTY C_STANDARD 99)
  target_include_directories(${target} PRIVATE src)
  if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /WX)
  else()
    target_compile_options(${target} PRIVATE -Wall -Wextra -pedantic -Werror)
  endif()
  add_sanitizers(${target})
endfunction(setup_target)

target_link_libraries(rescal    PRIVATE PNG::PNG Threads::Threads m)
target_link_libraries(genesis   PRIVATE PNG::PNG m)
target_link_libraries(regenesis PRIVATE m)
target_link_libraries(bin2png   PRIVATE PNG::PNG m)
target_link_libraries(csp2png   PRIVATE PNG::PNG m)

setup_target(rescal)
setup_target(genesis)
setup_target(regenesis)
setup_target(bin2png)
setup_target(csp2png)
setup_target(cspinfo)
setup_target(csp2bin)
