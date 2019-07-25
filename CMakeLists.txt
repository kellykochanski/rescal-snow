cmake_minimum_required (VERSION 3.0.2)

project (rescal-snow 
  VERSION 0.0.1
  DESCRIPTION "TODO"
  LANGUAGES C
)

# add_subdirectory(apps)
# add_subdirectory(tests)

add_subdirectory(lib/libpng-1.6.37)
add_subdirectory(lib/zlib-1.2.11)

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

target_link_libraries(rescal    PRIVATE png pthread)
target_link_libraries(genesis   PRIVATE png)  
target_link_libraries(bin2png   PRIVATE png)  
target_link_libraries(csp2png   PRIVATE png)  

target_include_directories(rescal    PRIVATE src)
target_include_directories(genesis   PRIVATE src)
target_include_directories(regenesis PRIVATE src)
target_include_directories(bin2png   PRIVATE src)
target_include_directories(csp2png   PRIVATE src)
target_include_directories(cspinfo   PRIVATE src)
target_include_directories(csp2bin   PRIVATE src)
