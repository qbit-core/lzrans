# liblzr + lzr CLI.  Equivalent to the Makefile:
#   cmake -B build/cmake -DCMAKE_BUILD_TYPE=Release && cmake --build build/cmake -j
#   ctest --test-dir build/cmake        (C API tests + CLI/tar round trips)
#   cmake --install build/cmake --prefix /usr/local
cmake_minimum_required(VERSION 3.13)

file(STRINGS include/lzr.h _v REGEX "^#define LZR_LIB_VERSION_(MAJOR|MINOR|PATCH) +[0-9]+")
string(REGEX REPLACE ".*MAJOR +([0-9]+).*MINOR +([0-9]+).*PATCH +([0-9]+).*" "\\1.\\2.\\3" LZR_VERSION "${_v}")
project(lzr VERSION ${LZR_VERSION} LANGUAGES C)

option(LZR_BUILD_SHARED "Build liblzr as a shared library" ON)
option(LZR_BUILD_STATIC "Build liblzr as a static library" ON)
option(LZR_BUILD_CLI "Build the lzr command-line tool" ON)
option(LZR_BUILD_TESTS "Build tests" ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_EXTENSIONS ON)          # gnu11, as the Python build
set(CMAKE_C_FLAGS_RELEASE "-O3")
find_package(Threads REQUIRED)
include(GNUInstallDirs)

set(LZR_SOURCES
  src/lzrans/api.c src/lzrans/enc.c src/lzrans/dec.c src/lzrans/opt.c
  src/lzrans/crc32.c src/lzrans/huf.c src/lzrans/fse.c)

function(lzr_lib target kind)
  add_library(${target} ${kind} ${LZR_SOURCES})
  target_include_directories(${target} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
  target_link_libraries(${target} PRIVATE Threads::Threads m)
  set_target_properties(${target} PROPERTIES OUTPUT_NAME lzr C_VISIBILITY_PRESET hidden
    VERSION ${PROJECT_VERSION} SOVERSION ${PROJECT_VERSION_MAJOR} PUBLIC_HEADER include/lzr.h)
  if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wno-unused-parameter -Wno-missing-field-initializers)
  endif()
  install(TARGETS ${target} EXPORT lzrTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
endfunction()

if(LZR_BUILD_STATIC)
  lzr_lib(lzr_static STATIC)
  target_link_libraries(lzr_static INTERFACE Threads::Threads m)
endif()
if(LZR_BUILD_SHARED)
  lzr_lib(lzr_shared SHARED)
endif()

if(LZR_BUILD_CLI)
  add_executable(lzr cli/lzr.c)
  if(LZR_BUILD_STATIC)
    target_link_libraries(lzr PRIVATE lzr_static)
  else()
    target_link_libraries(lzr PRIVATE lzr_shared)
  endif()
  install(TARGETS lzr RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(CODE "
    foreach(alias unlzr lzrcat)
      execute_process(COMMAND \${CMAKE_COMMAND} -E create_symlink lzr
        \"\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}/\${alias}\")
    endforeach()")
  install(FILES docs/lzr.1 DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
endif()

set(LZR_PC "${CMAKE_CURRENT_BINARY_DIR}/lzr.pc")
file(WRITE ${LZR_PC} "prefix=${CMAKE_INSTALL_PREFIX}
libdir=\${prefix}/${CMAKE_INSTALL_LIBDIR}
includedir=\${prefix}/${CMAKE_INSTALL_INCLUDEDIR}

Name: liblzr
Description: lzrans (LZR 1.0) compression library
Version: ${PROJECT_VERSION}
Libs: -L\${libdir} -llzr
Libs.private: -lm -pthread
Cflags: -I\${includedir}
")
install(FILES ${LZR_PC} DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
install(EXPORT lzrTargets NAMESPACE lzr:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/lzr FILE lzrConfig.cmake)

if(LZR_BUILD_TESTS)
  enable_testing()
  add_executable(test_api tests/c/test_api.c)
  target_link_libraries(test_api PRIVATE $<IF:$<BOOL:${LZR_BUILD_SHARED}>,lzr_shared,lzr_static>)
  add_test(NAME api COMMAND test_api)
  if(LZR_BUILD_CLI)
    add_test(NAME cli COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/tests/cli/run.sh $<TARGET_FILE:lzr>)
  endif()
endif()
