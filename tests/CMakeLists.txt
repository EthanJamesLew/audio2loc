add_executable(planeloc-tests
  unit/test_doppler.cpp
  unit/test_dsp.cpp
  unit/test_geo.cpp
  unit/test_util.cpp)
target_link_libraries(planeloc-tests PRIVATE planeloc_core Catch2::Catch2WithMain)
planeloc_warnings(planeloc-tests)

include(Catch)
catch_discover_tests(planeloc-tests)

# ---------------------------------------------------------------------------
# CLI contract: exact exit statuses, checked through a wrapper script because
# WILL_FAIL would pass on any nonzero status, including a crash.
#   planeloc_cli_test(NAME EXPECT code ARGS a b ... [STDOUT re] [STDERR re])
# ---------------------------------------------------------------------------
function(planeloc_cli_test name)
  cmake_parse_arguments(T "" "EXPECT;STDOUT;STDERR" "ARGS" ${ARGN})
  string(REPLACE ";" "|" _args "${T_ARGS}")
  set(_defs -DCLI=$<TARGET_FILE:planeloc> "-DARGS=${_args}" -DEXPECT=${T_EXPECT})
  if(DEFINED T_STDOUT)
    list(APPEND _defs "-DSTDOUT_MATCH=${T_STDOUT}")
  endif()
  if(DEFINED T_STDERR)
    list(APPEND _defs "-DSTDERR_MATCH=${T_STDERR}")
  endif()
  add_test(NAME cli.${name}
    COMMAND ${CMAKE_COMMAND} ${_defs} -P ${CMAKE_CURRENT_SOURCE_DIR}/cli/run_cli.cmake)
endfunction()

planeloc_cli_test(version EXPECT 0 ARGS --version STDOUT "^planeloc [0-9]+\\.[0-9]+\\.[0-9]+")
planeloc_cli_test(help EXPECT 0 ARGS --help STDOUT "usage: planeloc")
planeloc_cli_test(unknown_command EXPECT 1 ARGS frobnicate STDERR "unknown command")
planeloc_cli_test(bad_number EXPECT 1 ARGS simulate --seed abc STDERR "expects a number")
planeloc_cli_test(missing_value EXPECT 1 ARGS locate --audio STDERR "needs a value")

add_test(NAME version_from_tags
  COMMAND ${CMAKE_COMMAND} -DWORK=${CMAKE_CURRENT_BINARY_DIR}/version_from_tags
          -P ${CMAKE_CURRENT_SOURCE_DIR}/cli/version_from_tags.cmake)

# ---------------------------------------------------------------------------
# End to end on a simulated scenario: synthetic traffic and audio with a
# hidden observer, located by the suggested `locate` command.
# ---------------------------------------------------------------------------
add_test(NAME e2e.locate_synthetic_seed7
  COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/e2e/locate_synthetic.sh
          $<TARGET_FILE:planeloc> ${CMAKE_CURRENT_BINARY_DIR}/e2e-seed7 7 50)
set_tests_properties(e2e.locate_synthetic_seed7 PROPERTIES TIMEOUT 600 LABELS e2e)
