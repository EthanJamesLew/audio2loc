# Runs one CLI invocation and checks its exact exit status (WILL_FAIL
# would accept any nonzero status, including a crash).
#
#   cmake -DCLI=<planeloc> -DARGS="a|b|c" -DEXPECT=<code>
#         [-DSTDOUT_MATCH=<regex>] [-DSTDERR_MATCH=<regex>]
#         [-DOUTPUT=<file that must exist afterwards and is removed first>]
#         [-DOUTPUT_MATCH=<regex the OUTPUT file must contain>]
#         -P run_cli.cmake
#
# ARGS separates arguments with '|' so that a single -D value can carry them.
if(NOT DEFINED CLI OR NOT DEFINED EXPECT)
  message(FATAL_ERROR "run_cli.cmake needs CLI and EXPECT")
endif()
string(REPLACE "|" ";" _args "${ARGS}")
if(DEFINED OUTPUT)
  file(REMOVE "${OUTPUT}")
endif()
execute_process(COMMAND "${CLI}" ${_args}
  RESULT_VARIABLE _rc
  OUTPUT_VARIABLE _out
  ERROR_VARIABLE _err)
if(NOT "${_rc}" STREQUAL "${EXPECT}")
  message(FATAL_ERROR "exit status '${_rc}', expected ${EXPECT}\n"
    "--- stdout ---\n${_out}\n--- stderr ---\n${_err}")
endif()
if(DEFINED STDOUT_MATCH AND NOT _out MATCHES "${STDOUT_MATCH}")
  message(FATAL_ERROR "stdout does not match '${STDOUT_MATCH}':\n${_out}")
endif()
if(DEFINED STDERR_MATCH AND NOT _err MATCHES "${STDERR_MATCH}")
  message(FATAL_ERROR "stderr does not match '${STDERR_MATCH}':\n${_err}")
endif()
if(DEFINED OUTPUT)
  if(NOT EXISTS "${OUTPUT}")
    message(FATAL_ERROR "expected output file '${OUTPUT}' was not written")
  endif()
  if(DEFINED OUTPUT_MATCH)
    file(READ "${OUTPUT}" _content)
    if(NOT _content MATCHES "${OUTPUT_MATCH}")
      message(FATAL_ERROR "'${OUTPUT}' does not match '${OUTPUT_MATCH}'")
    endif()
  endif()
endif()
