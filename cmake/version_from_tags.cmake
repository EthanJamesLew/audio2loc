# The version of a source tree, from its release tags.
#
#   planeloc_version_from_tags(<source dir> <fallback> <out build> <out release>)
#
# Release tags are MAJOR.MINOR.PATCH or MAJOR.MINOR.PATCH-PRERELEASE (semver.org),
# with or without a leading v. <out release> is MAJOR.MINOR.PATCH of the newest tag behind the
# source, or <fallback> without one. <out build> is what a build reports:
#   on a clean tag            0.2.0                    or 1.0.0-rc.1
#   five commits after it     0.2.1-dev.5+gabc12345    or 1.0.0-rc.1.dev.5+gabc12345
#   uncommitted changes       adds +gabc12345.dirty, or .dirty after the hash
#   no tag in the history     <fallback>-dev+gabc12345
#   no git, no archived tag   <fallback>
# All of them sort after the tag and before the next release.
function(planeloc_version_from_tags source_dir fallback out_build out_release)
  set(patterns "v[0-9]*.[0-9]*.[0-9]*" "[0-9]*.[0-9]*.[0-9]*")
  set(core "v?([0-9]+)\\.([0-9]+)\\.([0-9]+)(-[0-9A-Za-z.-]+)?")
  set(tag "")
  set(ahead 0)
  set(hash "")
  set(dirty FALSE)

  # `git archive` fills this file in, so a source archive knows its tag: the tag
  # itself, or past it, as `git describe` puts it, <tag>-<ahead>-g<hash>.
  set(archived "$Format")
  if(EXISTS "${source_dir}/cmake/archive_version.txt")
    file(READ "${source_dir}/cmake/archive_version.txt" archived)
    string(STRIP "${archived}" archived)
  endif()
  if(NOT archived MATCHES "^\\$Format")
    if(archived MATCHES "^(.*)-([0-9]+)-g([0-9a-f]+)$")
      set(tag "${CMAKE_MATCH_1}")
      set(ahead "${CMAKE_MATCH_2}")
      set(hash "${CMAKE_MATCH_3}")
    else()
      set(tag "${archived}")
    endif()
  else()
    find_package(Git QUIET)
    if(GIT_FOUND)
      execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${source_dir}" rev-parse --short=8 HEAD
        OUTPUT_VARIABLE hash OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET RESULT_VARIABLE no_git)
      if(no_git)
        set(hash "")
      else()
        execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${source_dir}" diff --quiet HEAD --
          RESULT_VARIABLE dirty ERROR_QUIET)
        # The highest tag on this commit, else the highest behind it. A release
        # is above its pre-releases, and git orders those among themselves.
        # `git describe` is no use here: it picks arbitrarily among several tags
        # on one commit.
        foreach(filter --points-at --merged)
          execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${source_dir}"
              -c versionsort.suffix=- tag ${filter} HEAD --list ${patterns} --sort=-v:refname
            OUTPUT_VARIABLE tags OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
          string(REPLACE "\n" ";" tags "${tags}")
          set(best "")
          set(best_pre "")
          foreach(candidate IN LISTS tags)
            if(candidate MATCHES "^${core}$")
              set(number "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}.${CMAKE_MATCH_3}")
              set(pre "${CMAKE_MATCH_4}")
              if(tag STREQUAL "" OR number VERSION_GREATER best
                 OR (number VERSION_EQUAL best AND NOT best_pre STREQUAL "" AND pre STREQUAL ""))
                set(tag "${candidate}")
                set(best "${number}")
                set(best_pre "${pre}")
              endif()
            endif()
          endforeach()
          if(NOT tag STREQUAL "")
            break()
          endif()
        endforeach()
        if(NOT tag STREQUAL "")
          execute_process(COMMAND "${GIT_EXECUTABLE}" -C "${source_dir}" rev-list --count "${tag}..HEAD"
            OUTPUT_VARIABLE ahead OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
        endif()
      endif()
    endif()
  endif()

  set(build_meta "")
  if(NOT hash STREQUAL "")
    set(build_meta "g${hash}")
    if(dirty)
      string(APPEND build_meta ".dirty")
    endif()
  endif()

  if(tag MATCHES "^${core}$")
    set(major "${CMAKE_MATCH_1}")
    set(minor "${CMAKE_MATCH_2}")
    set(patch "${CMAKE_MATCH_3}")
    set(pre "${CMAKE_MATCH_4}")
    set(release "${major}.${minor}.${patch}")
    if(ahead EQUAL 0)
      set(build "${release}${pre}")
      if(dirty)
        string(APPEND build "+${build_meta}")
      endif()
    elseif(pre)
      set(build "${release}${pre}.dev.${ahead}+${build_meta}")
    else()
      math(EXPR next "${patch} + 1")
      set(build "${major}.${minor}.${next}-dev.${ahead}+${build_meta}")
    endif()
  elseif(NOT build_meta STREQUAL "")
    set(release "${fallback}")
    set(build "${fallback}-dev+${build_meta}")
  else()
    set(release "${fallback}")
    set(build "${fallback}")
  endif()
  set(${out_build} "${build}" PARENT_SCOPE)
  set(${out_release} "${release}" PARENT_SCOPE)
endfunction()
