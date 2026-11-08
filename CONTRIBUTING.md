# Contributing

## Build and test

```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DAUDIO2LOC_WERROR=ON
cmake --build build -j
ctest --test-dir build --output-on-failure --parallel
```

`ctest -LE e2e` skips the end-to-end locate, which simulates a scenario and requires the answer
within 50 m of the hidden observer. Sanitizer builds take `-DCMAKE_BUILD_TYPE=Debug
"-DAUDIO2LOC_SANITIZE=address;undefined"`.

CI (`.github/workflows/ci.yml`) checks formatting and builds with `-Werror`. It runs every test on
GCC 12 and 14, Clang 18, Linux arm64 and macOS. It also runs ASan/UBSan, builds CUDA (where
a CUDA build must fall back to the CPU and still pass), and builds the Docker image and runs a smoke
test.

## C++ style

The code follows the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html)
for layout and the [C++ Core Guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines)
for everything else. `clang-format` 18 decides layout and CI checks it:

```sh
pipx install clang-format==18.1.8 && clang-format -i src/*.{h,cpp,cu} tests/unit/*.cpp
pre-commit install   # optional: format and tidy whitespace on every commit
```

What the formatter cannot check:

- **Headers:** `#pragma once`, and a comment saying what the file is for.
- **Comments** say why, not what. Delete code rather than commenting it out.
- **Ownership:** values first, then `std::unique_ptr`; no owning raw pointers.
- **Errors:** bad user input (a file, a flag, an event box) produces a message and exit status 1,
  never a crash. `assert` is for this program's own invariants.
- **Dependencies:** the core needs only zlib. Add one only for a clear gain, and pin it to a
  release.

Hot paths (spectrogram scoring, the candidate search, GPU kernels) may trade these for speed. Say why
in a comment, and measure.

The CPU path is the reference. `src/gpu.cu` must reproduce it, and a change to one needs the same
change to the other. Compare the two on the same input (`--no-gpu` forces the CPU) and say in the
change which one was run on a GPU.

## Commits

One logical change per commit. The subject is an imperative sentence without a prefix ("Fit the
Doppler curve on emission time"). The body says why and what was measured. Work on a branch and
merge to `main` through a pull request once CI passes.

## Versions

Versions follow [Semantic Versioning](https://semver.org) and come from git tags, so nothing in the
source needs editing to release. `audio2loc --version` and the MCP server report the build exactly:

| Build | Reports |
| --- | --- |
| On a clean release tag `v0.2.0` | `0.2.0` |
| On a pre-release tag | `1.0.0-rc.1` |
| Five commits after `v0.2.0` | `0.2.1-dev.5+gabc12345` |
| Uncommitted changes | adds `+gabc12345.dirty`, or `.dirty` after the hash |
| No tag in the history, or a shallow clone | `0.1.0-dev+gabc12345` |
| GitHub's source archive of a tag | `0.2.0` |
| `-DAUDIO2LOC_VERSION=1.2.0` (builds without git, e.g. Docker) | `1.2.0` |

To release, go to **Releases > Draft a new release** on GitHub, type a new tag such as `v0.2.0`, and
publish. The Release workflow builds and tests the tag on Linux x86-64 and arm64 (glibc 2.35+) and
on macOS arm64. It checks that each binary reports the tag, then attaches the
`audio2loc-<version>-<target>.tar.gz` archives and their checksums to the release. It also pushes
`ghcr.io/<owner>/audio2loc:<version>` and `:<version>-cuda`, plus `:latest` and `:cuda` for
non-pre-releases. Running the workflow by hand (**Actions > Release > Run workflow**) is a dry run.
