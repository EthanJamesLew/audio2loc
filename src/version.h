// The version of this build, from the git release tags
// (cmake/version_from_tags.cmake).
#pragma once

namespace pl {
// "0.2.0" on a release tag, "0.2.1-dev.5+gabc12345.dirty" between releases.
const char* version_string() noexcept;
}  // namespace pl
