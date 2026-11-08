#pragma once
#include <string>

namespace pl {
// Runs an MCP stdio server until stdin closes. `self` is the audio2loc binary
// used for tool calls.
int run_mcp_server(const std::string& self, const std::string& work_dir,
                   const std::string& db_dir);
}  // namespace pl
