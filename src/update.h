// ADS-B database maintenance: per-UTC-day .plb caches built straight from
// adsb.lol's daily global history releases
// (github.com/adsblol/globe_history_YYYY), streamed - nothing is extracted.
#pragma once
#include <string>
#include <vector>

#include "adsb.h"

namespace pl {

struct UpdateOptions {
  std::string db = "planeloc_db";
  std::vector<std::string> dates;  // YYYY-MM-DD (UTC)
  LoadOptions filter;  // bbox / max altitude applied while ingesting
  bool force = false;
  std::string curl = "curl";
};

std::string db_day_path(const std::string& db, const std::string& date);
std::string utc_date(double epoch);  // YYYY-MM-DD
// Region filter a day file was built with (from its .meta sidecar). False if
// unknown.
bool db_day_bbox(const std::string& plb, double& w, double& s, double& e,
                 double& n);
// Dates (UTC) touched by [t0, t1].
std::vector<std::string> dates_between(double t0, double t1);
// Newest day with a published release (searches back from yesterday). Empty if
// none found.
std::string latest_available_date(const std::string& curl, int max_back = 5);
// Builds every requested day that is missing (or all with force). Returns
// number of failures.
int run_update(const UpdateOptions& o);
void print_db(const std::string& db);

}  // namespace pl
