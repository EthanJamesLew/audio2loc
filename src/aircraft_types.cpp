// ICAO type designator -> engine class. The video's approach: no clean
// "propeller vs jet" field in ADS-B, so filter out *known* jet types and keep
// everything else (props, helicopters, unknowns).
#include <fstream>
#include <mutex>
#include <unordered_map>

#include "adsb.h"
#include "util.h"

namespace pl {

static const char* kJets =
    // Airbus
    "A318 A319 A320 A321 A19N A20N A21N A306 A30B A310 A332 A333 A337 A338 "
    "A339 A342 A343 A345 A346 A359 A35K A388 "
    "A3ST BCS1 BCS3 "
    // Boeing / MDD
    "B701 B703 B712 B721 B722 B731 B732 B733 B734 B735 B736 B737 B738 B739 "
    "B37M B38M B39M B3XM B741 B742 B743 B744 "
    "B748 B74R B74S B752 B753 B762 B763 B764 B772 B773 B77L B77W B778 B779 "
    "B788 B789 B78X BLCF DC10 DC85 DC86 DC87 "
    "DC91 DC93 DC94 DC95 MD11 MD81 MD82 MD83 MD87 MD88 MD90 "
    // Regional jets
    "CRJ1 CRJ2 CRJ7 CRJ9 CRJX E135 E145 E35L E45X E170 E75L E75S E190 E195 "
    "E290 E295 E390 F70 F100 RJ70 RJ85 RJ1H "
    "B461 B462 B463 J328 SU95 A148 "
    // Business jets
    "BE40 C500 C501 C510 C525 C526 C550 C551 C55B C560 C56X C650 C680 C68A "
    "C700 C750 CL30 CL35 CL60 GL5T GL6T GL7T "
    "GLEX GLF2 GLF3 GLF4 GLF5 GLF6 GA5C GA6C GA7C GA8C G150 G200 G280 ASTR "
    "GALX E50P E55P E545 E550 LJ23 LJ24 LJ25 "
    "LJ28 LJ31 LJ35 LJ36 LJ40 LJ45 LJ55 LJ60 LJ70 LJ75 LJ85 F2TH F900 FA10 "
    "FA20 FA50 FA5X FA6X FA7X FA8X H25A H25B "
    "H25C HA4T HDJT PRM1 SF50 EA50 PC24 WW24 SBR1 SBR2 MU30 BE4W "
    // Military / other
    "C17 C5 C5M K35R K35E KC10 KC46 E3TF E3CF E6 B52 B1 B2 F15 F16 F18 F18S "
    "F22 F35 A10 T38 T6 C2 P8 E737 "
    "GLID1 ";

static std::unordered_map<std::string, ACClass>& table() {
  static std::unordered_map<std::string, ACClass> t;
  static std::once_flag once;
  std::call_once(once, [] {
    std::string cur;
    for (const char* p = kJets;; ++p) {
      if (*p == ' ' || *p == 0) {
        if (!cur.empty()) t[cur] = ACClass::Jet;
        cur.clear();
        if (!*p) break;
      } else {
        cur += *p;
      }
    }
    t.erase("T6");  // T-6 Texan II is a turboprop
    t.erase("GLID1");
    t["T6"] = ACClass::Turboprop;
    // A few common props/helis for nicer labelling (not required for
    // filtering).
    for (const char* s :
         {"C172", "C152", "C182", "C206", "C208", "C210", "P28A", "P28R",
          "PA32", "PA34", "PA44", "PA46", "SR20", "SR22", "BE36", "BE35",
          "BE58", "M20P", "M20T", "DA40", "DA42", "C150", "C177", "C180",
          "C185", "PA18", "PA28", "J3",   "RV7",  "RV8",  "RV10"})
      t[s] = ACClass::Piston;
    for (const char* s :
         {"AT43", "AT45", "AT72", "AT75", "AT76", "DH8A", "DH8B",
          "DH8C", "DH8D", "SF34", "B190", "BE20", "BE30", "BE9L",
          "PC12", "TBM7", "TBM8", "TBM9", "C30J", "C130", "E120",
          "JS41", "D328", "SW4",  "P46T", "KODI", "DHC6"})
      t[s] = ACClass::Turboprop;
    for (const char* s :
         {"R22",  "R44",  "R66",  "EC35", "EC45", "EC30", "EC55", "EC75",
          "AS50", "AS55", "AS65", "A109", "A119", "A139", "A169", "B06",
          "B407", "B412", "B429", "S76",  "S92",  "H60",  "UH1",  "MD50"})
      t[s] = ACClass::Heli;
  });
  return t;
}

ACClass classify_type(const std::string& code) {
  auto& t = table();
  auto it = t.find(code);
  return it == t.end() ? ACClass::Unknown : it->second;
}

void load_type_table(const std::string& path) {
  std::ifstream f(path);
  if (!f) die("cannot read type table %s", path.c_str());
  auto& t = table();
  std::string line;
  int n = 0;
  while (std::getline(f, line)) {
    auto c = split(line, ',');
    if (c.size() < 2) continue;
    std::string code = trim(c[0]), cl = lower(trim(c[1]));
    ACClass k = ACClass::Unknown;
    if (cl.find("jet") != std::string::npos)
      k = ACClass::Jet;
    else if (cl.find("turboprop") != std::string::npos ||
             cl.find("turboshaft") != std::string::npos)
      k = ACClass::Turboprop;
    else if (cl.find("piston") != std::string::npos)
      k = ACClass::Piston;
    else if (cl.find("heli") != std::string::npos ||
             cl.find("rotor") != std::string::npos)
      k = ACClass::Heli;
    else
      continue;
    t[code] = k;
    ++n;
  }
  log("types: loaded %d entries from %s", n, path.c_str());
}

const char* class_name(ACClass c) {
  switch (c) {
    case ACClass::Jet:
      return "jet";
    case ACClass::Turboprop:
      return "turboprop";
    case ACClass::Piston:
      return "piston";
    case ACClass::Heli:
      return "heli";
    default:
      return "unknown";
  }
}

}  // namespace pl
