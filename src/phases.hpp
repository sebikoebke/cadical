#ifndef _phases_hpp_INCLUDED
#define _phases_hpp_INCLUDED

#include <vector>

namespace CaDiCaL {

struct Phases {

  std::vector<signed char> best;     // The current largest trail phase.
  std::vector<signed char> forced;   // Forced through 'phase'.
  std::vector<signed char> prev;     // Previous during local search.
  std::vector<signed char> saved;    // The actual saved phase.
  std::vector<signed char> target;   // The current target phase.
  std::vector<signed char> conflicts; // Phases of literals involved in conflicts.
};

} // namespace CaDiCaL

#endif
