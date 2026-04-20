#ifndef _stats_hpp_INCLUDED
#define _stats_hpp_INCLUDED

#include "statistics.hpp"
#include <cstdint>
#include <cstdlib>
#include <vector>

namespace CaDiCaL {

struct Internal;

struct Stats {

  Internal *internal;

#define STATISTIC(NAME, VERBOSE, COMMAND, SYMBOL, OTHER) int64_t NAME = 0;

  CADICAL_STATISTICS

#undef STATISTIC

  struct {
    double real = 0;
    double process = 0;
  } time;

  uint64_t bump_used[2] = {0, 0};
  std::vector<uint64_t> used[2] = {{}, {}}; // used clauses in focused mode
  int64_t walk_minimum;

  Stats ();
  ~Stats () = default;

  void print (Internal *);
#ifndef QUIET
  void print_old (Internal *);
  void print_new (Internal *);
#endif
};

/*------------------------------------------------------------------------*/

} // namespace CaDiCaL

#endif
