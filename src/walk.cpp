#include "walk.hpp"
#include "internal.hpp"
#include "random.hpp"

namespace CaDiCaL {

/*------------------------------------------------------------------------*/

// Random walk local search based on 'ProbSAT' ideas.

// We (based on the Master project from Leah Hohl) tried to ticks
// local search similarly to the other parts of the solver with
// limited success however.
//
// On the problem `ncc_none_5047_6_3_3_3_0_435991723', the broken part
// of walk_flip is very cheap and should not be counted in ticks, but
// on various other problems `9pipe_k' it is very important to ticks
// this part too.

//  using ClauseOrBinary = std::variant <Clause*, TaggedBinary>;

struct Walker {

  Internal *internal;
  Random random;                 // local random number generator
  int64_t ticks;                 // ticks to approximate run time
  int64_t limit;                 // limit on number of propagations
  vector<ClauseOrBinary> broken; // currently unsatisfied clauses
  double epsilon;                // smallest considered score
  vector<double> table;          // break value to score table
  vector<double> scores;         // scores of candidate literals
  vector<pair<double, int>> scores_passat;   //maybe we can safe one loop with that structure
  std::vector<int> flips; // remember the flips compared to the last best saved model
                          // in walk_passat only contain the set of flipped variables: 
                          // repair_propagation_queue reads nothing but vidx() from
                          // the entry and takes the polarity fresh from vals[]
  vector<signed char> in_flips;   // shows if a variable already in flips
  int best_trail_pos;
  size_t minimum = (size_t)(-1);
  vector<int> propagation_queue;  // our replacement for the trail in walk_passat (with a difference): the list of all literals which are assigned or should. 
                                  // Changed if probSAT_repair is finished with LS_repair: Old propagations are not longer needed, therefore ther are cut. Just pending and flipped variables remain.
  size_t propagated = 0;          // how far propagation_queue has been processed (like Internal::propagated)
  vector<int> passat_trail;       // every variable we assigned via set_val during this walk
                                  // used at cleanup to reset exactly those vals in O(assigned)
  size_t activated = 0;           // counts assigned variables; up_expansion stops once all are
                                  // activated. Replaces ordering_O/ls_score: variable + polarity
                                  // selection is delegated to CaDiCaL's own decision heuristic.
  size_t pre_assigned = 0;        // vars already assigned at build time (root fixed/units);
                                  // counted in 'activated' but NOT activated by PASSAT
  size_t activatable = 0;         // active & unassigned vars: the universe PASSAT can decide/propagate
  vector<vector<int>> passat_lookup_table; // positions in `clauses` where v+/v- occurs
  vector<int> broken_clauses;     // all currently broken (conflict_counter == 0) clauses, used for probSAT_repair;
                                  // maintained incrementally by passat_assign and flip_and_repair
  vector<int> broken_pos;         // position of a clause inside broken_clauses (indexed by clause pos), -1 if not broken; enables fast removal (like in WalkerFO)
  vector<int> conflict_counter;   // counter which shows if there is a conflict inside a clause, if c_c == 0 => conflict, decreased if the opposite polarity is assigned to true
  vector<int> notfalse_xor;       // XOR of the literals conflict_counter counts.
                                  // This enable to finde the critical literal in a very fast way:
                                  // If the conflict_counter of clause c decrease to 1,
                                  // the number (xor value) in notfalse_xor[c] is equal 
                                  // to the last not false liter in the clause c
                                  // looking up is just a simple lookup and not a scanning proccess.
                                  // To be able to do that trick, a clause is not allowed to contain duplicates
                                  // (CaDiCal guarantee it)
  vector<int> bv;                 // incrementally exact break values, indexed by vlit:
                                  // bv[vlit(l)] = # clauses with conflict_counter == 1 whose unique
                                  // not-false (critical) literal is -l (=> flipping l to true breaks a clause)
  vector<int> broken_occ;         // occurence of a literal in a broken clause
  vector<int> clauses_critical_literal; // per clause: the critical literal (signed, like on the trail)
                                        // avoids a rescanning process of the clause in flip_and_repair
  vector<int> lsl;                // Last Satisfied Literal per literal, show how often the literal is the last satisfied literal per clause
  vector<int> sat_critical_lit;   // per clause: the unique true literal when satisfied_counter == 1, else 0
                                  // sat-side mirror of clauses_critical_literal, avoids rescanning the clause
  vector<int> flip_count;         // LS Hotspots, indexed by variables
  vector<signed char> mark;       // per-variable dedup flag, invariant 0 outside repair_propagation_queue
  vector<int> cache_queue;        // reusable cache to rebuild propagation_queue without allocating

  bool track_probSAT_repair = false;     // set to true, if we want to track the steps of the probSAT_repair (Local Search) modul
  std::vector<int> measure_start_assignment; // signed assignment snapshot taken at "Start Repair";
                                             // diffed against the "End Repair" assignment to report the
                                             // net (final) flips of one repair round
  bool cheap_break_value = true;   // if true, we calc a cheaper break value => O(1) instead of O(|clauses[-lit]|)
  bool track_break_value = false;  // set to true to look at the (real, cheap) break-value pair of every flippable
                                   // literal of each picked broken clause to break_value_measure.csv
                                   // Be aware, measuring a walk_passat version where the real break value is used make no sense,
                                   // because you should not find a difference
  size_t passat_expansion_barrier = 100; // upper barrier for the expansion the idea is to switch from time to time between expansion and repair,
                                         // because in some old cases we did only expansion and in some only repair, both had a poorer performance
                                         // note: if expansion == 0, the expansion is not limitted
                                         // good results with 50 and 100, really bad with 75
  bool use_up_expansion = false;  // if true, the main loop uses the original up_expansion
                                  // therefore assign until the first conflict arise instead of advanced_expansion
                                  // selected via --walkpassat=7

  bool dynamic_barrier = false;   // walkpassat=15 use a (dynamic) soft adaptive barrier toggling 1% <-> 10%
  size_t expansion_conflict_counter = 0;   // conflicts accumulated during the current advanced_expansion run
  size_t last_expansion_conflicts = 0;     // conflicts of the previous run, for the 20% comparison
  int dynamic_counter = 0;  // (unused) legacy counter from the earlier 50% ramp-up scheme, kept for ABI stability
  double avg_clause_size = 0.0; // average tracked clause length (set in passat_build)

  bool passat_track_improvement = false; // walkpassat=16 use the option to write the best found assignment (fewest broken clauses) in phases_best
  size_t last_start_broken = 0;          // broken clauses at the start of the last probSAT_repair
  size_t last_min_broken = 0;            // fewest broken clauses reached during the last repair
  std::vector<signed char> best_repair_model; 

  bool anti_stagnation = false;   // true if we want to expand further even if not all conflicts areresolved => goal: full assignment over complete conflict solving
                                  // on the one hand, with more variables we could find easier a solution via flipping for stagnating problems
                                  // on the other hand, with more variables the number of conflicts increase in the flipping process
  bool resolved_conflicts = true; // did the last probSAT_repair end with broken == 0? 
  size_t stagnation_counter = 0;  // flips since the last improvement of min_broken, reset at each repair start
  bool assumption_unsat = false;  // repair hit an assumption-only broken clause -> unrepairable, hard stop

  bool increased_passat_limit = false; // walkpassat=18..21: multiply the tick limit by 3

  bool autarky_mode = false;
  bool maintain_pick_stats = false; // if true we also count broken_occ[lit], lsl[lit] and
                                    // sat_critical_lit[clause] -- the extra bookkeeping that the
                                    // alternative pick scores (passat_score_mode != 0) read.
                                    // Set ONCE, right before passat_build, because the build
                                    // already has to initialise lsl/sat_critical_lit.
  bool advanced_pure_finding = false;  // if true we check also for pure literals after each repair round 
                                       // the idea is to find literals which become pure after finding and eliminating an autarky
  bool tuc_min_autarky_check = false;  // run build_autarky whenever |tuc_clauses| halves during repair
  int64_t tuc_autarky_ticks = 0;       // ticks spent in those in-repair build_autarky calls, reset per repair
  int passat_score_mode = 0;           // is used to determine which scoring mode is used for probSAT_pick_lit 
                                       // 0 = base^bv, 1 = base^lsl, 2 = broken_occ, 3 = broken_occ * base^bv
  bool advanced_picking_mode = false;  // if true, instead of picking just one clause, we pick clauses until we've found 10 flippable literals 
  bool autarky_check_expansion = true; // run build_autarky after a conflicting expansion
  bool autarky_check_repair = true;    // run build_autarky after a successful repair
  vector<int> satisfied_counter;  // satisfied_counter is initialized with 0 for every clause, is increased by one if one literals is assigned with true in the clause
  vector<int> tuc_clauses;        // tuc: touched unsatisfied clauses, a trail like broken clauses, with all clauses that are touched but not satisfied yet
                                  // if you want to find the broken clauses in tuc you could simply check the conflict_counter of the current clause
                                  // if we dont use broken clauses, we could also build a second array with pointers on the broken_clauses in tuc
  vector<int> tuc_pos;            // position of a clause in tuc_clauses for faster lookups, initialized with -1
  vector<int> autarky_set;        // set with all literals that build an autarky
                                  // the advantage is, that we could bould easily a set of clauses with passat_lookup_table if we want to generate the autarky clauses
  bool show_autarky = false;
  vector<signed char> unflippable;  // show if a literal is in an autarky => unflippable[lit] == 1 => lit is in an autarky
  vector<signed char> unvisitable;  // show if a clause is fullfilled by alsan autarky => unvisitable[lit] == 1 => clause is fullfilled by an autarky
  vector<signed char> not_autark;   // set of variables that cant be autark
  vector<signed char> pure_lits;   // pure_lits is a set of variables where passat_lookup_table[lit] == 0 
                                   // or passat_lookup_table[-lit] == 0 occure, this set can be fixed before the main walk_passat loop
  vector<int> autarky_worklist;     // clauses that lost their supporter and have to be processed
  bool autarky_elimination_mode = false; // true if we want to eliminate the autarky found by walk_passat
  bool frozen_autarky = false;          // true if the autarky contain a frozen literal => autarky cant be eliminated
  vector<int> autarky_trail;
  vector<signed char> autarky_val;
  vector<int> peel_counter;       // copy of the satisfied counter, for cheaper autarky peeling
  vector<unsigned> peel_stamp;    // timestamp which show if we already visited a clause during peeling
  unsigned peel_generation = 0;   // current peel generation

  std::vector<signed char> best_values; // best model stored so far
  double score (unsigned);              // compute score from break count
#ifndef NDEBUG
  std::vector<signed char> current_best_model; // best model found so far
#endif

  Walker (Internal *, int64_t limit);
  void populate_table (double size);


  // for efficiency, storing the model each time an improvement is
  // found is too costly. Instead we store some of the flips since
  // last time and the position of the best model found so far.
  //
  // Once that buffer is full, we save the last best assignment. This is
  // particularly important at the beginning where every or nearly every flip
  // improves the current assignment.

  // Push the literal on the buffer of flipped literals
  void push_flipped (int flipped);
  // save the best assignment found so far, if any was found
  void save_walker_trail (bool);
  // export the best model from walk to the main solver
  void save_final_minimum (size_t old_minimum);
};

// These are in essence the CB values from Adrian Balint's thesis.  They
// denote the inverse 'cb' of the base 'b' of the (probability) weight
// 'b^-i' for picking a literal with the break value 'i' (first column is
// the 'size', second the 'CB' value).

static double cbvals[][2] = {
    {0.0, 2.00}, {3.0, 2.50}, {4.0, 2.85}, {5.0, 3.70},
    {6.0, 5.10}, {7.0, 7.40}, // Adrian has '5.4', but '7.4' looks better.
};

static const int ncbvals = sizeof cbvals / sizeof cbvals[0];

// We interpolate the CB values for uniform random SAT formula to the non
// integer situation of average clause size by piecewise linear functions.
//
//   y2 - y1
//   ------- * (x - x1) + y1
//   x2 - x1
//
// where 'x' is the average size of clauses and 'y' the CB value.

inline static double fitcbval (double size) {
  int i = 0;
  while (i + 2 < ncbvals &&
         (cbvals[i][0] > size || cbvals[i + 1][0] < size))
    i++;
  const double x2 = cbvals[i + 1][0], x1 = cbvals[i][0];
  const double y2 = cbvals[i + 1][1], y1 = cbvals[i][1];
  const double dx = x2 - x1, dy = y2 - y1;
  assert (dx);
  const double res = dy * (size - x1) / dx + y1;
  assert (res > 0);
  return res;
}

// Initialize the data structures for one local search round.

Walker::Walker (Internal *i, int64_t l)
    : internal (i), random (internal->opts.seed), // global random seed
      ticks (0), limit (l), epsilon(-1), best_trail_pos (-1) {
  random += internal->stats.walk.count + internal->stats.walk.passat; // different seed every time
  flips.reserve (i->max_var / 4);
  mark.resize (i->max_var + 1, 0); // dedup flag for repair_propagation_queue, kept invariant 0
  flip_count.resize (i->max_var + 1, 0);
  best_values.resize (i->max_var + 1, 0);
#ifndef NDEBUG
  current_best_model.resize (i->max_var + 1, 0);
#endif
}

void Walker::populate_table (double size) {
  // This is the magic constant in ProbSAT (also called 'CB'), which we pick
  // according to the average size every second invocation and otherwise
  // just the default '2.0', which turns into the base '0.5'.
  //
  bool use_size_based_cb;
  if (internal->stats.walk.passat){
    use_size_based_cb = (internal->stats.walk.passat & 1);
  } else {
    use_size_based_cb = (internal->stats.walk.count & 1);
  }
  const double cb = use_size_based_cb ? fitcbval (size) : 2.0;
  assert (cb);
  const double base = 1 / cb; // scores are 'base^0,base^1,base^2,...

  double next = 1;
  for (epsilon = next; next; next = epsilon * base)
    table.push_back (epsilon = next);

  PHASE ("walk", internal->stats.walk.count,
         "CB %.2f with inverse %.2f as base and table size %zd", cb, base,
         table.size ());
}

// Add the literal to flip to the queue

void Walker::push_flipped (int flipped) {
  LOG ("push literal %s on the flips", LOGLIT (flipped));
  assert (flipped);
  if (best_trail_pos < 0) {
    LOG ("not pushing flipped %s to already invalid trail",
         LOGLIT (flipped));
    return;
  }

  const size_t size_trail = flips.size ();
  const size_t limit = internal->max_var / 4 + 1;
  if (size_trail < limit) {
    flips.push_back (flipped);
    LOG ("pushed flipped %s to trail which now has size %zd",
         LOGLIT (flipped), size_trail + 1);
    return;
  }

  if (best_trail_pos) {
    LOG ("trail reached limit %zd but has best position %d", limit,
         best_trail_pos);
    save_walker_trail (true);
    flips.push_back (flipped);
    LOG ("pushed flipped %s to trail which now has size %zu",
         LOGLIT (flipped), flips.size ());
    return;
  } else {
    LOG ("trail reached limit %zd without best position", limit);
    flips.clear ();
    LOG ("not pushing %s to invalidated trail", LOGLIT (flipped));
    best_trail_pos = -1;
    LOG ("best trail position becomes invalid");
  }
}

void Walker::save_walker_trail (bool keep) {
  assert (best_trail_pos != -1);
  assert ((size_t) best_trail_pos <= flips.size ());
//  assert (!keep || best_trail_pos == flips.size());
#ifdef LOGGING
  const size_t size_trail = flips.size ();
#endif
  const int kept = flips.size () - best_trail_pos;
  LOG ("saving %d values of flipped literals on trail of size %zd",
       best_trail_pos, flips.size ());

  const auto begin = flips.begin ();
  const auto best = flips.begin () + best_trail_pos;
  const auto end = flips.end ();

  auto it = begin;
  for (; it != best; ++it) {
    const int lit = *it;
    assert (lit);
    const signed char value = sign (lit);
    const int idx = std::abs (lit);
    best_values[idx] = value;
  }
  if (!keep) {
    LOG ("no need to shift and keep remaining %u literals", kept);
    return;
  }

#ifndef NDEBUG
  for (auto v : internal->vars) {
    if (internal->active (v))
      assert (best_values[v] == current_best_model[v]);
  }
#endif
  LOG ("flushed %u literals %.0f%% from trail", best_trail_pos,
       percent (best_trail_pos, size_trail));
  assert (it == best);
  auto jt = begin;
  for (; it != end; ++it, ++jt) {
    assert (jt <= it);
    assert (it < end);
    *jt = *it;
  }

  assert ((int) (end - jt) == best_trail_pos);
  assert ((int) (jt - begin) == kept);
  flips.resize (kept);
  LOG ("keeping %u literals %.0f%% on trail", kept,
       percent (kept, size_trail));
  LOG ("reset best trail position to 0");
  best_trail_pos = 0;
}

// finally export the final minimum
void Walker::save_final_minimum (size_t old_init_minimum) {
  assert (minimum <= old_init_minimum);
#ifdef NDEBUG
  (void) old_init_minimum;
#endif

  if (!best_trail_pos || best_trail_pos == -1)
    LOG ("minimum already saved");
  else
    save_walker_trail (false);

  ++internal->stats.walk.improved;
  for (auto v : internal->vars) {
    if (best_values[v])
      internal->phases.saved[v] = best_values[v];
    else
      assert (!internal->active (v));
  }
  internal->copy_phases (internal->phases.best);
  internal->copy_phases (internal->phases.prev);
}

// The scores are tabulated for faster computation (to avoid 'pow').

inline double Walker::score (unsigned i) {
  const double res = (i < table.size () ? table[i] : epsilon);
  LOG ("break %u mapped to score %g", i, res);
  return res;
}

/*------------------------------------------------------------------------*/

ClauseOrBinary Internal::walk_pick_clause (Walker &walker) {
  require_mode (WALK);
  assert (!walker.broken.empty ());
  size_t size = walker.broken.size ();
  if (size > INT_MAX)
    size = INT_MAX;
  int pos = walker.random.pick_int (0, (int)size - 1);
  ClauseOrBinary res = walker.broken[pos];
#ifdef LOGGING
  Clause *c;
  if (!res.is_binary ())
    c = res.clause ();
  else
    c = res.tagged_binary ().d;
  LOG (c, "picking random position %d", pos);
#endif
  return res;
}

/*------------------------------------------------------------------------*/

// Compute the number of clauses which would be become unsatisfied if 'lit'
// is flipped and set to false.  This is called the 'break-count' of 'lit'.

unsigned Internal::walk_break_value (int lit, int64_t &ticks) {
  require_mode (WALK);
  START (walkbreak);
  assert (val (lit) > 0);
  const int64_t oldticks = ticks;

  unsigned res = 0; // The computed break-count of 'lit'.
  ticks += (1 + cache_lines (watches (lit).size (), sizeof (Watch)));

  for (auto &w : watches (lit)) {
    assert (w.blit != lit);
    if (val (w.blit) > 0)
      continue;
    if (w.binary ()) {
      res++;
      continue;
    }

    Clause *c = w.clause;
#ifdef LOGGING
    assert (c != dummy_binary);
#endif
    ++ticks;

    assert (lit == c->literals[0]);

    // Now try to find a second satisfied literal starting at 'literals[1]'
    // shifting all the traversed literals to right by one position in order
    // to move such a second satisfying literal to 'literals[1]'.  This move
    // to front strategy improves the chances to find the second satisfying
    // literal earlier in subsequent break-count computations.
    //
    auto begin = c->begin () + 1;
    const auto end = c->end ();
    auto i = begin;
    int prev = 0;
    while (i != end) {
      const int other = *i;
      *i++ = prev;
      prev = other;
      if (val (other) < 0)
        continue;

      // Found 'other' as second satisfying literal.

      w.blit = other; // Update 'blit'
      *begin = other; // and move to front.

      break;
    }

    if (i != end)
      continue; // Double satisfied!

    // Otherwise restore literals (undo shift to the right).
    //
    while (i != begin) {
      const int other = *--i;
      *i = prev;
      prev = other;
    }
    res++; // Literal 'lit' single satisfies clause 'c'.
  }
  stats.ticks.walkbreak += (ticks - oldticks);
  STOP (walkbreak);

  return res;
}

/*------------------------------------------------------------------------*/

// Given an unsatisfied clause 'c', in which we want to flip a literal, we
// first determine the exponential score based on the break-count of its
// literals and then sample the literals based on these scores.  The CB
// value is smaller than one and thus the score is exponentially decreasing
// with the break-count increasing.  The sampling works as in 'ProbSAT' and
// 'YalSAT' by summing up the scores and then picking a random limit in the
// range of zero to the sum, then summing up the scores again and picking
// the first literal which reaches the limit.  Note, that during incremental
// SAT solving we can not flip assumed variables.  Those are assigned at
// decision level one, while the other variables are assigned at two.

int Internal::walk_pick_lit (Walker &walker, Clause *c) {
  LOG ("picking literal by break-count");
  assert (walker.scores.empty ());
  const int64_t old = walker.ticks;
  walker.ticks++;
  double sum = 0;
  int64_t propagations = 0;
  for (const auto lit : *c) {
    assert (active (lit));
    if (var (lit).level == 1) {
      LOG ("skipping assumption %d for scoring", -lit);
      continue;
    }
    assert (active (lit));
    propagations++;
    unsigned tmp = walk_break_value (-lit, walker.ticks);
    double score = walker.score (tmp);
    LOG ("literal %d break-count %u score %g", lit, tmp, score);
    walker.scores.push_back (score);
    sum += score;
  }
  (void) propagations; // TODO actually unused?
  LOG ("scored %zd literals", walker.scores.size ());
  assert (!walker.scores.empty ());
  assert (walker.scores.size () <= (size_t) c->size);
  const double lim = sum * walker.random.generate_double ();
  LOG ("score sum %g limit %g", sum, lim);
  const auto end = c->end ();
  auto i = c->begin ();
  auto j = walker.scores.begin ();
  int res;
  for (;;) {
    assert (i != end);
    res = *i++;
    if (var (res).level > 1)
      break;
    LOG ("skipping assumption %d without score", -res);
  }
  sum = *j++;
  while (sum <= lim && i != end) {
    res = *i++;
    if (var (res).level == 1) {
      LOG ("skipping assumption %d without score", -res);
      continue;
    }
    sum += *j++;
  }
  walker.scores.clear ();
  LOG ("picking literal %d by break-count", res);
  stats.ticks.walkpick += walker.ticks - old;
  return res;
}

int Internal::walk_pick_lit (Walker &walker, ClauseOrBinary c) {
  if (c.is_binary ())
    return walk_pick_lit (walker, c.tagged_binary ());
  return walk_pick_lit (walker, c.clause ());
}

int Internal::walk_pick_lit (Walker &walker, const TaggedBinary c) {
  LOG ("picking literal by break-count on binary clause [%" PRIu64 "]%s %s",
       c.d->id, LOGLIT (c.lit), LOGLIT (c.other));
  assert (walker.scores.empty ());
  const int64_t old = walker.ticks;
  double sum = 0;
  int64_t propagations = 0;
  const std::array<int, 2> clause = {c.lit, c.other};
  for (const auto lit : clause) {
    assert (active (lit));
    if (var (lit).level == 1) {
      LOG ("skipping assumption %d for scoring", -lit);
      continue;
    }
    assert (active (lit));
    assert (val (lit) < 0);
    propagations++;
    unsigned tmp = walk_break_value (-lit, walker.ticks);
    double score = walker.score (tmp);
    LOG ("literal %d break-count %u score %g", lit, tmp, score);
    walker.scores.push_back (score);
    sum += score;
  }
  (void) propagations; // TODO unused?
  LOG ("scored %zd literals", walker.scores.size ());
  assert (!walker.scores.empty ());
  assert (walker.scores.size () <= (size_t) 2);
  const double lim = sum * walker.random.generate_double ();
  LOG ("score sum %g limit %g", sum, lim);
  const auto end = clause.end ();
  auto i = clause.begin ();
  auto j = walker.scores.begin ();
  int res = 0;
  for (;;) {
    assert (i != end);
    res = *i++;
    if (var (res).level > 1)
      break;
    LOG ("skipping assumption %d without score", -res);
  }
  sum = *j++;
  while (sum <= lim && i != end) {
    res = *i++;
    if (var (res).level == 1) {
      LOG ("skipping assumption %d without score", -res);
      continue;
    }
    sum += *j++;
  }
  assert (res);
  walker.scores.clear ();
  LOG ("picking literal %d by break-count", res);
  stats.ticks.walkpick += walker.ticks - old;
  return res;
}

/*------------------------------------------------------------------------*/

// flips a literal unless we run out of ticks.
bool Internal::walk_flip_lit (Walker &walker, int lit) {
  START (walkflip);
  const int64_t old = walker.ticks;
  require_mode (WALK);
  LOG ("flipping assign %d", lit);
  assert (val (lit) < 0);

  // First flip the literal value.
  //
  const signed char tmp = sign (lit);
  const int idx = abs (lit);
  set_val (idx, tmp);
  assert (val (lit) > 0);

  // we are going to need it anyway and it probably still is in memory
  const Watches &ws = watches (-lit);
  if (!ws.empty ()) {
    const Watch &w = ws[0];
    __builtin_prefetch (&w, 0, 1);
  }

  // Then remove 'c' and all other now satisfied (made) clauses.
  {
    // Simply go over all unsatisfied (broken) clauses.

    LOG ("trying to make %zd broken clauses", walker.broken.size ());

    const auto eou = walker.broken.end ();
    auto j = walker.broken.begin (), i = j;
    // broken is in cache given how central it is... but not always (see the
    // ncc problems). Value was heuristically determined to give reasonnable
    // values.
    walker.ticks +=
        1 + cache_lines (walker.broken.size (), sizeof (*i));
#if defined(LOGGING) || !defined(NDEBUG)
    int64_t made = 0;
#endif

    while (i != eou) {

      ClauseOrBinary tagged = *j++ = *i++;

      if (tagged.is_binary ()) {
        const TaggedBinary &b = tagged.tagged_binary ();
        const int clit = b.lit;
        const int other = b.other;
        assert (val (clit) < 0 || val (other) < 0);
#if defined(LOGGING)
        assert (b.d->literals[0] == clit || b.d->literals[1] == clit);
        assert (b.d->literals[0] == other || b.d->literals[1] == other);
#endif
        if (clit == lit || other == lit) {
          LOG (b.d, "made");
          const int first_lit = lit;
          const int second_lit = clit ^ lit ^ other;
#ifdef LOGGING
          watch_binary_literal (first_lit, second_lit, b.d);
#else
          // placeholder for the clause, does not matter
          watch_binary_literal (first_lit, second_lit, dummy_binary);
#endif

          ++walker.ticks;
#if defined(LOGGING) || !defined(NDEBUG)
          made++;
#endif
          j--;
        } else {
          LOG (b.d, "still broken");
          assert (val (clit) < 0 && val (other) < 0);
        }
        continue;
      }

      // now the expansive part
      Clause *d = tagged.clause ();
      ++walker.ticks;
      int *literals = d->literals;
      LOG (d, "search for replacement");
      int prev = 0;
      // Find 'lit' in 'd'.
      //
      const int size = d->size;
      for (int i = 0; i < size; i++) {
        const int other = literals[i];
        assert (active (other));
        literals[i] = prev;
        prev = other;
        if (other == lit)
          break;
        assert (val (other) < 0);
      }
      // If 'lit' is in 'd' then move it to the front to watch it.
      //
      if (prev == lit) {
        literals[0] = lit;
        LOG (d, "made");
        watch_literal (literals[0], literals[1], d);
        ++walker.ticks;
#if defined(LOGGING) || !defined(NDEBUG)
        made++;
#endif
        j--;
      } else { // Otherwise the clause is not satisfied, undo shift.

        for (int i = size - 1; i >= 0; i--) {
          int other = literals[i];
          literals[i] = prev;
          prev = other;
        }
      }
      LOG (d, "clause after undoing shift");
    }
    assert ((int64_t) (j - walker.broken.begin ()) + made ==
            (int64_t) walker.broken.size ());
    walker.broken.resize (j - walker.broken.begin ());
    LOG ("made %" PRId64 " clauses by flipping %d, still %zu broken", made,
         lit, walker.broken.size ());
#ifndef NDEBUG
    for (auto d : walker.broken) {
      if (d.is_binary ()) {
        const TaggedBinary &b = d.tagged_binary ();
        assert (val (b.lit) < 0 && val (b.other) < 0);
      } else {
        for (auto lit : *d.clause ())
          assert (val (lit) < 0);
      }
    }
#endif
    if (walker.ticks > walker.limit) {
      STOP (walkflip);
      return false;
    }
  }

  stats.ticks.walkflipbroken += walker.ticks - old;

  const int64_t old_after_broken = walker.ticks;

  // Finally add all new unsatisfied (broken) clauses.
  {
#ifdef LOGGING
    int64_t broken = 0;
#endif
    Watches &ws = watches (-lit);
    // probably still in cache
    walker.ticks += 1 + cache_lines (ws.size (), sizeof (Watch));

    LOG ("trying to break %zd watched clauses", ws.size ());

    for (const auto &w : ws) {
      Clause *d = w.clause;
      const bool binary = w.binary ();
      if (binary) {
        const int other = w.blit;
        assert (w.blit != -lit);
        if (val (other) > 0) {
          LOG (d, "unwatch %d in", -lit);
          watch_binary_literal (other, -lit, d);
          ++walker.ticks;
          continue;
        }
        LOG (d, "broken");
#ifdef LOGGING
        assert (d != dummy_binary);
#endif
        walker.broken.push_back (TaggedBinary (d, -lit, other));
        ++walker.ticks;
#ifdef LOGGING
        broken++;
#endif
        continue;
      }

      if (walker.ticks > walker.limit) {
        STOP (walkflip);
        return false;
      }
      // now the expansive part
      assert (d->size != 2);
      ++walker.ticks;
      int *literals = d->literals, replacement = 0, prev = -lit;
      assert (d->size == w.size);
      const int size = d->size;
      assert (literals[0] == -lit);

      for (int i = 1; i < size; i++) {
        const int other = literals[i];
        assert (active (other));
        literals[i] = prev; // shift all to right
        prev = other;
        const signed char tmp = val (other);
        if (tmp < 0)
          continue;
        replacement = other; // satisfying literal
        break;
      }
      if (replacement) {
        assert (-lit != replacement);
        literals[1] = -lit;
        literals[0] = replacement;
        watch_literal (replacement, -lit, d);
        ++walker.ticks;
        LOG (d, "found replacement");
      } else {
        for (int i = size - 1; i > 0; i--) { // undo shift
          const int other = literals[i];
          literals[i] = prev;
          prev = other;
        }

        assert (literals[0] == -lit);
        LOG (d, "broken");
        walker.broken.push_back (d);
        ++walker.ticks;
#ifdef LOGGING
        broken++;
#endif
      }
    }
    LOG ("broken %" PRId64 " clauses by flipping %d", broken, lit);
    ws.clear ();
  }
  STOP (walkflip);
  stats.ticks.walkflipWL += walker.ticks - old_after_broken;
  stats.ticks.walkflip += walker.ticks - old;
  return true;
}

/*------------------------------------------------------------------------*/

// Check whether to save the current phases as new global minimum.

inline void Internal::walk_save_minimum (Walker &walker) {
  size_t broken = walker.broken.size ();
  if (broken >= walker.minimum)
    return;
  if (broken <= stats.walk.minimum) {
    stats.walk.minimum = broken;
    VERBOSE (3, "new global minimum %zd", broken);
  } else {
    VERBOSE (3, "new walk minimum %zd", broken);
  }

  walker.minimum = broken;

#ifndef NDEBUG
  for (auto i : vars) {
    const signed char tmp = vals[i];
    if (tmp)
      walker.current_best_model[i] = tmp;
  }
  if (walker.minimum == 0) {
    for (auto c : clauses) {
      if (c->garbage)
        continue;
      if (c->redundant)
        continue;
      int satisfied = 0;
      for (const auto &lit : *c) {
        const int tmp = internal->val (lit);
        if (tmp > 0) {
          LOG (c, "satisfied literal %d in", lit);
          satisfied++;
        }
      }
      assert (satisfied);
    }
  }
#endif
  if (walker.best_trail_pos == -1) {
    VERBOSE (3, "saving the new walk minimum %zd", broken);
    for (auto i : vars) {
      const signed char tmp = vals[i];
      if (tmp) {
        walker.best_values[i] = tmp;
#ifndef NDEBUG
        assert (tmp == walker.current_best_model[i]);
#endif
      } else {
        assert (!active (i));
      }
    }
    walker.best_trail_pos = 0;
  } else {
    walker.best_trail_pos = walker.flips.size ();
    LOG ("new best trail position %u", walker.best_trail_pos);
  }
}

/*------------------------------------------------------------------------*/

int Internal::walk_round (int64_t limit, bool prev) {

  stats.walk.count++;
  std::vector<int> propagated;
  bool failed = false; // Inconsistent assumptions?
  int res = decide_and_propagate_all_assumptions (propagated);
  if (res) {
    return res;
  }
  clear_watches ();

  // Remove all fixed variables first (assigned at decision level zero).
  //
  if (last.collect.fixed < stats.all.fixed)
    garbage_collection ();

#ifndef QUIET
  // We want to see more messages during initial local search.
  //
  if (localsearching) {
    assert (!force_phase_messages);
    force_phase_messages = true;
  }
#endif

  PHASE ("walk", stats.walk.count, "random walk limit of %" PRId64 " ticks",
         limit);

  // Instantiate data structures for this local search round.
  Walker walker (internal, limit);
#ifndef QUIET
  size_t old_global_minimum = stats.walk.minimum;
#endif

  level = 1; // Assumed variables assigned at level 1.

  if (assumptions.empty ()) {
    LOG ("no assumptions so assigning all variables to decision phase");
  } else {
    LOG ("assigning assumptions to their forced phase first");
    std::vector<int> propagated;
    for (auto lit : trail)
    propagated.push_back(lit);

    for (const auto lit : propagated) {
      signed char tmp = val (lit);
      if (tmp > 0)
        continue;
      assert (tmp == 0);
      if (!active (lit))
        continue;
      tmp = sign (lit);
      const int idx = abs (lit);
      LOG ("initial assign %d to assumption phase", tmp < 0 ? -idx : idx);
      set_val (idx, tmp);
      assert (level == 1);
      var (idx).level = 1;
    }
    if (!failed)
      LOG ("now assigning remaining variables to their decision phase");
  }

  level = 2; // All other non assumed variables assigned at level 2.

  if (!failed) {

    // warmup stores the result in phases, not in target
    const bool target = opts.warmup ? false : stable || opts.target == 2;
    for (auto idx : vars) {
      if (!active (idx)) {
        LOG ("skipping inactive variable %d", idx);
        continue;
      }
      if (vals[idx]) {
        assert (var (idx).level == 1);
        LOG ("skipping assumed variable %d", idx);
        continue;
      }
      int tmp = 0;
      if (prev)
        tmp = phases.prev[idx];
      if (!tmp)
        tmp = sign (decide_phase (idx, target));
      assert (tmp == 1 || tmp == -1);
      set_val (idx, tmp);
      assert (level == 2);
      var (idx).level = 2;
      LOG ("initial assign %d to decision phase", tmp < 0 ? -idx : idx);
    }

    LOG ("watching satisfied and registering broken clauses");
#ifdef LOGGING
    size_t watched = 0;
#endif

    double size = 0;
    int64_t n = 0;
    for (const auto c : clauses) {

      if (c->garbage)
        continue;
      if (c->redundant) {
        if (opts.walkredundant == 1 && c->size == 2 && !c->hyper)
	  LOG (c, "importing binary clause");
        else if (opts.walkredundant == 0)
          continue;
        if (!likely_to_be_kept_clause (c))
          continue;
      }

      bool satisfiable = false; // contains not only assumptions
      int satisfied = 0;        // clause satisfied?

      int *lits = c->literals;
      size += c->size;
      n++;
      const int size = c->size;

      // Move to front satisfied literals and determine whether there
      // is at least one (non-assumed) literal that can be flipped.
      //
      for (int i = 0; satisfied < 2 && i < size; i++) {
        const int lit = lits[i];
        assert (active (lit)); // Due to garbage collection.
        if (val (lit) > 0) {
          swap (lits[satisfied], lits[i]);
          if (!satisfied++)
            LOG ("first satisfying literal %d", lit);
        } else if (!satisfiable && var (lit).level > 1) {
          LOG ("non-assumption potentially satisfying literal %d", lit);
          satisfiable = true;
        }
      }

      if (!satisfied && !satisfiable) {
        LOG (c, "due to assumptions unsatisfiable");
        LOG ("stopping local search since assumptions falsify a clause");
        failed = true;
        break;
      }

      if (satisfied) {
        LOG (c, "pushing to satisfied");
        if (c->size == 2)
          watch_binary_literal (lits[0], lits[1], c);
        else
          watch_literal (lits[0], lits[1], c);
#ifdef LOGGING
        watched++;
#endif
      } else {
        assert (satisfiable); // at least one non-assumed variable ...
        LOG (c, "broken");
        assert (c->size == size);
        if (size == 2)
          walker.broken.push_back (TaggedBinary (c));
        else
          walker.broken.push_back (c);
      }
    }

    double average_size = relative (size, n);
    walker.populate_table (average_size);
    PHASE ("walk", stats.walk.count,
           "%" PRId64 " clauses average size %.2f over %d variables", n,
           average_size, active ());

#ifdef LOGGING
    if (!failed) {
      size_t broken = walker.broken.size ();
      size_t total = watched + broken;
      LOG ("watching %zd clauses %.0f%% "
           "out of %zd (watched and broken)",
           watched, percent (watched, total), total);
    }
#endif
  }

  assert (failed || walker.table.size ());

  if (!failed) {

    size_t broken = walker.broken.size ();
    size_t initial_minimum = broken;

    PHASE ("walk", stats.walk.count,
           "starting with %zd unsatisfied clauses "
           "(%.0f%% out of %" PRId64 ")",
           broken, percent (broken, stats.current.irredundant),
           stats.current.irredundant);

    walk_save_minimum (walker);
    assert (stats.walk.minimum <= walker.minimum);

    size_t minimum = broken;
#ifndef QUIET
    int64_t flips = 0;
#endif
    while (!terminated_asynchronously () && !walker.broken.empty () &&
           walker.ticks < walker.limit) {
#ifndef QUIET
      flips++;
#endif
      stats.walk.flips++;
      stats.walk.broken += (int64_t)broken;
      ClauseOrBinary c = walk_pick_clause (walker);
      const int lit = walk_pick_lit (walker, c);
      bool finished = walk_flip_lit (walker, lit);
      if (!finished)
        break;
      walker.push_flipped (lit);
      broken = walker.broken.size ();
      LOG ("now have %zd broken clauses in total", broken);
      if (broken >= minimum)
        continue;
      minimum = broken;
      VERBOSE (3, "new phase minimum %zd after %" PRId64 " flips",
               minimum, flips);
      walk_save_minimum (walker);
    }

    walker.save_final_minimum (initial_minimum);

#ifndef QUIET
    if (minimum == initial_minimum) {
      PHASE ("walk", internal->stats.walk.count,
             "%sno improvement %zd%s in %" PRId64 " flips and "
             "%" PRId64 " ticks",
             tout.bright_yellow_code (), minimum, tout.normal_code (),
             flips, walker.ticks);
    } else if (minimum < old_global_minimum)
      PHASE ("walk", stats.walk.count,
             "%snew global minimum %zd%s in %" PRId64 " flips and "
             "%" PRId64 " ticks",
             tout.bright_yellow_code (), minimum, tout.normal_code (),
             flips, walker.ticks);
    else
      PHASE ("walk", stats.walk.count,
             "best phase minimum %zd in %" PRId64 " flips and "
             "%" PRId64 " ticks",
             minimum, flips, walker.ticks);

    if (opts.profile >= 2) {
      PHASE ("walk", stats.walk.count, "%.2f million ticks per second",
             1e-6 *
                 relative (walker.ticks, time () - profiles.walk.started));

      PHASE ("walk", stats.walk.count, "%.2f thousand flips per second",
             relative (1e-3 * flips, time () - profiles.walk.started));

    } else {
      PHASE ("walk", stats.walk.count, "%.2f ticks", 1e-6 * walker.ticks);

      PHASE ("walk", stats.walk.count, "%.2f thousand flips", 1e-3 * flips);
    }
#endif

    if (minimum > 0) {
      LOG ("minimum %zd non-zero thus potentially continue",
           minimum);
      res = 0;
    } else {
      LOG ("minimum is zero thus stop local search");
      res = 10;
    }

  } else {

    res = 20;

    PHASE ("walk", stats.walk.count,
           "aborted due to inconsistent assumptions");
  }

  for (auto idx : vars)
    if (active (idx))
      set_val (idx, 0);

  assert (level == 2);
  level = 0;

  clear_watches ();
  connect_watches ();

#ifndef QUIET
  if (localsearching) {
    assert (force_phase_messages);
    force_phase_messages = false;
  }
#endif
  stats.ticks.walk += walker.ticks;
  return res;
}

void Internal::walk () {
  START_INNER_WALK ();

  backtrack ();
  if (propagated < trail.size () && !propagate ()) {
    LOG ("empty clause after root level propagation");
    learn_empty_clause ();
    STOP_INNER_WALK ();
    return;
  }

  int res = 0;
  if (opts.warmup)
    res = warmup ();
  if (res) {
    LOG ("stopping walk due to warmup");
    STOP_INNER_WALK ();
    return;
  }
  const int64_t ticks = stats.ticks.search[0] + stats.ticks.search[1];
  int64_t limit = ticks - last.walk.ticks;

  last.walk.ticks = ticks;
  limit *= 1e-3 * opts.walkeffort;
  if (limit < opts.walkmineff)
    limit = opts.walkmineff;
  // local search is very cache friendly, so we actually really go over a
  // lot of ticks
  if (limit > 1e3 * opts.walkmaxeff) {
    MSG ("reached maximum efficiency %" PRId64, limit);
    limit = 1e3 * opts.walkmaxeff;
  }
  VERBOSE (2,
           "walk scheduling: last %" PRId64 " current %" PRId64
           " delta %" PRId64,
           last.walk.ticks, ticks, limit);
  (void) walk_round (limit, false);
  STOP_INNER_WALK ();
}

/*----------------------------------------------------------------------------*/

// passat_build() prepares ...
// (a) passat_lookup_table
// (b) conflict_counter[pos] : number of not-false literals (true + unassigned)
//     of each clause; decreased when a literal inside the clause turns false.
//     conflict_counter == 0 => clause falsified (broken) and, since unassigned
//     literals count as not-false, automatically fully assigned
// (d) walker.activated : number of variables that already carry a value AND are activated, so
//     up_expansion knows when every variable that could been activated is activated (SAT case).
// (e) the ProbSAT score table, built from the average size of the tracked
//     clauses (like walk() does) and used by probSAT_pick_lit.
void Internal::passat_build (Walker &walker) {
  walker.passat_lookup_table.resize (2 * vsize);
  walker.conflict_counter.resize (clauses.size ());
  walker.notfalse_xor.resize (clauses.size (), 0);
  walker.broken_pos.resize (clauses.size (), -1);
  walker.bv.resize(2 * vsize, 0);
  walker.clauses_critical_literal.resize(clauses.size ());
  walker.tuc_pos.resize (clauses.size(), -1);
  walker.unvisitable.resize (clauses.size(), 0);
  walker.unflippable.resize (vsize, 0);
  walker.not_autark.resize(vsize, 0);
  walker.in_flips.resize(vsize, 0);
  walker.pure_lits.resize(vsize, 0);
  walker.satisfied_counter.resize (clauses.size(), 0);
  walker.autarky_val.resize(2 * max_var + 2, 0);

  // broken_occ, lsl and sat_critical_lit only feed the alternative pick scores
  // they only size if they are needed
  if (walker.maintain_pick_stats) {
    walker.broken_occ.resize(2 * vsize, 0);
    walker.lsl.resize(2 * vsize, 0);
    walker.sat_critical_lit.resize(clauses.size ());
  }


  // accumulate total literals and clause count over the tracked clauses to
  // derive the average clause size for the ProbSAT score table 
  double total_size = 0;
  int64_t counted = 0;

  for (size_t pos = 0; pos < clauses.size (); pos++) {
    Clause *c = clauses[pos];
    if (c->garbage)
      continue;
    if (c->redundant) {
      if (!opts.walkredundant)
        continue;
      if (!likely_to_be_kept_clause (c))
        continue;
    }
    // count this tracked clause for the average size (part (e))
    total_size += c->size;
    counted++;

    // because passat_build is cald first in walk_passat and 
    // we only work with complete unassigned clauses
    // we are able to set the conflict_counter for every clause equal to their size
    walker.conflict_counter[pos] = c->size;

    for (const auto lit : *c) {
      // No tracked clause can hold an assigned literal here:
      // walk_passat backtrack to level 0 and propagated to fixpoint,
      // => every assigned literal is root fixed 
      // => the garbage collection (runs before passat_build) drop the clauses containing 
      // the true polarity and delete the false polarity out of the remaining clauses
      // (walk_round use the same logic)
      // => no tracked clause should hold an assigned literal anymore
      assert (!val (lit));

      // xor every literal that starts not false.
      // Once the counter drops to 1 the xor is the critical literal
      // and no clause scan is needed to find it.
      walker.notfalse_xor[pos] ^= lit;

      // (a) if a clause contain the variable v, the clause-position in clauses is inserted
      // in the correct polarity of v in passat_lookup_table => passat_lookup_table[v] += [clause_position]
      walker.passat_lookup_table[vlit (lit)].push_back ((int) pos);
    }

    // Nothing is assigned yet, so no clause is satisfied
    // satisfied_counter, bv, clauses_critical_literal, lsl and sat_critical_lit
    // therefore should all be 0 (until the first passat_assign)
    assert (!walker.satisfied_counter[pos]);
  }

  // part (d)
  // an active variable is always unassigned here and an assigned variable is
  // always root fixed and therefore on the trail
  walker.pre_assigned = trail.size();
  // the pre-assigned variables are not activated by PASSAT, but they do count as activated:
  // the expansion loops stop at pre_assigned + activatable
  walker.activated = trail.size();
  walker.activatable = active();

#ifndef NDEBUG
  {
    size_t assigned = 0, activatable = 0;
    for (int idx = 1; idx <= max_var; idx++) {
      if (val (idx))
        assigned++;
      else if (active (idx))
        activatable++;
    }
    assert (assigned == walker.pre_assigned);
    assert (activatable == walker.activatable);
  }
#endif

  // build the ProbSAT score table from the average clause size, like walk()
  // average clause size is used in walkpassat=15
  const double avg_clause_size = relative (total_size, counted);
  walker.avg_clause_size = avg_clause_size; 
  walker.populate_table(avg_clause_size);
}

/*----------------------------------------------------------------------------*/

// passat_assign() assigns the literal lit to true and keeps all PASSAT
// bookkeeping consistent.
// passat_assign() should not be called on a inactive variable !
// It performs the following steps:
//   1. set lit to true (does nothing if it is already assigned)
//   2. push lit onto the propagation_queue so passat_up can propagate it
//   3. decrement the conflict_counter of every clause that contains -lit
//      (that literal just turned false); if it hits 0 the clause is falsified
//      (and thereby fully assigned), so it is appended to broken_clauses and
//      we report a conflict
//   4. increase walker.activated
//   5. autarky_mode only: keep satisfied_counter and the tuc structures consistent
// Returns true if no conflict arose (propagation may continue), false if
// assigning lit falsified at least one clause.
bool Internal::passat_assign(Walker &walker, int lit) {
  bool signal = true;
  // Only assign if 'lit' is currently unassigned; otherwise there is nothing
  // to do and we report that propagation may continue.
  // NOTE: the guard returns true for any val(lit) != 0. If lit were already
  // false (val(lit) == -1) this would silently hide a contradiction. In our
  // design that should never happens: passat_up only calls passat_assign on
  // unassigned unit literals, and a real conflict is caught earlier via
  // conflict_counter == 0 (just for debugging).
  if (val(lit) == 0){
    assert (active (lit));
    set_val(lit, 1);
    // record on the passat_trail so cleanup can reset exactly this assignment
    walker.passat_trail.push_back(vidx(lit));
    // (4) count this activation so up_expansion knows when all variables are assigned
    walker.activated++;
    // (2) enqueue for later propagation in passat_up
    walker.propagation_queue.push_back(lit);

    // (5) build the autarky bookkeeping
    if (walker.autarky_mode) {
      const auto &pos_clauses = walker.passat_lookup_table[vlit(lit)];
      // we have to increase ticks here because we load a whole line of passat_lookup_table => random Mem Access
      walker.ticks += 1 + cache_lines(pos_clauses.size(), sizeof(int));
      for (auto clause : pos_clauses){
        // we have to increase ticks here because we work on the counters of a clause => random Mem Acces
        walker.ticks++;
        // increase the clause because it is now fullfilled
        walker.satisfied_counter[clause]++;

        // lsl bookkeeping 
        if (walker.maintain_pick_stats) {
          if (walker.satisfied_counter[clause] == 1) {
            // clause´s unique true literal is found
            walker.lsl[vlit(-lit)]++;
            walker.sat_critical_lit[clause] = lit;
          } else if (walker.satisfied_counter[clause] == 2) {
            // clause lost its unique-satisfier status => undo the old critical entry
            const int m = walker.sat_critical_lit[clause];
            assert (m);
            walker.lsl[vlit(-m)]--;
            walker.sat_critical_lit[clause] = 0;
          }
        }

        // the clause just got its first true literal => it is satisfied now and has to leave tuc
        if (walker.satisfied_counter[clause] == 1 && walker.tuc_pos[clause] != -1){
          // swap - remove operation to remove clause from tuc
          int change_position = walker.tuc_pos[clause];
          int last_element = walker.tuc_clauses.back();
          walker.tuc_clauses[change_position] = last_element;
          walker.tuc_pos[last_element] = change_position;
          walker.tuc_pos[clause] = -1;
          walker.tuc_clauses.pop_back();
        }
      }
    }

    // (3) negative occurrences: -lit is now false, so the conflict_counter shrinks
    const auto &neg_clauses = walker.passat_lookup_table[vlit(-lit)];
    // we have to increase ticks here because we load a whole line of passat_lookup_table => random Mem Access
    walker.ticks += 1 + cache_lines(neg_clauses.size(), sizeof(int));
    for(auto clause : neg_clauses){
      // we have to increase ticks here because we work on the counters of a clause => random Mem Acces
      walker.ticks++;
      
      walker.conflict_counter[clause]--;
      // -lit just turned false, so it leaves the set the counter tracks
      walker.notfalse_xor[clause] ^= -lit;
      // a clause whose conflict_counter hit 0 is falsified
      // => it is added broken_clauses here
      // conflict_counter == 0 implies that all lits in the clause are assigned
      if (walker.conflict_counter[clause] == 0) {
        signal = false;
        walker.broken_pos[clause] = (int) walker.broken_clauses.size();
        walker.broken_clauses.push_back(clause);
        if (walker.dynamic_barrier) walker.expansion_conflict_counter++; 

        // bookeeping for the broken_occ list: increasing here because clause is now broken
        if (walker.maintain_pick_stats) {
          Clause *bc = clauses[clause];
          walker.ticks += cache_lines(bc->size, sizeof(int));
          for (auto i : *bc) {
            walker.broken_occ[vlit(i)]++;
          }
        }
        
        walker.ticks++;
        assert (walker.clauses_critical_literal[clause] == -lit);
        walker.bv[vlit(lit)]--;
        walker.clauses_critical_literal[clause] = 0;
      }

      // check if clause is in tuc, if not add it and add the postition in tuc_pos
      // tuc also include broken clauses (we could do the set also exclsuive => Optimization)
      if (walker.autarky_mode && walker.satisfied_counter[clause] == 0 && walker.tuc_pos[clause] == -1){
        walker.tuc_pos[clause] = (int) walker.tuc_clauses.size();
        walker.tuc_clauses.push_back(clause);
      }

      // adjust the break value
      if (walker.conflict_counter[clause] == 1) {
        // exactly one not-false literal is left, so notfalse_xor is equal to the critical literal
        assert (val (walker.notfalse_xor[clause]) >= 0);
        // other is the only non false literal => if it is flipped to fals the clause is broken
        // => bv of -other has to be incremented
        walker.bv[vlit(-walker.notfalse_xor[clause])]++;
        walker.clauses_critical_literal[clause] = walker.notfalse_xor[clause];
      }
    }
  }
  return signal;
}

/*----------------------------------------------------------------------------*/

// Unit propagation for the up_expansion module
// passat_up is working on walker.propagation_queue and with walker.propagated
bool Internal::passat_up(Walker &walker){
  // For every assigned literal still to be processed, look for clauses that
  // just became unit and propagate them through passat_assign
  while (walker.propagated < walker.propagation_queue.size()){
    const int lit = walker.propagation_queue[walker.propagated++];
    // 'lit' is true, so only clauses containing '-lit' can shrink: clauses
    // that contain 'lit' are already satisfied
    // Reading the occurrence row is charged like walk() charges reading a
    // watch list, so that the tick measurement stays comparable to walk()
    const auto &lit_clauses = walker.passat_lookup_table[vlit(-lit)];
    walker.ticks += 1 + cache_lines(lit_clauses.size(), sizeof(int));

    for (auto clause : lit_clauses){
      // one tick per clause we actually visit
      walker.ticks++; 
      // conflict_counter == 1: exactly one literal is still not false.
      // the searched literal is actually already in clauses_critical_literal
      if (walker.conflict_counter[clause] == 1){
        // one tick because clause_critical_literal[clause] may be random
        walker.ticks++;

        const int unit = walker.clauses_critical_literal[clause];

        assert (unit && val (unit) >= 0);

        if (val(unit) == 0) {
          if (!passat_assign(walker, unit)) return false;
        }
      }
    }
  }
  return true;
}

/*----------------------------------------------------------------------------*/

// In the PASSAT-Paper "the UP-guided Expansion module is responsible for
// enlarging the active variable set and construction the next subproblem."
// Therefore we use up_expansion to propagate as far as possible till we found 
// a new subproblem which can be sent to probSAT_repair.
// When UP cannot propagate further, we pick the next unnassigned variable to activate with
// CaDiCaL's own decision heuristic.
// Either next_decision_variable_with_best_score() or next_decision_variable_on_queue()
// is used for the variable picking and decide_phase() for the polarity.
// We use propagation_queue as replacement for the trail, because its easier to manipulate
// and dont break things when working on a local clone of the trail.
bool Internal::up_expansion(Walker &walker) {
  stats.walk.passatexpansion++;

  // First execute anything pending in propagation_queue
  // e.g. flips that probSAT_repair re-enqueued
  if (!passat_up(walker)) return false;

  // Loop until every variable is activated which could be activated.
  // Like advanced_expansion we only aim for pre_assigned + activatable:
  // counting up to max_var would include inactive (eliminated/substituted)
  // variables, which passat_assign must never see (assert (active (lit))).
  while (walker.activated < walker.pre_assigned + walker.activatable) {

    // check if we run out of ticks
    if (walker.ticks >= walker.limit) 
      return false;

    // pick a next unassigned variable to assign
    // because no propagation is left on the propagation_queue
    const int idx = use_scores () ? next_decision_variable_with_best_score ()
                                  : next_decision_variable_on_queue ();

    // next_decision_variable could pick an inactive variable
    // => we have to make shure we are only picking active variables
    if (!active (idx)) {
      set_val (idx, 1);
      walker.passat_trail.push_back (idx);
      continue;
    }

    const bool target = (stable || opts.target == 2);
    // chose the polarity for the choosen variable idx
    const int lit = decide_phase (idx, target);

    // activate the literal. On conflict hand over to probSAT_repair
    if (!passat_assign(walker, lit)) return false;
    // propagate the consequences. On conflict hand over to probSAT_repair
    if (!passat_up(walker)) return false;
  }
  // all activatable variables assigned, no conflict => SAT: no clause may be broken
  assert (walker.broken_clauses.empty ());
  return true;
}

/*----------------------------------------------------------------------------*/

// Idea: Propagate till the the propagation queue is empty
bool Internal::advanced_propagation(Walker &walker){
  bool no_conflict = true;

  while (walker.propagated < walker.propagation_queue.size()){

    const int lit = walker.propagation_queue[walker.propagated++];
    
    const auto &lit_clauses = walker.passat_lookup_table[vlit(-lit)];
    walker.ticks += 1 + cache_lines(lit_clauses.size(), sizeof(int));

    for (auto clause : lit_clauses){
      // one tick per clause we actually visit
      walker.ticks++; 

      // conflict_counter == 1: exactly one literal is still not false.
      // the searched literal is actually already in clauses_critical_literal
      if (walker.conflict_counter[clause] == 1){
        // one tick because clause_critical_literal[clause] may be random
        walker.ticks++;

        const int unit = walker.clauses_critical_literal[clause];

        assert (unit && val (unit) >= 0);

        if (val(unit) == 0) {
          if (!passat_assign(walker, unit)) no_conflict = false;
        }
      }
    }
  }
  return no_conflict;
}

/*----------------------------------------------------------------------------*/

// Idea: use advanced_propagation to force a fully assignment.
// probSAT_repair works quite good, but on a to small subset if I use up_expansion
// because up_expansion just assign till the first conflict arise.
// advanced_expansion assign all possible propagation on the propagation_queue 
// even if there is a conflict
bool Internal::advanced_expansion(Walker &walker) {
  stats.walk.passatexpansion++;

  bool no_conflict = true;

  // remember if probSAT_repair passed a solution with conflicts
  if (walker.anti_stagnation){
    no_conflict = walker.resolved_conflicts;
  }

  if (!passat_up(walker)) no_conflict = false;

  size_t start_activated = walker.activated;

  // Loop until every variable is activated which could be activated
  while (walker.activated < walker.pre_assigned + walker.activatable) {

    if (walker.ticks >= walker.limit)
      return false;

    // pick a next unassigned variable to assign
    // because no propagation is left on the propagation_queue
    const int idx = use_scores () ? next_decision_variable_with_best_score ()
                                  : next_decision_variable_on_queue ();

    // next_decision_variable could pick an inactive variable
    // => we have to make shure we are only picking active variables
    if (!active (idx)) {
      set_val (idx, 1);
      walker.passat_trail.push_back (idx);
      continue;
    }

    const bool target = (stable || opts.target == 2);
    // chose the polarity for the choosen variable idx
    const int lit = decide_phase (idx, target);

    // activate the literal. On conflict hand over to probSAT_repair
    if (!passat_assign(walker, lit)) no_conflict = false;
    // propagate till propagation_queue is empty, then hand over to probSAT_repair
    if (!advanced_propagation(walker)) no_conflict = false;

    // barrier active (!= 0) and reached (>= barrier newly activated vars)?
    if (walker.passat_expansion_barrier &&
        walker.passat_expansion_barrier <= walker.activated - start_activated) {
      // a conflict occurred in this run => hand over to probSAT_repair
      if (!no_conflict)
        return false;
      // nothing is broken, so there is nothing to repair. Reset the counter and start expanding again (no repair here)
      start_activated = walker.activated;
    }
  }
  // all activatable variables assigned => the subproblem is fully expanded
  // without a conflict no clause may be broken 
  assert (!no_conflict || walker.broken_clauses.empty ());
  return no_conflict;
}

/*----------------------------------------------------------------------------*/

// return a random clause-position (identifier) i from a list of claus-postions
// We get the clause than by clauses[i]
// Therefore we dont need the overhead from walk_pick_clause
// (I try to do it generic, so it could be used for other purposes than a random
// clause of broken)
int Internal::pick_random_clause(Walker &walker, const vector<int> &list_of_clauses) {
  assert (!list_of_clauses.empty ());
  int64_t size = list_of_clauses.size ();
  if (size > INT_MAX)
    size = INT_MAX;
  int picked = list_of_clauses[walker.random.pick_int (0, size - 1)];

  // an autarky clause should not be picked
  assert(!walker.unvisitable[picked]);

  return picked;
}

/*----------------------------------------------------------------------------*/

int Internal::advanced_picking(Walker &walker) {

  int picked_lits = 0;

  int64_t size = walker.broken_clauses.size ();
  if (size > INT_MAX) size = INT_MAX;
  int random_clause = 0;
  double sum = 0;

  // duplicated picks are possible, if a lit occurs in multiple broken clauses,
  // and we keep them on purpose: a literal in k drawn clauses gets k-fold roulette mass
  // => coverage bias toward literals that repair many broken clauses.
  // Cap the number of draws so we cannot spin forever if the drawn clauses only
  // hold unflippable/assumed literals; on an empty result the caller falls back
  // to a normal broken-clause pick.
  const int max_draws = 100;
  int draws = 0;
  while (picked_lits < 10 && draws < max_draws) {
    ++draws;
    random_clause = walker.broken_clauses[walker.random.pick_int (0, size - 1)];
    stats.walk.passatadvclauses++;
    Clause *c = clauses[random_clause];
    
    walker.ticks++;
    for(auto lit : *c){
      if (assumed(lit) || assumed(-lit) || !active(lit) || val(lit) == 0 || walker.unflippable[vidx(lit)])
      continue;

      unsigned bv = walker.bv[vlit(lit)];
      walker.ticks++;

      const double s = walker.score(bv);
      walker.scores_passat.push_back({s, lit});
      sum += s;
      picked_lits++;
    }
  }


  // every literal is assumed => nothing to flip, signal "not repairable"
  if (walker.scores_passat.empty()){
    return 0;
  }

  // roulette wheel: pick a random limit in [0, sum).
  const double limit = sum * walker.random.generate_double();

  // Phase 2: pick a random lit
  double current_value = 0;
  int final_pick = 0;
  for (const auto &pick : walker.scores_passat){
    final_pick = pick.second;
    current_value += pick.first;
    if (current_value > limit){
      break; 
    }
  }

  walker.scores_passat.clear();
  return final_pick;
}

/*----------------------------------------------------------------------------*/

// break-value of 'lit' for ProbSAT: the number of clauses that would become
// broken if 'lit' were flipped to true. Flipping 'lit' turns '-lit' false, so
// every clause that currently relies on '-lit' as its last satisfier
// (conflict_counter == 1) would break.
unsigned Internal::passat_break_value(Walker &walker, int lit){
  unsigned res = 0;
  const auto &row = walker.passat_lookup_table[vlit(-lit)];
  // reading the occurrence row is charged like walk() charges a watch list
  walker.ticks += 1 + cache_lines(row.size(), sizeof(int));
  for (int c : row){
    // one tick per clause we actually visit
    walker.ticks++;
    if (walker.conflict_counter[c] == 1)
      res++;
  }
  return res;
}

/*----------------------------------------------------------------------------*/
// exact broken-clause occurrence of 'lit' for debug validation of broken_occ:
// how many currently broken clauses (conflict_counter == 0) contain lit.
// tick-neutral (measurement only, wrapped in asserts).
unsigned Internal::passat_broken_occurence(Walker &walker, int lit){
  unsigned res = 0;
  for (int c : walker.passat_lookup_table[vlit(lit)])
    if (walker.conflict_counter[c] == 0)
      res++;
  return res;
}

/*----------------------------------------------------------------------------*/
// exact last-satisfied-literal value of 'lit' for debug validation of lsl:
// how many clauses have -lit as their unique true literal (satisfied_counter == 1),
// i.e. how many clauses would become touched-unsatisfied if lit were flipped true.
// tick-neutral (measurement only, wrapped in asserts).
unsigned Internal::passat_lsl_value(Walker &walker, int lit){
  // -lit can only be someone's sole satisfier if it is actually true
  if (val (-lit) <= 0) return 0;
  unsigned res = 0;
  for (int c : walker.passat_lookup_table[vlit(-lit)])
    if (walker.satisfied_counter[c] == 1)
      res++;
  return res;
}

/*----------------------------------------------------------------------------*/
// helper function which dont go through all clauses from vlit[-lit] to calc the break value,
// just count how many clauses contain vlit[-lit]
unsigned Internal::passat_fixed_occurence(Walker &walker, int lit){
  walker.ticks++;
  return (unsigned) walker.passat_lookup_table[vlit(-lit)].size();
}

/*-----------------------------------------------------------------------------*/

// pick via ProbSAT a random literal of the earlier picked clause
// The picked clause is broken (all its literals are false). Every flippable
// (non-assumed) literal is scored from its break-value (smaller break-value =>
// higher score) and one literal is sampled by a roulette wheel (like walk() does)
// Important: Returns 0 if every literal is an assumption: then the clause
// cannot be repaired without violating an assumption.
// In reality this conflict is catched earlier, but for safety its better to keep it
int Internal::probSAT_pick_lit(Walker &walker, int picked_clause){
  Clause *c = clauses[picked_clause];
  assert (walker.scores_passat.empty ());
  // one tick for entering the pick, like walk_pick_lit charges the framework
  walker.ticks++;
  double sum = 0;

  // start with a new broken clause => new id and row for break_value track file
  if (walker.track_break_value)
    ++break_value_pick;

  // Phase 1: score every flippable literal by its break-value.
  for (const auto lit : *c){
    // an assumed or inactive (fixed, eliminated or substituted) variable is not allowed to be flipped
    // otherwise flipping would mess the root level vals[] and those variables we cant cleanup here
    // a literal fixed should not be flipped => either autarky or pure lit
    if (assumed(lit) || assumed(-lit) || !active(lit) || val(lit) == 0 || walker.unflippable[vidx(lit)])
      continue;

    // measure the (real, cheap) break-value pair for this literal
    if (walker.track_break_value) {
      const int64_t saved_ticks = walker.ticks;
      const unsigned real_bv  = walker.bv[vlit(lit)];
      const unsigned cheap_bv = passat_fixed_occurence (walker, lit);
      walker.ticks = saved_ticks;
      write_log_file (walker, nullptr, picked_clause, lit, real_bv, cheap_bv);
    }

#ifndef NDEBUG
    // incrementally maintained break value must equal the exact break value
    {
      const int64_t saved_ticks = walker.ticks;
      assert (passat_break_value (walker, lit) == (unsigned) walker.bv[vlit (lit)]);
      
      if (walker.maintain_pick_stats) {
        assert (passat_broken_occurence (walker, lit) == (unsigned) walker.broken_occ[vlit(lit)]);
        assert (passat_lsl_value (walker, lit) == (unsigned) walker.lsl[vlit(lit)]);
      }
      walker.ticks = saved_ticks;
    }
#endif

    // different scoring values for probSAT:
    // 1. cheap_break_value[lit]
    // 2. lsl[lit]
    // 3. broken_occ[lit]
    // 4. broken_occ[lit] * base^bv[lit]
    double s;
    if (walker.cheap_break_value) {
      s = walker.score (passat_fixed_occurence (walker, lit));
    } 
    else if (walker.passat_score_mode == 1) {
      walker.ticks++;
      s = walker.score (walker.lsl[vlit(lit)]);
    } 
    else if (walker.passat_score_mode == 2) {
      walker.ticks++;
      s = walker.broken_occ[vlit(lit)];
    } 
    else if (walker.passat_score_mode == 3) {
      walker.ticks++;
      s = walker.broken_occ[vlit(lit)] * walker.score (walker.bv[vlit(lit)]);
    } 
    else {
      walker.ticks++;
      s = walker.score (walker.bv[vlit(lit)]);
    }

    walker.scores_passat.push_back({s, lit});
    sum += s;
  }

  // every literal is assumed => nothing to flip, signal "not repairable"
  if (walker.scores_passat.empty()){
    return 0;
  }

  // roulette wheel: pick a random limit in [0, sum).
  const double limit = sum * walker.random.generate_double();

  // Phase 2: pick a random lit
  double current_value = 0;
  int res = 0;
  for (const auto &pick : walker.scores_passat){
    res = pick.second;
    current_value += pick.first;
    if (current_value > limit){
      break; 
    }
  }

  walker.scores_passat.clear();
  return res;
}

/*----------------------------------------------------------------------------*/

// function flipps and repair based on the given literal:
// 1. flip the literal
// 2. increase the conflict_counter in all clauses containing lit.
//    clause whose counter rises from 0 to 1 was broken and is removed from broken_clauses
//    also adjust the clauses_critical_literal
// 3. decrease the conflict_counter in all clauses containing -lit
//    clause whose counter drops to 0 becomes broken and is appended to broken_clauses
// 4. add the flipped variable to walker.flips
void Internal::flip_and_repair(Walker &walker, int lit){
  // thrashing: was this var already flipped in this walk_passat run?
  const int fidx = vidx(lit);
  if (walker.flip_count[fidx]++ > 0) stats.walk.passatreflips++;
  // 1.
  set_val(lit, 1);

  // 2.
  const auto &row = walker.passat_lookup_table[vlit(lit)];
  walker.ticks += 1 + cache_lines(row.size(), sizeof(int));
  for (int c : row){
    walker.ticks++;

    // first adjust the old critical literal break value => it is not longer the only positive variable
    if (walker.conflict_counter[c] == 1){
      walker.ticks++;
      // the entry that was incremented for critical literal m is bv[vlit(-m)]
      const int m = walker.clauses_critical_literal[c];
      assert (m);
      walker.bv[vlit(-m)]--;
      walker.clauses_critical_literal[c] = 0;
    }

    walker.conflict_counter[c]++;
    // lit was false and just turned true, so it re-enters the set of not false literals per clause
    walker.notfalse_xor[c] ^= lit;
    // conflict_counter == 1 means the clause was broken and is now
    // satisfied => remove it from broken_clauses
    if (walker.conflict_counter[c] == 1){
      
      // bookeeping for the broken_occ list: decreasing here because clause is now satisfied
      if (walker.maintain_pick_stats) {
        Clause *bc = clauses[c];
        walker.ticks += cache_lines(bc->size, sizeof(int));
        for (auto i : *bc) {
            walker.broken_occ[vlit(i)]--;
        }
      }

      int change_position = walker.broken_pos[c];
      int last_element = walker.broken_clauses.back();
      // move the last element into the new free slot
      walker.broken_clauses[change_position] = last_element;
      // update broken_pos
      walker.broken_pos[last_element] = change_position;
      walker.broken_pos[c] = -1;
      // cut off duplicate of last element
      walker.broken_clauses.pop_back();

      //adjust the break_value of lit
      walker.ticks++;
      walker.bv[vlit(-lit)]++;
      // the clause is critical now and lit is its only true literal
      walker.clauses_critical_literal[c] = lit;
    }

    // autarky bookkeeping: lit is now a true literal in c, so c is satisfied
    if (walker.autarky_mode) {

      walker.satisfied_counter[c]++;

      // lsl bookkeeping
      if (walker.maintain_pick_stats) {
        if (walker.satisfied_counter[c] == 1) {
          walker.lsl[vlit(-lit)]++;
          walker.sat_critical_lit[c] = lit;
        } else if (walker.satisfied_counter[c] == 2) {
          const int m = walker.sat_critical_lit[c];
          assert (m);
          walker.lsl[vlit(-m)]--;
          walker.sat_critical_lit[c] = 0;
        }
      }

      if (walker.satisfied_counter[c] == 1 && walker.tuc_pos[c] != -1) {
        int change_position = walker.tuc_pos[c];
        int last_element = walker.tuc_clauses.back();
        walker.tuc_clauses[change_position] = last_element;
        walker.tuc_pos[last_element] = change_position;
        walker.tuc_pos[c] = -1;
        walker.tuc_clauses.pop_back();
      }
    }
  }

  // 3.
  const auto &neg_row = walker.passat_lookup_table[vlit(-lit)];
  walker.ticks += 1 + cache_lines(neg_row.size(), sizeof(int));
  for (int c : neg_row){
    walker.ticks++;
    walker.conflict_counter[c]--;
    // -lit was true before the flip and is now false => it leaves the not false literals of a clause
    walker.notfalse_xor[c] ^= -lit;

    // autarky bookkeeping: -lit was true before the flip and is now false, so c
    // loses a true literal
    if (walker.autarky_mode) {
      walker.satisfied_counter[c]--;

      // lsl bookkeeping: -lit was true before the flip and is now false
      if (walker.maintain_pick_stats) {
        if (walker.satisfied_counter[c] == 0) {
          // clause lost its last satisfier => undo the critical entry
          const int m = walker.sat_critical_lit[c];
          assert (m == -lit);
          walker.lsl[vlit(-m)]--;
          walker.sat_critical_lit[c] = 0;
        } else if (walker.satisfied_counter[c] == 1) {
          // clause dropped back to a unique satisfier
          Clause *cl = clauses[c];
          walker.ticks += cache_lines(cl->size, sizeof(int));
          for (auto t : *cl) {
            if (val (t) > 0) {
              walker.lsl[vlit(-t)]++;
              walker.sat_critical_lit[c] = t;
              break;
              }
          }
        }
      }

      if (walker.satisfied_counter[c] == 0 && walker.tuc_pos[c] == -1) {
        walker.tuc_pos[c] = (int) walker.tuc_clauses.size();
        walker.tuc_clauses.push_back(c);
      }
    }

    // conflict_counter == 0 => clause is now broken, append it
    if (walker.conflict_counter[c] == 0) {
      walker.broken_pos[c] = (int) walker.broken_clauses.size();
      walker.broken_clauses.push_back(c);

      // the critical literal was -lit, which just turned false, 
      // so the entry bv[vlit(lit)] that counted this clause has to go down
      walker.ticks++;
      assert (walker.clauses_critical_literal[c] == -lit);
      walker.bv[vlit(lit)]--;
      walker.clauses_critical_literal[c] = 0;

      // bookeeping for the broken_occ list: increasing here because clause is now broken
      if (walker.maintain_pick_stats) {
        Clause *bc = clauses[c];
        walker.ticks += cache_lines(bc->size, sizeof(int));
        for (auto i : *bc) {
          walker.broken_occ[vlit(i)]++;
        }
      }
    }

    // adjust the break value
    if (walker.conflict_counter[c] == 1) {
      // exactly one not-false literal is left, so notfalse_xor is equal to the critical literal
      assert (val (walker.notfalse_xor[c]) >= 0);
      walker.bv[vlit(-walker.notfalse_xor[c])]++;
      walker.clauses_critical_literal[c] = walker.notfalse_xor[c];
    }
  }

  // 4. Only the first flip of a variable is pushed and then marked
  if (!walker.in_flips[fidx]) {
    walker.in_flips[fidx] = 1;
    walker.flips.push_back(lit);
  }
  stats.walk.passatflips++;
  // accumulate the number of broken clauses still present after this flip,
  // analogous to walk()'s stats.walk.broken (printed as "per flip")
  stats.walk.passatbroken += walker.broken_clauses.size();

}

/*----------------------------------------------------------------------------*/

// Rebuilds the propagation_queue after a probSAT_repair run so that up_expansion
// can resume propagation
void Internal::repair_propagation_queue(Walker &walker){
  // collect pending variables, mark is empty here, so no check is needed
  // we dont need to think about the correct polarity insertion of v in cache_queue
  // because vals[] already contain the correct polarity after probSAT_repair
  // therefore cache_queue is up to date here
  walker.ticks += 1 + cache_lines(walker.propagation_queue.size() - walker.propagated, sizeof(int));
  for (size_t k = walker.propagated; k < walker.propagation_queue.size(); k++){
    walker.ticks++;
    const int v = vidx(walker.propagation_queue[k]);
    walker.mark[v] = 1;
    walker.cache_queue.push_back(val(v) > 0 ? v : -v);
  }
  // collect flipped variables, mark is not empty here, so we need to check if
  // variable is already in the cache_queue
  walker.ticks += 1 + cache_lines(walker.flips.size(), sizeof(int));
  for (const int lit : walker.flips){
    walker.ticks++;
    const int v = vidx(lit);
    if (!walker.mark[v]){
      walker.mark[v] = 1;
      walker.cache_queue.push_back(val(v) > 0 ? v : -v);
    }
  }
  // reset mark: We consider only the instances of chache_queue
  // which should be cheaper than reset each of the |max_var| - positions by 0
  // not charged: cache_queue and every mark[] entry touched here were just
  // written by the two loops above, so they are still in cache
  for (const int lit : walker.cache_queue)
    walker.mark[vidx(lit)] = 0;

  // make cache_queue the new queue in O(1) and discard the old one
  std::swap(walker.propagation_queue, walker.cache_queue);
  walker.cache_queue.clear();
  walker.propagated = 0;
}

/*----------------------------------------------------------------------------*/

// Unified measurement-log writer for walk_passat
// Depending on which tracking flag is set it appends to the corresponding file
// The log files are opened in Internal, so they are opened once and stay open across all
// walk_passat runs of this solver run
void Internal::write_log_file (Walker &walker, const char *label, int picked_clause, int lit, unsigned real_bv, unsigned cheap_bv) {
  if (walker.track_probSAT_repair && label) {
    // open once; never truncated again for the rest of the solver run
    if (!measure_file) {
      measure_file = fopen ("local_search_modul_measure.log", "w");
      if (!measure_file)
        return;
    }
    FILE *f = measure_file;
    // bump the round on each Start so the matching Start/End share one index
    if (!strcmp (label, "Start Repair"))
      ++measure_round;
    fprintf (f, "[passat] Round %zu - %s:\n", measure_round, label);
    fprintf (f, "[passat] |conflicts| = %zu\n", walker.broken_clauses.size ());
    // number of currently activated variables; repair only flips, never activates,
    // so this must be identical before and after probSAT_repair
    fprintf (f, "[passat] |variables| = %zu\n", walker.activated);
    // collect the variables PASSAT activated (recorded on passat_trail), in
    // index order, signed by their current assignment in vals[]. Root-fixed
    // variables are not listed; repair never flips those anyway, so the
    // Start/End diff below is unaffected.
    std::vector<char> seen (max_var + 1, 0);
    for (const int idx : walker.passat_trail)
      seen[vidx (idx)] = 1;
    std::vector<int> assignment;
    for (int idx = 1; idx <= max_var; idx++)
      if (seen[idx])
        assignment.push_back (val (idx) > 0 ? idx : -idx);

    fputs ("[passat] variables assignment = [", f);
    for (size_t i = 0; i < assignment.size (); i++)
      fprintf (f, "%s%d", i ? ", " : "", assignment[i]);
    fputs ("]\n", f);

    const bool is_start = !strcmp (label, "Start Repair");
    if (is_start) {
      // snapshot the assignment so the matching "End Repair" can diff against it
      walker.measure_start_assignment = assignment;
    } else {
      // net (final) flips of this round: same variable set in the same order
      // (repair never activates new vars), so compare element-wise and report the
      // entries whose sign actually changed between Start and End.
      const std::vector<int> &before = walker.measure_start_assignment;
      if (before.size () == assignment.size ()) {
        // number of variables whose net assignment changed over this round
        size_t final_flips = 0;
        for (size_t i = 0; i < assignment.size (); i++)
          if (before[i] != assignment[i])
            ++final_flips;
        fprintf (f, "[passat] Flips: %zu\n", final_flips);
        fputs ("[passat] From: [", f);
        bool first = true;
        for (size_t i = 0; i < assignment.size (); i++)
          if (before[i] != assignment[i])
            fprintf (f, "%s%d", first ? (first = false, "") : ", ", before[i]);
        fputs ("]\n[passat] To: [", f);
        first = true;
        for (size_t i = 0; i < assignment.size (); i++)
          if (before[i] != assignment[i])
            fprintf (f, "%s%d", first ? (first = false, "") : ", ", assignment[i]);
        fputs ("]\n", f);
      }
    }
    // flush so the file can be read live while the solver is still running
    fflush (f);
  }

  if (walker.track_break_value && !label) {
    // open once and write the header on creation; persists across walk_passat calls
    if (!break_value_file) {
      break_value_file = fopen ("break_value_measure.csv", "w");
      if (!break_value_file)
        return;
      fprintf (break_value_file,
               "pick_id,clause_idx,lit,real_break_value,cheap_break_value\n");
    }
    // One picked broken clause = one pick_id, its flippable literals are the rows
    fprintf (break_value_file, "%zu,%d,%d,%u,%u\n", break_value_pick,
             picked_clause, lit, real_bv, cheap_bv);
    fflush (break_value_file);
  }
}

/*----------------------------------------------------------------------------*/

// write a log to show the autarkies and pure literals
void Internal::write_autarky_log(Walker &walker, bool pure) {
  if (!autarky_file)
    autarky_file = fopen ("autarky.log", "w");
  if (!autarky_file)
    return;
  FILE *f = autarky_file;

  if (pure)
    fprintf (f, "%zu. Pure literals (pre-loop):\n", ++pure_log_count);
  else
    fprintf (f, "%zu. Autarky:\n", ++autarky_log_count);

  // which trail variables belong to this block
  // both predicates are vidx-based, so they take the trail entry directly
  auto selected = [&] (int idx) {
    return pure ? (bool) walker.pure_lits[vidx (idx)]
                : (!walker.not_autark[vidx (idx)] && !walker.pure_lits[vidx (idx)]);
  };

  // the literal a trail variable currently stands for: passat_trail only records
  // which variables PASSAT assigned, the polarity lives in vals[] because
  // flip_and_repair changes it without touching the trail
  auto trail_lit = [&] (int idx) { return val (idx) > 0 ? idx : -idx; };

  // V: the literals of this block (pure lits, or the peeled autarky literals)
  fprintf (f, "%s = {", pure ? "V(pure)" : "V(autarky)");
  bool first = true;
  for (const auto idx : walker.passat_trail)
    if (selected (idx))
      fprintf (f, "%s%d", first ? (first = false, "") : ", ", trail_lit (idx));
  fputs ("}\n", f);

  // F: clauses satisfied by V = union of occ(lit) over V, de-duplicated
  std::vector<int> fset;
  for (const auto idx : walker.passat_trail)
    if (selected (idx))
      for (const auto c : walker.passat_lookup_table[vlit (trail_lit (idx))])
        fset.push_back (c);
  std::sort (fset.begin (), fset.end ());
  fset.erase (std::unique (fset.begin (), fset.end ()), fset.end ());
  fprintf (f, "%s = {", pure ? "F(pure)" : "F(autarky)");
  first = true;
  for (const auto c : fset)
    fprintf (f, "%sc%d", first ? (first = false, "") : ", ", c);
  fputs ("}\n\n", f);

  // per-literal table: for each autarky literal, one row for the literal itself
  // occ(lit): the clauses it satisfies and directly below one row for the
  // opposite polarity occ(-lit): clauses it touches negatively, which the
  // autarky property guarantees are satisfied by some other autarky literal
  // Not meaningful for pure literals so skipped there
  if (!pure) {
    fputs ("lit\t| clauses\n", f);
    fputs ("--------------------------------------------------\n", f);
    for (const auto idx : walker.passat_trail)
      if (selected (idx)) {
        const int lit = trail_lit (idx);
        fprintf (f, "%d\t| ", lit);
        bool fc = true;
        for (const auto c : walker.passat_lookup_table[vlit (lit)])
          fprintf (f, "%s%d", fc ? (fc = false, "") : ", ", c);
        fputs ("\n", f);
        // opposite polarity -lit directly below
        fprintf (f, "%d\t| ", -lit);
        fc = true;
        for (const auto c : walker.passat_lookup_table[vlit (-lit)])
          fprintf (f, "%s%d", fc ? (fc = false, "") : ", ", c);
        fputs ("\n", f);
      }
    fputs ("--------------------------------------------------\n\n", f);
  }
  fflush (f);
}

/*----------------------------------------------------------------------------*/

void Internal::passat_assign_pure_literals(Walker &walker) {
  const int64_t pure_before = stats.walk.passatpureliterals;

  for (int i = 1; i <= max_var; i++) {
    // only care about active and unassigned variables
    if (!active(i) || val(i)) continue;

    // check if the positive polarity of variable i occur in clauses => if so occ_lit == 1
    const bool occ_lit = !walker.passat_lookup_table[vlit(i)].empty();
    walker.ticks++;

    // check if the negative polarity of variable i occur in clauses => if so occ_not_lit == 1
    const bool occ_not_lit = !walker.passat_lookup_table[vlit(-i)].empty();
    walker.ticks++;

    // a literal is pure if the opposite polarity never occurs
    // if both polarities are absent we found a free variable 
    // if both polarities are present, literal cant be pure
    int pure_lit;
    if (occ_lit && !occ_not_lit) {
      pure_lit = i;
    }
    else if (occ_not_lit && !occ_lit) {
      pure_lit = -i;
    } else {
      continue;
    }

    // frozen check: if a variable is marked as frozen, it is not allowed to be removed
    if (frozen(pure_lit)) walker.frozen_autarky = true;

    bool check = passat_assign (walker, pure_lit);
    assert(check);
    (void) check;

    // fix the literal (writes are not charged, only memory accesses are)
    walker.unflippable[vidx(pure_lit)] = 1;
    walker.pure_lits[vidx(pure_lit)] = 1;

    // bookeeping if we want to eliminate the founded autarky
    walker.autarky_trail.push_back(pure_lit);
    walker.autarky_val[vlit(pure_lit)]  =  1;
    walker.autarky_val[vlit(-pure_lit)] = -1;

    stats.walk.passatpureliterals++;
    // every clause of pure literal is satisfied and we dont need to visit them in the future
    const auto &row = walker.passat_lookup_table[vlit(pure_lit)];
    walker.ticks += cache_lines (row.size (), sizeof (int));
    stats.walk.passatpureclauses += (int64_t) row.size ();
    if (walker.autarky_mode) {
      for (const auto c : row) {
        walker.unvisitable[c] = 1;
        walker.ticks++;
      }
    }
  }

  if (walker.show_autarky && stats.walk.passatpureliterals > pure_before)
    write_autarky_log (walker, /*pure=*/true);
}

/*----------------------------------------------------------------------------*/

// passat_assign_pure_literals runs once before the main loop and 
// requires true purity: the opposite polarity cant occur at all. 
// this test is weakened: 
// lit can be added to the current autark set A, if every clause containing -lit is already satisfied by A.
// (Kullmann's autarky closure)
//
// Only unassigned variables are candidates: an assigned variable that survived
// the peeling is already in A, and a peeled one has an open clause containing
// its negation, so it fails the test by construction.
//
// The literals found here count as AUTARKY literals, not as pure literals:
// they exist only because A already covers the clauses of their negation,
// without the autarky they would not be addable
void Internal::find_pure_literals(Walker &walker){
  assert (walker.autarky_mode);

  bool changed = true;
  
  while (changed) {
    changed = false;

    for (int i = 1; i <= max_var; i++) {
      // try to find active but not assigned literal
      if (!active(i) || val(i) != 0) continue;

      const auto &pos = walker.passat_lookup_table[vlit(i)];
      const auto &neg = walker.passat_lookup_table[vlit(-i)];
      walker.ticks += 2;

      if (pos.empty() && neg.empty()) continue;

      // all clauses containing +i covered => -i is addable,
      // because for lit = -i the clauses containing -lit are exactly those of +i
      int pure_lit = 0;
      bool addable = true;
      // both scans can break early, so the sequential part is charged by the
      // number of entries actually read
      size_t visited = 0;
      for (int c : pos) {
        visited++;
        walker.ticks++;
        if (!walker.unvisitable[c]) { addable = false; break; }
      }
      walker.ticks += cache_lines (visited, sizeof (int));

      if (addable) {
        pure_lit = -i;
      } else {
        // mirrored test: all clauses containing -i covered => +i is addable
        addable = true;
        visited = 0;
        for (int c : neg) {
          visited++;
          walker.ticks++;
          if (!walker.unvisitable[c]) { addable = false; break; }
        }
        walker.ticks += cache_lines (visited, sizeof (int));
        if (addable) pure_lit = i;
      }

      // both polarities still occur in an open clause
      if (!pure_lit) continue;

      // frozen check: if a variable is marked as frozen, it is not allowed to be removed
      if (frozen(pure_lit)) walker.frozen_autarky = true;

      // cannot fail: every clause containing -pure_lit is satisfied by an
      // unflippable literal of A, so none of them can break
      bool check = passat_assign (walker, pure_lit);
      assert(check);
      (void) check;

      // fix the literal
      walker.unflippable[vidx(pure_lit)] = 1;

      // bookeeping if we want to eliminate the founded autarky
      walker.autarky_trail.push_back(pure_lit);
      walker.autarky_val[vlit(pure_lit)]  =  1;
      walker.autarky_val[vlit(-pure_lit)] = -1;

      stats.walk.passatautarkylits++;
      stats.walk.passatextrapure++;
      changed = true;

      // mark the clauses of the new (pure) autarky lit as unvisitable
      const auto &row = walker.passat_lookup_table[vlit(pure_lit)];
      walker.ticks += cache_lines (row.size (), sizeof (int));
      for (const auto c : row) {
        walker.ticks++;
        if (!walker.unvisitable[c]) {
          walker.unvisitable[c] = 1;
          stats.walk.passatautarkyclauses++;
        }
      }
    }

    if (changed) stats.walk.passatextrapurerounds++;
  }
}

/*----------------------------------------------------------------------------*/

// here we build the autarky_set
void Internal::build_autarky(Walker &walker) {
  assert (walker.autarky_mode);

  // make shure not_autark array is every zero at the start
  std::fill (walker.not_autark.begin (), walker.not_autark.begin () + max_var + 1, 0);
  walker.ticks += 1 + cache_lines ((size_t) max_var + 1, sizeof (signed char));

  walker.ticks += 1 + cache_lines(walker.tuc_clauses.size(), sizeof(int));
  walker.autarky_worklist.assign(walker.tuc_clauses.begin(), walker.tuc_clauses.end());


  // lazy sizing if we need peel_counter and peel_stamp
  if(walker.peel_stamp.empty()){
    walker.peel_counter.resize (clauses.size ());
    walker.peel_stamp.resize (clauses.size (), 0);
  }
  walker.peel_generation++;

  // Phase 1:
  // remove all assigned literals lit that cannot be in an autarky set,
  // because clauses with -lit are not fullfilled
  size_t i = 0;
  while (i < walker.autarky_worklist.size()) {
    walker.ticks++;
    Clause *clause = clauses[walker.autarky_worklist[i]];
    i++;

    walker.ticks += cache_lines(clause->size, sizeof(int));
    for (auto lit : *clause){
      assert(val(lit) != 1 || walker.not_autark[vidx(lit)]);

      if (val(lit) == -1) {
        if (walker.not_autark[vidx(-lit)]) continue;

        walker.not_autark[vidx(-lit)] = 1;

        // -lit is assigned true and the clause we are processing is touched by it
        // => -lit cannot be part of the autarky candidate
        // every clause with -lit loses a satisifier
        // therefore we have to decrease the temporary satisfied_counter (peel_counter) by one
        const auto &row = walker.passat_lookup_table[vlit(-lit)];
        walker.ticks += 1 + cache_lines(row.size(), sizeof(int));
        for (auto c : row){
          walker.ticks++;
          // check if clause is visited before in that call
          // if not we have to assign the correct peel_counter
          if (walker.peel_stamp[c] != walker.peel_generation){
            walker.peel_stamp[c] = walker.peel_generation;
            walker.peel_counter[c] = walker.satisfied_counter[c];
          }
          // c is not visited yet for -lit, so satisfied_counter[c] still should be positiv
          assert(walker.peel_counter[c] > 0);
          // decrement the satisfied_counter because -lit cant be the reason for autarky
          walker.peel_counter[c]--;
          if (walker.peel_counter[c] == 0){
            // we need to take care of the rest of c without -lit
            walker.autarky_worklist.push_back(c);
          }
        }
      }
    }
  }

  // Phase 2:
  // Build the actual autarky set 
  // and update unflippable and unvisitable
  bool found_autarky = false;
  // did the autark set actually grow since the last call of Repair (but in the same walk_passat run)
  bool autarky_grew = false;

  walker.ticks += 1 + cache_lines(walker.passat_trail.size(), sizeof(int));
  for (auto idx : walker.passat_trail){
    walker.ticks++;
    // make shure that inactive (parked) literals and free variables are not in the set of an autarky
    // because they dont cause any clause
    if (!active(idx)) continue;

    assert (val(idx));
    const int lit = val(idx) > 0 ? idx : -idx;
    
    walker.ticks++;
    if (walker.passat_lookup_table[vlit(lit)].empty ()) {
      walker.ticks++;
      if (walker.passat_lookup_table[vlit(-lit)].empty ()) continue;
    }

    // take care of a new founded autarky lit which is not pure.
    // Pure literals are assigned and tracked before the main loop because they are trivially autark.
    // Together with the !unflippable guard below this keeps the invariant that
    // every literal enters autarky_trail exactly once. 
    // Additional pure literals in later loop runs are processed by find_pure_literals().
    if (walker.not_autark[vidx(lit)] != 1 && !walker.pure_lits[vidx(lit)]) {
      found_autarky = true;

      // frozen check: if a variable is marked as frozen, it is not allowed to be removed
      if (frozen(lit)) walker.frozen_autarky = true;

      if (!walker.unflippable[vidx (lit)]) {
        walker.unflippable[vidx (lit)] = 1;
        // bookeeping if we want to eliminate the founded autarky
        walker.autarky_trail.push_back(lit);
        walker.autarky_val[vlit(lit)]  =  1;
        walker.autarky_val[vlit(-lit)] = -1;

        stats.walk.passatautarkylits++;
        autarky_grew = true;

        // take care of all the clause we are not allowed to touch anymore 
        // because of the autarky
        const auto &row = walker.passat_lookup_table[vlit(lit)];
        walker.ticks += cache_lines(row.size(), sizeof(int));
        for (auto c : row){
          walker.ticks++;
          if (!walker.unvisitable[c]) {
            walker.unvisitable[c] = 1;
            stats.walk.passatautarkyclauses++;
          }
        }
      }
      #ifndef NDEBUG
        else {
          // proof of the skip above: everything must already be covered
          for (auto c : walker.passat_lookup_table[vlit (lit)])
            assert (walker.unvisitable[c]);
        }
      #endif
    }
  }

  if (found_autarky) stats.walk.passatautarky++;

  // optional flag for showing the grwing autarky
  if (walker.show_autarky && autarky_grew)
    write_autarky_log (walker);
}

/*----------------------------------------------------------------------------*/

// Returns true if the conflict was fully resolved (broken == 0) so that
// up_expansion can resume; false if it could not be repaired (=> UNSAT).
bool Internal::probSAT_repair(Walker &walker) {
  /*
  1. operate on walker.broken_clauses, which passat_assign and flip_and_repair
     maintain incrementally (a broken clause is always fully assigned, since
     conflict_counter counts not-false = true + unassigned literals)
  2. Flip decision of Literal l (wie auch in walk() vernwendet):
     Wir wählen zufällig eine clause c aus broken und dann mit probsat ein literal l aus c.
  3. Flip Literal l (which is decide in 2.).
  4. Make all adjustment and book keeping stuff.
  5. Look if there is still broken clauses, then return to (2), otherwise SAT and break the loop
  6. If a solution is found, write the solution in vals[] (maybe it is already written), wirte the new values on the propagation_queue (we need to do this correct?)
     and pass to up_expansion
  */

  stats.walk.passatrepair++;

  //Clear all earlier made flips
  // resetting via the list itself costs O(distinct variables), not O(max_var)
  for (const int l : walker.flips)
    walker.in_flips[vidx (l)] = 0;
  walker.flips.clear();

  walker.stagnation_counter = 0;

  // measurement of the input in the  LS step right after expansion
  if (walker.track_probSAT_repair) write_log_file (walker, "Start Repair", 0, 0, 0, 0);

  // convergence: track broken at start and the best (lowest) broken ever reached
  const size_t start_broken = walker.broken_clauses.size();
  // safe the broken size of the last expansion run
  walker.last_start_broken = start_broken;
  size_t min_broken = start_broken;
  walker.last_min_broken = start_broken;

  // tuc min is used in v41 and try to find a bigger autarky whenever |TUC| has halved.
  // (does not work good)
  size_t min_tuc = walker.tuc_clauses.size ();
  walker.tuc_autarky_ticks = 0;

  // the stagnation limit must be determined experimentally
  // Keep in mind that the larger the limit is the longer probSAT_repair runs and can resolve conflicts,
  // but less variables become active
  // We have to find a good balance between number of active variables and correct solved conflcits by probSAT_repair
  // 1.000 was better than 10.000
  const size_t stagnation_limit = 1000 * start_broken;

  // remember the last assignment of the walk if improvement option is activ
  if (walker.passat_track_improvement) {
    walker.best_repair_model.resize (walker.passat_trail.size ());
    for (size_t i = 0; i < walker.passat_trail.size (); i++)
      walker.best_repair_model[i] = val (walker.passat_trail[i]);
  }

  while(!walker.broken_clauses.empty() && walker.ticks < walker.limit){
    // pick random clause, then a literal of it via ProbSAT
    int lit = 0;
    if (walker.advanced_picking_mode) {
      lit = advanced_picking(walker);
    }
    if (lit == 0) {
      lit = probSAT_pick_lit(walker, pick_random_clause(walker, walker.broken_clauses));
    }

    // lit == 0 from a broken clause: it consists only of assumptions, so it cannot
    // be repaired without violating an assumption => unsatisfiable under the
    // assumptions, stop and report failure.
    if (lit == 0) {
      walker.assumption_unsat = true;
      walker.resolved_conflicts = false;
      return false;
    }

    // make the actual LS_repair
    flip_and_repair(walker, lit);

    // check if we found an improvement for the current solution.
    // track the minimum ALWAYS (it is O(1) and feeds the broken-min statistic for every
    // version); only snapshot the assignment when the improvement option (v16) is active.
    if (walker.broken_clauses.size() < min_broken) {
      min_broken = walker.broken_clauses.size();
      walker.last_min_broken = min_broken;

      // save the current best solution (values are in val, active literals are on passat_trail)
      // in best_repair_model
      if (walker.passat_track_improvement) {
        for (size_t i = 0; i < walker.passat_trail.size (); i++) {
          walker.best_repair_model[i] = val(walker.passat_trail[i]);
        }
      }
      
      // reset the stagnation counter because we found a better solution
      walker.stagnation_counter = 0;
    } else {
      // increase the stagnation counter because the last found solution was not better than the min solution
      walker.stagnation_counter++;
    }

    // check if there is a new tuc minimum
    if (walker.tuc_min_autarky_check) {
      const size_t cur_tuc = walker.tuc_clauses.size ();

      if (cur_tuc * 2 <= min_tuc) {
        const int64_t tuc_ticks_before = walker.ticks;
        const int64_t tuc_lits_before = stats.walk.passatautarkylits;
        const int64_t tuc_clauses_before = stats.walk.passatautarkyclauses;
        build_autarky (walker);
        walker.tuc_autarky_ticks += walker.ticks - tuc_ticks_before;

        min_tuc = walker.tuc_clauses.size ();

        stats.walk.passatautarkylitstuc +=
            stats.walk.passatautarkylits - tuc_lits_before;
        stats.walk.passatautarkyclausestuc +=
            stats.walk.passatautarkyclauses - tuc_clauses_before;
        stats.walk.passatautarkychecktuc++;
      }
    }

    if (walker.anti_stagnation && walker.stagnation_counter >= stagnation_limit) {
      // repair is stagnating and only burns ticks without improving
      // => stop earlier and invest the remaining ticks in further expansion
      stats.walk.passatstagnationbreaks++;
      break;
    }
  }

  stats.walk.passatstagnation += walker.stagnation_counter;

  // measure the output of the LS step after probSAT_repair
  if (walker.track_probSAT_repair) write_log_file (walker, "End Repair", 0, 0, 0, 0);

  // record broken at start for every repair (assumption-only exit excluded)
  stats.walk.passatbrokenstart += start_broken;

  // tick limit reached while clauses are still broken: the conflict could not be
  // resolved within the budget
  if (!walker.broken_clauses.empty()) {
    // convergence: for FAILED repairs only, how close did we get to broken==0?
    stats.walk.passatbrokenmin += min_broken;
    walker.resolved_conflicts = false;
    
    if (walker.anti_stagnation && !walker.assumption_unsat)
      repair_propagation_queue(walker);
    return false;
  }

  // broken == 0: the conflict is fully resolved. Rebuild the propagation_queue so
  // up_expansion can resume on the repaired partial assignment, then report success.
  stats.walk.passatrepairsuccess++;
  repair_propagation_queue(walker);

  // the next expansion starts with an empty broken list
  assert (walker.broken_clauses.empty ());

  walker.resolved_conflicts = true;

  return true;
}

/*---------------------------------------------------------------------------*/

// PASSAT-Algorithm
// procedure:
// 1. backtrack() to empty the trail
// 2. propagate() unpropagated literals in the trail
// 3. limit calculation
// 4. walker instantiation
// 5. build all necessary structures for passat
// 6. make assumptions (if an assumption is failed break!)
// 7. PASSAT-Loop: up_expansion and probSAT_repair
// 8. clean up and return the best found phase
void Internal::walk_passat() {
  START_INNER_WALK ();

  // increase the statistic counter for passat
  stats.walk.passat++;

  // wall-clock of this call, accumulated into stats.walk.passatseconds at every
  // exit below; only feeds the flips-per-second line in the statistics
  const double passat_start_time = time ();

    backtrack ();
  
  //propagate() is called if unpropagated literals are still present in the trail after backtrack()
  if (propagated < trail.size () && !propagate ()) {
    LOG ("empty clause after root level propagation");
    learn_empty_clause ();
    stats.walk.passatseconds += time () - passat_start_time;
    STOP_INNER_WALK ();
    return;
  }

  // make shure we do not work with clauses that contain a root level fixed literal
  if (last.collect.fixed < stats.all.fixed) 
    garbage_collection ();

  //calc limit, identically as in walk()
  const int64_t start_ticks = stats.ticks.search[0] + stats.ticks.search[1];
  int64_t limit = start_ticks - last.walk.ticks;
  last.walk.ticks = start_ticks;
  limit *= 1e-3 * opts.walkeffort;
  if (limit < opts.walkmineff) limit = opts.walkmineff;
  if (limit > 1e3 * opts.walkmaxeff) limit = 1e3 * opts.walkmaxeff;

  Walker walker (internal, limit);

  walker.maintain_pick_stats = (opts.walkpassat >= 36 && opts.walkpassat <= 39);

  // build occurrence lists, counters and the activation count for this run
  passat_build (walker);

  /* ---------------------------------------------------------------------------------------------------------------------------------------

  // select the configuration from --walkpassat=n:
  //
  // --------- versions of the bachelorproject: build PASSAT Algorithm into walk ------------------
  // versions 1 to 7 use the exact break value, versions 8 to 14 the cheap break value
  // version 7 and 14 use up_expansion as described in the pap
  // First check if we use the cheap break value, second we decide which barrier size we use
  // version 15 uses the exact break value with a (dynamic) soft-adaptive barrier + imrpovement tracking (like v16)
  // version 16 uses the exact break value, a 0.1% barrier and improvement tracking : keep the better of the
  // post-expansion and the post-repair assignment when writing phases.saved
  // version 17 = version 5 + improvement tracking
  // version 18 = version 5 + 3x tick limit
  // version 19 = version 15 + 3x tick limit
  // version 20 = version 16 + 3x tick limit
  // version 21 = version 17 + 3x tick limit
  // version 22 = version 5 + anti-stagnation
  // version 23 = version 22 + 3x tick limit
  //
  // ------------ versions of the bachelor thesis: find and eliminate autarkies during local search --------------------------------
  // ------------ autarky check with classic PASSAT Algorithm ------------------
  // version 24 = version 7 (classic PASSAT, up_expansion) + autarky check after repair + elimination
  // version 25 = version 24 + advanced pure lits finding
  // version 26 = version 25 + 3x tick limit
  // ------------ autarky check with best version of previous tests => dynamic barrier (does not wrk good here) ------------------
  // version 27 = dynamic barrier + 3x tick limit (=v19) + autarky check after repair + elimination
  // ------------ autarky check with anti stagnation, works best for autarkies because reaches full coverage ------------------
  // version 28 = version 22 (v5 + anti-stagnation) + autarky check (only) after expansion + elimination
  // version 29 = version 22 (v5 + anti-stagnation) + autarky check (only) after repair + elimination
  // version 30 = version 29 + autarky check after expansion as well
  // ------------ version 29 is the reference base for all further optimations and tests ------------------
  // version 31 = version 29 + advanced pure lits finding
  // version 32 = version 29 + 3x tick limit
  // version 33 = version 29 + advanced pure lits finding + 3x tick limit
  // version 34 = version 29 + build autarky at a TUC minimum (50% fewer tuc clauses than the last minimum)
  // version 35 = version 34 + advanced pure lits finding + 3x tick limit
  // version 36 = version 29, pick scored by base^lsl instead of the break value
  // version 37 = version 29, pick scored by the occurrence in broken clauses (broken_occ)
  // version 38 = version 29, pick scored by broken_occ[lit] * base^bv[lit]
  // version 39 = version 38 + advanced pure lits finding + 3x tick limit
  // version 40 = version 29 + advanced_picking (multi-clause pool)
  // version 41 = version 40 + advanced pure lits finding + 3x tick limit

  ---------------------------------------------------------------------------------------------------------------------------------------- */

  if (opts.walkpassat == 24 || opts.walkpassat == 25 || opts.walkpassat == 26) {
    // up_expansion base (classic PASSAT)
    walker.cheap_break_value = false;
    walker.use_up_expansion = true;
    walker.autarky_mode = true;
    walker.increased_passat_limit = (opts.walkpassat == 26);
    walker.autarky_check_expansion = false;
    walker.autarky_check_repair = true;
    walker.autarky_elimination_mode = true;
    walker.advanced_pure_finding = (opts.walkpassat == 25 || opts.walkpassat == 26);
  }
  else if (opts.walkpassat == 27) {
    walker.cheap_break_value = false;
    walker.dynamic_barrier = true;
    walker.passat_expansion_barrier = (walker.avg_clause_size > 3.5)
        ? std::max ((size_t) 1, walker.activatable / 100)   // 1%
        : std::max ((size_t) 1, walker.activatable / 10);    // 10%
    walker.passat_track_improvement = false;
    walker.increased_passat_limit = true;
    walker.autarky_mode = true;
    walker.autarky_check_expansion = false;
    walker.autarky_check_repair = true;
    walker.autarky_elimination_mode = true;
  }
  else if (opts.walkpassat >= 28 && opts.walkpassat <= 41) {
    // anti-stagnation base; 29 is the reference every later version builds on
    walker.cheap_break_value = false;
    walker.passat_expansion_barrier = std::max ((size_t) 1, walker.activatable / 10); // 10% like v5
    walker.anti_stagnation = true;
    walker.autarky_mode = true;
    walker.autarky_elimination_mode = true;

    walker.increased_passat_limit = (opts.walkpassat == 32 || opts.walkpassat == 33 ||
                                     opts.walkpassat == 35 || opts.walkpassat == 39 ||
                                     opts.walkpassat == 41);
    walker.advanced_pure_finding = (opts.walkpassat == 31 || opts.walkpassat == 33 ||
                                    opts.walkpassat == 35 || opts.walkpassat == 39 ||
                                    opts.walkpassat == 41);
    walker.tuc_min_autarky_check = (opts.walkpassat == 34 || opts.walkpassat == 35);
    walker.passat_score_mode = (opts.walkpassat == 36) ? 1   // base^lsl
                             : (opts.walkpassat == 37) ? 2   // broken_occ
                             : (opts.walkpassat == 38 ||
                                opts.walkpassat == 39) ? 3   // broken_occ * base^bv
                             : 0;
    walker.advanced_picking_mode = (opts.walkpassat == 40 || opts.walkpassat == 41);

    // 28 checks only after expansion, 30 checks after both, everything else
    // only after repair.
    walker.autarky_check_expansion = (opts.walkpassat == 28 || opts.walkpassat == 30);
    walker.autarky_check_repair = (opts.walkpassat != 28);
  }
  else if (opts.walkpassat == 15 || opts.walkpassat == 19) {
    // 1. Pick the starting barrier from the average clause length:
    //    avg clause-length > 3.5 => start at 1% (static 1% works better on long clauses), else 10%.
    // 2a. After each advanced_expansion run: if this run had >20% MORE conflicts than the previous
    //     run, drop the barrier straight down to 1% (repair more, expand less).
    // 2b. If this run had >20% FEWER conflicts, raise the barrier to 10% to expand faster.
    walker.cheap_break_value = false;
    walker.dynamic_barrier = true;
    // start from the average clause length: wide clauses (many conflicts) -> 1%, binary-dominated -> 10%
    walker.passat_expansion_barrier = (walker.avg_clause_size > 3.5)
        ? std::max ((size_t) 1, walker.activatable / 100)   // 1%
        : std::max ((size_t) 1, walker.activatable / 10);    // 10%
    walker.passat_track_improvement = true;
    walker.increased_passat_limit = (opts.walkpassat == 19);
  } 
  else if (opts.walkpassat == 16 || opts.walkpassat == 20) {
    walker.cheap_break_value = false;
    walker.passat_expansion_barrier = std::max ((size_t) 1, walker.activatable / 1000); // 0.1%
    walker.passat_track_improvement = true;
    walker.increased_passat_limit = (opts.walkpassat == 20); // v20 = v16 + 3x
  } 
  else if (opts.walkpassat == 17 || opts.walkpassat == 18 || opts.walkpassat == 21) {
    walker.cheap_break_value = false;
    walker.passat_expansion_barrier = std::max ((size_t) 1, walker.activatable / 10); // 10%
    walker.passat_track_improvement = (opts.walkpassat == 17 || opts.walkpassat == 21);
    walker.increased_passat_limit = (opts.walkpassat == 18 || opts.walkpassat == 21);
  }
  else if (opts.walkpassat == 22 || opts.walkpassat == 23) {
    walker.cheap_break_value = false;
    walker.passat_expansion_barrier = std::max ((size_t) 1, walker.activatable / 10); // 10% like v5
    walker.anti_stagnation = true;
    walker.increased_passat_limit = (opts.walkpassat == 23); // v23 = v22 + 3x
  }
  else {
    walker.cheap_break_value = (opts.walkpassat > 7);
    switch (((opts.walkpassat - 1) % 7) + 1) {
    case 1: walker.passat_expansion_barrier = 10; break; // s = 10
    case 2: walker.passat_expansion_barrier = 100; break; // s = 100
    case 3: walker.passat_expansion_barrier = 0; break; // unlimited
    case 4: walker.passat_expansion_barrier = std::max ((size_t) 1, walker.activatable / 100); break; // s = 1% of active variables
    case 5: walker.passat_expansion_barrier = std::max ((size_t) 1, walker.activatable / 10); break; // s = 10% of active variables
    case 6: walker.passat_expansion_barrier = std::max ((size_t) 1, walker.activatable / 2); break; // s = 50% of active variables
    case 7: walker.use_up_expansion = true; break; // use the original up_expansion instead of advanced_expansion
    }
  }

  // if the increasing limit flag is true => 3x larger tick limit
  if (walker.increased_passat_limit) walker.limit *= 3;

  // care about the assumptions
  bool consistent_with_assumptions = true;
  for (int lit : assumptions){
    // get the value of the lit in vals[]
    signed char assumptions_value = val(lit);
    // already correctly assigned assumptions
    if (assumptions_value > 0) {
      continue;
    }
    // wrongly assigned assumption => inconsistent => UNSAT (???)
    else if (assumptions_value < 0){
      consistent_with_assumptions = false;
      break;
    }
    // force an assignment immediately
    else{
      // eliminated/substituted assumption: leave it to reconstruction, don't force
      if (!active(lit)) continue;
      // try to assign the assumption
      if(!passat_assign(walker, lit)){
        consistent_with_assumptions = false;
        break;
      }
    }
  }

  // PASSAT main loop on an (initially) empty assignment:
  // up_expansion activates variables via CaDiCaL's decision heuristic and
  // propagates (UP) until
  // SAT (all activated, no conflict) or a conflict; probSAT_repair then repairs
  // the fully-activated subproblem. Resume only if the conflict was resolved.
  bool no_conflict = false;
  
  const int64_t autarky_lits_at_start = stats.walk.passatautarkylits;
  const int64_t autarky_clauses_at_start = stats.walk.passatautarkyclauses;
  // flips done by this walk_passat call, for the flips-per-second report below
  const int64_t flips_at_start = stats.walk.passatflips;
  if (consistent_with_assumptions){
    no_conflict = true;

    // detect and fix pure literals before the main loop
    int64_t pure_ticks_before = walker.ticks;
    passat_assign_pure_literals (walker);
    stats.walk.passatpureticks += walker.ticks - pure_ticks_before;

    // for dynamic barrier (walkpassat=15): no previous run to compare against the first expansion
    bool first_run = true;
    while (walker.ticks < walker.limit) {
      int64_t ticks_before = walker.ticks;
      // count the conflicts of this advanced_expansion run
      if (walker.dynamic_barrier) walker.expansion_conflict_counter = 0;
      no_conflict = walker.use_up_expansion ? up_expansion(walker)
                                            : advanced_expansion(walker);
      stats.walk.passatexpansionticks += walker.ticks - ticks_before;

      // SAT over the activated set
      if (no_conflict)
        break;

      // check for autarkies before Repair, because maybe we found a autarky
      // if so it would make no sense to flipp on a existing autarky,
      // because the conflict is not inside the autarky
      int64_t autarky_ticks_before = walker.ticks;
      int64_t autarky_lits_before = stats.walk.passatautarkylits;
      int64_t autarky_clauses_before = stats.walk.passatautarkyclauses;

      if (walker.autarky_mode && walker.autarky_check_expansion)
        build_autarky (walker);
      stats.walk.passatautarkyticksexp += walker.ticks - autarky_ticks_before;
      stats.walk.passatautarkylitsexp += stats.walk.passatautarkylits - autarky_lits_before;
      stats.walk.passatautarkyclausesexp += stats.walk.passatautarkyclauses - autarky_clauses_before;

      // update the dynamic barrier: after each run compare this runs conflicts to the previous runs conflicts
      if (walker.dynamic_barrier) {
        const size_t cur = walker.expansion_conflict_counter;
        const size_t prev = walker.last_expansion_conflicts;
        if (!first_run) {
          if (cur > prev * 1.2) {
            // 20% more conflicts than the run before => barrier down to 1% immediately
            const size_t one = std::max ((size_t) 1, walker.activatable / 100);
            if (walker.passat_expansion_barrier > one) stats.walk.passatbarrierdown++;
            walker.passat_expansion_barrier = one;
          } else if (cur < prev * 0.8) {
            // 20% fewer conflicts => barrier up to 10% => expand faster
            const size_t ten = std::max ((size_t) 1, walker.activatable / 10);
            if (walker.passat_expansion_barrier < ten) stats.walk.passatbarrierup++;
            walker.passat_expansion_barrier = ten;
          }
        }
        walker.last_expansion_conflicts = cur;
        first_run = false;
      }
      
      ticks_before = walker.ticks;
      const bool repaired = probSAT_repair(walker);
      
      stats.walk.passatrepairticks += walker.ticks - ticks_before - walker.tuc_autarky_ticks;
      stats.walk.passatautarkytickstuc += walker.tuc_autarky_ticks;

      // check for autarkies after Repair
      // we could find autarkies if Repair luckily flipped one
      autarky_ticks_before = walker.ticks;
      autarky_lits_before = stats.walk.passatautarkylits;
      autarky_clauses_before = stats.walk.passatautarkyclauses;

      if (walker.autarky_mode && walker.autarky_check_repair) 
        build_autarky (walker);

      stats.walk.passatautarkyticksrep += walker.ticks - autarky_ticks_before;
      stats.walk.passatautarkylitsrep += stats.walk.passatautarkylits - autarky_lits_before;
      stats.walk.passatautarkyclausesrep += stats.walk.passatautarkyclauses - autarky_clauses_before;

      // try to find new pure literals
      if (walker.autarky_mode && walker.advanced_pure_finding) {
        const int64_t extra_pure_ticks_before = walker.ticks;
        find_pure_literals (walker);
        const int64_t spent = walker.ticks - extra_pure_ticks_before;
        stats.walk.passatpureticks += spent;
        stats.walk.passatextrapureticks += spent;
      }
      
      // if improvement version (walkpassat=16) is used and a reapir failed, 
      // we check if the current assignment is the best assignment the walk passed through
      if (walker.passat_track_improvement && !repaired && !walker.anti_stagnation) {
        // restore only the activated vars
        // everything else was never touched by probSAT
        // iterate over the snapshot, not over passat_trail: the trail is only
        // appended to, but find_pure_literals above can have grown it since
        // probSAT_repair sized best_repair_model, and the extra entries have no
        // snapshot entry to restore from
        assert (walker.best_repair_model.size () <= walker.passat_trail.size ());
        for (size_t i = 0; i < walker.best_repair_model.size (); i++)
          set_val (walker.passat_trail[i], walker.best_repair_model[i]);
        if (walker.last_min_broken >= walker.last_start_broken)
          stats.walk.passatexpkept++;
      }

      // conflict not resolvable -> UNSAT
      // only anti_stagnation keeps going with expansion if tick limit is not exceeded
      if (!repaired && (!walker.anti_stagnation || walker.assumption_unsat))
        break;
    }
  }

  stats.walk.passatactivations += walker.activated - walker.pre_assigned;
  stats.walk.passatactivatable += walker.activatable;

  // count the size of the autarky of the current walk_passat run
  {
    const int64_t run_lits =
        stats.walk.passatautarkylits - autarky_lits_at_start;
    const int64_t run_clauses =
        stats.walk.passatautarkyclauses - autarky_clauses_at_start;
    if (run_lits) {
      stats.walk.passatautarkyruns++;
      if (run_lits > stats.walk.passatautarkylitsmax)
        stats.walk.passatautarkylitsmax = run_lits;
      if (run_clauses > stats.walk.passatautarkyclausesmax)
        stats.walk.passatautarkyclausesmax = run_clauses;
    }
  }

  LOG("walk_passat: %s", no_conflict ? "SAT" : "limit reached");

  // PHASE to show the progress
  PHASE ("walk_passat", stats.walk.passat, "%s after %" PRId64 " ticks",
         no_conflict ? "satisfied activated set" : "limit reached",
         walker.ticks);

#ifndef QUIET
  if (opts.profile >= 2) {
    const double seconds = time () - profiles.walk.started;
    const int64_t run_flips = stats.walk.passatflips - flips_at_start;
    PHASE ("walk_passat", stats.walk.passat, "%.2f million ticks per second",
           1e-6 * relative (walker.ticks, seconds));
    // the interesting number for data-structure work: the tick budget is fixed,
    // so a faster implementation shows up here and nowhere else
    PHASE ("walk_passat", stats.walk.passat, "%.2f million flips per second",
           1e-6 * relative (run_flips, seconds));
  }
#endif

  // Save the result: the whole point of walk_passat is to leave a better
  // polarity assignment in phases.saved, which the following CDCL search uses
  // as decision phases
  for (int id = 1; id <= max_var; id++)
    if (val(id))
      phases.saved[id] = val(id);

  // In walk_passat we do not push on the trail if we assign, 
  // therefore we have to restore values of variables we assigned during walk_passat to 0.
  // Then we have to reset the decision level to the root,
  // otherwise the next CDCL search runs on a corrupted state (like in walk_round). 
  // Fixed vars are never on the trail, so their real root-level vals stay untouched.
  for (const auto idx : walker.passat_trail) {

    if (walker.autarky_mode && walker.autarky_val[vlit(idx)]){
      // check if the correct literal of an autarky variable is assigned to true and 
      // not incorrect flipped
      assert(val(idx) == walker.autarky_val[vlit(idx)]);
    }

    set_val(idx, 0);
    if (!scores.contains(idx)) scores.push_back (idx);
    if (queue.bumped < btab[idx]) update_queue_unassigned (idx);
  }
  
  level = 0;

  stats.walk.passatseconds += time () - passat_start_time;
  STOP_INNER_WALK();

  // we do the same as a call from autarky() in autarky.cpp does
  // 1. build a lookup array for the autarky set
  // 1. claer all watches
  // 2. call autarky_apply
  // 3. eiminate all literals
  // 4. reconnect all watches
#ifndef NDEBUG
  // Check that everey lit is at most once on the autarky_trail
  {
    vector<signed char> seen (vsize, 0);
    for (const auto lit : walker.autarky_trail) {
      assert (!seen[vidx (lit)]);
      seen[vidx (lit)] = 1;
    }
  }
#endif

  if (!unsat && walker.autarky_mode && walker.autarky_elimination_mode) {

    // Counted calls of autarky_applay() as autarky() in autarky.cpp does
    ++stats.autarkies.tries;

    // if an autarky contain a frozen literal, the autarky cant be removed
    if (!walker.frozen_autarky && !walker.autarky_trail.empty ()) {
      clear_watches ();
      autarky_apply (walker.autarky_val, walker.autarky_trail);

      for (auto idx : vars) {
        if (!walker.autarky_val[vlit(idx)]) continue;
        assert(active(idx));
        mark_eliminated(idx);
      }

      connect_watches ();

      // Count all successfull elimination
      ++stats.autarkies.successful;
      stats.autarkies.eliminated += (int64_t) walker.autarky_trail.size ();

      PHASE ("walk_passat", stats.walk.passat,
             "eliminated autarky of %zu literals",
             walker.autarky_trail.size ());
    } 
  } 
}

} // namespace CaDiCaL
