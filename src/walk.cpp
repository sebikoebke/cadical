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

  // for efficiency, storing the model each time an improvement is
  // found is too costly. Instead we store some of the flips since
  // last time and the position of the best model found so far.
  Random random;                 // local random number generator
  int64_t ticks;                 // ticks to approximate run time
  int64_t limit;                 // limit on number of propagations
  vector<ClauseOrBinary> broken; // currently unsatisfied clauses
  double epsilon;                // smallest considered score
  vector<double> table;          // break value to score table
  vector<double> scores;         // scores of candidate literals
  vector<pair<double, int>> scores_passat;   //maybe we can safe one loop with that structure
  std::vector<int> flips; // remember the flips compared to the last best saved model
  int best_trail_pos;
  int64_t minimum = INT64_MAX;
  vector<int> propagation_queue;  // our replacement for the trail in walk_passat (with a difference): the list of all literals which are assigned or should. 
                                  // Changed if probSAT_repair is finished with LS_repair: Old propagations are not longer needed, therefore ther are cut. Just pending and flipped variables remain.
  size_t propagated = 0;          // how far propagation_queue has been processed (like Internal::propagated)
  size_t activated = 0;           // counts assigned variables; up_expansion stops once all are
                                  // activated. Replaces ordering_O/ls_score: variable + polarity
                                  // selection is delegated to CaDiCaL's own decision heuristic.
  vector<vector<int>> passat_lookup_table; // positions in `clauses` where v+/v- occurs
  vector<int> passat_clauses;     // all clauses where all variables inside are activated
  vector<int> broken_clauses;     // all broken clauses out of passat_clauses, used for probSAT_repair
  vector<int> broken_pos;         // position of a clause inside broken_clauses (indexed by clause pos), -1 if not broken; enables fast removal (like in WalkerFO)
  vector<int> conflict_counter;   // counter which shows if there is a conflict inside a clause, if c_c == 0 => conflict, decreased if the opposite polarity is assigned to true
  vector<int> activation_counter; // counter which shows when a clauses contain only activated clauses. a_c is decreased, if one variable in a clause is activated
  vector<int> flip_count;         // LS Hotspots, indexed by variables
  vector<signed char> mark;       // per-variable dedup flag, invariant 0 outside repair_propagation_queue
  vector<int> cache_queue;            // reusable cache to rebuild propagation_queue without allocating
  std::vector<signed char> best_values; // best model stored so far
  double score (unsigned);              // compute score from break count
#ifndef NDEBUG
  std::vector<signed char> current_best_model; // best model found so far
  size_t tracked_clauses = 0;   // # clauses passat_build tracks (non-garbage,
                                // non-skipped-redundant); used to assert in
                                // up_expansion that no clause was forgotten
#endif
  Walker (Internal *, int64_t limit);
  void populate_table (double size);
  void push_flipped (int flipped);
  void save_walker_trail (bool);
  void save_final_minimum (int64_t old_minimum);
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
      ticks (0), limit (l), best_trail_pos (-1) {
  random += internal->stats.walk.count; // different seed every time
  flips.reserve (i->max_var / 4);
  mark.resize (i->max_var + 1, 0); // dedup flag for repair_propagation_queue, kept invariant 0
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
  const bool use_size_based_cb = (internal->stats.walk.count & 1);
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
void Walker::save_final_minimum (int64_t old_init_minimum) {
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
  int64_t size = walker.broken.size ();
  if (size > INT_MAX)
    size = INT_MAX;
  int pos = walker.random.pick_int (0, size - 1);
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
  ticks += (1 + cache_lines (watches (lit).size (), sizeof (Clause *)));

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
  walker.ticks += 1;
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
  const int tmp = sign (lit);
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
    // broken is in cache given how central it is... but not always (see the
    // ncc problems). Value was heuristically determined to give reasonnable
    // values.
    walker.ticks +=
        1 + cache_lines (walker.broken.size (), sizeof (Clause *));
    auto j = walker.broken.begin (), i = j;
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
    walker.ticks += 1 + cache_lines (ws.size (), sizeof (Clause *));

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
  int64_t broken = walker.broken.size ();
  if (broken >= walker.minimum)
    return;
  if (broken <= stats.walk.minimum) {
    stats.walk.minimum = broken;
    VERBOSE (3, "new global minimum %" PRId64 "", broken);
  } else {
    VERBOSE (3, "new walk minimum %" PRId64 "", broken);
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
    VERBOSE (3, "saving the new walk minimum %" PRId64 "", broken);
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
  int old_global_minimum = stats.walk.minimum;
#endif

  bool failed = false; // Inconsistent assumptions?

  level = 1; // Assumed variables assigned at level 1.

  if (assumptions.empty ()) {
    LOG ("no assumptions so assigning all variables to decision phase");
  } else {
    LOG ("assigning assumptions to their forced phase first");
    for (const auto lit : assumptions) {
      signed char tmp = val (lit);
      if (tmp > 0)
        continue;
      if (tmp < 0) {
        LOG ("inconsistent assumption %d", lit);
        failed = true;
        break;
      }
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
    int64_t watched = 0;
#endif

    double size = 0;
    int64_t n = 0;
    for (const auto c : clauses) {

      if (c->garbage)
        continue;
      if (c->redundant) {
        if (!opts.walkredundant)
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
      int64_t broken = walker.broken.size ();
      int64_t total = watched + broken;
      LOG ("watching %" PRId64 " clauses %.0f%% "
           "out of %" PRId64 " (watched and broken)",
           watched, percent (watched, total), total);
    }
#endif
  }

  assert (failed || walker.table.size ());

  int res; // Tells caller to continue with local search.

  if (!failed) {

    int64_t broken = walker.broken.size ();
    int64_t initial_minimum = broken;

    PHASE ("walk", stats.walk.count,
           "starting with %" PRId64 " unsatisfied clauses "
           "(%.0f%% out of %" PRId64 ")",
           broken, percent (broken, stats.current.irredundant),
           stats.current.irredundant);

    walk_save_minimum (walker);
    assert (stats.walk.minimum <= walker.minimum);

    int64_t minimum = broken;
#ifndef QUIET
    int64_t flips = 0;
#endif
    while (!terminated_asynchronously () && !walker.broken.empty () &&
           walker.ticks < walker.limit) {
#ifndef QUIET
      flips++;
#endif
      stats.walk.flips++;
      stats.walk.broken += broken;
      ClauseOrBinary c = walk_pick_clause (walker);
      const int lit = walk_pick_lit (walker, c);
      bool finished = walk_flip_lit (walker, lit);
      if (!finished)
        break;
      walker.push_flipped (lit);
      broken = walker.broken.size ();
      LOG ("now have %" PRId64 " broken clauses in total", broken);
      if (broken >= minimum)
        continue;
      minimum = broken;
      VERBOSE (3, "new phase minimum %" PRId64 " after %" PRId64 " flips",
               minimum, flips);
      walk_save_minimum (walker);
    }

    walker.save_final_minimum (initial_minimum);

#ifndef QUIET
    if (minimum == initial_minimum) {
      PHASE ("walk", internal->stats.walk.count,
             "%sno improvement %" PRId64 "%s in %" PRId64 " flips and "
             "%" PRId64 " ticks",
             tout.bright_yellow_code (), minimum, tout.normal_code (),
             flips, walker.ticks);
    } else if (minimum < old_global_minimum)
      PHASE ("walk", stats.walk.count,
             "%snew global minimum %" PRId64 "%s in %" PRId64 " flips and "
             "%" PRId64 " ticks",
             tout.bright_yellow_code (), minimum, tout.normal_code (),
             flips, walker.ticks);
    else
      PHASE ("walk", stats.walk.count,
             "best phase minimum %" PRId64 " in %" PRId64 " flips and "
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
      LOG ("minimum %" PRId64 " non-zero thus potentially continue",
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
  VERBOSE (2,
           "walk scheduling: last %" PRId64 " current %" PRId64
           " delta %" PRId64,
           last.walk.ticks, ticks, limit);
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
  (void) walk_round (limit, false);
  STOP_INNER_WALK ();
  assert (!unsat);
}

/*----------------------------------------------------------------------------*/

// passat_build() prepares ...
// (a) passat_lookup_table
// (b) conflict_counter[pos] : counter for each clause, starts with 
//     conflict_counter[pos] = clauses[pos].size() and decrease, if a literal 
//     inside the clause is wrong. If conflict_counter ==0 => empty clause (conflict)
// (c) activation_counter[pos] : counter for the repair modul to show which clauses
//     are fully activated
// (d) walker.activated : number of variables that already carry a value, so
//     up_expansion knows when every variable has been activated (SAT case).
// (e) the ProbSAT score table, built from the average size of the tracked
//     clauses (like walk() does) and used by probSAT_pick_lit.
void Internal::passat_build (Walker &walker) {
  walker.passat_lookup_table.resize (2 * vsize);
  walker.conflict_counter.resize (clauses.size ());
  walker.activation_counter.resize (clauses.size ());
  walker.broken_pos.resize (clauses.size (), -1);

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
#ifndef NDEBUG
    // this clause is tracked (passes the same filter as up_expansion's check)
    walker.tracked_clauses++;
#endif

    // count this tracked clause for the average size (part (e))
    total_size += c->size;
    counted++;

    // literals with val != -1
    int not_false = 0;
    // literals with val == 0
    int unassigned = 0;
    for (const auto lit : *c) {
      const signed char v = val (lit);
      // if a literal is not false, we have to increase the conflict_counter of the clause
      if (v >= 0)
        not_false++;
      // if a literal is unassigned, we have to increase the activation_counter of the clause
      if (v == 0)
        unassigned++;
      // part of (a)
      walker.passat_lookup_table[vlit (lit)].push_back ((int) pos);
    }
    // (b)
    walker.conflict_counter[pos] = not_false;
    // (c)
    walker.activation_counter[pos] = unassigned;
    // already fully activated clause
    if (unassigned == 0)
      walker.passat_clauses.push_back ((int) pos);
  }

  // part (d)
  walker.activated = 0;
  for (int idx = 1; idx <= max_var; idx++)
    if (val (idx))
      walker.activated++;

  // build the ProbSAT score table from the average clause size, like walk()
  walker.populate_table (relative (total_size, counted));
}

/*----------------------------------------------------------------------------*/

// passat_assign() assigns the literal 'lit' to true and keeps all PASSAT
// bookkeeping consistent. It performs the following steps:
//   1. set 'lit' to true (does nothing if it is already assigned)
//   2. push 'lit' onto the propagation_queue so passat_up can propagate it
//   3. decrement the activation_counter of every clause that contains the
//      variable; once it hits 0 the clause is fully activated and gets
//      appended to passat_clauses
//   4. decrement the conflict_counter of every clause that contains '-lit'
//      (that literal just turned false); if it hits 0 the clause is falsified
//      and we report a conflict
// Returns true if no conflict arose (propagation may continue), false if
// assigning 'lit' falsified at least one clause.
bool Internal::passat_assign(Walker &walker, int lit) {
  bool signal = true;
  // Only assign if 'lit' is currently unassigned; otherwise there is nothing
  // to do and we report that propagation may continue.
  // NOTE: the guard returns true for any val(lit) != 0. If 'lit' were already
  // false (val(lit) == -1) this would silently hide a contradiction. In our
  // design that should never happens: passat_up only calls passat_assign on
  // unassigned unit literals, and a real conflict is caught earlier via
  // conflict_counter == 0 (just for debugging).
  if (val(lit) == 0){
    set_val(lit, 1);
    // count this activation so up_expansion knows when all variables are assigned
    walker.activated++;
    // (2) enqueue for later propagation in passat_up
    walker.propagation_queue.push_back(lit);
    // (3) positive occurrences: the variable becomes active in these clauses
    for(auto clause : walker.passat_lookup_table[vlit(lit)]) {
      walker.activation_counter[clause] -= 1;
      if (walker.activation_counter[clause] == 0) {
        walker.passat_clauses.push_back(clause);
      }
    }
    // (3) negative occurrences: the variable becomes active here too, and
    // (4) '-lit' is now false, so the conflict_counter shrinks as well
    for(auto clause : walker.passat_lookup_table[vlit(-lit)]){
      walker.activation_counter[clause] -= 1;
      if (walker.activation_counter[clause] == 0) {
        walker.passat_clauses.push_back(clause);
      }

      walker.conflict_counter[clause] -= 1;
      // a clause whose conflict_counter hit 0 is falsified => report a conflict.
      // probSAT_repair rebuilds all broken clauses from passat_clauses anyway, so
      // we only need to signal that some conflict occurred, not which clause.
      if (walker.conflict_counter[clause] == 0)
        signal = false;
    }
  }
  return signal;
}

/*----------------------------------------------------------------------------*/

// Unit propagation for the up_expansion module.
// passat_up is working on walker.propagation_queue and with walker.propagated
bool Internal::passat_up(Walker &walker){
  // For every assigned literal still to be processed, look for clauses that
  // just became unit and propagate them through passat_assign.
  while (walker.propagated < walker.propagation_queue.size()){
    const int lit = walker.propagation_queue[walker.propagated++];
    // 'lit' is true, so only clauses containing '-lit' can shrink: clauses
    // that contain 'lit' are already satisfied.
    // Reading the occurrence row is charged like walk() charges reading a
    // watch list, so that the tick measurement stays comparable to walk().
    const auto &lit_clauses = walker.passat_lookup_table[vlit(-lit)];
    walker.ticks += 1 + cache_lines(lit_clauses.size(), sizeof(int));

    for (auto clause : lit_clauses){
      // one tick per clause we actually visit
      walker.ticks++; 
      // conflict_counter == 1: exactly one literal is still not false. Find it:
      // if it is unassigned the clause is unit and must be propagated; if it is
      // already true the clause is satisfied and we skip it.
      if (walker.conflict_counter[clause] == 1){
        Clause *c = clauses[clause];
        int unit = 0;
        for (auto clause_lit : *c){
          // scanning the clause touches each literal
          walker.ticks++;
          if (val(clause_lit) >= 0){
            unit = clause_lit;
            break;
          }
        }
        if (val(unit) == 0){
          // Propagating the unit may itself falsify a clause => passat_assign
          // then returns false.
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
  // First execute anything pending in propagation_queue
  // e.g. flips that probSAT_repair re-enqueued
  if (!passat_up(walker)) return false;

  // Loop until every variable is activated
  while (walker.activated < (size_t) max_var) {
    // pick a next unassigned variable to assign 
    // because no propagation is left on the propagation_queue
    const int idx = use_scores () ? next_decision_variable_with_best_score ()
                                  : next_decision_variable_on_queue ();
    const bool target = (stable || opts.target == 2);
    // chose the polarity for the choosen variable idx
    const int lit = decide_phase (idx, target);

    // activate the literal. On conflict hand over to probSAT_repair
    if (!passat_assign(walker, lit)) return false;
    // propagate the consequences. On conflict hand over to probSAT_repair
    if (!passat_up(walker)) return false;
  }
  // all variables assigned, no conflict => SAT
#ifndef NDEBUG
  // Completeness: with every variable assigned, every tracked clause must be
  // fully activated, i.e. present in passat_clauses. Comparing the count rules
  // out a "forgotten" clause that the satisfaction loop below would not see.
  assert (walker.passat_clauses.size () == walker.tracked_clauses);
  // Satisfaction: no tracked clause is broken. conflict_counter only drops in
  // passat_assign (which returns false on 0), so this should always hold here.
  for (const int pos : walker.passat_clauses)
    assert (walker.conflict_counter[pos] > 0);
#endif
  return true;
}

/*----------------------------------------------------------------------------*/

// Function which build the list of broken clauses
void Internal::build_broken(Walker &walker){
  walker.broken_clauses.clear();
  for (int i : walker.passat_clauses){
    // satisfied => not broken
    if (walker.conflict_counter[i] > 0) {
      walker.broken_pos[i] = -1;
      continue;
    }
    // remember the position of the broken clause inside broken_clauses
    // it is later cheaper to look up and flip the variable and remove the
    // corresponding clause from broken_clauses
    walker.broken_pos[i] = (int) walker.broken_clauses.size();
    walker.broken_clauses.push_back(i);
  }
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
  return list_of_clauses[walker.random.pick_int (0, size - 1)];
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
// - Problem: Currently, our bottleneck is in probSAT_pick_lit,
// because probSAT_break_value performs a large number of loop-iterations
// - Approach: Instead of calculating the break value, we simply use the size of
// |occ(-lit)|
// - Loss: Weighted randomness (as before), flips that cause little
// damage are not favored.
// - Gain: One query (O(1)) instead of having to go through the list occ(-lit)
// - Side effect: We have weighted randomness even though we do not read
// the conflict_counter, but this does not need to exist and is smaller.
// Example: Two variables a and b. If |occ(a)| = 4 and |occ(b)| = 100,
// then break_value(a) < break_value(b) is more likely
unsigned Internal::passat_fixed_occurence(Walker &walker, int lit){
  walker.ticks += 1;
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
  walker.ticks += 1;
  double sum = 0;

  // Phase 1: score every flippable literal by its break-value.
  for (const auto lit : *c){
    // an assumed variable is not allowed to be flipped
    if (assumed(lit) || assumed(-lit))
      continue;
    // opts.passatfixedocc: trade the exact (state-aware) break value for the
    // static occurrence count of -lit (O(1) instead of O(|occ(-lit)|)).
    const unsigned bv = opts.passatfixedocc ? passat_fixed_occurence(walker, lit)
                                            : passat_break_value(walker, lit);
    const double s = walker.score(bv);
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
// 2. delete all clauses c´ from broken_clauses
// 3. increase the conflict_counter in all clauses of c´
// 4. decrease the conflict_counter in all clauses which dont occur in broken_clauses
//    but in passat_clauses. => If conflict_counter goes down to 0, add clause to
//    broken clauses
// 5. add the flipped variable to walker.flips 
void Internal::flip_and_repair(Walker &walker, int lit){
  // 1.
  set_val(lit, 1);

  // 2. and 3. 
  const auto &row = walker.passat_lookup_table[vlit(lit)];
  walker.ticks += 1 + cache_lines(row.size(), sizeof(int));
  for (int c : row){
    walker.ticks++;
    walker.conflict_counter[c] += 1;
    // conflict_counter == 1 means the clause was broken and is now
    // satisfied => remove it from broken_clauses
    if (walker.conflict_counter[c] == 1){
      int change_position = walker.broken_pos[c];
      int last_element = walker.broken_clauses.back();
      // move the last element into the new free slot
      walker.broken_clauses[change_position] = last_element;
      // update broken_pos
      walker.broken_pos[last_element] = change_position;
      walker.broken_pos[c] = -1;
      // cut off duplicate of last element
      walker.broken_clauses.pop_back();
    }
  }

  // 4.
  const auto &neg_row = walker.passat_lookup_table[vlit(-lit)];
  walker.ticks += 1 + cache_lines(neg_row.size(), sizeof(int));
  for (int c : neg_row){
    walker.ticks++;
    walker.conflict_counter[c] -= 1;
    // conflict_counter == 0 => clause is now broken, append it
    if (walker.conflict_counter[c] == 0) {
      walker.broken_pos[c] = (int) walker.broken_clauses.size();
      walker.broken_clauses.push_back(c);
    }
  }

  // 5.
  walker.flips.push_back(lit);

}

/*----------------------------------------------------------------------------*/

// Rebuilds the propagation_queue after a probSAT_repair run so that up_expansion
// can resume propagation. We do NOT keep the old queue: since walk_passat never
// backtracks, already-processed entries are history that is never read
// again (their value lives in vals[] and their consequences were drawn earlier).
void Internal::repair_propagation_queue(Walker &walker){
  // collect pending variables, mark is empty here, so no check is needed
  // we dont need to think about the correct polarity insertion of v in cache_queue
  // because vals[] already contain the correct polarity after probSAT_repair
  // therefore cache_queue is up to date here
  for (size_t k = walker.propagated; k < walker.propagation_queue.size(); k++){
    const int v = abs(walker.propagation_queue[k]);
    walker.mark[v] = 1;
    walker.cache_queue.push_back(val(v) > 0 ? v : -v);
  }
  // collect flipped variables, mark is not empty here, so we need to check if 
  // variable is already in the cache_queue
  for (const int lit : walker.flips){
    const int v = abs(lit);
    if (!walker.mark[v]){
      walker.mark[v] = 1;
      walker.cache_queue.push_back(val(v) > 0 ? v : -v);
    }
  }
  // reset mark. We consider only the instances of chache_queue 
  // which should be cheaper than reset each of the |max_var| - positions by 0
  for (const int lit : walker.cache_queue)
    walker.mark[abs(lit)] = 0;

  // make cache_queue the new queue in O(1) and discard the old one
  std::swap(walker.propagation_queue, walker.cache_queue);
  walker.cache_queue.clear();
  walker.propagated = 0;
}

// Returns true if the conflict was fully resolved (broken == 0) so that
// up_expansion can resume; false if it could not be repaired (=> UNSAT).
bool Internal::probSAT_repair(Walker &walker) {
  /*
  1. build a list of broken clauses where ProbSAT should operate on. The list with 
     all clauses where ProbSAT is able to look at already exist => PASSAT_CLAUSES
  2. Flip decision of Literal l (wie auch in walk() vernwendet):
     Wir wählen zufällig eine clause c aus broken und dann mit probsat ein literal l aus c.
  3. Flip Literal l (which is decide in 2.). 
  4. Make all adjustment and book keeping stuff.
  5. Look if there is still broken clauses, then return to (2), otherwise SAT and break the loop
  6. If a solution is found, write the solution in vals[] (maybe it is already written), wirte the new values on the propagation_queue (we need to do this correct?) 
     and pass to up_expansion
  */

  build_broken(walker);
  //Clear all earlier made flips
  walker.flips.clear();
  while(!walker.broken_clauses.empty() && walker.ticks < walker.limit){
    // pick random clause, then a literal of it via ProbSAT
    int lit = probSAT_pick_lit(walker, pick_random_clause(walker, walker.broken_clauses));

    // lit == 0: the clause consists only of assumptions, so it cannot be
    // repaired without violating an assumption => unsatisfiable under the
    // assumptions, stop and report failure.
    if (lit == 0)
      return false;

    // make the actual LS_repair
    flip_and_repair(walker, lit);
  }
  // Tick limit reached while clauses are still broken: the conflict could not be
  // resolved within the budget
  if (!walker.broken_clauses.empty())
    return false;

  // broken == 0: the conflict is fully resolved. Rebuild the propagation_queue so
  // up_expansion can resume on the repaired partial assignment, then report success.
  repair_propagation_queue(walker);
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

    backtrack ();
  
  //propagate() is called if unpropagated literals are still present in the trail after backtrack()
  if (propagated < trail.size () && !propagate ()) {
    LOG ("empty clause after root level propagation");
    learn_empty_clause ();
    STOP_INNER_WALK ();
    return;
  }

  //calc limit, identically as in walk()
  const int64_t start_ticks = stats.ticks.search[0] + stats.ticks.search[1];
  int64_t limit = start_ticks - last.walk.ticks;
  last.walk.ticks = start_ticks;
  limit *= 1e-3 * opts.walkeffort;
  if (limit < opts.walkmineff) limit = opts.walkmineff;

  Walker walker (internal, limit);

  // build occurrence lists, counters and the activation count for this run
  passat_build (walker);

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
  if (consistent_with_assumptions){
    no_conflict = true;
    while (walker.ticks < walker.limit) {
      no_conflict = up_expansion (walker);
      if (no_conflict)
        break;                        // SAT over the activated set
      if (!probSAT_repair (walker))
        break;                        // conflict not resolvable -> UNSAT
    }
  }

  LOG("walk_passat: %s", no_conflict ? "SAT" : "limit reached");

  // Save the result: the whole point of walk_passat is to leave a better
  // polarity assignment in phases.saved, which the following CDCL search uses
  // as decision phases
  for (int id = 1; id <= max_var; id++)
    if (val (id))
      phases.saved[id] = val (id);

  // vals[] is the solver's single global assignment
  // table with the invariant "vals[v] != 0 <=> v is on the trail" 
  // We broke that invariant by assigning via passat_assign() without pushing to the trail
  // Restore it by clearing vals[] for all active variables and reset
  // the decision level to the root, otherwise the next CDCL search runs on a
  // corrupted state (like in walk_round)
  for (int id = 1; id <= max_var; id++)
    if (active (id))
      set_val (id, 0);
  level = 0;

  STOP_INNER_WALK();
}

} // namespace CaDiCaL
