#include <cassert>
#include <deque>
#include <vector>

#include "../../src/cadical.hpp"

static CaDiCaL::Solver solver;

class TheorySolver : CaDiCaL::ExternalPropagator,
                     CaDiCaL::FixedAssignmentListener {
  std::deque<std::vector<int>> current_trail;
  unsigned dec_count = 0;
  unsigned model_count = 0;
  unsigned clause_count = 0;
  unsigned lit_count = 0;
  bool ready_to_add = false;
  std::vector<int> dec_queue{6, 3, -4, 0, 0, 0, 0};
  std::vector<int> original_clauses{
      1,   0,  -2, 0,   3,   4,  0,  8,   3,  0,  8, -7, 0,   -3, 7,
      -8,  0,  -9, 4,   0,   -9, 8,  0,   -4, -8, 9, 0,  -11, 4,  0,
      -11, 10, 0,  -11, 8,   0,  -4, -10, -8, 11, 0, 12, -9,  0,  12,
      -11, 0,  9,  11,  -12, 0,  -6, -12, 0,  -3, 5, 0};
  std::vector<int> clause_queue{
      -13, 4,  0,   -13, -3, 0,  -4,  3,   13, 0,  -14, 4,  0,
      -14, 10, 0,   -14, -3, 0,  -4,  -10, 3,  14, 0,   15, -13,
      0,   15, -14, 0,   13, 14, -15, 0,   -5, 15, 0};

  void check_model () {
    bool falsified = true;
    int propagating = 0;
    int idx = 0;
    for (auto &t : current_trail) {
      printf ("trail[%d]: ", idx++);
      for (auto &lit : t) {
        printf ("%d ", lit);
      }
      printf ("\n");
    }
    for (auto &original : original_clauses) {
      printf ("%d", original);
      if (!original) {
        printf ("\n");
        assert (!falsified);
        assert (propagating != 1);
        falsified = true;
        propagating = 0;
        continue;
      }
      bool found = false;
      int lvl_idx = -1;
      for (auto &t : current_trail) {
        lvl_idx++;
        if (found)
          break;
        for (auto &lit : t) {
          if (lit == original) {
            printf ("=1@%d", lvl_idx);
            propagating = -1;
            found = true;
            falsified = false;
            break;
          } else if (lit == -original) {
            printf ("=-1@%d", lvl_idx);
            found = true;
            break;
          }
        }
      }
      if (!found) {
        printf ("=0");
        falsified = false;
        if (!propagating)
          propagating = 1;
        else
          propagating = -1;
      }
      printf (" ");
    }
    if (!lit_count) {
      printf ("not yet added additional clauses\n");
      return;
    }
    printf ("added additional clauses\n");
    int cq_idx = 0;
    for (auto &original : clause_queue) {
      if (lit_count <= cq_idx++)
        break;
      printf ("%d", original);
      if (!original) {
        printf ("\n");
        assert (!falsified);
        assert (propagating != 1);
        falsified = true;
        propagating = 0;
        continue;
      }
      bool found = false;
      int lvl_idx = -1;
      for (auto &t : current_trail) {
        lvl_idx++;
        if (found)
          break;
        for (auto &lit : t) {
          if (lit == original) {
            printf ("=1@%d", lvl_idx);
            found = true;
            propagating = -1;
            falsified = false;
            break;
          } else if (lit == -original) {
            printf ("=-1@%d", lvl_idx);
            found = true;
            break;
          }
        }
      }
      if (!found) {
        printf ("=0");
        falsified = false;
        if (!propagating)
          propagating = 1;
        else
          propagating = -1;
      }
      printf (" ");
    }
  }

public:
  TheorySolver () {
    solver.connect_fixed_listener (this);
    solver.connect_external_propagator (this);

    solver.add_observed_var (1);
    solver.add_observed_var (2);

    current_trail.push_back (std::vector<int> ());
  }

  ~TheorySolver () { solver.disconnect_external_propagator (); }

  int cb_decide () override {
    check_model ();
    int res = dec_queue[dec_count];
    if (res)
      dec_count++;
    else if (dec_count == 3) {
      solver.add_observed_var (13);
      solver.add_observed_var (14);
      solver.add_observed_var (15);
      ready_to_add = true;
      dec_count++;
    }
    /*
    else if (dec_count == 4) {
      dec_count++;
      return -4;
    } else if (dec_count == 5) {
      dec_count++;
      return 10;
    }
    */
    return res;
  }

  void notify_fixed_assignment (int) override {};

  int cb_propagate () override { return 0; };

  void notify_assignment (const std::vector<int> &lits) override {
    for (auto const lit : lits) {
      current_trail.back ().push_back (lit);
    }
  };

  void notify_new_decision_level () override {
    current_trail.push_back (std::vector<int> ());
  };

  void notify_backtrack (size_t new_level) override {
    while (current_trail.size () > new_level + 1) {
      current_trail.pop_back ();
    }
  };

  bool cb_check_found_model (const std::vector<int> &model) override {
    return true;
  }

  bool cb_has_external_clause (bool &is_forgettable) override {
    return ready_to_add;
  };

  int cb_add_external_clause_lit () override {
    if (lit_count == clause_queue.size () - 1)
      ready_to_add = false;
    return clause_queue[lit_count++];
  }
};

int main () {
  // Some solver specific configurations
  solver.set ("log", 1);
  solver.set ("check", 0);
  solver.set ("walk", 0);
  solver.set ("lucky", 0);
  solver.set ("ilb", 0);
  solver.set ("ilbassumptions", 0);

  solver.set ("chronoadd", -1);

  TheorySolver tsolver;

  solver.add (1);
  solver.add (0);
  solver.add (-2);
  solver.add (0);

  solver.add_observed_var (3);
  solver.add_observed_var (4);
  solver.add (3);
  solver.add (4);
  solver.add (0);

  solver.add_observed_var (5);
  solver.add_observed_var (6);
  solver.add_observed_var (7);
  solver.add_observed_var (8);

  solver.add (8);
  solver.add (3);
  solver.add (0);
  solver.add (8);
  solver.add (-7);
  solver.add (0);
  solver.add (-3);
  solver.add (7);
  solver.add (-8);
  solver.add (0);

  solver.add_observed_var (9);
  solver.add (-9);
  solver.add (4);
  solver.add (0);
  solver.add (-9);
  solver.add (8);
  solver.add (0);
  solver.add (-4);
  solver.add (-8);
  solver.add (9);
  solver.add (0);

  solver.add_observed_var (10);
  solver.add_observed_var (11);
  solver.add (-11);
  solver.add (4);
  solver.add (0);
  solver.add (-11);
  solver.add (10);
  solver.add (0);
  solver.add (-11);
  solver.add (8);
  solver.add (0);
  solver.add (-4);
  solver.add (-10);
  solver.add (-8);
  solver.add (11);
  solver.add (0);

  solver.add_observed_var (12);
  solver.add (12);
  solver.add (-9);
  solver.add (0);
  solver.add (12);
  solver.add (-11);
  solver.add (0);
  solver.add (9);
  solver.add (11);
  solver.add (-12);
  solver.add (0);

  solver.add (-6);
  solver.add (-12);
  solver.add (0);

  solver.phase (6);

  solver.add (-3);
  solver.add (5);
  solver.add (0);

  solver.solve ();
  return 0;
}
