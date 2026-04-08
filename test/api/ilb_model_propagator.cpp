#include "../../src/cadical.hpp"

#ifdef NDEBUG
#undef NDEBUG
#endif

#include <cassert>
#include <iostream>
#include <vector>

#define BIG_NUM 100000

//--------------------------------------------------------------------------//
//
class ILBPropagator : CaDiCaL::ExternalPropagator {
  CaDiCaL::Solver *solver;

  int last = 0;
  bool adding = 0;

public:
  ILBPropagator (CaDiCaL::Solver *s) : solver (s) {
    is_lazy = true;
    are_reasons_forgettable = false;

    solver->connect_external_propagator (this);

    // Fatal error: can not connect more than one propagator
    // solver->connect_external_propagator (this);

    std::cout << "[ilb-prop] ILBPropagator is created and connected "
              << std::endl;

    for (int i = 1; i < BIG_NUM; i++)
      solver->add_observed_var (i);
  }

  ~ILBPropagator () {
    std::cout << "[ilb-prop][~ILBPropagator starts]" << std::endl;

    solver->disconnect_external_propagator ();
  }

  void notify_assignment (const std::vector<int> &lits) { (void) lits; }

  void notify_new_decision_level () {}

  void notify_backtrack (size_t new_level) { (void) new_level; }

  bool cb_check_found_model (const std::vector<int> &model) {
    (void) model;
    if (__builtin_popcount (last) == 1)
      std::cout << "[ilb-prop][cb_check_model][" << last << "]"
                << std::endl;
    adding = true;
    if (++last == BIG_NUM)
      return true;
    return false;
  }

  int cb_decide () { return 0; };

  int cb_propagate () { return 0; };

  int cb_add_reason_clause_lit (int propagated_lit) {
    (void) propagated_lit;
    // In this example we have no used external propagations, so this
    // function will not be called by the SAT solver.

    // For the possible error scenarios see the ActivePropagator example.

    assert (false);
    return 0;
  };

  bool cb_has_external_clause (bool &is_forgettable) {
    (void) is_forgettable;

    // Fatal error: not allowed to force backtrack in that state of the
    // solver
    // solver->force_backtrack(0);

    return adding;
  }

  int cb_add_external_clause_lit () {
    if (adding) {
      adding = false;
      return last;
    }
    return 0;
  }
};

int main () {
  std::cout << "-----------------------------------------------------------"
            << std::endl;
  CaDiCaL::Solver *solver = new CaDiCaL::Solver;
  ILBPropagator *ilb = new ILBPropagator (solver);
  for (int i = 1; i < BIG_NUM - 1; i++) {
    solver->add (-i);
    solver->add (i + 1);
    solver->add (0);
  }
  solver->solve ();
  delete ilb;
  delete solver;
  return 0;
}
