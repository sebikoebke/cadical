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
  int adding = 0;
  bool done = 0;

public:
  ILBPropagator (CaDiCaL::Solver *s) : solver (s) {
    is_lazy = false;
    are_reasons_forgettable = false;

    solver->connect_external_propagator (this);

    // Fatal error: can not connect more than one propagator
    // solver->connect_external_propagator (this);

    std::cout << "[ilb-prop] ILBPropagator is created and connected "
              << std::endl;

    for (int i = 1; i < 2 * BIG_NUM; i++) {
      /*
      if (__builtin_popcount (i) == 1)
        std::cout << "[ilb-prop][add_observed_var][" << i << "]"
                  << std::endl;
                  */

      solver->add_observed_var (i);
    }
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
    std::cout << "[ilb-prop][cb_check_model][" << done << "]" << std::endl;
    if (done)
      return true;
    adding = 1;
    done = true;
    return false;
  }

  int cb_decide () {
    if (last < 2 * BIG_NUM - 1)
      return ++last;
    return 0;
  };

  int cb_propagate () {
    if (last > 1 && last % 2 == 0 && last < 2 * BIG_NUM - 1) {
      if (__builtin_popcount (last) == 1)
        std::cout << "[ilb-prop][cb_propagate][" << last + 1 << "]"
                  << std::endl;
      return ++last;
    }
    return 0;
  };

  int cb_add_reason_clause_lit (int propagated_lit) {
    if (!adding) {
      adding = 1;
      return propagated_lit;
    } else if (adding) {
      adding = -1;
      return -(propagated_lit - 2);
    }
    adding = 0;
    return 0;
  };

  bool cb_has_external_clause (bool &is_forgettable) {
    (void) is_forgettable;

    return adding;
  }

  int cb_add_external_clause_lit () {
    if (adding > 0) {
      std::cout << "[ilb-prop][cb_add_external_clause_lit][-1]"
                << std::endl;
      adding = -1;
      return -1;
    } else if (adding < 0) {
      std::cout << "[ilb-prop][cb_add_external_clause_lit]["
                << -(2 * BIG_NUM - 1) << "]" << std::endl;
      adding = 0;
      return -(2 * BIG_NUM - 1);
    }
    std::cout << "[ilb-prop][cb_add_external_clause_lit][0]" << std::endl;
    return 0;
  }
};

int main () {
  std::cout << "-----------------------------------------------------------"
            << std::endl;
  CaDiCaL::Solver *solver = new CaDiCaL::Solver;
  solver->set ("exteagerrecalc", 0);
  solver->set ("exteagerreasons", 0);
  ILBPropagator *ilb = new ILBPropagator (solver);
  for (int i = 1; i < 2 * BIG_NUM; i++) {
    solver->add (i);
  }
  solver->add (0);
  solver->solve ();
  delete ilb;
  delete solver;
  return 0;
}
