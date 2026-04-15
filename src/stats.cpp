// vim: set tw=300: set VIM text width to 300 characters for this file.

#include "internal.hpp"

namespace CaDiCaL {

/*------------------------------------------------------------------------*/

Stats::Stats () {
  time.real = absolute_real_time ();
  time.process = absolute_process_time ();
  walk_minimum = LONG_MAX;
  used[0].resize (127);
  used[1].resize (127);
}

/*------------------------------------------------------------------------*/

#define PRT(FMT, ...) \
  do { \
    if (FMT[0] == ' ' && !all) \
      break; \
    MSG (FMT, __VA_ARGS__); \
  } while (0)

/*------------------------------------------------------------------------*/

void Stats::print (Internal *internal) {

#ifdef QUIET
  (void) internal;
#else

  Stats &stats = internal->stats;

  int all = internal->opts.verbose > 0 || internal->opts.stats;
#ifdef LOGGING
  if (internal->opts.log)
    all = true;
#endif // ifdef LOGGING

  if (internal->opts.profile)
    internal->print_profile ();

  double t = internal->solve_time ();

  int64_t propagations = 0;
  propagations += stats.propagations_cover;
  propagations += stats.propagations_probe;
  propagations += stats.propagations_search;
  propagations += stats.propagations_transred;
  propagations += stats.propagations_vivify;

  int64_t vivified = stats.vivify_subsumed + stats.vivify_strengthened +
                     stats.vivify_implied;
  int64_t searchticks =
      stats.ticks_search_stable + stats.ticks_search_unstable;
  int64_t inprobeticks = stats.ticks_vivify + stats.ticks_probe +
                         stats.ticks_factor + stats.ticks_ternary +
                         stats.ticks_sweep + stats.ticks_backbone;
  int64_t totalticks = searchticks + inprobeticks;

  size_t extendbytes = internal->external->extension.size ();
  extendbytes *= sizeof (int);

  SECTION ("statistics");

  if (all || stats.blocked) {
    PRT ("blocked:         %15" PRId64
         "   %10.2f %%  of irredundant clauses",
         stats.blocked,
         percent (stats.blocked, stats.clauses_added_irredundant));
    PRT ("  blockings:     %15" PRId64 "   %10.2f    internal",
         stats.blockings, relative (stats.conflicts, stats.blockings));
    PRT ("  candidates:    %15" PRId64 "   %10.2f    per blocking ",
         stats.blocked_candidates,
         relative (stats.blocked_candidates, stats.blockings));
    PRT ("  blockres:      %15" PRId64 "   %10.2f    per candidate",
         stats.blocked_resolutions,
         relative (stats.blocked_resolutions, stats.blocked_candidates));
    PRT ("  pure:          %15" PRId64 "   %10.2f %%  of all variables",
         stats.vars_all_pure, percent (stats.vars_all_pure, stats.vars));
    PRT ("  pureclauses:   %15" PRId64 "   %10.2f    per pure literal",
         stats.blocked_pure,
         relative (stats.blocked_pure, stats.vars_all_pure));
  }
  if (all || stats.conflicts_chrono)
    PRT ("chronological:   %15" PRId64 "   %10.2f %%  of conflicts",
         stats.conflicts_chrono,
         percent (stats.conflicts_chrono, stats.conflicts));
  if (all)
    PRT ("compacts:        %15" PRId64 "   %10.2f    interval",
         stats.compacts, relative (stats.conflicts, stats.compacts));
  if (all || stats.conflicts) {
    PRT ("conflicts:       %15" PRId64 "   %10.2f    per second",
         stats.conflicts, relative (stats.conflicts, t));
    PRT ("  backtracked:   %15" PRId64 "   %10.2f %%  of conflicts",
         stats.backtracked, percent (stats.backtracked, stats.conflicts));
  }
  if (all || stats.incremental_decay) {
    PRT ("inc-decay:       %15" PRId64 "   %10.2f %%   per search",
         stats.incremental_decay,
         percent (stats.incremental_decay, stats.searches));
  }
  if (all || stats.backbone_rounds) {
    PRT ("backbone:        %15" PRId64 "   %10.2f %%  of vars",
         stats.backbone_probes,
         percent (stats.backbone_probes, stats.vars));
    PRT ("   rounds:       %15" PRId64 "   %10.2f    per phase",
         stats.backbone_rounds,
         relative (stats.backbone_rounds, stats.backbone_phases));
    PRT ("   phases:       %15" PRId64 "   %10.2f    interval",
         stats.backbone_phases,
         relative (stats.conflicts, stats.backbone_phases));
    PRT ("   units:        %15" PRId64 "   %10.2f    per phase",
         stats.backbone_units,
         percent (stats.backbone_units, stats.backbone_phases));
  }
  if (all || stats.deduplicate_init) {
    PRT ("dedup-init-rnds: %15" PRId64 "   %10.2f %%  of interval",
         stats.deduplicate_init_rounds,
         percent (stats.deduplicate_init_rounds, stats.conflicts));
    PRT ("dedup-init:      %15" PRId64 "   %10.2f %%  of subsumed",
         stats.deduplicate_init,
         percent (stats.deduplicate_init, stats.deduplicate_init));
  }
  if (all || stats.conditioned) {
    PRT ("conditioned:     %15" PRId64
         "   %10.2f %%  of irredundant clauses",
         stats.conditioned,
         percent (stats.conditioned, stats.clauses_added_irredundant));
    PRT ("  conditionings: %15" PRId64 "   %10.2f    interval",
         stats.conditionings,
         relative (stats.conflicts, stats.conditionings));
    PRT ("  condcands:     %15" PRId64 "   %10.2f    candidate clauses",
         stats.condition_candidates,
         relative (stats.condition_candidates, stats.conditionings));
    PRT ("  condassinit:   %17.1f  %9.2f %%  initial assigned",
         relative (stats.condition_init_assigned, stats.conditionings),
         percent (stats.condition_init_assigned,
                  stats.condition_final_assigned));
    PRT ("  condcondinit:  %17.1f  %9.2f %%  initial condition",
         relative (stats.condition_init_conditional, stats.conditionings),
         percent (stats.condition_init_conditional,
                  stats.condition_init_assigned));
    PRT ("  condautinit:   %17.1f  %9.2f %%  initial autarky",
         relative (stats.condition_init_autarky, stats.conditionings),
         percent (stats.condition_init_autarky,
                  stats.condition_init_assigned));
    PRT ("  condassrem:    %17.1f  %9.2f %%  final assigned",
         relative (stats.condition_final_assigned, stats.conditioned),
         percent (stats.condition_final_assigned,
                  stats.condition_final_assigned));
    PRT ("  condcondrem:   %19.3f  %7.2f %%  final conditional",
         relative (stats.condition_final_conditional, stats.conditioned),
         percent (stats.condition_final_conditional,
                  stats.condition_final_assigned));
    PRT ("  condautrem:    %19.3f  %7.2f %%  final autarky",
         relative (stats.condition_final_autarky, stats.conditioned),
         percent (stats.condition_final_autarky,
                  stats.condition_final_assigned));
    PRT ("  condprops:     %15" PRId64 "   %10.2f    per candidate",
         stats.condition_propagated,
         relative (stats.condition_propagated, stats.condition_candidates));
  }
  if (all || stats.cover_total) {
    PRT ("covered:         %15" PRId64
         "   %10.2f %%  of irredundant clauses",
         stats.cover_total,
         percent (stats.cover_total, stats.clauses_added_irredundant));
    PRT ("  coverings:     %15" PRId64 "   %10.2f    interval",
         stats.coverings, relative (stats.conflicts, stats.coverings));
    PRT ("  asymmetric:    %15" PRId64 "   %10.2f %%  of covered clauses",
         stats.cover_asymmetric,
         percent (stats.cover_asymmetric, stats.cover_total));
    PRT ("  blocked:       %15" PRId64 "   %10.2f %%  of covered clauses",
         stats.cover_blocked,
         percent (stats.cover_blocked, stats.cover_total));
  }
  if (all || stats.decisions) {
    PRT ("decisions:       %15" PRId64 "   %10.2f    per second",
         stats.decisions, relative (stats.decisions, t));
    PRT ("  searched:      %15" PRId64 "   %10.2f    per decision",
         stats.searches, relative (stats.searches, stats.decisions));
  }
  if (all || stats.decisions_random) {
    PRT ("rand. dec phase: %15" PRId64 "   %10.2f    per interval",
         stats.decisions_random,
         relative (stats.decisions_random_phases, stats.decisions));
    PRT ("random decs:     %15" PRId64 "   %10.2f    per phase",
         stats.decisions_random,
         relative (stats.decisions_random, stats.decisions_random_phases));
  }
  if (all || stats.vars_all_eliminated) {
    PRT ("eliminated:      %15" PRId64 "   %10.2f %%  of all variables",
         stats.vars_all_eliminated,
         percent (stats.vars_all_eliminated, stats.vars));
    PRT ("  fastelim:      %15" PRId64 "   %10.2f %%  of eliminated",
         stats.vars_all_eliminated_fast,
         percent (stats.vars_all_eliminated_fast,
                  stats.vars_all_eliminated));
    PRT ("  elimphases:    %15" PRId64 "   %10.2f    interval",
         stats.eliminate_phases,
         relative (stats.conflicts, stats.eliminate_phases));
    PRT ("  elimrounds:    %15" PRId64 "   %10.2f    per phase",
         stats.eliminations,
         relative (stats.eliminations, stats.eliminate_phases));
    PRT ("  elimtried:     %15" PRId64 "   %10.2f %%  eliminated",
         stats.eliminate_tried,
         percent (stats.vars_all_eliminated, stats.eliminate_tried));
    PRT ("  elimgates:     %15" PRId64 "   %10.2f %%  gates per tried",
         stats.eliminate_gates,
         percent (stats.eliminate_gates, stats.eliminate_tried));
    PRT ("  elimequivs:    %15" PRId64 "   %10.2f %%  equivalence gates",
         stats.eliminate_equivalence,
         percent (stats.eliminate_equivalence, stats.eliminate_gates));
    PRT ("  elimands:      %15" PRId64 "   %10.2f %%  and gates",
         stats.eliminate_and,
         percent (stats.eliminate_and, stats.eliminate_gates));
    PRT ("  elimites:      %15" PRId64 "   %10.2f %%  if-then-else gates",
         stats.eliminate_ite,
         percent (stats.eliminate_ite, stats.eliminate_gates));
    PRT ("  elimxors:      %15" PRId64 "   %10.2f %%  xor gates",
         stats.eliminate_xor,
         percent (stats.eliminate_xor, stats.eliminate_gates));
    PRT ("  elimdefs:      %15" PRId64 "   %10.2f %%  definitions",
         stats.eliminate_defs_extracted,
         percent (stats.eliminate_defs_extracted, stats.eliminate_gates));
    PRT ("  elimsubst:     %15" PRId64 "   %10.2f %%  substituted",
         stats.eliminate_substituted,
         percent (stats.eliminate_substituted, stats.vars_all_eliminated));
    PRT (
        "  elimsubstequi: %15" PRId64 "   %10.2f %%  equivalence gates",
        stats.eliminate_equivalence,
        percent (stats.eliminate_equivalence, stats.eliminate_substituted));
    PRT ("  elimsubstands: %15" PRId64 "   %10.2f %%  and gates",
         stats.eliminated_and,
         percent (stats.eliminated_and, stats.eliminate_substituted));
    PRT ("  elimsubstites: %15" PRId64 "   %10.2f %%  if-then-else gates",
         stats.eliminated_ite,
         percent (stats.eliminated_ite, stats.eliminate_substituted));
    PRT ("  elimsubstxors: %15" PRId64 "   %10.2f %%  xor gates",
         stats.eliminated_xor,
         percent (stats.eliminated_xor, stats.eliminate_substituted));
    PRT ("  elimsubstdefs: %15" PRId64 "   %10.2f %%  definitions",
         stats.eliminated_defs,
         percent (stats.eliminated_defs, stats.eliminate_substituted));
    PRT ("  elimres:       %15" PRId64 "   %10.2f    per eliminated",
         stats.eliminate_resolved,
         relative (stats.eliminate_resolved, stats.vars_all_eliminated));
    PRT ("  elimrestried:  %15" PRId64 "   %10.2f %%  per resolution",
         stats.eliminate_resolve_tried,
         percent (stats.eliminate_resolve_tried, stats.eliminate_resolved));
    PRT ("  def checked:   %15" PRId64 "   %10.2f    per phase",
         stats.eliminate_defs_checked,
         relative (stats.eliminate_defs_checked, stats.eliminations));
    PRT ("  def extracted: %15" PRId64 "   %10.2f %%  per checked",
         stats.eliminate_defs_extracted,
         percent (stats.eliminate_defs_extracted,
                  stats.eliminate_defs_checked));
    PRT ("  def units:     %15" PRId64 "   %10.2f %%  per checked",
         stats.eliminate_defs_unit,
         percent (stats.eliminate_defs_unit, stats.eliminate_defs_checked));
  }
  if (all || stats.propagator_cb) {
    PRT ("ext.prop. calls: %15" PRId64 "   %10.2f %%  of queries",
         stats.propagator_cb_propagate,
         percent (stats.propagator_cb_propagate, stats.propagator_cb));
    PRT ("  propagating:   %15" PRId64 "   %10.2f %%  per eprop-call",
         stats.propagator_cb_propagate_assign,
         percent (stats.propagator_cb_propagate_assign,
                  stats.propagator_cb_propagate));
    PRT ("  explained:     %15" PRId64 "   %10.2f %%  per eprop-call",
         stats.propagator_cb_propagate_explain,
         percent (stats.propagator_cb_propagate_explain,
                  stats.propagator_cb_propagate));
    PRT ("  falsified:     %15" PRId64 "   %10.2f %%  per eprop-call",
         stats.propagator_cb_propagate_clash,
         percent (stats.propagator_cb_propagate_clash,
                  stats.propagator_cb_propagate_clash));
    PRT ("ext.clause calls:%15" PRId64 "   %10.2f %%  of queries",
         stats.propagator_cb_add,
         percent (stats.propagator_cb_add, stats.propagator_cb));
    PRT ("  learned:       %15" PRId64 "   %10.2f %%  per called",
         stats.propagator_learned,
         percent (stats.propagator_learned, stats.propagator_cb_add));
    PRT ("  conflicting:   %15" PRId64 "   %10.2f %%  per learned",
         stats.propagator_learned_conflict,
         percent (stats.propagator_learned_conflict,
                  stats.propagator_learned));
    PRT ("  propagating:   %15" PRId64 "   %10.2f %%  per learned",
         stats.propagator_learned_propagating,
         percent (stats.propagator_learned_propagating,
                  stats.propagator_learned));
    PRT ("ext.final check: %15" PRId64 "   %10.2f %%  of queries",
         stats.propagator_cb_check_model,
         percent (stats.propagator_cb_check_model, stats.propagator_cb));
  }
  if (all || stats.factored) {
    PRT ("factored:        %15" PRId64 "   %10.2f %%  of variables",
         stats.factored, percent (stats.factored, internal->max_var));
    PRT ("  ands:          %15" PRId64 "   %10.2f %%  of factored",
         stats.factored_and, percent (stats.factored_and, stats.factored));
    PRT ("  xors:          %15" PRId64 "   %10.2f %%  of factored",
         stats.factored_xor, percent (stats.factored_xor, stats.factored));
    PRT ("  ites:          %15" PRId64 "   %10.2f %%  of factored",
         stats.factored_ite, percent (stats.factored_ite, stats.factored));
    PRT ("  eliminated:    %15" PRId64 "   %10.2f %%  of factored",
         stats.factored_eliminated,
         percent (stats.factored_eliminated, stats.factored));
    PRT ("  factor:        %15" PRId64 "   %10.2f    conflict interval",
         stats.factorings, relative (stats.conflicts, stats.factorings));
    PRT ("  cls factored:  %15" PRId64 "   %10.2f    per factored",
         stats.factor_added_clauses,
         relative (stats.factor_added_clauses, factored));
    PRT ("  lits factored: %15" PRId64 "   %10.2f    per factored",
         stats.literals_factored,
         relative (stats.literals_factored, factored));
    PRT ("  cls unfactored:%15" PRId64 "   %10.2f    per factored",
         stats.clauses_unfactored,
         relative (stats.clauses_unfactored, factored));
    PRT ("  cls redundant: %15" PRId64 "   %10.2f %%  of unfactored",
         stats.clauses_unfactored_redundant,
         percent (stats.clauses_unfactored_redundant,
                  stats.clauses_unfactored));
    PRT ("  lits unfactored:%14" PRId64 "   %10.2f    per factored",
         stats.literals_unfactored,
         relative (stats.literals_unfactored, factored));
  }
  if (all || stats.vars_all_fixed) {
    PRT ("fixed:           %15" PRId64 "   %10.2f %%  of all variables",
         stats.vars_all_fixed, percent (stats.vars_all_fixed, stats.vars));
    PRT ("  failed:        %15" PRId64 "   %10.2f %%  of all variables",
         stats.failed, percent (stats.failed, stats.vars));
    PRT ("  probefailed:   %15" PRId64 "   %10.2f %%  per failed",
         stats.probefailed, percent (stats.probefailed, stats.failed));
    PRT ("  transredunits: %15" PRId64 "   %10.2f %%  per failed",
         stats.transredunits, percent (stats.transredunits, stats.failed));
    PRT ("  inprobephases: %15" PRId64 "   %10.2f    interval",
         stats.inprobingphases,
         relative (stats.conflicts, stats.inprobingphases));
    PRT ("  inprobesuccess:%15" PRId64 "   %10.2f %%  phases",
         stats.inprobesuccess,
         percent (stats.inprobesuccess, stats.inprobingphases));
    PRT ("  probingrounds: %15" PRId64 "   %10.2f    per phase",
         stats.probingrounds,
         relative (stats.probingrounds, stats.inprobingphases));
    PRT ("  probed:        %15" PRId64 "   %10.2f    per failed",
         stats.probed, relative (stats.probed, stats.failed));
    PRT ("  hbrs:          %15" PRId64 "   %10.2f    per probed",
         stats.hbrs, relative (stats.hbrs, stats.probed));
    PRT ("  hbrsizes:      %15" PRId64 "   %10.2f    per hbr",
         stats.hbrsizes, relative (stats.hbrsizes, stats.hbrs));
    PRT ("  hbreds:        %15" PRId64 "   %10.2f %%  per hbr",
         stats.hbreds, percent (stats.hbreds, stats.hbrs));
    PRT ("  hbrsubs:       %15" PRId64 "   %10.2f %%  per hbr",
         stats.hbrsubs, percent (stats.hbrsubs, stats.hbrs));
  }
  PRT ("  units:         %15" PRId64 "   %10.2f    interval", stats.learned_units,
       relative (stats.conflicts, stats.learned_units));
  PRT ("  binaries:      %15" PRId64 "   %10.2f    interval",
       stats.learned_binaries, relative (stats.conflicts, stats.learned_binaries));
  if (all || stats.flush.learned) {
    PRT ("flushed:         %15" PRId64 "   %10.2f %%  per conflict",
         stats.flush.learned,
         percent (stats.flush.learned, stats.conflicts));
    PRT ("  hyper:         %15" PRId64 "   %10.2f %%  per conflict",
         stats.flush.hyper, relative (stats.flush.hyper, stats.conflicts));
    PRT ("  flushings:     %15" PRId64 "   %10.2f    interval",
         stats.flush.count, relative (stats.conflicts, stats.flush.count));
  }
  if (all || stats.instantiated) {
    PRT ("instantiated:    %15" PRId64 "   %10.2f %%  of tried",
         stats.instantiated, percent (stats.instantiated, stats.instried));
    PRT ("  instrounds:    %15" PRId64 "   %10.2f %%  of elimrounds",
         stats.instrounds, percent (stats.instrounds, stats.elimrounds));
  }
  if (all || stats.conflicts) {
    PRT ("learned:         %15" PRId64 "   %10.2f %%  per conflict",
         stats.learned_clauses,
         percent (stats.learned_clauses, stats.conflicts));
    PRT ("  bumped:        %15" PRId64 "   %10.2f    per learned",
         stats.vars_bumped, relative (stats.vars_bumped, stats.learned_clauses));
    PRT ("  recomputed:    %15" PRId64 "   %10.2f %%  per learned",
         stats.recomputed,
         percent (stats.recomputed, stats.learned_clauses));
    PRT ("  promoted1:     %15" PRId64 "   %10.2f %%  per learned",
         stats.promoted1, percent (stats.promoted1, stats.learned_clauses));
    PRT ("  promoted2:     %15" PRId64 "   %10.2f %%  per learned",
         stats.promoted2, percent (stats.promoted2, stats.learned_clauses));
    PRT ("  improvedglue:  %15" PRId64 "   %10.2f %%  per learned",
         stats.improvedglue,
         percent (stats.improvedglue, stats.learned_clauses));
  }
  if (all || stats.lucky.succeeded) {
    PRT ("lucky:           %15" PRId64 "   %10.2f %%  of tried",
         stats.lucky.succeeded,
         percent (stats.lucky.succeeded, stats.lucky.tried));
    PRT ("  constantzero   %15" PRId64 "   %10.2f %%  of tried",
         stats.lucky.constant.zero,
         percent (stats.lucky.constant.zero, stats.lucky.tried));
    PRT ("  constantone    %15" PRId64 "   %10.2f %%  of tried",
         stats.lucky.constant.one,
         percent (stats.lucky.constant.one, stats.lucky.tried));
    PRT ("  backwardone    %15" PRId64 "   %10.2f %%  of tried",
         stats.lucky.backward.one,
         percent (stats.lucky.backward.one, stats.lucky.tried));
    PRT ("  backwardzero   %15" PRId64 "   %10.2f %%  of tried",
         stats.lucky.backward.zero,
         percent (stats.lucky.backward.zero, stats.lucky.tried));
    PRT ("  forwardone     %15" PRId64 "   %10.2f %%  of tried",
         stats.lucky.forward.one,
         percent (stats.lucky.forward.one, stats.lucky.tried));
    PRT ("  forwardzero    %15" PRId64 "   %10.2f %%  of tried",
         stats.lucky.forward.zero,
         percent (stats.lucky.forward.zero, stats.lucky.tried));
    PRT ("  positivehorn   %15" PRId64 "   %10.2f %%  of tried",
         stats.lucky.horn.positive,
         percent (stats.lucky.horn.positive, stats.lucky.tried));
    PRT ("  negativehorn   %15" PRId64 "   %10.2f %%  of tried",
         stats.lucky.horn.negative,
         percent (stats.lucky.horn.negative, stats.lucky.tried));
    PRT ("  units:         %15" PRId64 "   %10.2f  of tried",
         stats.lucky.units,
         relative (stats.lucky.units, stats.lucky.tried));
  }
  PRT ("  extendbytes:   %15zd   %10.2f    bytes and MB", extendbytes,
       extendbytes / (double) (1l << 20));
  if (all || stats.learned_clauses)
    PRT ("learned_lits:    %15" PRId64 "   %10.2f %%  learned literals",
         stats.learned_literals,
         percent (stats.learned_literals, stats.learned_literals));
  PRT ("minimized:       %15" PRId64 "   %10.2f %%  learned literals",
       stats.minimized, percent (stats.minimized, stats.learned_literals));
  PRT ("shrunken:        %15" PRId64 "   %10.2f %%  learned literals",
       stats.shrunken, percent (stats.shrunken, stats.learned_literals));
  PRT ("minishrunken:    %15" PRId64 "   %10.2f %%  learned literals",
       stats.minishrunken,
       percent (stats.minishrunken, stats.learned_literals));

  if (all || stats.conflicts) {
    PRT ("otfs:            %15" PRId64 "   %10.2f %%  of conflict",
         stats.otfs_subsumed + stats.otfs.strengthened,
         percent (stats.otfs_subsumed + stats.otfs.strengthened,
                  stats.conflicts));
    PRT ("  subsumed       %15" PRId64 "   %10.2f %%  of conflict",
         stats.otfs_subsumed,
         percent (stats.otfs_subsumed, stats.conflicts));
    PRT ("  strengthened   %15" PRId64 "   %10.2f %%  of conflict",
         stats.otfs.strengthened,
         percent (stats.otfs.strengthened, stats.conflicts));
  }

  PRT ("propagations:    %15" PRId64 "   %10.2f M  per second",
       propagations, relative (propagations / 1e6, t));
  PRT ("  coverprops:    %15" PRId64 "   %10.2f %%  of propagations",
       stats.propagations_cover,
       percent (stats.propagations_cover, propagations));
  PRT ("  probeprops:    %15" PRId64 "   %10.2f %%  of propagations",
       stats.propagations_probe,
       percent (stats.propagations_probe, propagations));
  PRT ("  searchprops:   %15" PRId64 "   %10.2f %%  of propagations",
       stats.propagations_search,
       percent (stats.propagations_search, propagations));
  PRT ("  transredprops: %15" PRId64 "   %10.2f %%  of propagations",
       stats.propagations_transred,
       percent (stats.propagations_transred, propagations));
  PRT ("  vivifyprops:   %15" PRId64 "   %10.2f %%  of propagations",
       stats.propagations_vivify,
       percent (stats.propagations_vivify, propagations));
  if (all || stats.reactivated) {
    PRT ("reactivated:     %15" PRId64 "   %10.2f %%  of all variables",
         stats.reactivated, percent (stats.reactivated, stats.vars));
  }
  if (all || stats.reduced) {
    PRT ("reduced:         %15" PRId64 "   %10.2f %%  per conflict",
         stats.reduced, percent (stats.reduced, stats.conflicts));
    PRT ("  reductions:    %15" PRId64 "   %10.2f    interval",
         stats.reductions, relative (stats.conflicts, stats.reductions));
    PRT ("  sqrt scheme:   %15" PRId64 "   %10.2f %%  reductions",
         stats.reduced_sqrt,
         relative (stats.reduced_sqrt, stats.reductions));
    PRT ("  prct scheme:   %15" PRId64 "   %10.2f %%  reductions",
         stats.reduced_prct,
         relative (stats.reduced_prct, stats.reductions));
    PRT ("  collections:   %15" PRId64 "   %10.2f    interval",
         stats.collections, relative (stats.conflicts, stats.collections));
  }
  if (all || stats.rephased.total) {
    PRT ("rephased:        %15" PRId64 "   %10.2f    interval",
         stats.rephased.total,
         relative (stats.conflicts, stats.rephased.total));
    PRT ("  rephasedbest:  %15" PRId64 "   %10.2f %%  rephased best",
         stats.rephased.best,
         percent (stats.rephased.best, stats.rephased.total));
    PRT ("  rephasedflip:  %15" PRId64 "   %10.2f %%  rephased flipping",
         stats.rephased.flipped,
         percent (stats.rephased.flipped, stats.rephased.total));
    PRT ("  rephasedinv:   %15" PRId64 "   %10.2f %%  rephased inverted",
         stats.rephased.inverted,
         percent (stats.rephased.inverted, stats.rephased.total));
    PRT ("  rephasedorig:  %15" PRId64 "   %10.2f %%  rephased original",
         stats.rephased.original,
         percent (stats.rephased.original, stats.rephased.total));
    PRT ("  rephasedrand:  %15" PRId64 "   %10.2f %%  rephased random",
         stats.rephased.random,
         percent (stats.rephased.random, stats.rephased.total));
    PRT ("  rephasedwalk:  %15" PRId64 "   %10.2f %%  rephased walk",
         stats.rephased.walk,
         percent (stats.rephased.walk, stats.rephased.total));
  }
  if (all)
    PRT ("rescored:        %15" PRId64 "   %10.2f    interval",
         stats.scores_rescored, relative (stats.conflicts, stats.scores_rescored));
  if (all || stats.restarts) {
    PRT ("restarts:        %15" PRId64 "   %10.2f    interval",
         stats.restarts, relative (stats.conflicts, stats.restarts));
    PRT ("  reused:        %15" PRId64 "   %10.2f %%  per restart",
         stats.reused, percent (stats.reused, stats.restarts));
    PRT ("  reusedlevels:  %15" PRId64 "   %10.2f %%  per restart levels",
         stats.reusedlevels,
         percent (stats.reusedlevels, stats.restartlevels));
  }
  if (all || stats.restored) {
    PRT ("restored:        %15" PRId64 "   %10.2f %%  per weakened",
         stats.restored, percent (stats.restored, stats.weakened));
    PRT ("  restorations:  %15" PRId64 "   %10.2f %%  per extension",
         stats.restorations,
         percent (stats.restorations, stats.extensions));
    PRT ("  literals:      %15" PRId64 "   %10.2f    per restored clause",
         stats.restoredlits, relative (stats.restoredlits, stats.restored));
  }
  if (all || stats.stabphases) {
    PRT ("stabilizing:     %15" PRId64 "   %10.2f %%  of conflicts",
         stats.stabphases, percent (stats.stabconflicts, stats.conflicts));
    PRT ("  restartstab:   %15" PRId64 "   %10.2f %%  of all restarts",
         stats.restartstable,
         percent (stats.restartstable, stats.restarts));
    PRT ("  reusedstab:    %15" PRId64 "   %10.2f %%  per stable restarts",
         stats.reusedstable,
         percent (stats.reusedstable, stats.restartstable));
  }
  if (all || stats.vars_all_substituted) {
    PRT ("substituted:     %15" PRId64 "   %10.2f %%  of all variables",
         stats.vars_all_substituted,
         percent (stats.vars_all_substituted, stats.vars));
    PRT ("  decompositions:%15" PRId64 "   %10.2f    per phase",
         stats.decompositions,
         relative (stats.decompositions, stats.inprobingphases));
  }
  if (all || stats.sweep_equivalences) {
    PRT ("sweep equivs:    %15" PRId64 "   %10.2f %%  of swept variables",
         stats.sweep_equivalences,
         percent (stats.sweep_equivalences, stats.sweep_variables));
    PRT ("  sweepings:     %15" PRId64 "   %10.2f    vars per sweeping",
         stats.sweep, relative (stats.sweep_variables, stats.sweep));
    PRT ("  swept vars:    %15" PRId64 "   %10.2f %%  of all variables",
         stats.sweep_variables,
         percent (stats.sweep_variables, stats.vars));
    PRT ("  sweep units:   %15" PRId64 "   %10.2f %%  of all variables",
         stats.sweep_units, percent (stats.sweep_units, stats.vars));
    PRT ("  solved:        %15" PRId64 "   %10.2f    per swept variable",
         stats.sweep_solved,
         relative (stats.sweep_solved, stats.sweep_variables));
    PRT ("  sat:           %15" PRId64 "   %10.2f %%  solved",
         stats.sweep_sat, percent (stats.sweep_sat, stats.sweep_solved));
    PRT ("  unsat:         %15" PRId64 "   %10.2f %%  solved",
         stats.sweep_unsat,
         percent (stats.sweep_unsat, stats.sweep_solved));
    PRT ("  backbone solved:%14" PRId64 "   %10.2f %%  solved",
         stats.sweep_solved_backbone,
         percent (stats.sweep_solved_backbone, stats.sweep_solved));
    PRT ("    sat:         %15" PRId64 "   %10.2f %%  backbone solved",
         stats.sweep_sat_backbone,
         percent (stats.sweep_sat_backbone, stats.sweep_solved_backbone));
    PRT ("    unsat:       %15" PRId64 "   %10.2f %%  backbone solved",
         stats.sweep_unsat_backbone,
         percent (stats.sweep_unsat_backbone, stats.sweep_solved_backbone));
    PRT ("    unknown:     %15" PRId64 "   %10.2f %%  backbone solved",
         stats.sweep_unknown_backbone,
         percent (stats.sweep_unknown_backbone,
                  stats.sweep_solved_backbone));
    PRT ("    fixed:       %15" PRId64 "   %10.2f   per swept variable",
         stats.sweep_fixed_backbone,
         relative (stats.sweep_fixed_backbone, stats.sweep_variables));
    PRT ("    flip:        %15" PRId64 "   %10.2f   per swept variable",
         stats.sweep_flip_backbone,
         relative (stats.sweep_flip_backbone, stats.sweep_variables));
    PRT ("    flipped:     %15" PRId64 "   %10.2f %%  of backbone flip",
         stats.sweep_flipped_backbone,
         percent (stats.sweep_flipped_backbone, stats.sweep_flip_backbone));
    PRT ("  equiv solved:  %15" PRId64 "   %10.2f %%  solved",
         stats.sweep_solved_equivalences,
         percent (stats.sweep_solved_equivalences, stats.sweep_solved));
    PRT ("    sat:         %15" PRId64 "   %10.2f %%  equiv solved",
         stats.sweep_sat_equivalences,
         percent (stats.sweep_sat_equivalences,
                  stats.sweep_solved_equivalences));
    PRT ("    unsat:       %15" PRId64 "   %10.2f %%  equiv solved",
         stats.sweep_unsat_equivalences,
         percent (stats.sweep_unsat_equivalences,
                  stats.sweep_solved_equivalences));
    PRT ("    unknown:     %15" PRId64 "   %10.2f %%  equiv solved",
         stats.sweep_unknown_equivalences,
         percent (stats.sweep_unknown_equivalences,
                  stats.sweep_solved_equivalences));
    PRT ("    flip:        %15" PRId64 "   %10.2f    per swept variable",
         stats.sweep_flip_equivalences,
         relative (stats.sweep_flip_equivalences, stats.sweep_variables));
    PRT ("    flipped:     %15" PRId64 "   %10.2f %%  of equiv flip",
         stats.sweep_flipped_equivalences,
         percent (stats.sweep_flipped_equivalences,
                  stats.sweep_flip_equivalences));
    PRT ("  depth:         %15" PRId64 "   %10.2f    per swept variable",
         stats.sweep_depth,
         relative (stats.sweep_depth, stats.sweep_variables));
    PRT ("  environment:   %15" PRId64 "   %10.2f    per swept variable",
         stats.sweep_environment,
         relative (stats.sweep_environment, stats.sweep_variables));
    PRT ("  clauses:       %15" PRId64 "   %10.2f    per swept variable",
         stats.sweep_clauses,
         relative (stats.sweep_clauses, stats.sweep_variables));
    PRT ("  completed:     %15" PRId64 "   %10.2f    sweeps to complete",
         stats.sweep_completed,
         relative (stats.sweep, stats.sweep_completed));
  }
  if (all || stats.subsumed) {
    PRT ("subsumed:        %15" PRId64 "   %10.2f %%  of all clauses",
         stats.subsumed, percent (stats.subsumed, stats.added.total));
    PRT ("  subsumephases: %15" PRId64 "   %10.2f    interval",
         stats.subsumephases,
         relative (stats.conflicts, stats.subsumephases));
    PRT ("  subsumerounds: %15" PRId64 "   %10.2f    per phase",
         stats.subsumerounds,
         relative (stats.subsumerounds, stats.subsumephases));
    PRT ("  deduplicated:  %15" PRId64 "   %10.2f %%  per subsumed",
         stats.deduplicated, percent (stats.deduplicated, stats.subsumed));
    PRT ("  transreds:     %15" PRId64 "   %10.2f    interval",
         stats.transreds, relative (stats.conflicts, stats.transreds));
    PRT ("  transitive:    %15" PRId64 "   %10.2f %%  per subsumed",
         stats.transitive, percent (stats.transitive, stats.subsumed));
    PRT ("  subirr:        %15" PRId64 "   %10.2f %%  of subsumed",
         stats.subsumed_irredundant, percent (stats.subsumed_irredundant, stats.subsumed));
    PRT ("  subred:        %15" PRId64 "   %10.2f %%  of subsumed",
         stats.subsumed_redundant, percent (stats.subsumed_redundant, stats.subsumed));
    PRT ("  subtried:      %15" PRId64 "   %10.2f    tried per subsumed",
         stats.subtried, relative (stats.subtried, stats.subsumed));
    PRT ("  subchecks:     %15" PRId64 "   %10.2f    per tried",
         stats.subchecks, relative (stats.subchecks, stats.subtried));
    PRT ("  subchecks2:    %15" PRId64 "   %10.2f %%  per subcheck",
         stats.subchecks2, percent (stats.subchecks2, stats.subchecks));
    PRT ("  elimotfsub:    %15" PRId64 "   %10.2f %%  of subsumed",
         stats.elimotfsub, percent (stats.elimotfsub, stats.subsumed));
    PRT ("  elimbwsub:     %15" PRId64 "   %10.2f %%  of subsumed",
         stats.elimbwsub, percent (stats.elimbwsub, stats.subsumed));
    PRT ("  eagersub:      %15" PRId64 "   %10.2f %%  of subsumed",
         stats.eager_subsumed, percent (stats.eager_subsumed, stats.subsumed));
    PRT ("  eagertried:    %15" PRId64 "   %10.2f    tried per eagersub",
         stats.eager_subsumtions, relative (stats.eager_subsumtions, stats.eager_subsumed));
  }
  if (all || stats.strengthened) {
    PRT ("strengthened:    %15" PRId64 "   %10.2f %%  of all clauses",
         stats.strengthened,
         percent (stats.strengthened, stats.added.total));
    PRT ("  elimotfstr:    %15" PRId64 "   %10.2f %%  of strengthened",
         stats.elimotfstr, percent (stats.elimotfstr, stats.strengthened));
    PRT ("  elimbwstr:     %15" PRId64 "   %10.2f %%  of strengthened",
         stats.elimbwstr, percent (stats.elimbwstr, stats.strengthened));
  }
  if (all || stats.htrs) {
    PRT ("ternary:         %15" PRId64 "   %10.2f %%  of resolved",
         stats.htrs, percent (stats.htrs, stats.ternres));
    PRT ("  phases:        %15" PRId64 "   %10.2f    interval",
         stats.ternary, relative (stats.conflicts, stats.ternary));
    PRT ("  htr3:          %15" PRId64
         "   %10.2f %%  ternary hyper ternres",
         stats.htrs3, percent (stats.htrs3, stats.htrs));
    PRT ("  htr2:          %15" PRId64 "   %10.2f %%  binary hyper ternres",
         stats.htrs2, percent (stats.htrs2, stats.htrs));
  }
  PRT ("ticks:           %15" PRId64 "   %10.2f    propagation", totalticks,
       relative (totalticks, stats.propagations_search));
  PRT (" searchticks:    %15" PRId64 "   %10.2f %%  totalticks",
       searchticks, percent (searchticks, totalticks));
  PRT ("   stableticks:  %15" PRId64 "   %10.2f %%  searchticks",
       stats.ticks_search[1], percent (stats.ticks_search[1], searchticks));
  PRT ("   unstableticks:%15" PRId64 "   %10.2f %%  searchticks",
       stats.ticks_search[0], percent (stats.ticks_search[0], searchticks));
  PRT (" inprobeticks:   %15" PRId64 "   %10.2f %%  totalticks",
       inprobeticks, percent (inprobeticks, totalticks));
  PRT ("   backboneticks:%15" PRId64 "   %10.2f %%  searchticks",
       stats.ticks_backbone, percent (stats.ticks_backbone, searchticks));
  PRT ("   factorticks:  %15" PRId64 "   %10.2f %%  searchticks",
       stats.ticks_factor, percent (stats.ticks_factor, searchticks));
  PRT ("   probeticks:   %15" PRId64 "   %10.2f %%  searchticks",
       stats.ticks_probe, percent (stats.ticks_probe, searchticks));
  PRT ("   sweepticks:   %15" PRId64 "   %10.2f %%  searchticks",
       stats.ticks_sweep, percent (stats.ticks_sweep, searchticks));
  PRT ("   ternaryticks: %15" PRId64 "   %10.2f %%  searchticks",
       stats.ticks_ternary, percent (stats.ticks_ternary, searchticks));
  PRT ("   vivifyticks:  %15" PRId64 "   %10.2f %%  searchticks",
       stats.ticks_vivify, percent (stats.ticks_vivify, searchticks));
  PRT ("   walkticks:    %15" PRId64 "   %10.2f %%  searchticks",
       stats.ticks_walk, percent (stats.ticks_walk, searchticks));
  PRT ("   walkflipticks:%15" PRId64 "   %10.2f %%  searchticks",
       stats.ticks_walkflip, percent (stats.ticks_walkflip, searchticks));
  PRT ("   walkflipbrk:  %15" PRId64 "   %10.2f %%  searchticks",
       stats.ticks_walkflipbroken,
       percent (stats.ticks_walkflipbroken, searchticks));
  PRT ("   walkflipWL:   %15" PRId64 "   %10.2f %%  searchticks",
       stats.ticks_walkflipWL,
       percent (stats.ticks_walkflipWL, searchticks));
  PRT ("   walkpickticks:%15" PRId64 "   %10.2f %%  searchticks",
       stats.ticks_walkpick, percent (stats.ticks_walkpick, searchticks));
  PRT ("   walkbreak:    %15" PRId64 "   %10.2f %%  searchticks",
       stats.ticks_walkbreak, percent (stats.ticks_walkbreak, searchticks));
  if (all) {
    PRT ("tier recomputed: %15" PRId64 "   %10.2f    interval",
         stats.tierecomputed,
         relative (stats.conflicts, stats.tierecomputed));
  }
  if (all || stats.ilbtriggers) {
    PRT ("trail reuses:    %15" PRId64 "   %10.2f %%  of incremental calls",
         stats.ilbsuccess, percent (stats.ilbsuccess, stats.ilbtriggers));
    PRT ("  levels:        %15" PRId64 "   %10.2f    per reuse",
         stats.levelsreused,
         relative (stats.levelsreused, stats.ilbsuccess));
    PRT ("  literals:      %15" PRId64 "   %10.2f    per reuse",
         stats.literalsreused,
         relative (stats.literalsreused, stats.ilbsuccess));
    PRT ("  assumptions:   %15" PRId64 "   %10.2f    per reuse",
         stats.ilb_reused_assumptions,
         relative (stats.ilb_reused_assumptions, stats.ilbsuccess));
  }
  if (all || vivified) {
    PRT ("vivified:        %15" PRId64 "   %10.2f %%  of all clauses",
         vivified, percent (vivified, stats.added.total));
    PRT ("  vivifications: %15" PRId64 "   %10.2f    interval",
         stats.vivifications,
         relative (stats.conflicts, stats.vivifications));
    PRT ("  vivifychecks:  %15" PRId64 "   %10.2f %%  per conflict",
         stats.vivifychecks, percent (stats.vivifychecks, stats.conflicts));
    const int64_t vivified = stats.vivifiedtier1 + stats.vivifiedtier2 +
                             stats.vivifiedtier3 + stats.vivifiedirred;
    PRT ("  vivified:      %15" PRId64 "   %10.2f %%  per check", vivified,
         percent (vivified, stats.vivifychecks));
    PRT ("  vified-irred:  %15" PRId64 "   %10.2f %%  per vivified",
         stats.vivifiedirred, percent (stats.vivifiedirred, vivified));
    PRT ("  vified-tier1:  %15" PRId64 "   %10.2f %%  per vivified",
         stats.vivifiedtier1, percent (stats.vivifiedtier1, vivified));
    PRT ("  vified-tier2:  %15" PRId64 "   %10.2f %%  per vivified",
         stats.vivifiedtier2, percent (stats.vivifiedtier2, vivified));
    PRT ("  vified-tier3:  %15" PRId64 "   %10.2f %%  per vivified",
         stats.vivifiedtier3, percent (stats.vivifiedtier3, vivified));
    PRT ("  vivifysched:   %15" PRId64 "   %10.2f %%  checks per scheduled",
         stats.vivifysched,
         percent (stats.vivifychecks, stats.vivifysched));
    PRT ("  vivifyunits:   %15" PRId64 "   %10.2f %%  per vivify check",
         stats.vivifyunits,
         percent (stats.vivifyunits, stats.vivifychecks));
    PRT ("  vivifyinst:    %15" PRId64 "   %10.2f %%  per vivify check",
         stats.vivifyinst, percent (stats.vivifyinst, stats.vivifychecks));
    PRT ("  vivifysubs:    %15" PRId64 "   %10.2f %%  per subsumed",
         stats.vivifysubs, percent (stats.vivifysubs, stats.subsumed));
    PRT ("  vivifyflushed: %15" PRId64 "   %10.2f %%  per subsumed",
         stats.vivifyflushed,
         percent (stats.vivifyflushed, stats.subsumed));
    PRT ("  vivifysubred:  %15" PRId64 "   %10.2f %%  per subs",
         stats.vivifysubred,
         percent (stats.vivifysubred, stats.vivifysubs));
    PRT ("  vivifysubirr:  %15" PRId64 "   %10.2f %%  per subs",
         stats.vivifysubirr,
         percent (stats.vivifysubirr, stats.vivifysubs));
    PRT ("  vivifystrs:    %15" PRId64 "   %10.2f %%  per strengthened",
         stats.vivifystrs, percent (stats.vivifystrs, stats.strengthened));
    PRT ("  vivifystrirr:  %15" PRId64 "   %10.2f %%  per vivifystrs",
         stats.vivifystrirr,
         percent (stats.vivifystrirr, stats.vivifystrs));
    PRT ("  vivifystred1:  %15" PRId64 "   %10.2f %%  per vivifystrs",
         stats.vivifystred1,
         percent (stats.vivifystred1, stats.vivifystrs));
    PRT ("  vivifystred2:  %15" PRId64 "   %10.2f %%  per viviyfstrs",
         stats.vivifystred2,
         percent (stats.vivifystred2, stats.vivifystrs));
    PRT ("  vivifystred3:  %15" PRId64 "   %10.2f %%  per vivifystrs",
         stats.vivifystred3,
         percent (stats.vivifystred3, stats.vivifystrs));
    PRT ("  vivifydemote:  %15" PRId64 "   %10.2f %%  per vivifystrs",
         stats.vivifydemote,
         percent (stats.vivifydemote, stats.vivifystrs));
    PRT ("  vivifydecs:    %15" PRId64 "   %10.2f    per checks",
         stats.vivifydecs, relative (stats.vivifydecs, stats.vivifychecks));
    PRT ("  vivifyreused:  %15" PRId64
         "   %10.2f %%  per non-reused decision",
         stats.vivifyreused,
         percent (stats.vivifyreused, stats.vivifydecs));
  }
  if (all || stats.walk.count) {
    PRT ("walked:          %15" PRId64 "   %10.2f    interval",
         stats.walk.count, relative (stats.conflicts, stats.walk.count));
    if (all || stats.warmup.count) {
      PRT ("  prop-warmup:   %15" PRId64 "   %10.2f    per warmup",
           stats.warmup.propagated,
           relative (stats.warmup.propagated, stats.warmup.count));
      PRT ("  dec-warmup:    %15" PRId64 "   %10.2f    per warmup",
           stats.warmup.decision,
           relative (stats.warmup.decision, stats.warmup.count));
      PRT ("  dummydec-w:    %15" PRId64 "   %10.2f    per warmup",
           stats.warmup.dummydecision,
           relative (stats.warmup.dummydecision, stats.warmup.count));
      PRT ("  conflicts:     %15" PRId64 "   %10.2f    per warmup",
           stats.warmup.conflicts,
           relative (stats.warmup.conflicts, stats.warmup.count));
      PRT ("  warmup:        %15" PRId64 "   %10.2f    per walk",
           stats.warmup.count,
           relative (stats.warmup.count, stats.walk.count));
    }
#ifndef QUIET
    if (internal->profiles.walk.value > 0)
      PRT ("  flips:         %15" PRId64 "   %10.2f M  per second",
           stats.walk.flips,
           relative (1e-6 * stats.walk.flips,
                     internal->profiles.walk.value));
    else
#endif
      PRT ("  flips:         %15" PRId64 "   %10.2f    per walk",
           stats.walk.flips, relative (stats.walk.flips, stats.walk.count));
    if (stats.walk.minimum < LONG_MAX)
      PRT ("  minimum:       %15" PRId64 "   %10.2f %%  clauses",
           (int64_t) stats.walk.minimum,
           percent (stats.walk.minimum, stats.added.irredundant));
    PRT ("  broken:        %15" PRId64 "   %10.2f    per flip",
         stats.walk.broken, relative (stats.walk.broken, stats.walk.flips));
    PRT ("  improved:      %15" PRId64 "   %10.2f    per walk",
         stats.walk.improved,
         relative (stats.walk.improved, stats.walk.count));
    PRT ("  wrv-flip:      %15" PRId64 "   %10.2f %% flip",
         stats.walk.weight_reducing_var,
         percent (stats.walk.weight_reducing_var, stats.walk.flips));
    PRT ("  side-flip:     %15" PRId64 "   %10.2f %% flip",
         stats.walk.sideways,
         percent (stats.walk.sideways, stats.walk.flips));
    PRT ("  wght-transfer: %15" PRId64 "   %10.2f %% flip",
         stats.walk.weight_transfer,
         percent (stats.walk.weight_transfer, stats.walk.flips));
  }
  if (all || stats.weakened) {
    PRT ("weakened:        %15" PRId64 "   %10.2f    average size",
         stats.weakened, relative (stats.weakenedlen, stats.weakened));
    PRT ("  extensions:    %15" PRId64 "   %10.2f    interval",
         stats.extensions, relative (stats.conflicts, stats.extensions));
    PRT ("  flipped:       %15" PRId64 "   %10.2f    per weakened",
         stats.extended, relative (stats.extended, stats.weakened));
  }

  if (all || stats.congruence.gates) {
    PRT ("congruence:      %15" PRId64 "   %10.2f    interval",
         stats.congruence.rounds,
         relative (stats.conflicts, stats.congruence.rounds));
    PRT ("   units:        %15" PRId64 "   %10.2f    per congruent",
         stats.congruence.units,
         relative (stats.congruence.units, stats.congruence.congruent));
    PRT ("   cong-and:     %15" PRId64 "   %10.2f    per found gates",
         stats.congruence.ands,
         relative (stats.congruence.ands, stats.congruence.gates));
    PRT ("   cong-ite:     %15" PRId64 "   %10.2f    per found gates",
         stats.congruence.ites,
         relative (stats.congruence.ites, stats.congruence.gates));
    PRT ("   cong-xor:     %15" PRId64 "   %10.2f    per found gates",
         stats.congruence.xors,
         relative (stats.congruence.xors, stats.congruence.gates));
    PRT ("   congruent:    %15" PRId64 "   %10.2f    per round",
         stats.congruence.congruent,
         relative (stats.congruence.rounds, stats.congruence.congruent));
    PRT ("   unaries:      %15" PRId64 "   %10.2f    per round",
         stats.congruence.unaries,
         relative (stats.congruence.rounds, stats.congruence.unaries));
    int64_t rewritten = stats.congruence.rewritten_ands +
                        stats.congruence.rewritten_xors +
                        stats.congruence.rewritten_ites;
    PRT ("   rewritten:    %15" PRId64 "   %10.2f    per round", rewritten,
         percent (rewritten, stats.congruence.rounds));
    PRT ("   rewri.-ands:  %15" PRId64 "   %10.2f    per rewritten",
         stats.congruence.rewritten_ands,
         percent (stats.congruence.rewritten_ands, rewritten));
    PRT ("   rewri.-xors:  %15" PRId64 "   %10.2f%%  per rewritten",
         stats.congruence.rewritten_xors,
         percent (stats.congruence.rewritten_xors, rewritten));
    PRT ("   rewri.-ites:  %15" PRId64 "   %10.2f%%  per rewritten",
         stats.congruence.rewritten_ites,
         percent (stats.congruence.rewritten_ites, rewritten));
    PRT ("   subsumed:     %15" PRId64 "   %10.2f%%  per round",
         stats.congruence.subsumed,
         relative (stats.congruence.subsumed, stats.congruence.rounds));
    PRT ("   dummy-ands:   %15" PRId64 "   %10.2f%%  per round",
         stats.congruence.congruent_dummy_ands,
         relative (stats.congruence.congruent_dummy_ands,
                   stats.congruence.rounds));
  }

  LINE ();
  MSG ("%sseconds are measured in %s time for solving%s",
       tout.magenta_code (), internal->opts.realtime ? "real" : "process",
       tout.normal_code ());

  SECTION ("glue usage");

  internal->print_tier_usage_statistics ();

#endif // ifndef QUIET
}

void Internal::print_resource_usage () {
#ifndef QUIET
  SECTION ("resources");
  uint64_t m = maximum_resident_set_size ();
  MSG ("total process time since initialization: %12.2f    seconds",
       internal->process_time ());
  MSG ("total real time since initialization:    %12.2f    seconds",
       internal->real_time ());
  MSG ("maximum resident set size of process:    %12.2f    MB",
       m / (double) (1l << 20));
#endif
}

/*------------------------------------------------------------------------*/

void Checker::print_stats () {

  if (!stats.added && !stats.deleted)
    return;

  SECTION ("checker statistics");

  MSG ("checks:          %15" PRId64 "", stats.checks);
  MSG ("assumptions:     %15" PRId64 "   %10.2f    per check",
       stats.assumptions, relative (stats.assumptions, stats.checks));
  MSG ("propagations:    %15" PRId64 "   %10.2f    per check",
       stats.propagations, relative (stats.propagations, stats.checks));
  MSG ("original:        %15" PRId64 "   %10.2f %%  of all clauses",
       stats.original, percent (stats.original, stats.added));
  MSG ("derived:         %15" PRId64 "   %10.2f %%  of all clauses",
       stats.derived, percent (stats.derived, stats.added));
  MSG ("deleted:         %15" PRId64 "   %10.2f %%  of all clauses",
       stats.deleted, percent (stats.deleted, stats.added));
  MSG ("insertions:      %15" PRId64 "   %10.2f %%  of all clauses",
       stats.insertions, percent (stats.insertions, stats.added));
  MSG ("collections:     %15" PRId64 "   %10.2f    deleted per collection",
       stats.collections, relative (stats.collections, stats.deleted));
  MSG ("collisions:      %15" PRId64 "   %10.2f    per search",
       stats.collisions, relative (stats.collisions, stats.searches));
  MSG ("searches:        %15" PRId64 "", stats.searches);
  MSG ("units:           %15" PRId64 "", stats.learned_units);
}

void LratChecker::print_stats () {

  if (!stats.added && !stats.deleted)
    return;

  SECTION ("lrat checker statistics");

  MSG ("checks:          %15" PRId64 "", stats.checks);
  MSG ("insertions:      %15" PRId64 "   %10.2f %%  of all clauses",
       stats.insertions, percent (stats.insertions, stats.added));
  MSG ("original:        %15" PRId64 "   %10.2f %%  of all clauses",
       stats.original, percent (stats.original, stats.added));
  MSG ("derived:         %15" PRId64 "   %10.2f %%  of all clauses",
       stats.derived, percent (stats.derived, stats.added));
  MSG ("rat:             %15" PRId64 "   %10.2f %%  of derived clauses",
       stats.rat, percent (stats.rat, stats.derived));
  MSG ("deleted:         %15" PRId64 "   %10.2f %%  of all clauses",
       stats.deleted, percent (stats.deleted, stats.added));
  MSG ("finalized:       %15" PRId64 "   %10.2f %%  of all clauses",
       stats.finalized, percent (stats.finalized, stats.added));
  MSG ("collections:     %15" PRId64 "   %10.2f    deleted per collection",
       stats.collections, relative (stats.collections, stats.deleted));
  MSG ("collisions:      %15" PRId64 "   %10.2f    per search",
       stats.collisions, relative (stats.collisions, stats.searches));
  MSG ("searches:        %15" PRId64 "", stats.searches);
}

} // namespace CaDiCaL
