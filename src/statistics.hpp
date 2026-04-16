#ifndef _statistics_h_INCLUDED
#define _statistics_h_INCLUDED

// clang-format off

//           NAME,    VERBOSITY (0-3), PRINTING
#define CADICAL_STATISTICS \
\
  STATISTIC (backbone_phases,                      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (backbone_probes,                      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (backbone_rounds,                      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (backbone_units,                       1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (backtracked,                          3, NO_SECONDARY,       "",                    "") \
  STATISTIC (blocked,                              1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (blocked_candidates,                   1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (blocked_pure,                         1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (blocked_pure_literals,                1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (blocked_resolutions,                  1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (blockings,                            1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (clauses_added_irredundant,            3, NO_SECONDARY,        0,                     0) \
  STATISTIC (clauses_added_redundant,              3, NO_SECONDARY,        0,                     0) \
  STATISTIC (clauses_added_total,                  3, NO_SECONDARY,        0,                     0) \
  STATISTIC (clauses_current_irredundant,          3, NO_SECONDARY,        0,                     0) \
  STATISTIC (clauses_current_redundant,            3, NO_SECONDARY,        0,                     0) \
  STATISTIC (clauses_current_total,                3, NO_SECONDARY,        0,                     0) \
  STATISTIC (clauses_improved_glue,                3, NO_SECONDARY,       "",                    "") \
  STATISTIC (clauses_promoted_tier1,               3, NO_SECONDARY,       "",                    "") \
  STATISTIC (clauses_promoted_tier2,               3, NO_SECONDARY,       "",                    "") \
  STATISTIC (clauses_recomputed_glue,              3, NO_SECONDARY,       "",                    "") \
  STATISTIC (collected,                            1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (collections,                          1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (compacts,                             1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (condition_active,                     2, PER_CONDITIONING,   "",    "per conditioning") \
  STATISTIC (condition_candidates,                 2, PER_CONDITIONING,   "",    "per conditioning") \
  STATISTIC (conditioned,                          1, PER_CONDITIONING,   "",    "per conditioning") \
  STATISTIC (condition_final_assigned_init,        2, PER_CONDITIONING,   "",    "per conditioning") \
  STATISTIC (condition_final_assigned,             2, PER_CONDITIONING,   "",    "per conditioning") \
  STATISTIC (condition_final_autarky,              2, PER_CONDITIONING,   "",    "per conditioning") \
  STATISTIC (condition_final_conditional,          2, PER_CONDITIONING,   "",    "per conditioning") \
  STATISTIC (conditionings,                        2, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (condition_init_assigned,              2, PER_CONDITIONING,   "",    "per conditioning") \
  STATISTIC (condition_init_autarky,               2, PER_CONDITIONING,   "",    "per conditioning") \
  STATISTIC (condition_init_conditional,           2, PER_CONDITIONING,   "",    "per conditioning") \
  STATISTIC (condition_propagated,                 2, PER_CONDITIONING,   "",    "per conditioning") \
  STATISTIC (conflicts,                            0, PER_SECOND,         "",          "per second") \
  STATISTIC (conflicts_chrono,                     3, NO_SECONDARY,       "",                    "") \
  STATISTIC (congruence_ands,                      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_dummy_ands,                1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_gates,                     1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_gates_and,                 1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_gates_ite,                 1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_gates_xor,                 1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_ites,                      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_rewritten_ands,            1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_rewritten_ites,            1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_rewritten_xors,            1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_rounds,                    1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_simplified,                1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_simplified_ands,           1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_simplified_ites,           1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_simplified_xors,           1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_subsumed,                  1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_trivial_ite,               1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_unary,                     1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_unary_and,                 1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_unary_ite,                 1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_units,                     1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruence_xors,                      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (congruent,                            1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (cover_asymmetric,                     2, PER_COVERING,       "",        "per covering") \
  STATISTIC (cover_blocked,                        2, PER_COVERING,       "",        "per covering") \
  STATISTIC (coverings,                            2, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (cover_total,                          2, PER_COVERING,       "",        "per covering") \
  STATISTIC (decisions,                            0, PER_SECOND,         "",          "per second") \
  STATISTIC (decisions_random,                     1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (decisions_random_phases,              1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (decisions_searched,                   0, PER_SECOND,         "",          "per second") \
  STATISTIC (decompositions,                       1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (deduplicated,                         1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (deduplicate_init,                     1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (deduplicate_init_rounds,              1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (deduplicate_failed,                   1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (deduplicate_hyper_unary,              1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (deduplications,                       1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eager_subsumed,                       1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eager_subsumtions,                    1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_and,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_complete,                   1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminated,                           1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminated_and,                       1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminated_defs,                      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_defs_checked,               1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_defs_extracted,             1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_defs_ticks,                 1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_defs_unit,                  1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminated_equivalence,               1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminated_ite,                       1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminated_xor,                       1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_equivalence,                1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_fast_phases,                1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_fast_rounds,                1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_gates,                      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_ite,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_otf_strengthened,           1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_otf_subsumed,               1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_phases,                     1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_resolved,                   1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_resolve_tried,              1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_strengthened_bw,            1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_substituted,                1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_subsumed_bw,                1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_tried,                      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminate_xor,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (eliminations,                         1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (extended,                             1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (extensions,                           1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (factor_added_clauses,                 1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (factor_added_literals,                1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (factored,                             1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (factored_and,                         1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (factored_eliminated,                  1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (factored_ite,                         1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (factored_xor,                         1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (factorings,                           1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (factor_removed_clauses,               1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (factor_removed_literals,              1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (factor_removed_redundant,             1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (failed_literals,                      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (flushed,                              1, PCNT_CONFLICTS,    "%",           "conflicts") \
  STATISTIC (flush_hyper,                          2, PCNT_CONFLICTS,    "%",           "conflicts") \
  STATISTIC (flushings,                            2, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (flush_learned,                        2, PCNT_CONFLICTS,    "%",           "conflicts") \
  STATISTIC (garbage_bytes,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (garbage_clauses,                      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (garbage_literals,                     1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (hbr_redundant,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (hbrs,                                 1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (hbr_sizes,                            1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (hbr_subsuming,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (ilb_reused_assumptions,               1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (ilb_reused_levels,                    1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (ilb_reused_literals,                  1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (ilb_success,                          1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (ilb_triggers,                         1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (incremental_decay,                    1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (inprobingphases,                      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (inprobingsuccess,                     1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (instantiated,                         1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (instantiate_tried,                    1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (instantiations,                       1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (irredundant_literals,                 1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (learned_binaries,                     1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (learned_clauses,                      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (learned_literals,                     1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (learned_units,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (lucky,                                1, PCNT_SEARCH,       "%",              "search") \
  STATISTIC (lucky_backward_one,                   2, CONFLICT_INTERVAL, "%",              "search") \
  STATISTIC (lucky_backward_zero,                  2, CONFLICT_INTERVAL, "%",              "search") \
  STATISTIC (lucky_constant_one,                   2, CONFLICT_INTERVAL, "%",              "search") \
  STATISTIC (lucky_constant_zero,                  2, CONFLICT_INTERVAL, "%",              "search") \
  STATISTIC (lucky_forward_one,                    2, CONFLICT_INTERVAL, "%",              "search") \
  STATISTIC (lucky_forward_zero,                   2, CONFLICT_INTERVAL, "%",              "search") \
  STATISTIC (lucky_horn_negative,                  2, CONFLICT_INTERVAL, "%",              "search") \
  STATISTIC (lucky_horn_positive,                  2, CONFLICT_INTERVAL, "%",              "search") \
  STATISTIC (lucky_random,                         2, CONFLICT_INTERVAL, "%",              "search") \
  STATISTIC (lucky_tried,                          2, CONFLICT_INTERVAL, "%",              "search") \
  STATISTIC (lucky_units,                          2, CONFLICT_INTERVAL, "%",              "search") \
  STATISTIC (mark_block,                           3, NO_SECONDARY,        0,                     0) \
  STATISTIC (mark_elim,                            3, NO_SECONDARY,        0,                     0) \
  STATISTIC (mark_factor,                          3, NO_SECONDARY,        0,                     0) \
  STATISTIC (mark_subsume,                         3, NO_SECONDARY,        0,                     0) \
  STATISTIC (mark_ternary,                         3, NO_SECONDARY,        0,                     0) \
  STATISTIC (minimized,                            1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (otfs_strengthened,                    1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (otfs_subsumed,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (preprocessings,                       1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (probed,                               1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (probe_failed_literals,                1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (probingrounds,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (propagations,                         0, PER_SECOND,         "",          "per second") \
  STATISTIC (propagations_backbone,                2, PCNT_PROPAGATIONS, "%",        "propagations") \
  STATISTIC (propagations_cover,                   2, PCNT_PROPAGATIONS, "%",        "propagations") \
  STATISTIC (propagations_instantiate,             2, PCNT_PROPAGATIONS, "%",        "propagations") \
  STATISTIC (propagations_probe,                   2, PCNT_PROPAGATIONS, "%",        "propagations") \
  STATISTIC (propagations_search,                  2, PCNT_PROPAGATIONS, "%",        "propagations") \
  STATISTIC (propagations_transred,                2, PCNT_PROPAGATIONS, "%",        "propagations") \
  STATISTIC (propagations_vivify,                  2, PCNT_PROPAGATIONS, "%",        "propagations") \
  STATISTIC (propagator_cb,                        1, NO_SECONDARY,        0,                     0) \
  STATISTIC (propagator_cb_add,                    1, PCNT_EXT_CB,       "%",           "callbacks") \
  STATISTIC (propagator_cb_check_model,            1, PCNT_EXT_CB,       "%",           "callbacks") \
  STATISTIC (propagator_cb_propagate,              1, PCNT_EXT_CB,       "%",           "callbacks") \
  STATISTIC (propagator_cb_propagate_assign,       2, PCNT_EXT_CB_PROP,  "%", "propagate callbacks") \
  STATISTIC (propagator_cb_propagate_clash,        2, PCNT_EXT_CB_PROP,  "%", "propagate callbacks") \
  STATISTIC (propagator_cb_propagate_explain,      2, PCNT_EXT_CB_PROP,  "%", "propagate callbacks") \
  STATISTIC (propagator_learned,                   1, NO_SECONDARY,        0,                     0) \
  STATISTIC (propagator_learned_conflict,          2, PCNT_EXT_LEARNED,  "%",  "propagator clauses") \
  STATISTIC (propagator_learned_elevating,         2, PCNT_EXT_LEARNED,  "%",  "propagator clauses") \
  STATISTIC (propagator_learned_propagating,       2, PCNT_EXT_LEARNED,  "%",  "propagator clauses") \
  STATISTIC (propagator_learned_unit,              2, PCNT_EXT_LEARNED,  "%",  "propagator clauses") \
  STATISTIC (recomputed_tiers,                     1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (reduced,                              1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (reductions,                           1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (rephased,                             1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (rephased_best,                        2, PCNT_REPHASED,     "%",            "rephased") \
  STATISTIC (rephased_flipped,                     2, PCNT_REPHASED,     "%",            "rephased") \
  STATISTIC (rephased_inverted,                    2, PCNT_REPHASED,     "%",            "rephased") \
  STATISTIC (rephased_original,                    2, PCNT_REPHASED,     "%",            "rephased") \
  STATISTIC (rephased_random,                      2, PCNT_REPHASED,     "%",            "rephased") \
  STATISTIC (rephased_walk,                        2, PCNT_REPHASED,     "%",            "rephased") \
  STATISTIC (restart,                              1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (restart_levels,                       3, CONFLICT_INTERVAL,  "",                    "") \
  STATISTIC (restart_stable,                       3, CONFLICT_INTERVAL,  "",                    "") \
  STATISTIC (restorations,                         1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (restored_clauses,                     1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (restored_literals,                    1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (reused,                               1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (reused_levels,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (reused_stable,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (scores_rescored,                      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (scores_shuffled,                      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (searches,                             2, NO_SECONDARY,        0,                     0) \
  STATISTIC (sections,                             3, NO_SECONDARY,       "",                    "") \
  STATISTIC (shrunken,                             1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (shrunken_minimize,                    1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (stable_conflicts,                     2, PCNT_CONFLICTS,    "%",           "conflicts") \
  STATISTIC (stable_phases_incremental,            2, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (stable_phases_total,                  1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (strengthened,                         1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (subsume_checks,                       1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (subsume_checks_binary,                1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (subsumed,                             1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (subsumed_irredundant,                 1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (subsumed_redundant,                   1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (subsume_phases,                       1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (subsume_rounds,                       1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (subsume_tried,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep,                                1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_backbone_fixed,                 1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_backbone_flip,                  1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_backbone_flipped,               1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_backbone_solved,                1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_backbone_solved_sat,            1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_backbone_solved_unknown,        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_backbone_solved_unsat,          1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_clauses,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_completed,                      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_depth,                          1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_environment,                    1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_equivalences,                   1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_equivalences_flip,              1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_equivalences_flipped,           1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_equivalences_solved,            1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_equivalences_solved_sat,        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_equivalences_solved_unknown,    1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_equivalences_solved_unsat,      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_solved,                         1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_solved_sat,                     1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_solved_unsat,                   1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_units,                          1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (sweep_variables,                      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (ternary,                              1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (ternary_htrs,                         1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (ternary_htrs_binary,                  1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (ternary_htrs_ternary,                 1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (ternary_resolutions,                  1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (ticks,                                0, PER_PROPAGATION,    "",         "propagation") \
  STATISTIC (ticks_backbone,                       2, PCNT_TICKS,        "%",               "ticks") \
  STATISTIC (ticks_factor,                         2, PCNT_TICKS,        "%",               "ticks") \
  STATISTIC (ticks_probe,                          2, PCNT_TICKS,        "%",               "ticks") \
  STATISTIC (ticks_search_stable,                  2, PCNT_TICKS,        "%",               "ticks") \
  STATISTIC (ticks_search_unstable,                2, PCNT_TICKS,        "%",               "ticks") \
  STATISTIC (ticks_sweep,                          2, PCNT_TICKS,        "%",               "ticks") \
  STATISTIC (ticks_ternary,                        2, PCNT_TICKS,        "%",               "ticks") \
  STATISTIC (ticks_vivify,                         2, PCNT_TICKS,        "%",               "ticks") \
  STATISTIC (ticks_walk,                           2, PCNT_TICKS,        "%",               "ticks") \
  STATISTIC (ticks_walk_break,                     2, PCNT_TICKS_WALK,   "%",          "walk ticks") \
  STATISTIC (ticks_walk_flip,                      2, PCNT_TICKS_WALK,   "%",          "walk ticks") \
  STATISTIC (ticks_walk_flip_broken,               2, PCNT_TICKS_WALK,   "%",          "walk ticks") \
  STATISTIC (ticks_walk_flip_wl,                   2, PCNT_TICKS_WALK,   "%",          "walk ticks") \
  STATISTIC (ticks_walk_pick,                      2, PCNT_TICKS_WALK,   "%",          "walk ticks") \
  STATISTIC (transitive_rounds,                    1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (transitive_clauses,                   1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (transitive_units,                     1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (variables_extension,                  1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (variables_original,                   1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vars,                                 2, NO_SECONDARY,        0,                     0) \
  STATISTIC (vars_active,                          1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vars_all_eliminated,                  1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vars_all_eliminated_fast,             1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vars_all_fixed,                       1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vars_all_pure,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vars_all_substituted,                 1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vars_bumped,                          3, NO_SECONDARY,       "",                    "") \
  STATISTIC (vars_declared,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vars_inactive,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vars_now_eliminated,                  1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vars_now_fixed,                       1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vars_now_pure,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vars_now_substituted,                 1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vars_reactivated,                     1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vars_unused,                          1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivifications,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivified,                             1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivified_irredundant,                 1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivified_redundant_tier1,             1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivified_redundant_tier2,             1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivified_redundant_tier3,             1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivify_checks,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivify_decisions,                     1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivify_demote,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivify_flushed,                       1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivify_implied,                       1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivify_instantiated,                  1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivify_reused,                        1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivify_scheduled,                     1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivify_strengthened,                  1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivify_strengthened_irredundant,      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivify_strengthened_redundant_tier1,  1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivify_strengthened_redundant_tier2,  1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivify_strengthened_redundant_tier3,  1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivify_subsumed,                      1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivify_subsumed_irredundant,          1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivify_subsumed_redundant,            1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (vivify_units,                         1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (walk,                                 1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (walk_broken,                          2, CONFLICT_INTERVAL,  "",                    "") \
  STATISTIC (walk_flips,                           2, CONFLICT_INTERVAL,  "",                    "") \
  STATISTIC (walk_improved,                        2, CONFLICT_INTERVAL,  "",                    "") \
  STATISTIC (walk_minimum,                         2, CONFLICT_INTERVAL,  "",                    "") \
  STATISTIC (walk_sideways,                        2, CONFLICT_INTERVAL,  "",                    "") \
  STATISTIC (walk_warmup,                          2, CONFLICT_INTERVAL,  "",                    "") \
  STATISTIC (walk_warmup_conflicts,                2, CONFLICT_INTERVAL,  "",                    "") \
  STATISTIC (walk_warmup_decision,                 2, CONFLICT_INTERVAL,  "",                    "") \
  STATISTIC (walk_warmup_decision_dummy,           2, CONFLICT_INTERVAL,  "",                    "") \
  STATISTIC (walk_warmup_propagated,               2, CONFLICT_INTERVAL,  "",                    "") \
  STATISTIC (walk_weight_reducing,                 2, CONFLICT_INTERVAL,  "",                    "") \
  STATISTIC (walk_weight_transfer,                 2, CONFLICT_INTERVAL,  "",                    "") \
  STATISTIC (weakened,                             1, CONFLICT_INTERVAL,  "",            "interval") \
  STATISTIC (weakened_lengths,                     1, CONFLICT_INTERVAL,  "",            "interval")

// clang-format on

#endif
