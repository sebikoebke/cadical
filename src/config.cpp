#include "internal.hpp"

namespace CaDiCaL {

/*------------------------------------------------------------------------*/

struct NameVal {
  const char *name;
  int val;
};

/*------------------------------------------------------------------------*/

// These are dummy configurations, which require additional code.

static NameVal default_config[1]; // With '-pedantic' just '[]' or
static NameVal plain_config[1];   // '[0]' gave a warning.

/*------------------------------------------------------------------------*/

// Here we have the pre-defined default configurations.

static NameVal sat_config[] = {
    {"elimeffort", 10},
    {"stabilizeonly", 1},
    {"subsumeeffort", 60},
};

static NameVal unsat_config[] = {
    {"stabilize", 0},
    {"walk", 0},
};

/*------------------------------------------------------------------------*/

// Configurations for runs with diferent identification points 
// (after expansion, after repair, both, after PALSAT)

static NameVal palsatids1_config[] = {
    {"autarkies", 0},           {"walkpalsat", 1},
    {"walkpalsatautarky", 1},   {"walkpalsatautarkypure", 0},
    {"walkpalsatautarkytuc", 0}, {"walkpalsatiwtl", 0},
    {"walkpalsatpick", 0},
};

static NameVal palsatids2_config[] = {
    {"autarkies", 0},           {"walkpalsat", 1},
    {"walkpalsatautarky", 2},   {"walkpalsatautarkypure", 0},
    {"walkpalsatautarkytuc", 0}, {"walkpalsatiwtl", 0},
    {"walkpalsatpick", 0},
};

static NameVal palsatids3_config[] = {
    {"autarkies", 0},           {"walkpalsat", 1},
    {"walkpalsatautarky", 3},   {"walkpalsatautarkypure", 0},
    {"walkpalsatautarkytuc", 0}, {"walkpalsatiwtl", 0},
    {"walkpalsatpick", 0},
};

static NameVal palsatids4_config[] = {
    {"autarkies", 0},           {"walkpalsat", 1},
    {"walkpalsatautarky", 4},   {"walkpalsatautarkypure", 0},
    {"walkpalsatautarkytuc", 0}, {"walkpalsatiwtl", 0},
    {"walkpalsatpick", 0},
};

// Configurations for different optimization techniques fpr autarkies

static NameVal palsatbase_config[] = {
    {"autarkies", 0},           {"walkpalsat", 5},
    {"walkpalsatautarky", 4},   {"walkpalsatautarkypure", 0},
    {"walkpalsatautarkytuc", 0}, {"walkpalsatiwtl", 0},
    {"walkpalsatpick", 0},
};

static NameVal palsatiwtl_config[] = {
    {"autarkies", 0},           {"walkpalsat", 5},
    {"walkpalsatautarky", 4},   {"walkpalsatautarkypure", 0},
    {"walkpalsatautarkytuc", 0}, {"walkpalsatiwtl", 1},
    {"walkpalsatpick", 0},
};

static NameVal palsatpure_config[] = {
    {"autarkies", 0},           {"walkpalsat", 5},
    {"walkpalsatautarky", 4},   {"walkpalsatautarkypure", 1},
    {"walkpalsatautarkytuc", 0}, {"walkpalsatiwtl", 0},
    {"walkpalsatpick", 0},
};

static NameVal palsatpureiwtl_config[] = {
    {"autarkies", 0},           {"walkpalsat", 5},
    {"walkpalsatautarky", 4},   {"walkpalsatautarkypure", 1},
    {"walkpalsatautarkytuc", 0}, {"walkpalsatiwtl", 1},
    {"walkpalsatpick", 0},
};

// Configurations for different scores for ProbSAT to
// check influence of Local Search on autarkies

static NameVal palsatlsl_config[] = {
    {"autarkies", 0},           {"walkpalsat", 5},
    {"walkpalsatautarky", 4},   {"walkpalsatautarkypure", 1},
    {"walkpalsatautarkytuc", 0}, {"walkpalsatiwtl", 1},
    {"walkpalsatpick", 1},
};

static NameVal palsatmakebreak_config[] = {
    {"autarkies", 0},           {"walkpalsat", 5},
    {"walkpalsatautarky", 4},   {"walkpalsatautarkypure", 1},
    {"walkpalsatautarkytuc", 0}, {"walkpalsatiwtl", 1},
    {"walkpalsatpick", 2},
};

/*------------------------------------------------------------------------*/

#define CONFIGS \
\
  CONFIG (default, "set default advanced internal options") \
  CONFIG (plain, "disable all internal preprocessing options") \
  CONFIG (sat, "set internal options to target satisfiable instances") \
  CONFIG (unsat, "set internal options to target unsatisfiable instances") \
\
  CONFIG (palsatids1, "classic palsat, autarky check after repair") \
  CONFIG (palsatids2, "classic palsat, autarky check after expansion") \
  CONFIG (palsatids3, "classic palsat, autarky check after both") \
  CONFIG (palsatids4, "classic palsat, autarky check at the end") \
\
  CONFIG (palsatbase, "anti-stagnation palsat, autarky check at the end") \
  CONFIG (palsatiwtl, "as '--palsatbase' plus increased walk tick limit") \
  CONFIG (palsatpure, "as '--palsatbase' plus pure literal fixpoint") \
  CONFIG (palsatpureiwtl, "as '--palsatbase' plus both of them") \
\
  CONFIG (palsatlsl, "as '--palsatpureiwtl' but picking on base^lsl") \
  CONFIG (palsatmakebreak, "as '--palsatpureiwtl' but picking on make-break")

static const char *configs[] = {
#define CONFIG(N, D) #N,
    CONFIGS
#undef CONFIG
};

static size_t num_configs = sizeof configs / sizeof *configs;

/*------------------------------------------------------------------------*/

bool Config::has (const char *name) {
#define CONFIG(N, D) \
  if (!strcmp (name, #N)) \
    return true;
  CONFIGS
#undef CONFIG
  return false;
}

bool Config::set (Options &opts, const char *name) {
  if (!strcmp (name, "default")) {
    opts.reset_default_values ();
    return true;
  }
  if (!strcmp (name, "plain")) {
    opts.disable_preprocessing ();
    return true;
  }
#define CONFIG(N, D) \
  do { \
    if (strcmp (name, #N)) \
      break; \
    const NameVal *BEGIN = N##_config; \
    const NameVal *END = BEGIN + sizeof N##_config / sizeof (NameVal); \
    for (const NameVal *P = BEGIN; P != END; P++) { \
      assert (Options::has (P->name)); \
      opts.set (P->name, P->val); \
    } \
    return true; \
  } while (0);
  CONFIGS
#undef CONFIG
  return false;
}

/*------------------------------------------------------------------------*/

void Config::usage () {
#define CONFIG(N, D) printf ("  %-18s " D "\n", "--" #N);
  CONFIGS
#undef CONFIG
}

/*------------------------------------------------------------------------*/

const char **Config::begin () { return configs; }
const char **Config::end () { return &configs[num_configs]; }

} // namespace CaDiCaL
