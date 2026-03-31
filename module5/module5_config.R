# =============================================================================
# Module 5 — Configuration (v2: denser tables)
# =============================================================================
#
# Ground truth parameters for the boundary validation experiment.
#
# v2 changes: K reduced from 7 to 5, N increased to ~1800 per gender.
# This gives ~1800 unions / 25 cells = ~72 per cell average, vs the
# previous ~450 / 49 = ~9. Eliminates the sparse-table problem that
# dominated v1 results (zero cells, clogit separation, z-constant artifacts).
#
# =============================================================================

# --- Simulation parameters ---
N_MC         <- 100       # Monte Carlo repetitions
N_NULL_SIMS  <- 30        # Null simulation iterations per MC rep (Paper 2 uses 50)
J_CHOICE     <- 50        # Choice set size for clogit
SEED         <- 42

# --- Population ---
K <- 5  # number of ancestry groups

GROUP_NAMES <- c("Nordic", "European", "LatAm", "MENA", "HornAfr")

GROUP_SIZES <- c(2200, 1100, 700, 700, 300)  # per gender (total = 5000)
# Design rationale (v3: tripled from v2):
#   Nordic = large majority (44%)
#   European = second largest (22%)
#   LatAm = medium (14%)
#   MENA = medium (14%)
#   HornAfr = small (6%) — tests small-group behavior
#
# Expected cell counts at ~100% pairing: 5000 / 25 = 200 per cell average
#   Nordic-Nordic: ~2200 (very well-estimated)
#   MENA-MENA: ~300+ (solid)
#   HornAfr-HornAfr: ~100+ (good)
#   Nordic-HornAfr: ~30-80 (adequate — was 0-4 in v2)
#   HornAfr-LatAm: ~30-60 (adequate)

# --- Geography (1D for tractability) ---
GEO_CENTERS <- c(
  Nordic   = 0,
  European = 4,
  LatAm    = 10,
  MENA     = 16,
  HornAfr  = 18
)
GEO_SPREAD <- 4.0  # within-group SD — slightly wider to reduce geographic determinism

# --- Distance weight in utility ---
DELTA <- 0.12   # slightly lower than v1 — let preferences matter more vs geography

# --- Ground truth: K×K preference matrix alpha_gh ---
# Row = ego's group (men), Column = partner's group (women)
# Diagonal = endogamy strength
# Off-diagonal = pairwise cross-group affinity/aversion
# ASYMMETRIC: alpha_gh != alpha_hg

ALPHA <- matrix(c(
  #  Nordic  Euro   LatAm   MENA   HAfr
     0.8,    0.3,   -0.1,   -0.6,  -0.4,   # Nordic men: weak endogamy, likes European, dislikes MENA/HAfr
     0.2,    1.5,    0.0,   -0.2,  -0.1,   # European men: moderate endogamy, mild aversion to MENA
    -0.1,    0.0,    2.2,    0.1,   0.2,   # LatAm men: strong endogamy, slight affinity to MENA/HAfr
    -0.8,   -0.3,    0.1,    2.8,   0.5,   # MENA men: strong endogamy, strong aversion to Nordic, affinity to HAfr
    -0.5,   -0.1,    0.2,    0.3,   3.5    # HornAfr men: very strong endogamy, aversion to Nordic, affinity to MENA
), nrow = K, ncol = K, byrow = TRUE)

rownames(ALPHA) <- GROUP_NAMES
colnames(ALPHA) <- GROUP_NAMES

# Design notes:
#   - Endogamy gradient: Nordic(0.8) < European(1.5) < LatAm(2.2) < MENA(2.8) < HornAfr(3.5)
#   - Key asymmetry: Nordic→MENA = -0.6, MENA→Nordic = -0.8 (0.2 difference)
#     This is the main test of whether methods detect directional asymmetry
#   - Cross-block affinity: MENA↔HornAfr (+0.5/+0.3)
#   - Expected block structure: {Nordic, European}, {LatAm}, {MENA, HornAfr}
#     → should produce 2-3 clusters depending on cut
#   - Sign reversal test: LatAm→Nordic = -0.1, Nordic→LatAm = -0.1 (symmetric mild aversion)
#     vs MENA→Nordic = -0.8, Nordic→MENA = -0.6 (asymmetric strong aversion)

# --- Bilateral matching parameters ---
ACCEPTANCE_RATES <- c(0.25, 0.50, 0.75)  # core experimental conditions
N_ROUNDS         <- 150                    # more rounds — 5000 agents need more encounters

# --- Output ---
RESULTS_DIR <- if (exists("script_dir")) {
  file.path(script_dir, "results")
} else {
  file.path(getwd(), "results")
}
if (!dir.exists(RESULTS_DIR)) dir.create(RESULTS_DIR, recursive = TRUE)

cat("Module 5 config loaded (v2: denser tables).\n")
cat(sprintf("  K = %d groups, N = %d per gender, %d MC reps\n",
            K, sum(GROUP_SIZES), N_MC))
cat(sprintf("  Avg cell count: ~%.0f unions/cell\n", sum(GROUP_SIZES) / K^2))
cat(sprintf("  Acceptance rates: %s\n",
            paste(ACCEPTANCE_RATES, collapse = ", ")))
