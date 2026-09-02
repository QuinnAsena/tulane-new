> ## SUPERSEDED — 2 September 2026
>
> This document reports an analysis that weighted observation error **per bin**, using the
> uncertainties in `ocfs_uncertainty.csv`. The current pipeline
> (`kalman-smoother/ocfs_smoother.r`) instead fixes observation error to a **single number**,
> pending confirmation of how those uncertainties were derived. That one change moves several
> results, so read the numbers below with the following in mind.
>
> **Superseded:**
>
> | | this document (per-bin `ME`) | current (constant `ME`) |
> |---|---|---|
> | ABRUPT (`heinrich`) importance | 0.812 | **0.391** |
> | HUMAN importance | 0.564 | 0.517 |
> | CLIM / TIME / FIRE importance | 0.969 / 0.917 / 0.233 | 0.992 / 0.877 / 0.279 |
> | omnibus *P*, all five covariates | 0.0054 | 0.0025 |
> | `Σb` (persistence) | 0.85 | 0.62 |
> | long-run multiplier `1/(1−Σb)` | 6.5 | 2.6 |
> | `se` as a fraction of sd(response) | 0.31 | 0.58 |
> | model ranking and weights | as tabulated below | see `results/ocfs/ocfs_model_ranking.csv` |
>
> The headline claim below that **Heinrich stadials carry the signal does not survive** the
> change: ABRUPT drops from clearly supported to below the 0.5 no-information line. CLIM and
> TIME remain supported and FIRE remains unsupported under both. Resolving the measurement-error
> question therefore decides the ABRUPT result specifically.
>
> **Still valid, and not affected:** the λ = 0.25 derivation (§ on the transform); `c` beating
> `d` by 10.3 AICc; no support for time-varying coefficients; the PrDens/time confound and the
> reason TIME must compete with HUMAN; individual lag coefficients being unidentified while
> their sum is; and the methodological finding that no individual coefficient survives model
> averaging at |z| > 2. The Kalman-smoother bug and its fix are unaffected.
>
> Kept for the record because the model-averaging, standard-error and sensitivity work
> documented here has not yet been redone in the current pipeline.

# OCFS drivers: results from the rebuilt TVARSS analysis

**19 August 2026.** Produced by `kalman-smoother/ocfs_analysis.r` (Tier 2) and
`kalman-smoother/ocfs_model_selection.r` (Tier 3). Companion to
`2026-08-12_tvarss-ocfs-diagnostic.md`, which explains why the original analysis had to be
rebuilt.

Response: `OCFS_conc^0.25`, 83 observations over 293 bins of 200 yr (59.4 ka to present).
Predictors enter through **c** (the process equation). AR(2), coefficients fixed in time.

---

## The short version

1. **The predictors are jointly supported.** Omnibus LRT `dev = 16.58, df = 5, P = 0.0054`;
   AICc 727.1 → 722.5. This holds across every choice tested.
2. **Climate and abrupt climate carry the signal.** Summed Akaike weights: CLIM 0.97 / 0.73,
   ABRUPT 0.81 / 0.92 (the two figures are the two TIME definitions).
3. **Fire does not.** `char_acc` has summed weight 0.23–0.29, below the 0.5 no-information
   line, and a model-averaged coefficient indistinguishable from zero.
4. **Human population density does not survive a temporal trend.** `PrDens` looks significant
   until a plain linear trend in time is allowed to compete; then it drops to *P* = 0.081 and
   summed weight 0.32–0.56. The single strongest term in the best model is **time itself**.
5. **No individual coefficient survives model averaging at |z| > 2.** With unconditional SEs
   that fold in selection uncertainty, the closest is `mean_co2` at z = −1.98. Report the
   model ranking and variable importance, not individual *P*-values.
6. **No support for time-varying coefficients.** Freeing `sb0` costs 5.1 AICc; `sb` costs 6.9;
   both cost 12.1.

---

## 1. The model set

Five a priori hypothesis blocks, 2⁵ = 32 models, each block in exactly 16 of them so that
summed weights are comparable between blocks (Burnham & Anderson 2002 §6.9.6).

| block | terms | why grouped this way |
|---|---|---|
| CLIM | `d18O` + `mean_co2` | r = 0.89; not separately identifiable, so tested as a pair |
| ABRUPT | `heinrich` | 0/1 stadial indicator |
| FIRE | `char_acc` | |
| HUMAN | `PrDens` | |
| TIME | late-record term | the competitor that makes HUMAN interpretable |

**TIME exists because `PrDens` is confounded with time.** It is exactly zero for bins 62–309
and nonzero only in the contiguous block 1–61 (1.1–13.1 ka); just 18 of the 83 OCFS
observations fall in its support. Without a competitor, any late-record change in OCFS would
be credited to human population. The set was therefore fitted twice:

- `time_lin` — a linear trend in time (standardised row index)
- `time_hol` — an indicator for bins 1–61, i.e. exactly `PrDens`'s own support

## 2. HUMAN versus TIME, head to head

| comparison | vs linear trend | vs Holocene indicator |
|---|---|---|
| adding HUMAN to TIME | dev 6.42, **P = 0.011** | dev 1.58, P = 0.209 |
| adding TIME to HUMAN | dev 2.67, P = 0.102 | dev 0.63, P = 0.428 |

Against a linear trend, `PrDens` adds something. Against the Holocene indicator, neither term
adds anything over the other — on this record they are close to the same variable. Note also
that TIME *alone* is worse than the null under `time_lin` (AICc 729.4 vs 727.1), so the
pattern is not a bare monotone trend either.

## 3. Model ranking

`time_lin` set — top 7 of 32, and the better-supported of the two sets:

| model | k | logLik | AICc | ΔAICc | w | evidence ratio |
|---|---|---|---|---|---|---|
| CLIM+ABRUPT+HUMAN+TIME | 9 | −346.8 | 714.1 | 0.00 | 0.315 | 1.0 |
| CLIM+ABRUPT+TIME | 8 | −348.3 | 714.5 | 0.49 | 0.246 | 1.3 |
| CLIM+HUMAN+TIME | 8 | −349.3 | 716.5 | 2.41 | 0.095 | 3.3 |
| CLIM+ABRUPT+FIRE+HUMAN+TIME | 10 | −346.8 | 716.6 | 2.56 | 0.088 | 3.6 |
| CLIM+ABRUPT+FIRE+TIME | 9 | −348.2 | 716.9 | 2.84 | 0.076 | 4.1 |
| CLIM+TIME | 7 | −351.3 | 718.0 | 3.97 | 0.043 | 7.3 |
| CLIM+ABRUPT | 7 | −351.6 | 718.6 | 4.57 | 0.032 | 9.8 |

Only **2 models fall within 2 AICc** and 6 within 4, so there is real selection uncertainty —
but every model in the top five contains CLIM, and four of five contain ABRUPT.

`time_hol` set — best model is CLIM+ABRUPT (AICc 718.6, w = 0.247), 3 models within 2 AICc.
Its best AICc is 4.5 worse than the `time_lin` set's, so a linear trend describes the temporal
structure better than a Holocene step.

## 4. Variable importance (summed Akaike weights)

| block | time_lin | time_hol | reading |
|---|---|---|---|
| CLIM | **0.969** | **0.729** | supported under both |
| ABRUPT | **0.812** | **0.922** | supported under both |
| TIME | 0.917 | 0.361 | depends entirely on how it is defined |
| HUMAN | 0.564 | 0.320 | ambiguous at best; never exceeds TIME |
| FIRE | 0.233 | 0.285 | not supported |

Interpretable as importance only because the set is balanced. ~0.5 is the no-information
value.

## 5. The best model in detail

`CLIM + ABRUPT + HUMAN + TIME`, AICc 714.0, dominant eigenvalue 0.792, `1 − Σb` = 0.247.

| term | c | Wald SE | drop-one dev | *P* |
|---|---|---|---|---|
| time | 0.329 | 0.132 | 9.36 | **0.0022** |
| mean_co2 | −0.486 | 0.218 | 6.96 | **0.0084** |
| heinrich | 0.466 | 0.246 | 4.96 | **0.026** |
| PrDens | −0.318 | 0.195 | 3.04 | 0.081 |
| d18O | 0.156 | 0.190 | 0.74 | 0.391 |

`c` is the effect on the latent series **per 200-yr bin**, per SD of predictor (`heinrich`
per stadial). Effects accumulate through the autoregression while the covariate persists.

## 6. Model-averaged coefficients (`time_lin`)

Unconditional SEs per Burnham & Anderson eq. 4.9, i.e. within-model variance plus
between-model spread. `full` substitutes zero where a term is absent, so it is a shrinkage
estimate that folds in selection uncertainty — this is the column to quote.

| term | summed w | conditional | ± | full | ± | z |
|---|---|---|---|---|---|---|
| mean_co2 | 0.969 | −0.477 | 0.222 | −0.462 | 0.233 | −1.98 |
| time | 0.917 | 0.307 | 0.137 | 0.281 | 0.156 | 1.80 |
| heinrich | 0.812 | 0.493 | 0.254 | 0.401 | 0.299 | 1.34 |
| PrDens | 0.564 | −0.328 | 0.200 | −0.185 | 0.222 | −0.84 |
| d18O | 0.969 | 0.144 | 0.200 | 0.139 | 0.198 | 0.70 |
| char_acc | 0.233 | 0.006 | 0.114 | 0.001 | 0.055 | 0.03 |

**Nothing reaches |z| = 2.** That is the honest state of the evidence with 83 observations,
two collinear climate proxies and genuine selection uncertainty. The joint result and the
variable-importance ordering are on much firmer ground than any single coefficient.

## 7. Sensitivity

| test | result |
|---|---|
| λ = 0.20 / 0.25 / 0.30 | dev 24.96 / 24.99 / 24.96, P ≈ 1.4×10⁻⁴ throughout — not knife-edge |
| p = 1 / 2 / 3 | AICc 713.6 / 714.1 / 713.5; coefficients barely move; dominant eigenvalue 0.77–0.84 |
| without the Holocene (bins ≥ 54; 70 obs) | dev 23.97, df 4, **P = 8.1×10⁻⁵**; all signs unchanged |
| dropping `na.interp`-filled predictor bins | dev 25.13, df 5, P = 1.3×10⁻⁴ |
| observation error τ = 0 / 0.5 / 1.0 | predictor rank order stable to within one place |

The pre-Holocene test excludes HUMAN by design — `PrDens` is identically zero before ~13 ka,
so it cannot be estimated on that subset. That the climate and stadial signal *strengthens*
there, with no human term available, is the cleanest evidence in the analysis that the result
is not an artefact of the late record.

## 8. What to be careful about in writing

- **Quote `c` per 200-yr bin, not the long-run effect.** `Σb` = 0.85 in the five-predictor
  model, so `1/(1 − Σb)` has a 95% range of 3.8–23.4. Individual lag coefficients are
  unidentified (SEs ≈ 0.43 while their *sum* has SE 0.057); only persistence is estimable,
  never which lag carries it. That is expected when the median gap between observations is
  three bins.
- **`mean_co2`'s Wald interval touches zero while its LRT gives P = 0.008.** Prefer the LRT.
- **Do not report `PrDens` without saying TIME was in the model.** It is significant when TIME
  is absent and marginal when it is present, and that difference is the finding.
- **The observation model is deliberately conservative.** τ = 0 uses only counting error from
  `ocfs_uncertainty.csv`; allowing extra observation error strengthens every result
  (omnibus dev 16.6 → 24.5), so nothing is hidden by the conservative choice. See section 2 of
  `ocfs_analysis.r` for why τ is not estimated.
- **7 bins have zero spore counts**, handled with a Poisson upper-bound SD on the λ scale. A
  binomial model of `OCFS_total` against the `Lycopodium` tracer would handle them properly and
  remains the natural robustness check if a reviewer presses on the count structure.

## Outputs

| file | contents |
|---|---|
| `results/ocfs/model_ranking_time_lin.csv`, `_time_hol.csv` | all 32 models, AICc, weights, evidence ratios |
| `results/ocfs/variable_importance.csv` | summed weights per block |
| `results/ocfs/model_averaged_time_lin.csv`, `_time_hol.csv` | averaged coefficients with unconditional SEs |
| `results/ocfs/coefficients_with_ci.csv` | five-predictor model, c and long-run with CIs |
| `results/ocfs/drop_one_lambda025_c.csv` | drop-one LRTs |
| `results/ocfs/smoothed_lambda025.csv` | smoothed series, SE, and back-transformed band |
| `results/ocfs/ocfs_smoothed.png`, `tau_profile.png` | figures |
| `results/ocfs/tier2_fits.rds`, `tier3_selection.rds` | all fitted objects |
