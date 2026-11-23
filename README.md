# Different instruments, the same chord: the comb of resonances and the first-order phase for Dirichlet L-functions with complex characters, the cyclotomic union, and the robustness of the main resonance of study VI (study VII)

**David Hess (MAYAN ALFA V1008 LAB v.34.116) · 2026** · [Česká verze](README.cs.md)

Seventh study of the series, after *Berry saturation and the statistics of the Riemann zeros* (study I, [doi:10.5281/zenodo.23815328](https://doi.org/10.5281/zenodo.23815328)), *Beyond the diagonal approximation* (study II, [doi:10.5281/zenodo.23928017](https://doi.org/10.5281/zenodo.23928017)), *Low-lying zeros of 795,381 quadratic Dirichlet L-functions* (study III, [doi:10.5281/zenodo.23991165](https://doi.org/10.5281/zenodo.23991165)), *Three-point correlation of the Riemann zeros at finite height* (study IV, [doi:10.5281/zenodo.24094318](https://doi.org/10.5281/zenodo.24094318)), *Central values of quadratic Dirichlet L-functions* (study V, [doi:10.5281/zenodo.24256378](https://doi.org/10.5281/zenodo.24256378)) and *The fractional parts of γ ln p/2π beyond the first order* (study VI, [doi:10.5281/zenodo.24381946](https://doi.org/10.5281/zenodo.24381946)). Study VI measured a comb of resonances in the variance of the windowed sums Σ g x^{iγ} over the zeros of ζ and of Dirichlet L-functions with real characters. This study extends the measurement to complex characters and asks what depends on the character and what only on its modulus. Every rule was recorded before the zeros were computed or, for the later additions, before the quantity was evaluated.

| result | value |
|---|---|
| data · zeros of 16 Dirichlet L-functions modulo 5, 7, 13, 15, 21 between heights 36,000 and 80,000, computed for this study with PARI/GP | 1,296,658 zeros (twelve L-functions with complex characters, 958,685 zeros; four that complete the families modulo 15 and 21); 178,844 zeros of four of them between 400,000 and 420,000 |
| comb of a single L-function against the two-point function with the Euler factors of p \| q removed | mean z² over 201 points 0.13–0.25 for the twelve; GUE without arithmetic 11.5–49; the 18 pairs of L-functions of one modulus have the same comb (0.13–0.29) |
| first order | the first Fourier coefficient of {γ ln p/2π} equals −w χ(p) ln p/√p to 5.2·10⁻⁵ in 56 combinations; the conjugate phase is excluded (p ≤ 7·10⁻¹⁵⁴) |
| the union of all characters of a modulus (zeros of the Dedekind zeta function of Q(ζ_q), q = 5, 7, 15, 21, up to 947,573 zeros) | a prime sounds only in the harmonics that are multiples of its order modulo q; coefficients to 7.7·10⁻⁵ |
| cross-covariance of two L-functions of one modulus | equals the prime trace with the relative phase χ_a χ̄_b (19 errors at ln x = 3, relative difference 4·10⁻⁵) |
| number variance | saturates at a level set by the modulus (Q = 0.32–1.72 against the two-point function) |
| robustness of the main resonance of study VI (zeros of ζ at t ≈ 10⁸ and 3·10¹⁰) | holds for the Hann window, six error estimates and the full covariance of the scan; GUE between 12 and 69 errors away depending on the error estimate |
| three blind tests | modulus and character, primes dividing the modulus, and the orders of 16 primes in four cyclotomic fields identified from the zeros alone, 16 of 16 in each test |
| random matrices | pair correlation of each L-function up to three mean spacings: χ² per bin 0.78–1.86 against the two-point function of Bogomolny and Keating, 6.3–16.5 against GUE alone; spacings of the unions = superposition of independent GUE spectra to 0.010; between the two bands the deviation from GUE falls to 0.719 ± 0.021 of its value, predicted 0.727 |
| completeness | a count test of Turing type found two pairs of close zeros missed by the search grid; repaired before the evaluation |
| verdicts of the 18 hypotheses recorded beforehand | all SUPPORTED |

The statements behind the first order, the union and the cross-covariance are theorems or direct consequences of them (Landau's formula, the explicit formula, orthogonality of characters); the two-point function behind the comb, the number variance and the pair correlation is a conjecture. What is added is data for complex characters with null models and errors, zeros that did not exist when the rules were written, and the blind tests. Agreement is numerical support, not a proof; nothing here bears on the Riemann hypothesis.

How the rules recorded beforehand, the laboratory and an AI assistant were used is described in Section 3.10 of the article.

## Contents

| folder | what |
|---|---|
| `paper/` | the article `comb_complex_characters_en` (PDF, DOCX, Markdown, LaTeX `main.tex`); the Czech version in `paper/cs/` |
| `figures/` | figures 1–9 in English (`_en`) and Czech (`_cz`), PNG 300 dpi and PDF |
| `data/results/` | the results files written by the evaluation scripts, unchanged (every number of the article; see `data/results/README.md` for which file holds which section) |
| `data/manifests/` | manifests of the zeros (SHA-256 of every window, environment, script checksum): 5,280 windows of set hr7, 1,760 of set hr7b with the two repaired windows (`oprava_hr7b_dodatek7.json`), 800 of set hr7c |
| `data/rules/` | the rules recorded before the computation with their ten dated addenda: `PRAVIDLA_STUDIE_VII.md` (rules and addenda 1–8), `PRAVIDLA_STUDIE_VII_DODATEK_9_RMT.md`, `PRAVIDLA_STUDIE_VII_DODATEK_10_RMT_VYSKA.md` (Czech, as written) |
| `code/` | `independent_check_vii.py` with its log (recomputes the statistics, the closed-form predictions and the 18 verdicts from the results files with new code), `check_text_vii.py` (every number of the English article against the results files), `check_numbers_en_cs.py` (numbers of the English and Czech texts), `figures_vii.py` (figures from the results files); `provenance/` holds the scripts as they were run (see `code/provenance/README.md`) |

## Reproducing

```
pip install -r requirements.txt
python3 -I code/independent_check_vii.py .            # 25 checks: statistics, predictions and verdicts from the results files
python3 -I code/check_text_vii.py .                   # every number of the English article against the results files
python3 -I code/check_numbers_en_cs.py paper/comb_complex_characters_en.md paper/cs/comb_complex_characters_cs.md
python3 -I code/figures_vii.py . figures              # figures 1-9, English and Czech, PNG and PDF
```

The zeros were computed with PARI/GP 2.17.2 (cypari2) by `code/provenance/zeros/fz_nuly_chi_v3.py` (set hr7) and `fz_nuly_chi_v3b.py` (sets hr7b and hr7c) on a desktop computer; the zeros themselves are not included, but the manifests with the SHA-256 of every window are, so a recomputation can be compared window by window. The zeros of ζ are available from A. Odlyzko and from the LMFDB (zeros of D. Platt). The scripts in `code/provenance/` contain paths of the computers on which they ran and are given as a record; their SHA-256 fingerprints are those recorded in the rules before the data were evaluated.

Text, figures and data: CC BY 4.0 (`LICENSE-CC-BY-4.0.md`); code: MIT (`LICENSE`). Czech summary of the licences: `LICENCE.cs.md`.
