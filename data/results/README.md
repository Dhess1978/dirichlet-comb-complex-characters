# Results files of study VII

The files are the outputs of the evaluation scripts in `code/provenance/`, unchanged (keys in Czech, as written by the scripts). Sections refer to the article.

| file | written by | sections |
|---|---|---|
| `hr7_teorie.json` | `evaluation/hr7_teorie.py` (fingerprint recorded in addendum 1) | predictions of the comb (two-point function without the Euler factors of p \| q, naive ζ, GUE) |
| `hr7_meze_teorie.json` | `evaluation/hr7_meze_teorie.py` | predictions of the comb for every subset S of {2, 3, 5, 7} (test 3.12) |
| `hr7_eval_hr7.json` | `evaluation/hr7_eval.py` | 4.1 data checks, 4.2 comb and pairs, teeth, 4.3 first-order phase, explicit formula, diagonal regime |
| `hr7_eval2_hr7.json` | `evaluation/hr7_eval2.py` | 4.4 union modulo 5 and 7, 4.5 cross-covariance, 4.6 number variance, 4.7 Hann window for the L-functions |
| `hr7_eval3_final.json` | `evaluation/hr7_eval3.py` | 4.4 and 5.6 orders in four fields (test 3.13), 4.8 and 5.4 instrument (3.11), 5.5 tuning (3.12), 5.2 comb by prime limit |
| `hr7_r_P8.json`, `hr7_r_P3e10.json`, `hr7_r_kontrola.json`, `hr7_r_verdikty.json` | `evaluation/hr7_r.py` | 4.7 robustness of the main resonance (zeros of ζ) |
| `hr7_dvojce_hr7.json` | `evaluation/hr7_dvojce.py` | 4.3 independent second implementation of the comb and of the coefficients |
| `hr7_kontroly_K1K4.json`, `hr7_k23_vysledek_mac.json`, `hr7b_k23_vysledek_mac.json` | `evaluation/hr7_kontroly.py`, `zeros/hr7_k23.py` | 4.9 certification (K1, K4) and recomputation (K2, K3) |
| `hr7_lmfdb.json`, `hr7b_lmfdb.json`, `hr7b_pari25.json` | comparison with the LMFDB (first zeros) | 2 checks |
| `hr7_uplnost.json`, `hr7_skoky.json`, `hr7_skoky_po_oprave.json` | `evaluation/hr7_uplnost.py`, `evaluation/hr7_skoky.py` | 4.10 completeness, count test before and after the repair |
| `hr7_druhy_vypocet_porovnani.json` | comparison of the second computation on the desktop computer | 4.1 (4,397 windows identical bit for bit) |
| `hr7_soused11.json` | `evaluation/hr7_soused11.py` | 4.11 resolution near 11^5 and 11^6 |
| `hr7_rmt_teorie.json`, `hr7_rmt_mereni.json`, `hr7_rmt_vysledek.json` | `random_matrices/hr7_rmt_teorie.py`, `hr7_rmt.py` | 4.12 pair correlation and spacings, 4.13 unions (addendum 9) |
| `hr7_rmt_teorie_vyska_036.json`, `hr7_rmt_teorie_vyska_400.json`, `hr7_rmt_vyska_vysledek.json`, `hr7_rmt_vyska_obr.json` | `random_matrices/hr7_rmt_teorie_vyska.py`, `hr7_rmt_vyska.py`, `hr7_rmt_vyska_obr.py` | 4.14 height (addendum 10) |
| `hr7c_skoky.json`, `hr7c_vzorek.json`, `hr7c_mp_*.json` | `evaluation/hr7_skoky.py`, `random_matrices/hr7c_vzorek.py`, `hr7c_mpmath.py` | 4.14 checks of the zeros of the second band (count test, 24 zeros in mpmath) |

The check files `hr7_lmfdb.json`, `hr7b_lmfdb.json`, `hr7b_pari25.json`, `hr7c_vzorek.json`, `hr7c_mp_*.json` and `../manifests/oprava_hr7b_dodatek7.json` contain the few individual zeros that were compared (the first 25 zeros of each L-function with the LMFDB, 24 zeros of the second band in mpmath, the four added zeros of the repair); the zeros themselves are not part of the repository.
