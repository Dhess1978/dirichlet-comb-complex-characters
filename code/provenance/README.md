# Scripts as they were run

These are the scripts of study VII as they were run, with the paths of the computers they ran on. Their SHA-256 fingerprints are those recorded in the rules (`data/rules/`) before the data they evaluate were opened; the logs of the runs are in `logs/`.

- `zeros/` — computation of the zeros with PARI/GP 2.17.2 through cypari2: `fz_nuly_chi_v3.py` (set hr7, fingerprint cf90cceb…), `fz_nuly_chi_v3b.py` (sets hr7b and hr7c, f350356b…), `doplnek9_mac.py` and `doplnek9.py` (the nine windows interrupted by a keyboard interrupt, recomputed on the desktop computer and independently in the cloud), `hr7b_oprava.py` (the repair of the two windows of addendum 7), `hr7_k23.py` (checks K2 and K3 on the desktop computer), `manifest_z_oken.py` and `hr7b_slouc.py` (manifests from the window files).
- `evaluation/` — `hr7_teorie.py`, `hr7_meze_teorie.py`, `hr7_sigma_bk.py`, `hr7_berry.py` (predictions), `hr7_eval.py`, `hr7_eval2.py`, `hr7_eval3.py` (sections 3.1–3.13 of the rules), `hr7_r.py` (part R), `hr7_kontroly.py` (K1, K4 with python-flint), `hr7_dvojce.py` (independent second implementation), `hr7_uplnost.py`, `hr7_skoky.py` (count test), `hr7b_pari25.py`, `hr7_soused11.py` (addendum 8), `graf_hr7.py` (working figures).
- `random_matrices/` — `hr7_rmt_teorie.py`, `hr7_rmt.py` (addendum 9), `hr7_rmt_teorie_vyska.py`, `hr7_rmt_vyska.py`, `hr7_rmt_vyska_obr.py`, `hr7c_vzorek.py`, `hr7c_mpmath.py` (addendum 10).

The zeros of ζ used for part R and for the unions (the bands P8 and P3e10 and the band 36,000–80,000) are those of studies I and VI (tables of D. Platt in the LMFDB and of A. Odlyzko).
