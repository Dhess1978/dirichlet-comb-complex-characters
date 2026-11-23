#!/usr/bin/env python3
"""Studie VII · dodatek 10 (M3): data pro obrazek - parova korelace a rozestupy horniho pasma s bloky (pro chybove usecky).
Zadny verdikt, zadne nove kriterium: tytez funkce (hladky, parova, rozestupy) jako zmrazene vyhodnoceni; R2 se musi shodovat s R2_nahore ve vysledku.
Pouziti: python3 -I hr7_rmt_vyska_obr.py SLOZKA_SE_SKRIPTEM_hr7_rmt SLOZKA_HR7C VYSLEDEK_VYSKA.json VYSTUP.json
Vstupy jen cte."""
import gzip, hashlib, json, os, sys, numpy as np
SK, V, FV, OUT = sys.argv[1:5]; sys.path.insert(0, SK); import hr7_rmt as R
TA, TB = 400000.0, 420000.0; K0, K1 = 4000, 4200; LF = ['q15_c2-1', 'q21_c2-1', 'q21_c3-1', 'q21_c4-1']
VY = json.load(open(FV)); assert VY['uplnost_v_poradku']
O = dict(cas_UTC=R.cas(), skript_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest(), vysledek_sha256=hashlib.sha256(open(FV, 'rb').read()).hexdigest(), TA=TA, TB=TB, L={})
for jm in LF:
    g = R.okna(os.path.join(V, jm), K0, K1); g = g[(g > TA) & (g <= TB)]; q = R.modul(jm); w = R.hladky(g, q); P = R.parova(w); S = R.rozestupy(w, R.HS)
    shoda = bool(P['R2'] == VY['L'][jm]['R2_nahore'] and S['p'] == VY['L'][jm]['p_nahore'] and len(g) == VY['L'][jm]['nul_nahore'])
    O['L'][jm] = dict(q=q, nul=int(len(g)), parova=P, rozestupy=S, shoda_s_vyhodnocenim=shoda); print(jm, len(g), 'shoda s vyhodnocenim:', shoda)
O['shoda_vse'] = bool(all(x['shoda_s_vyhodnocenim'] for x in O['L'].values()))
json.dump(O, open(OUT, 'x'), ensure_ascii=False); print('shoda vse:', O['shoda_vse'])
