#!/usr/bin/env python3
"""Studie VII · dodatek 10: vzorek nul horniho pasma pro kontrolu druhou knihovnou (mpmath, cloud). Jen cte.
Pro kazdou L-funkci 6 nahodnych nul (pevne seminko) a k nim sousedni nula (stred mezi nimi = kontrolni bod, kde |L| nema byt male).
Pouziti: python3 -I hr7c_vzorek.py SLOZKA_HR7C VYSTUP.json"""
import gzip, hashlib, json, os, sys, numpy as np
V, OUT = sys.argv[1:3]; LF = ['q15_c2-1', 'q21_c2-1', 'q21_c3-1', 'q21_c4-1']; rng = np.random.default_rng(34110); O = dict(skript_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest(), L={})
for jm in LF:
    g = []
    for k in range(4000, 4200): b = json.loads(gzip.decompress(open(os.path.join(V, jm, 'okno_%06d.json.gz' % k), 'rb').read())); g.extend(b['g'])
    g = np.array(g); i = np.sort(rng.choice(len(g) - 1, 6, replace=False)); O['L'][jm] = dict(n=int(len(g)), index=[int(x) for x in i], gamma=[repr(float(g[x])) for x in i], dalsi=[repr(float(g[x + 1])) for x in i])
json.dump(O, open(OUT, 'x'), indent=1); print({j: O['L'][j]['gamma'][:2] for j in LF})
