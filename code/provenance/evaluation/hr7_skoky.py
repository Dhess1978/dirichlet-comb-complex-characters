"""Studie VII · kontrola úplnosti nul (Turingova typu, dodatek 7): r_n = θ_χ(γ_n)/π − n; chybějící nula = trvalý skok r o +1.
Skok v indexu i: průměr r přes 400 následujících nul minus průměr přes 400 předchozích; hlásí se |skok| > 0,6. Jen čte data.
Použití: python3 hr7_skoky.py VYSTUP.json SLOZKA [SLOZKA ...]"""
import json, os, gzip, math, sys, hashlib, numpy as np
from scipy.special import loggamma
out = dict(skript_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest(), L={}); W = 400
for V in sys.argv[2:]:
    M = json.load(open(os.path.join(V, 'manifest.json')))
    for jm, z in M['L'].items():
        g = []; ok = []
        for k in sorted(int(x) for x in z['bloky']):
            b = json.loads(gzip.decompress(open(os.path.join(V, jm, z['bloky'][str(k)]['soubor']), 'rb').read())); g += b['g']; ok += [k] * b['n']
        g = np.array(g); q = int(z['popis']['q']); odd = int(z['popis']['lichy'])
        r = (0.5 * g * math.log(q / math.pi) + np.imag(loggamma((0.5 + odd) / 2 + 0.5j * g))) / math.pi - np.arange(len(g))
        c = np.cumsum(np.r_[0, r]); d = np.array([(c[i + W] - c[i]) / W - (c[i] - c[i - W]) / W for i in range(W, len(r) - W)])
        sk = []; i = 0
        while i < len(d):
            if abs(d[i]) > 0.6: j = i + int(np.argmax(np.abs(d[i:i + 2 * W]))); sk.append(dict(okno=int(ok[j + W]), gamma=float(g[j + W]), skok=round(float(d[j]), 3))); i = j + 2 * W
            else: i += 1
        out['L'][jm] = dict(slozka=os.path.basename(os.path.normpath(V)), n=int(len(g)), nejvetsi_abs_mimo_skoky=float(np.max(np.abs(d))) if not sk else None, skoky=sk)
        print('%-10s n=%d  %s' % (jm, len(g), ('skoky: ' + ', '.join('okno %d γ≈%.2f %+.2f' % (s['okno'], s['gamma'], s['skok']) for s in sk)) if sk else 'bez skoku (největší |Δ| %.2f)' % np.max(np.abs(d))), flush=True)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
