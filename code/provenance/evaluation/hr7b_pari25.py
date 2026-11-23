"""Studie VII · 3.9 pro sadu hr7b: prvních 25 kladných nul každé L-funkce v PARI (k porovnání s LMFDB) + Conreyho označení."""
import json, math, cypari2
P = cypari2.Pari(); P.allocatemem(10 ** 9, silent=True); P.default('realprecision', 38)
S = {'q15_c2-1': (15, [2, 1]), 'q21_c3-1': (21, [3, 1]), 'q21_c2-1': (21, [2, 1]), 'q21_c4-1': (21, [4, 1])}
out = {}
for jm, (q, chi) in S.items():
    G = 'znstar(%d,1)' % q; c = '[' + ','.join(map(str, chi)) + ']'
    m = [m for m in range(1, q) if math.gcd(m, q) == 1 and str(P('znconreychar(%s,%d)' % (G, m))).replace(' ', '') == c]
    z = sorted(float(v) for v in P('lfunzeros(lfuncreate([%s,%s]), 60)' % (G, c)))[:25]
    out[jm] = dict(q=q, chi=chi, conrey='%d.%s' % (q, m[0] if len(m) == 1 else m), pari=z); print(jm, out[jm]['conrey'], len(z), z[:3])
json.dump(out, open('data/hr7b_pari25.json', 'w'), indent=1)
