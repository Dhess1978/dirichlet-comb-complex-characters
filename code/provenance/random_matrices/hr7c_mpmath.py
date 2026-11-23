"""Studie VII · dodatek 10: kontrola nul horniho pasma (400 000-420 000) druhou knihovnou - |L(1/2 + i*gamma, chi)| v mpmath pres Hurwitzovu zetu.
Charakter z PARI (chareval) jako v hr7_eval.py; zadna cast vypoctu nul se nepouziva. Kontrolni bod: stred mezi nulou a jejim sousedem.
Pouziti: python3 hr7c_mpmath.py VZOREK.json VYSTUP.json JMENO_L_FUNKCE"""
import json, math, sys, time, numpy as np, mpmath as mp, cypari2
PARI = cypari2.Pari(); mp.mp.dps = 20
CH = {'q15_c2-1': (15, [2, 1]), 'q21_c3-1': (21, [3, 1]), 'q21_c2-1': (21, [2, 1]), 'q21_c4-1': (21, [4, 1])}
def chi_tab(q, chi):
    G = 'znstar(%d,1)' % q; c = '[' + ','.join(str(int(x)) for x in chi) + ']'; T = np.zeros(q, complex)
    for a in range(q):
        if math.gcd(a, q) == 1: T[a] = complex(mp.expjpi(2 * mp.mpf(str(PARI('chareval(%s,%s,%d)' % (G, c, a))))))
    for a in range(q):
        for b in range(q): assert abs(T[a * b % q] - T[a] * T[b]) < 1e-12
    assert int(PARI('zncharconductor(%s,%s)' % (G, c))) == q
    return T
def Lhod(q, Tab, t):
    s = mp.mpc(0.5, t); return float(abs(mp.mpf(q) ** (-s) * mp.fsum(mp.mpc(Tab[a].real, Tab[a].imag) * mp.zeta(s, mp.mpf(a) / q) for a in range(1, q) if abs(Tab[a]) > 0)))
Vz = json.load(open(sys.argv[1])); jm = sys.argv[3]; q, chi = CH[jm]; Tab = chi_tab(q, chi); d = Vz['L'][jm]; t0 = time.time(); R = dict(q=q, chi=chi, nuly=[], stredy=[])
for g, g2 in zip(d['gamma'], d['dalsi']):
    g = float(g); g2 = float(g2); R['nuly'].append(dict(gamma=g, absL=Lhod(q, Tab, g))); R['stredy'].append(dict(t=(g + g2) / 2, absL=Lhod(q, Tab, (g + g2) / 2)))
    print(jm, 'gamma %.6f |L| %.3e | stred |L| %.3e | %.0f s' % (g, R['nuly'][-1]['absL'], R['stredy'][-1]['absL'], time.time() - t0), flush=True)
R['max_absL_v_nulach'] = max(x['absL'] for x in R['nuly']); R['min_absL_ve_stredech'] = min(x['absL'] for x in R['stredy'])
json.dump(R, open(sys.argv[2], 'w'), indent=1)
