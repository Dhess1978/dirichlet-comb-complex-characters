"""Studie VII · grafy hudebního čtení (3.12 slyšet ladění, 3.13 slyšet těleso) z výstupu hr7_eval3.py — jen čte JSON, nic nepočítá znovu.
Použití: python3 graf_hr7_hudba.py VYSTUP_hr7_eval3.json PREFIX_PNG"""
import json, math, sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
E = json.load(open(sys.argv[1])); PRE = sys.argv[2]
INK, MUTED, GRID = '#22252e', '#8c91a6', '#d9dbe3'; BLUE, GREEN, ORANGE = '#3987e5', '#199e70', '#d95926'
# 1) 3.12: Q pro 16 kandidátních „ladění“ S (která prvočísla chybějí v hřebeni), panel na modul
L = E['ladeni']; qs = sorted(set(v['q'] for v in L.values())); S = list(next(iter(L.values()))['Q'].keys())
fig, axs = plt.subplots(len(qs), 1, figsize=(12, 2.35 * len(qs)), sharex=True, squeeze=False)
for a, q in zip(axs[:, 0], qs):
    js = [j for j, v in L.items() if v['q'] == q]; pr = L[js[0]]['pravda']; ip = S.index(pr)
    a.axvspan(ip - 0.42, ip + 0.42, color=GREEN, alpha=0.14, lw=0)
    for i, j in enumerate(js):
        y = [L[j]['Q'][s] for s in S]; x = np.arange(len(S)) + (i - (len(js) - 1) / 2) * 0.14
        a.plot(x, y, 'o', ms=4.5, color=[BLUE, ORANGE, '#9b59b6', '#e3b341', '#5dcaa5'][i % 5], label='%s: chosen %s (Q %.2f)' % (j, L[j]['vybrano'], L[j]['Q'][L[j]['vybrano']]))
    a.set_yscale('log'); a.set_ylim(0.05, 2000); a.axhline(1, color=MUTED, lw=0.8, ls=':'); a.grid(axis='y', color=GRID, lw=0.6)
    a.set_ylabel('Q (log)', color=INK); a.set_title('q = %d · actually missing %s' % (q, pr if pr != '[]' else 'none of the candidates 2, 3, 5, 7 (13 is not a candidate)'), fontsize=9, loc='left', color=INK)
    a.legend(fontsize=7, loc='lower right', frameon=False, ncol=2 if len(js) > 3 else 1)
    for s_ in ('top', 'right'): a.spines[s_].set_visible(False)
axs[-1, 0].set_xticks(range(len(S))); axs[-1, 0].set_xticklabels([s.replace('[]', '∅') for s in S], fontsize=8, rotation=30)
axs[-1, 0].set_xlabel('candidate tuning: primes removed from the comb (S)', color=INK)
ok = sum(v['spravne'] for v in L.values())
fig.suptitle('3.12 hearing the tuning: the correct S chosen from the zeros alone for %d / %d L-functions (green band = truth)' % (ok, len(L)), fontsize=11, color=INK, x=0.01, ha='left')
plt.tight_layout(); plt.savefig(PRE + '_ladeni.png', dpi=130); plt.close()
# 2) 3.13: normované koeficienty c(p, k)/(w·φ·ln p/p^{k/2}) sjednocení všech charakterů modulu; předpověď −1 právě když p^k ≡ 1 (mod q)
T = E['teleso2']; qs = sorted(T, key=int); NC = 2; NR = (len(qs) + 1) // 2; fig, axs = plt.subplots(NR, NC, figsize=(6.6 * NC, 5.3 * NR), squeeze=False)
def rad(p, q): return next(k for k in range(1, 50) if pow(p, k, q) == 1)
for a, q in zip(axs.ravel(), qs):
    r = T[q]; w = r['w']; phi = len(r['clenove']); ps = sorted(set(int(k.split(',')[0]) for k in r['c'])); qq = int(q)
    for i, p in enumerate(ps):
        for k in range(1, 7):
            v = r['c']['%d,%d' % (p, k)]; nrm = w * phi * math.log(p) / p ** (k / 2); c = complex(*v['c']) / nrm; pred = v['c_pred'] / nrm
            x = k + (i - (len(ps) - 1) / 2) * 0.11
            col = BLUE if pred == 0 else GREEN; e = r['sigma'] / nrm; big = v['x'] > 200
            a.errorbar(x, max(min(c.real, 0.55), -1.45), yerr=min(e, 3), fmt='none', ecolor=col, elinewidth=0.8, alpha=0.5, zorder=2)
            a.plot(x, max(min(c.real, 0.55), -1.45), 'o', ms=5, color=col if not big else 'white', mec=col, mew=1.0 if big else 0.6, zorder=3)
            if pred != 0: a.plot(x, pred, '_', ms=12, mew=2, color=INK, zorder=2)
    a.axhline(0, color=MUTED, lw=0.8); a.axhline(-1, color=INK, lw=0.8, ls=':'); a.set_ylim(-1.5, 0.6); a.set_xticks(range(1, 7))
    a.set_xlabel('harmonic k (x = p^k)', color=INK); a.grid(axis='y', color=GRID, lw=0.6)
    fo = r['f_odhad']; fp = r['f_pred']; ok_ = all(fo[p] == f for p, f in fp.items()); nav = [p for p in fo if p not in fp]
    pop = ('\ndescriptive, outside the test: ' + ', '.join('%s→%s (%s)' % (p, fo[p], 'correct' if fo[p] == rad(int(p), qq) else 'wrong, the order is %d; x = %s^%s > 200' % (rad(int(p), qq), p, fo[p])) for p in nav)) if nav else ''
    a.set_title('Q(ζ%s): %d sectors, N = %d zeros\norder from the data: %s · prediction %s%s' % (q, phi, r['N'], ', '.join('%s→%s' % (p, fo[p]) for p in fp), 'met' if ok_ else 'NOT MET', pop), fontsize=8.5, loc='left', color=INK)
    for s_ in ('top', 'right'): a.spines[s_].set_visible(False)
[axs[i, 0].set_ylabel('Re c(p, k) / (w·φ(q)·ln p / p^{k/2})', color=INK) for i in range(NR)]; [a.set_visible(False) for a in axs.ravel()[len(qs):]]
axs[0, 0].plot([], [], 'o', color=GREEN, label='p^k ≡ 1 (mod q): tooth predicted (−1)'); axs[0, 0].plot([], [], 'o', color=BLUE, label='p^k ≢ 1: prediction 0 (the sectors cancel)'); axs[0, 0].plot([], [], 'o', color='white', mec=MUTED, label='open circle: x = p^k > 200 (outside the reliable range)'); axs[0, 0].legend(fontsize=7.5, loc='lower left', frameon=False)
nr = sum(len(T[q]['f_pred']) for q in qs); nok = sum(T[q]['f_odhad'][p] == f for q in qs for p, f in T[q]['f_pred'].items())
fig.suptitle('3.13 hearing the body: a prime p sounds in the whole orchestra only in the harmonic equal to its order mod q — %d / %d orders correct\n(bars ±1 σ; points off the axis pressed to the edge)' % (nok, nr), fontsize=10.5, color=INK, x=0.01, ha='left')
plt.tight_layout(); plt.savefig(PRE + '_teleso.png', dpi=130); plt.close(); print('ok')
