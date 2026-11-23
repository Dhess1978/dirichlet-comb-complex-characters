"""Průběžné grafy studie VII z výstupu hr7_eval.py (jen čte JSON; nic nepočítá znovu).
Použití: python3 graf_hr7.py VYSTUP_hr7_eval.json PREFIX_PNG"""
import json, math, sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
E = json.load(open(sys.argv[1])); PRE = sys.argv[2]; d = np.array(json.load(open('/home/claude/br107/data/hr7_teorie.json'))['delta'])
BAR = ['#3987e5', '#e3b341', '#9b59b6', '#d95926', '#5dcaa5']
# 1) hřeben: pro každý modul data všech L-funkcí téhož modulu + předpověď
qs = sorted(set(v['q'] for v in E['hreben'].values())); fig, axs = plt.subplots(len(qs), 1, figsize=(12, 3.0 * len(qs)), squeeze=False)
for a, q in zip(axs[:, 0], qs):
    js = [j for j, v in E['hreben'].items() if v['q'] == q]; v0 = E['hreben'][js[0]]
    a.plot(d, v0['V_GUE'], '--', color='#d95926', lw=1.2, label='GUE without primes'); a.plot(d, v0['V_naiv'], ':', color='#8c91a6', lw=1.2, label='naive ζ (p | q not removed)')
    a.plot(d, v0['V_BK_q'], color='#199e70', lw=2, label='prediction with primes, modulus %d' % q)
    for i, j in enumerate(js + [k for k in E.get('reference_VI', {}) if abs(E['reference_VI'][k]['D']) == q]):
        src = E['hreben'].get(j) or E['reference_VI'][j]; lab = j if j in E['hreben'] else '%s (real, study VI)' % j
        a.errorbar(d, src['V'], yerr=src['err'], fmt='o', ms=2.2, elinewidth=0.5, color=BAR[i % len(BAR)], alpha=0.85, label=lab)
    Qs = ', '.join('%s %.2f' % (j, E['hreben'][j]['Q']) for j in js); a.set_title('q = %d · Q against the prediction: %s' % (q, Qs), fontsize=9); a.legend(fontsize=7, ncol=3); a.grid(alpha=.3); a.set_ylabel('V')
axs[-1, 0].set_xlabel('δ = ln x − ln(qT/2π)'); plt.tight_layout(); plt.savefig(PRE + '_hreben.png', dpi=120); plt.close()
# 2) fázové kolo: naměřené c proti předpovědi −w·χ(p)·ln p/√p (normováno na w·ln p/√p)
fig, ax = plt.subplots(1, 2, figsize=(12, 5.6))
a = ax[0]; a.plot([-180, 180], [-180, 180], color='#199e70', lw=1.5, label='prediction: arg c = arg(−χ(p))'); a.plot([-180, 180], [180, -180], '--', color='#d95926', lw=1.2, label='conjugate phase: arg c = −arg(−χ(p))')
for k, r in E['faze'].items():
    c = complex(*r['c']); cp = complex(*r['c_pred']); ang = lambda z: math.degrees(math.atan2(z.imag, z.real))
    a.plot(ang(cp), ang(c), 'o', color='#3987e5', ms=4 + 2 * (r['p'] == 2))
a.set_xlim(-185, 185); a.set_ylim(-185, 185); a.set_xticks(range(-180, 181, 60)); a.set_yticks(range(-180, 181, 60)); a.grid(alpha=.3); a.legend(fontsize=8, loc='upper left')
a.set_xlabel('predicted phase of the first coefficient (°)'); a.set_ylabel('measured phase (°)')
mx = max(r['rozdil'] for r in E['faze'].values())
a.set_title('the melody carries the phase: %d combinations L-function × prime\nlargest |c − c_pred| = %.1e' % (len(E['faze']), mx), fontsize=9)
# 3) histogram jedné kombinace s komplexním χ(p)
k = max(E['faze'], key=lambda k: abs(E['faze'][k]['chi_p'][1]) * (E['faze'][k]['p'] == 2)); r = E['faze'][k]; x = (np.arange(50) + 0.5) / 50; N = r['N']
b = ax[1]; b.bar(x, np.array(r['O']) / (N / 50), width=0.018, color='#3987e5', alpha=.7, label='data')
b.plot(x, np.array(r['E']) / (N / 50), color='#199e70', lw=2, label='prediction with the phase χ(p)'); b.plot(x, np.array(r['E_sdruzena']) / (N / 50), '--', color='#d95926', lw=1.5, label='conjugate phase χ̄(p)')
b.axhline(1, color='#8c91a6', lw=.8); b.set_ylim(0.4, 1.5); b.set_xlabel('{γ·ln p/2π}'); b.set_ylabel('density'); b.legend(fontsize=8); b.grid(alpha=.3)
b.set_title('%s, p = %d, χ(p) = %+.3f%+.3fi: χ² p = %.2g, conjugate p = %.1g' % (k.split(',')[0], r['p'], r['chi_p'][0], r['chi_p'][1], r['p_hodnota'], r['p_sdruzena']), fontsize=9)
plt.tight_layout(); plt.savefig(PRE + '_faze.png', dpi=120); plt.close(); print('ok')
