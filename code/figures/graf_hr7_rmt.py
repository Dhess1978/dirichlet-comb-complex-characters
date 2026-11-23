"""Studie VII · grafy navazani na RMT (dodatek 9): M1 parova korelace jedne L-funkce, M2 rozestupy sjednoceni.
Jen cte JSON (mereni, vysledek, teorie); modely GUE a Poisson pro M2 pocita funkcemi merici ho skriptu hr7_rmt.py.
Pouziti: python3 -I graf_hr7_rmt.py MERENI.json VYSLEDEK.json TEORIE.json SLOZKA_SE_SKRIPTEM_hr7_rmt PREFIX_PNG en|cs"""
import json, math, sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
FM, FV, FT, SK, PRE, LG = sys.argv[1:7]; sys.path.insert(0, SK); import hr7_rmt as R
M = json.load(open(FM)); V = json.load(open(FV)); T = json.load(open(FT)); EN = LG == 'en'
INK, MUTED, GRID = '#22252e', '#8c91a6', '#d9dbe3'; BLUE, GREEN, ORANGE = '#3987e5', '#199e70', '#d95926'
def cz(x, nd):                                    # desetinna carka v ceske verzi
    s = ('%.' + str(nd) + 'f') % x; return s if EN else s.replace('.', ',')
def tx(en, cs): return en if EN else cs
def osy(a):
    for s_ in ('top', 'right'): a.spines[s_].set_visible(False)
    a.grid(color=GRID, lw=0.6); a.tick_params(colors=INK, labelsize=8)
sem = lambda B: np.array(B).std(0, ddof=1) / math.sqrt(len(B))
# Conreyho znaceni: u sady hr7 urceno z hodnot chi(p) tabulky pravidel (oddil 2.1) podle definice chi_q(m, n), u sady hr7b podle data/hr7b_lmfdb.json
CONREY = {'q5_c1': '5.2', 'q5_c3': '5.3', 'q7_c2': '7.2', 'q7_c4': '7.4', 'q7_c1': '7.3', 'q7_c5': '7.5', 'q13_c1': '13.2', 'q13_c11': '13.7', 'q15_c1-1': '15.2', 'q15_c3-1': '15.8',
          'q21_c1-1': '21.17', 'q21_c5-1': '21.5', 'q15_c2-1': '15.14', 'q21_c3-1': '21.20', 'q21_c2-1': '21.2', 'q21_c4-1': '21.11'}

# ---------- M1: R2 - GUE po modulech + chi2 na kos pro 16 L-funkci
hr = np.array(T['hranice']); d = (hr[:-1] + hr[1:]) / 2; gue = np.array(T['GUE']); qs = sorted(T['moduly'], key=int)
fig, axs = plt.subplots(2, 3, figsize=(13.2, 7.4)); axs = axs.ravel()
for a, q in zip(axs, qs):
    js = [j for j, x in M['M1'].items() if str(x['q']) == q]; bk = np.array(T['moduly'][q]['BK'])
    r = np.array([M['M1'][j]['parova']['R2'] for j in js]); e = np.array([sem(M['M1'][j]['parova']['bloky']) for j in js])
    for rr in r: a.plot(d, rr - gue, '-', color=BLUE, lw=0.6, alpha=0.35, zorder=2)
    a.errorbar(d, r.mean(0) - gue, yerr=np.sqrt((e ** 2).sum(0)) / len(js), fmt='o', ms=3.2, color=BLUE, ecolor=BLUE, elinewidth=0.8, zorder=4,
               label=tx('zeros: mean of %d L-functions', 'nuly: průměr %d L-funkcí') % len(js))
    a.plot(d, bk - gue, '-', color=INK, lw=1.6, zorder=3, label=tx('Bogomolny–Keating [8] − GUE', 'Bogomolny–Keating [8] − GUE'))
    a.axhline(0, color=ORANGE, lw=1.4, ls='--', zorder=2.4, label=tx('GUE alone', 'samotné GUE'))
    cb = [V['M1'][j]['chi2_BK'] for j in js]; cg = [V['M1'][j]['chi2_GUE'] for j in js]
    a.set_title(tx('q = %s · $\\chi^2$ per bin: [8] %s–%s, GUE %s–%s', 'q = %s · $\\chi^2$ na koš: [8] %s–%s, GUE %s–%s') % (q, cz(min(cb), 2), cz(max(cb), 2), cz(min(cg), 1), cz(max(cg), 1)),
                fontsize=9.5, loc='left', color=INK)
    a.set_ylim(-0.10, 0.10); a.set_xlim(0, 3); osy(a)
    a.set_xlabel(tx('distance d (mean spacings)', 'vzdálenost d (střední rozestupy)'), color=INK, fontsize=9)
    a.set_ylabel(tx('$R_2(d)$ − GUE', '$R_2(d)$ − GUE'), color=INK, fontsize=9)
    if q == qs[0]: a.legend(fontsize=7.5, frameon=False, loc='lower right')
a = axs[5]; js = list(V['M1']); x = np.arange(len(js))
a.plot(x, [V['M1'][j]['chi2_GUE'] for j in js], 's', ms=5.5, color=ORANGE, label=tx('against GUE alone', 'proti samotnému GUE'))
a.plot(x, [V['M1'][j]['chi2_BK'] for j in js], 'o', ms=5.5, color=BLUE, label=tx('against Bogomolny–Keating [8]', 'proti Bogomolnému–Keatingovi [8]'))
a.axhline(2.0, color=INK, lw=0.9, ls=':'); a.text(len(js) - 0.5, 2.15, tx('pass mark 2.0', 'mez 2,0'), fontsize=7.5, color=INK, ha='right', va='bottom')
a.set_yscale('log'); a.set_ylim(0.5, 30); a.set_xticks(x); a.set_xticklabels([CONREY[j] for j in js], rotation=60, fontsize=7.5); osy(a)
a.set_ylabel(tx('$\\chi^2$ per bin (60 bins, log scale)', '$\\chi^2$ na koš (60 košů, log. měřítko)'), color=INK, fontsize=9)
a.set_xlabel(tx('L-function (Conrey label of the character)', 'L-funkce (Conreyho značení charakteru)'), color=INK, fontsize=9)
a.set_title(tx('all 16 L-functions', 'všech 16 L-funkcí'), fontsize=9.5, loc='left', color=INK); a.legend(fontsize=7.5, frameon=False, loc='center right')
fig.suptitle(tx('Pair correlation of the zeros of one L-function at short range, 36,000 < γ ≤ 80,000: measured minus GUE (thin lines: single L-functions)',
                'Párová korelace nul jedné L-funkce na krátké vzdálenosti, 36 000 < γ ≤ 80 000: měření minus GUE (tenké čáry: jednotlivé L-funkce)'), fontsize=11, color=INK, x=0.01, ha='left')
plt.tight_layout(rect=(0, 0, 1, 0.965)); plt.savefig(PRE + '_parova_%s.png' % LG, dpi=150); plt.close()

# ---------- M2: rozestupy sjednoceni proti slozeni nezavislych spekter GUE
HU = np.linspace(0, 4, 41); s = (HU[:-1] + HU[1:]) / 2; pg = R.kose(R.E_gue, HU); po = -np.diff(np.exp(-HU)) / np.diff(HU)
jm = {'U5': '$\\mathbb{Q}(\\zeta_5)$', 'U7': '$\\mathbb{Q}(\\zeta_7)$', 'U15': '$\\mathbb{Q}(\\zeta_{15})$', 'U21': '$\\mathbb{Q}(\\zeta_{21})$'}
fig, axs = plt.subplots(2, 4, figsize=(14.4, 6.6), gridspec_kw=dict(height_ratios=[2.0, 1.0]), sharex='col')
for i, U in enumerate(V['M2']):
    x = V['M2'][U]; p = np.array(x['p']); pS = np.array(x['p_S']); e = np.array(x['chyba']); a = axs[0, i]
    a.plot(s, po, ':', color=MUTED, lw=1.5, label=tx('Poisson', 'Poisson'))
    a.plot(s, pg, '--', color=ORANGE, lw=1.4, label=tx('one GUE spectrum', 'jedno spektrum GUE'))
    a.plot(s, pS, '-', color=INK, lw=1.6, label=tx('%d independent GUE spectra', '%d nezávislých spekter GUE') % x['sektoru'])
    a.plot(s, p, 'o', ms=3.4, color=BLUE, label=tx('zeros of the union', 'nuly sjednocení'))
    a.set_title(tx('%s · %d sectors · %s zeros', '%s · %d sektorů · %s nul') % (jm[U], x['sektoru'], format(x['nul'], ',d').replace(',', ',' if EN else ' ')), fontsize=9.5, loc='left', color=INK)
    a.set_ylim(0, 1.0); a.set_xlim(0, 4); osy(a)
    if i == 0: a.set_ylabel(tx('spacing density p(s)', 'hustota rozestupů p(s)'), color=INK, fontsize=9); a.legend(fontsize=7.5, frameon=False, loc='upper right')
    b = axs[1, i]; b.axhline(0, color=INK, lw=1.0); b.errorbar(s, p - pS, yerr=e, fmt='o', ms=3.0, color=BLUE, ecolor=BLUE, elinewidth=0.8)
    b.set_ylim(-0.015, 0.015); osy(b); b.set_xlabel(tx('spacing s (mean spacings of the union)', 'rozestup s (střední rozestupy sjednocení)'), color=INK, fontsize=9)
    b.text(0.98, 0.06, tx('largest $|p - p_S|$ %s · $\\chi^2$ per bin %s', 'největší $|p - p_S|$ %s · $\\chi^2$ na koš %s') % (cz(x['max_rozdil_od_S'], 4), cz(x['chi2_S'], 1)), transform=b.transAxes, fontsize=7.5, color=INK, ha='right', va='bottom')
    if i == 0: b.set_ylabel(tx('measured − independent', 'měření − nezávislá'), color=INK, fontsize=9)
fig.suptitle(tx('Nearest-neighbour spacings of the zeros of all characters of a modulus (cyclotomic Dedekind zeta), 36,000 < γ ≤ 80,000',
                'Rozestupy sousedních nul všech charakterů modulu (Dedekindova zeta kruhového tělesa), 36 000 < γ ≤ 80 000'), fontsize=11, color=INK, x=0.01, ha='left')
plt.tight_layout(rect=(0, 0, 1, 0.965)); plt.savefig(PRE + '_orchestr_%s.png' % LG, dpi=150); plt.close()
print('ok', PRE, LG)
