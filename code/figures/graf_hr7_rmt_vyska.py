"""Studie VII · graf M3 (dodatek 10): parova korelace minus GUE ve dvou vyskach, amplituda A a velikost V proti teorii.
Jen cte JSON: mereni M1 (dolni pasmo), data horniho pasma (hr7_rmt_vyska_obr.py), teorie obou pasem, vysledek zmrazeneho vyhodnoceni.
Pouziti: python3 -I graf_hr7_rmt_vyska.py MERENI_M1.json OBR_HORNI.json TEORIE_036.json TEORIE_400.json VYSLEDEK_VYSKA.json PREFIX_PNG en|cs"""
import json, math, sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
FM, FO, FT0, FT1, FV, PRE, LG = sys.argv[1:8]; EN = LG == 'en'
M = json.load(open(FM)); O = json.load(open(FO)); T0 = json.load(open(FT0)); T1 = json.load(open(FT1)); V = json.load(open(FV)); assert O['shoda_vse']
INK, MUTED, GRID = '#22252e', '#8c91a6', '#d9dbe3'; BLUE, GREEN, ORANGE = '#3987e5', '#199e70', '#d95926'
CONREY = {'q15_c2-1': '15.14', 'q21_c3-1': '21.20', 'q21_c2-1': '21.2', 'q21_c4-1': '21.11'}; LF = list(V['L'])
def cz(x, nd):
    s = ('%.' + str(nd) + 'f') % x; return s if EN else s.replace('.', ',')
def tx(en, cs): return en if EN else cs
def osy(a):
    for s_ in ('top', 'right'): a.spines[s_].set_visible(False)
    a.grid(color=GRID, lw=0.6); a.tick_params(colors=INK, labelsize=8)
sem = lambda B: np.array(B).std(0, ddof=1) / math.sqrt(len(B))
hr = np.array(T0['hranice']); d = (hr[:-1] + hr[1:]) / 2; gue = np.array(T0['GUE'])
PAS = [(tx('36,000 < γ ≤ 80,000', '36 000 < γ ≤ 80 000'), lambda j: M['M1'][j]['parova'], T0, BLUE, 'dole'), (tx('400,000 < γ ≤ 420,000', '400 000 < γ ≤ 420 000'), lambda j: O['L'][j]['parova'], T1, GREEN, 'nahore')]
fig, axs = plt.subplots(2, 3, figsize=(13.6, 7.4))
for i, q in enumerate(['15', '21']):
    js = [j for j in LF if str(V['L'][j]['q']) == q]
    for k, (nazev, zdroj, T, barva, kl) in enumerate(PAS):
        a = axs[i, k]; bk = np.array(T['moduly'][q]['BK']); r = np.array([zdroj(j)['R2'] for j in js]); e = np.array([sem(zdroj(j)['bloky']) for j in js])
        if len(js) > 1:
            for rr in r: a.plot(d, rr - gue, '-', color=barva, lw=0.6, alpha=0.35, zorder=2)
        a.axhline(0, color=ORANGE, lw=1.4, ls='--', zorder=2.4, label=tx('GUE alone', 'samotné GUE'))
        a.plot(d, bk - gue, '-', color=INK, lw=1.6, zorder=3, label=tx('Bogomolny–Keating [8] − GUE', 'Bogomolny–Keating [8] − GUE'))
        a.errorbar(d, r.mean(0) - gue, yerr=np.sqrt((e ** 2).sum(0)) / len(js), fmt='o', ms=3.2, color=barva, ecolor=barva, elinewidth=0.8, zorder=4,
                   label=(tx('zeros: mean of %d L-functions', 'nuly: průměr %d L-funkcí') % len(js)) if len(js) > 1 else tx('zeros of L(s, χ), χ = %s', 'nuly L(s, χ), χ = %s') % CONREY[js[0]])
        cb = [V['L'][j]['chi2_BK_' + kl] for j in js]; rm = V['L'][js[0]]['rms_teorie_' + kl]
        a.set_title('q = %s · %s' % (q, nazev), fontsize=9.5, loc='left', color=INK)
        a.text(0.02, 0.97, tx('rms of the prediction %s · $\\chi^2$ per bin against [8]: %s', 'rms předpovědi %s · $\\chi^2$ na koš proti [8]: %s') % (cz(rm, 4), cz(cb[0], 2) if len(cb) == 1 else cz(min(cb), 2) + '–' + cz(max(cb), 2)),
               transform=a.transAxes, fontsize=7.6, color=INK, ha='left', va='top')
        a.set_ylim(-0.10, 0.10); a.set_xlim(0, 3); osy(a); a.set_xlabel(tx('distance d (mean spacings)', 'vzdálenost d (střední rozestupy)'), color=INK, fontsize=9); a.set_ylabel('$R_2(d)$ − GUE', color=INK, fontsize=9)
        if k == 0: a.legend(fontsize=7.2, frameon=False, loc='lower right')
S = V['sdruzeno']
a = axs[0, 2]; x = np.arange(len(LF) + 1); sir = 0.16; lo, hi = 0.5, 1.5
for kl, barva, pos, naz in (('dole', BLUE, -sir, tx('36,000–80,000', '36 000–80 000')), ('nahore', GREEN, sir, tx('400,000–420,000', '400 000–420 000'))):
    y = [V['L'][j]['A_' + kl] for j in LF] + [S['A_' + kl]]; e = [V['L'][j]['chyba_A_' + kl] for j in LF] + [S['chyba_A_' + kl]]
    a.errorbar(x + pos, y, yerr=e, fmt='o', ms=5, color=barva, ecolor=barva, elinewidth=1.0, label=naz); lo = min(lo, min(np.array(y) - np.array(e)) - 0.25); hi = max(hi, max(np.array(y) + np.array(e)) + 0.1)
a.axhline(1, color=INK, lw=1.2); a.axvline(len(LF) - 0.5, color=GRID, lw=1.0)
a.set_xticks(x); a.set_xticklabels([CONREY[j] for j in LF] + [tx('pooled', 'sdruženo')], fontsize=8); osy(a); a.set_ylim(lo, hi)
a.set_ylabel(tx('amplitude A of the measured deviation\nalong the predicted shape', 'amplituda A měřené odchylky\nve tvaru předpovědi'), color=INK, fontsize=9)
a.set_title(tx('amplitude (1 = prediction) · pooled, upper band %s ± %s', 'amplituda (1 = předpověď) · sdruženo nahoře %s ± %s') % (cz(S['A_nahore'], 3), cz(S['chyba_A_nahore'], 3)), fontsize=8.6, loc='left', color=INK)
a.legend(fontsize=7.5, frameon=False, loc='lower left', ncol=2, title=tx('height γ', 'výška γ'), title_fontsize=7.5)
a = axs[1, 2]
for xx, kl, barva in ((0, 'dole', BLUE), (1, 'nahore', GREEN)):
    a.hlines(S['teorie_V_' + kl], xx - 0.24, xx + 0.24, color=ORANGE, lw=2.4, zorder=2)
    a.errorbar([xx], [S['V_' + kl]], yerr=[S['chyba_V_' + kl]], fmt='o', ms=6.5, color=barva, ecolor=barva, elinewidth=1.2, zorder=4)
    a.text(xx + 0.28, S['V_' + kl], cz(S['V_' + kl], 4) + ' ± ' + cz(S['chyba_V_' + kl], 4), fontsize=7.6, color=INK, va='center', ha='left')
a.set_xticks([0, 1]); a.set_xticklabels([tx('36,000–80,000', '36 000–80 000'), tx('400,000–420,000', '400 000–420 000')], fontsize=8); a.set_xlim(-0.6, 1.9); osy(a)
a.set_ylim(0, 1.3 * max(S['V_dole'], S['teorie_V_dole'], S['V_nahore'], S['teorie_V_nahore'])); a.set_ylabel(tx('size V of the deviation from GUE (rms)', 'velikost V odchylky od GUE (rms)'), color=INK, fontsize=9); a.set_xlabel(tx('height γ', 'výška γ'), color=INK, fontsize=9)
a.set_title(tx('ratio upper / lower: %s ± %s · predicted %s', 'poměr nahoře / dole: %s ± %s · předpověď %s') % (cz(S['pomer_V'], 3), cz(S['chyba_pomeru'], 3), cz(S['pomer_teorie'], 3)), fontsize=8.6, loc='left', color=INK)
a.legend(handles=[Line2D([], [], marker='o', ls='', ms=6, color=INK, label=tx('measured (pooled over 4 L-functions)', 'měření (sdruženo přes 4 L-funkce)')), Line2D([], [], color=ORANGE, lw=2.4, label=tx('prediction [8]', 'předpověď [8]'))], fontsize=7.5, frameon=False, loc='lower left')
fig.suptitle(tx('Does the deviation of the pair correlation from GUE fade with height as predicted? Four L-functions at two heights',
                'Slábne odchylka párové korelace od GUE s výškou tak, jak říká předpověď? Čtyři L-funkce ve dvou výškách'), fontsize=11, color=INK, x=0.01, ha='left')
plt.tight_layout(rect=(0, 0, 1, 0.965)); plt.savefig(PRE + '_vyska_%s.png' % LG, dpi=150); plt.close(); print('ok', PRE, LG)
