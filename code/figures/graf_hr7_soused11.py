"""Studie VII · dodatek 8: graf „co orchestr slyší kolem 11^5 a 11^6“ z výstupu hr7_soused11.py (jen čte JSON).
Použití: python3 graf_hr7_soused11.py VYSTUP_hr7_soused11.json VYSTUP.png"""
import json, sys, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
E = json.load(open(sys.argv[1])); INK, MUTED, GRID = '#22252e', '#8c91a6', '#d9dbe3'; BLUE, GREEN, ORANGE = '#3987e5', '#199e70', '#d95926'
fig, axs = plt.subplots(2, 1, figsize=(12, 8.6))
for a, k in zip(axs, ('5', '6')):
    m = E['mrizka'][k]; d = np.array(m['d']); A = m['A']; Z = np.array([complex(*v) for v in m['Z']]); P = np.array([complex(*v) for v in m['P']]); Pb = np.array([complex(*v) for v in m['P_bez_x0']]); R = Z - Pb
    krok = 8 if k == '5' else 8
    a.plot(d, np.abs(P) / A, color=GREEN, lw=2, label='předpověď z prvočísel (explicitní vzorec, bez volného parametru)', zorder=2)
    a.plot(d[::krok], np.abs(Z[::krok]) / A, 'o', ms=3.6, color=BLUE, mec='white', mew=0.4, label='měření z nul (12 L-funkcí, 947 573 nul)', zorder=3)
    a.plot(d, np.abs(R) / A, color=ORANGE, lw=1.8, label='měření po odečtení sousedních prvočísel', zorder=2)
    ymax = float(np.max(np.abs(Z)) / A) * 1.18; a.set_ylim(0, ymax)
    pr = [p for p in m['prvocisla'] if d[0] <= p['d'] <= d[-1]]
    for p in pr:
        a.plot([p['d']] * 2, [0, ymax * 0.045], color=INK, lw=1.1, zorder=4)
        if k == '5': a.annotate('%d' % p['n'], (min(max(p['d'], d[0] + 8), d[-1] - 8), ymax * 0.06), ha='center', va='bottom', fontsize=7.5, color=INK, rotation=90)
    a.axvline(0, color=MUTED, lw=0.9, ls=':'); a.axhline(1, color=MUTED, lw=0.8, ls=':')
    b = E['body'][k]; zb = b['zbytek_po_odecteni_sousedu']
    if k == '5':
        a.annotate('11⁵ = 161 051\nsama mlčí (zbytek %.0e)' % abs(complex(*zb)), (0, 0.0), xytext=(-150, ymax * 0.78), fontsize=9, color=INK, arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.8))
        a.annotate('soused 161 071\n(o 20 dál, ≡ 1 mod 21)', (20, float(np.abs(Z[d == 20][0]) / A)), xytext=(110, ymax * 0.80), fontsize=9, color=INK, arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.8))
        a.set_title('kolem 11⁵: jedenáctka tu neleží v tělese (11⁵ ≡ 2 mod 21); naměřený „zub“ je ozvěna sousedního prvočísla · šířka rozlišení ±%.0f' % b['sirka_jadra'], fontsize=9.5, loc='left', color=INK)
    else:
        a.annotate('zub jedenáctky po odečtení sousedů:\n%.4f plné síly (čeká se 1)' % abs(complex(*zb)), (0, 1.0), xytext=(250, 1.9), fontsize=9, color=INK, arrowprops=dict(arrowstyle='-', color=MUTED, lw=0.8))
        a.set_title('kolem 11⁶ = 1 771 561 (≡ 1 mod 21): zub jedenáctky je schovaný pod %d sousedními prvočísly (čárky dole) · šířka rozlišení ±%.0f' % (len(pr), b['sirka_jadra']), fontsize=9.5, loc='left', color=INK)
    a.set_xlim(d[0], d[-1]); a.set_ylabel('|Z(x)| / plná síla jedenáctky', color=INK); a.grid(axis='y', color=GRID, lw=0.6); a.set_xlabel('x − 11^%s' % k, color=INK)
    for s_ in ('top', 'right'): a.spines[s_].set_visible(False)
axs[0].legend(fontsize=8, loc='upper left', frameon=False)
fig.suptitle('Proč čtení řádu 11 v Q(ζ₂₁) selhalo: sousední prvočísla a mez rozlišení (shoda měření s předpovědí %.0e plné síly)' % max(E['mrizka'][k]['max_rozdil_k_plne_sile'] for k in ('5', '6')), fontsize=11, color=INK, x=0.01, ha='left')
plt.tight_layout(); plt.savefig(sys.argv[2], dpi=130); print('ok')
