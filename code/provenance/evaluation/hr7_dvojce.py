"""Studie VII · NEZÁVISLÉ DVOJČE hlavních měření (bez společného kódu s hr7_eval.py): hřeben V(δ) a první Fourierův koeficient c pro každou L-funkci,
porovnání s výstupem hr7_eval.py. Jiná cesta: okno po okně (výběr nul searchsorted místo bincount), fáze γ·(ln x) přes rozdíl γ − T₀ a fázi T₀·ln x v mpmath,
c přes ⟨p^{iγ}⟩ s fází počítanou v mpmath po blocích (ne přes zlomkovou část), chyba se nepočítá (jen V a c).
Použití: python3 hr7_dvojce.py SLOZKA_HR7 VYSTUP_hr7_eval.json VYSTUP.json [--zkouska T1]"""
import gzip, json, math, os, sys, time, numpy as np, mpmath as mp
t0 = time.time(); VN, EV, OUT = sys.argv[1:4]; T1 = float(sys.argv[sys.argv.index('--zkouska') + 1]) if '--zkouska' in sys.argv else 36000.0
E = json.load(open(EV)); T2 = float(E['T2']); H = 20.0; DEL = np.round(np.arange(-2.5, 2.5 + 1e-9, 0.025), 3); mp.mp.dps = 30
MN = json.load(open(os.path.join(VN, 'manifest.json'))); OK = float(MN['okno']); k0 = int(round(T1 / OK)); k1 = int(round(T2 / OK)); out = dict(T1=T1, T2=T2, L={})
for jm in E['L']:
    B = MN['L'][jm]['bloky']; q = int(MN['L'][jm]['popis']['q']); g = []
    for k in range(k0, k1): g.extend(json.loads(gzip.decompress(open(os.path.join(VN, jm, B[str(k)]['soubor']), 'rb').read()))['g'])
    g = np.array(sorted(g)); TT = T1 + 9 * H; Mw = int((T2 - 9 * H - TT) // H) + 1; S = np.zeros((len(DEL), Mw), complex); N2 = np.zeros(Mw)
    for i in range(Mw):
        T0 = TT + i * H; a = np.searchsorted(g, T0 - 9.5 * H); b = np.searchsorted(g, T0 + 9.5 * H); d = g[a:b] - T0
        # stejný výběr nul jako zaokrouhlení na nejbližší střed v hr7_eval.py: nula patří oknu, je-li její nejbližší střed vzdálen nejvýš 9 kroků
        mm = np.rint((g[a:b] - TT) / H).astype(np.int64); sel = np.abs(mm - i) <= 9; d = d[sel]; w = np.exp(-0.5 * (d / H) ** 2); N2[i] = float(np.sum(w * w))
        lx0 = mp.log(mp.mpf(q) * mp.mpf(T0) / (2 * mp.pi)); gg = g[a:b][sel]
        for j, dl in enumerate(DEL):
            # γ·(lx0 + δ) = (γ − T0)(lx0 + δ) + T0(lx0 + δ); druhý člen v mpmath mod 2π
            f0 = float(mp.fmod(mp.mpf(T0) * (lx0 + mp.mpf(float(dl))), 2 * mp.pi)); ph = d * float(lx0 + mp.mpf(float(dl))) + f0
            S[j, i] = complex(np.sum(w * np.cos(ph)), np.sum(w * np.sin(ph)))
    V = (np.abs(S - S.mean(1, keepdims=True)) ** 2).mean(1) / N2.mean(); Ve = np.array(E['hreben'][jm]['V'])
    rc = {}
    for p in (2, 3, 5, 7, 11, 13):
        if q % p == 0: continue
        lp = mp.log(p); acc = 0j
        for a in range(0, len(g), 20000):
            gb = g[a:a + 20000]; ref = mp.mpf(float(gb[0])); f0 = float(mp.fmod(ref * lp, 2 * mp.pi)); ph = (gb - float(gb[0])) * float(lp) + f0; acc += complex(np.sum(np.cos(ph)), np.sum(np.sin(ph)))
        c = acc / len(g); ce = complex(*E['faze']['%s,%d' % (jm, p)]['c']); rc[str(p)] = dict(c=[c.real, c.imag], rozdil_proti_eval=abs(c - ce))
    out['L'][jm] = dict(nul=len(g), max_rel_rozdil_V=float(np.max(np.abs(V - Ve) / Ve)), c=rc, max_rozdil_c=max(r['rozdil_proti_eval'] for r in rc.values()))
    print('%-10s V: největší relativní rozdíl proti hr7_eval %.1e | c: největší rozdíl %.1e | %.0f s' % (jm, out['L'][jm]['max_rel_rozdil_V'], out['L'][jm]['max_rozdil_c'], time.time() - t0), flush=True)
json.dump(out, open(OUT, 'w'))
