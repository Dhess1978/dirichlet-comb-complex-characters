#!/usr/bin/env python3
"""Studie VII · dodatek 10 (M3): slabne odchylka parove korelace od GUE s vyskou podle teorie? (PRAVIDLA_STUDIE_VII_DODATEK_10_RMT_VYSKA)
Pouziti: python3 -I hr7_rmt_vyska.py HR7_RMT.py SLOZKA_HR7C MERENI_M1.json TEORIE_036.json TEORIE_400.json VYSTUP.json
Vstupy jen cte. Otisk tohoto skriptu se uklada drive, nez se nuly vyssiho pasma otevrou."""
import gzip, hashlib, json, math, os, sys, numpy as np
PYR, V, FM, FT0, FT1, OUT = sys.argv[1:7]; argv0 = sys.argv; sys.argv = ['x']; exec(open(PYR).read().split("if __name__ == '__main__':")[0]); sys.argv = argv0
TA, TB = 400000.0, 420000.0; K0, K1 = 4000, 4200; LF = ['q15_c2-1', 'q21_c2-1', 'q21_c3-1', 'q21_c4-1']
M1 = json.load(open(FM)); T0 = json.load(open(FT0)); Th = json.load(open(FT1)); assert T0['T1'] == 36000.0 and Th['T1'] == TA and Th['T2'] == TB
sem = lambda B: np.array(B).std(0, ddof=1) / math.sqrt(len(B)); gue = np.array(T0['GUE']); R = dict(cas_UTC=cas(), skript_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest(), L={}, uplnost={})
# ---------------------------------------------------------------- nacteni a uplnost
G = {}; ok_u = True
for jm in LF:
    g = []; kum = 0.0; rada = []
    for k in range(K0, K1):
        b = json.loads(gzip.decompress(open(os.path.join(V, jm, 'okno_%06d.json.gz' % k), 'rb').read())); assert b['n'] == len(b['g']) and abs(b['a'] - k * 100.0) < 1e-9
        g.extend(b['g']); kum += b['n'] - b['ocek']; rada.append(kum)
    g = np.array(g, dtype=np.float64); assert np.all(np.diff(g) > 0); g = g[(g > TA) & (g <= TB)]; q = modul(jm); hl = float(hladky(np.array([TB]), q)[0] - hladky(np.array([TA]), q)[0])
    u = dict(nul=int(len(g)), hladky_pocet=hl, rozdil=float(len(g) - hl), kumulativni_konec=float(rada[-1]), kumulativni_rozpeti=float(max(rada) - min(rada)), nejmensi_rozestup=float(np.diff(hladky(g, q)).min()))
    u['v_poradku'] = bool(abs(u['kumulativni_konec']) <= 2.5 and u['kumulativni_rozpeti'] <= 4.0 and abs(u['rozdil']) <= 2.5); ok_u = ok_u and u['v_poradku']; R['uplnost'][jm] = u; G[jm] = g
R['uplnost_v_poradku'] = bool(ok_u)
if not ok_u:
    json.dump(R, open(OUT, 'x'), ensure_ascii=False, indent=1); print('UPLNOST NENI V PORADKU - vyhodnoceni se nespousti:', json.dumps(R['uplnost'])); sys.exit(3)
# ---------------------------------------------------------------- amplitudy
def amplituda(r, bloky, bk, n):
    e = sem(bloky); e = np.where(e > 0, e, np.sqrt(np.maximum(bk, 1e-4) / (n / len(bloky) * 0.05) / len(bloky)))   # pojistka: kos s nulovym rozptylem bloku dostane Poissonuv odhad
    x = bk - gue; w = x / e ** 2; jm_ = float((x * w).sum()); A = float(((r - gue) * w).sum() / jm_); Ab = [float(((np.array(b) - gue) * w).sum() / jm_) for b in bloky]
    return A, float(np.std(Ab, ddof=1) / math.sqrt(len(Ab))), float(np.sqrt((x * x).mean())), float((((r - bk) / e) ** 2).mean()), float((((r - gue) / e) ** 2).mean())
for jm in LF:
    q = str(modul(jm)); w = hladky(G[jm], modul(jm)); P = parova(w); S = rozestupy(w, HS)
    Ah, sAh, rmh, cbh, cgh = amplituda(np.array(P['R2']), P['bloky'], np.array(Th['moduly'][q]['BK']), P['n'])
    d = M1['M1'][jm]; Al, sAl, rml, cbl, cgl = amplituda(np.array(d['parova']['R2']), d['parova']['bloky'], np.array(T0['moduly'][q]['BK']), d['parova']['n'])
    pg = kose(E_gue, HS)
    R['L'][jm] = dict(q=int(q), nul_nahore=int(len(G[jm])), nul_dole=d['nul'], A_nahore=Ah, chyba_A_nahore=sAh, A_dole=Al, chyba_A_dole=sAl, rms_teorie_nahore=rmh, rms_teorie_dole=rml, V_nahore=Ah * rmh, V_dole=Al * rml,
                       chi2_BK_nahore=cbh, chi2_GUE_nahore=cgh, chi2_BK_dole=cbl, chi2_GUE_dole=cgl, rozestupy_max_rozdil_od_GUE_nahore=float(np.abs(np.array(S['p']) - pg).max()), rozestupy_max_rozdil_od_GUE_dole=float(np.abs(np.array(d['rozestupy']['p']) - pg).max()),
                       R2_nahore=P['R2'], p_nahore=S['p'])
def sdruz(klic_h, klic_s, nasob=None):
    x = np.array([R['L'][j][klic_h] * (R['L'][j][nasob] if nasob else 1.0) for j in LF]); s = np.array([R['L'][j][klic_s] * (R['L'][j][nasob] if nasob else 1.0) for j in LF]); w = 1 / s ** 2
    return float((x * w).sum() / w.sum()), float(1 / math.sqrt(w.sum())), w / w.sum()
Ah, sAh, _ = sdruz('A_nahore', 'chyba_A_nahore'); Al, sAl, _ = sdruz('A_dole', 'chyba_A_dole')
Vh, sVh, wh = sdruz('A_nahore', 'chyba_A_nahore', 'rms_teorie_nahore'); Vl, sVl, wl = sdruz('A_dole', 'chyba_A_dole', 'rms_teorie_dole')
th = float(sum(w * R['L'][j]['rms_teorie_nahore'] for w, j in zip(wh, LF))); tl = float(sum(w * R['L'][j]['rms_teorie_dole'] for w, j in zip(wl, LF)))
pom = Vh / Vl; spom = pom * math.hypot(sVh / Vh, sVl / Vl); pomt = th / tl
a = sum(R['L'][j]['chi2_BK_nahore'] <= 2.0 for j in LF) >= 3; b = abs(Ah - 1) <= 3 * sAh; c = (Vl - Vh) >= 3 * math.hypot(sVh, sVl); d_ = abs(pom - pomt) <= 3 * spom
verd = 'O' if (a and b and c and d_) else ('Z' if (not a or not b or not d_) else 'N')
R['sdruzeno'] = dict(A_nahore=Ah, chyba_A_nahore=sAh, A_dole=Al, chyba_A_dole=sAl, V_nahore=Vh, chyba_V_nahore=sVh, V_dole=Vl, chyba_V_dole=sVl, rozdil_V=Vl - Vh, chyba_rozdilu=math.hypot(sVh, sVl), pomer_V=pom, chyba_pomeru=spom,
                     pomer_teorie=pomt, teorie_V_nahore=th, teorie_V_dole=tl, kriteria=dict(a=bool(a), b=bool(b), c=bool(c), d=bool(d_)), verdikt=verd)
json.dump(R, open(OUT, 'x'), ensure_ascii=False, indent=1)
print('uplnost', {j: (R['uplnost'][j]['nul'], round(R['uplnost'][j]['rozdil'], 2), round(R['uplnost'][j]['kumulativni_rozpeti'], 2)) for j in LF})
for j in LF: x = R['L'][j]; print(j, 'A nahore %.3f +- %.3f, dole %.3f +- %.3f | chi2 BK nahore %.2f GUE %.1f | rozestupy max od GUE nahore %.3f dole %.3f' % (x['A_nahore'], x['chyba_A_nahore'], x['A_dole'], x['chyba_A_dole'], x['chi2_BK_nahore'], x['chi2_GUE_nahore'], x['rozestupy_max_rozdil_od_GUE_nahore'], x['rozestupy_max_rozdil_od_GUE_dole']))
s = R['sdruzeno']; print('sdruzeno: A nahore %.3f +- %.3f | V nahore %.5f +- %.5f, dole %.5f +- %.5f | pomer %.3f +- %.3f (teorie %.3f) | kriteria %s | VERDIKT %s' % (s['A_nahore'], s['chyba_A_nahore'], s['V_nahore'], s['chyba_V_nahore'], s['V_dole'], s['chyba_V_dole'], s['pomer_V'], s['chyba_pomeru'], s['pomer_teorie'], s['kriteria'], s['verdikt']))
