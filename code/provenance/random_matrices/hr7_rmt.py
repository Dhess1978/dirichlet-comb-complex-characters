#!/usr/bin/env python3
"""Studie VII · dodatek 9 (RMT): M1 parova korelace a rozestupy jedne L-funkce, M2 rozestupy sjednoceni (PRAVIDLA_STUDIE_VII_DODATEK_9_RMT).
Pouziti:
  python3 -I hr7_rmt.py mereni KS_NULY LMFDB_PLATT VYSTUP.json [--zkouska POCET_OKEN]
  python3 -I hr7_rmt.py vyhodnot MERENI.json TEORIE.json VYSLEDEK.json
Vstupy jen cte. Zadna sit."""
import gzip, hashlib, json, math, os, struct, sys, datetime, numpy as np
T1, T2 = 36000.0, 80000.0; BL = 10
HR7 = ['q5_c1', 'q5_c3', 'q7_c1', 'q7_c2', 'q7_c4', 'q7_c5', 'q13_c1', 'q13_c11', 'q15_c1-1', 'q15_c3-1', 'q21_c1-1', 'q21_c5-1']; HR7B = ['q15_c2-1', 'q21_c2-1', 'q21_c3-1', 'q21_c4-1']
CEKANO = {'q5_c1': 75032, 'q5_c3': 75033, 'q7_c2': 77389, 'q7_c4': 77388, 'q7_c1': 77389, 'q7_c5': 77389, 'q13_c1': 81725, 'q13_c11': 81724, 'q15_c1-1': 82726, 'q15_c3-1': 82726, 'q21_c1-1': 85082, 'q21_c5-1': 85082,
          'q15_c2-1': 82726, 'q21_c3-1': 85082, 'q21_c2-1': 85083, 'q21_c4-1': 85082, 'VI_D-3': 71455, 'VI_D5': 75033, 'VI_D-7': 77389, 'zeta': 63763}     # pocty nul v pasmu z vyhodnoceni 3.11-3.13 (hr7_eval3_final.json)
SJ = {'U5': ['zeta', 'VI_D5', 'q5_c1', 'q5_c3'], 'U7': ['zeta', 'VI_D-7', 'q7_c1', 'q7_c2', 'q7_c4', 'q7_c5'],
      'U15': ['zeta', 'VI_D-3', 'VI_D5', 'q5_c1', 'q5_c3', 'q15_c1-1', 'q15_c3-1', 'q15_c2-1'],
      'U21': ['zeta', 'VI_D-3', 'VI_D-7', 'q7_c1', 'q7_c2', 'q7_c4', 'q7_c5', 'q21_c1-1', 'q21_c5-1', 'q21_c2-1', 'q21_c3-1', 'q21_c4-1']}
def modul(jm): return 1 if jm == 'zeta' else (abs(int(jm[4:])) if jm.startswith('VI_D') else int(jm[1:].split('_')[0]))
cas = lambda: datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S')
HP = np.round(np.arange(0, 3.0001, 0.05), 6); HS = np.round(np.arange(0, 3.0001, 0.1), 6); HU = np.round(np.arange(0, 4.0001, 0.1), 6)

# ---------------------------------------------------------------- nacteni
def okna(slozka, k0, k1, nahrady=None):
    g = []
    for k in range(k0, k1):
        p = (nahrady or {}).get(k) or os.path.join(slozka, 'okno_%06d.json.gz' % k); b = json.loads(gzip.decompress(open(p, 'rb').read())); assert b['n'] == len(b['g']); g.extend(b['g'])
    g = np.array(g, dtype=np.float64); assert np.all(np.diff(g) > 0); return g
def platt_do(path, tmax):
    out = []
    with open(path, 'rb') as f:
        nb = struct.unpack('<Q', f.read(8))[0]; pos = 8
        for i in range(nb):
            f.seek(pos); t0, t1, N0, N1 = struct.unpack('<ddQQ', f.read(32)); n = N1 - N0
            if t0 > tmax: break
            raw = np.frombuffer(f.read(13 * n), dtype=np.dtype([('a', '<u8'), ('b', '<u4'), ('c', 'u1')])); L_ = np.longdouble
            inc = raw['a'].astype(L_) * L_(2.0) ** -101 + raw['b'].astype(L_) * L_(2.0) ** -37 + raw['c'].astype(L_) * L_(2.0) ** -5
            out.append((L_(t0) + np.cumsum(inc)).astype(np.float64)); pos += 32 + 13 * n
    return np.concatenate(out)
def hladky(g, q): return g / (2 * math.pi) * (np.log(q * g / (2 * math.pi)) - 1)

# ---------------------------------------------------------------- statistiky
def parova(w):
    n = int(np.searchsorted(w, w[-1] - 3.0, side='right')); idx = np.arange(n); C = np.array([np.searchsorted(w, w[:n] + e, side='right') - (idx + 1) for e in HP]); d = np.diff(C, axis=0).astype(np.float64)   # (60, n)
    hr = [(n * b) // BL for b in range(BL + 1)]
    return dict(n=n, R2=(d.sum(1) / (n * 0.05)).tolist(), bloky=[(d[:, hr[b]:hr[b + 1]].sum(1) / ((hr[b + 1] - hr[b]) * 0.05)).tolist() for b in range(BL)])
def rozestupy(w, hrany):
    s = np.diff(w); n = len(s); dl = hrany[1] - hrany[0]; hr = [(n * b) // BL for b in range(BL + 1)]
    H = lambda x: (np.histogram(x, bins=hrany)[0] / (len(x) * dl)).tolist()
    return dict(n=n, stred=float(s.mean()), p=H(s), nad=float((s > hrany[-1]).mean()), bloky=[H(s[hr[b]:hr[b + 1]]) for b in range(BL)], nejmensi=float(s.min()))

def mereni(KS, PL, OUT, zk=None):
    k0, k1 = 360, (360 + zk if zk else 800); t2 = k1 * 100.0; G = {}
    VA = os.path.join(KS, 'HR7_nuly_2026_11_20'); VB = os.path.join(KS, 'HR7B_nuly_2026_11_20'); VO = os.path.join(KS, 'HR7B_oprava_2026_11_21'); VR = os.path.join(KS, 'FZ_L_nuly_2026_11_16')
    O = json.load(open(os.path.join(VO, 'oprava.json'))); assert O['v_poradku'] and not O['zkouska']; NAH = {}
    for wk in O['okna']:
        p = os.path.join(VO, wk['soubor']); assert hashlib.sha256(open(p, 'rb').read()).hexdigest() == wk['sha256']; NAH.setdefault(wk['L'], {})[int(wk['k'])] = p
    for jm in HR7: G[jm] = okna(os.path.join(VA, jm), k0, k1)
    for jm in HR7B: G[jm] = okna(os.path.join(VB, jm), k0, k1, NAH.get(jm))
    for D, sl in ((-3, 'D_m3'), (5, 'D_5'), (-7, 'D_m7')): G['VI_D%d' % D] = okna(os.path.join(VR, sl), k0, k1)
    z = platt_do(os.path.join(PL, 'zeros_26000.dat'), t2 + 200); G['zeta'] = z
    for jm in G: G[jm] = G[jm][(G[jm] > T1) & (G[jm] <= t2)]
    pocty = {jm: int(len(G[jm])) for jm in G}
    if not zk: assert pocty == CEKANO, {jm: (pocty[jm], CEKANO[jm]) for jm in G if pocty[jm] != CEKANO[jm]}
    R = dict(cas_UTC=cas(), T1=T1, T2=t2, zkouska=bool(zk), pocty=pocty, pocty_souhlasi=(pocty == CEKANO), opravena_okna=sum(len(v) for v in NAH.values()), skript_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest(), M1={}, M2={})
    for jm in HR7 + HR7B:
        q = modul(jm); w = hladky(G[jm], q); R['M1'][jm] = dict(q=q, nul=pocty[jm], parova=parova(w), rozestupy=rozestupy(w, HS))
    for U, sek in SJ.items():
        g = np.sort(np.concatenate([G[j] for j in sek])); assert np.all(np.diff(g) > 0); W = sum(hladky(g, modul(j)) for j in sek)
        f = lambda T: [math.log(modul(j) * T / (2 * math.pi)) for j in sek]; fm = f((T1 + t2) / 2)
        R['M2'][U] = dict(sektory=sek, moduly=[modul(j) for j in sek], nul=int(len(g)), podily=[x / sum(fm) for x in fm], podily_dole=[x / sum(f(T1)) for x in f(T1)], podily_nahore=[x / sum(f(t2)) for x in f(t2)], rozestupy=rozestupy(W, HU))
    json.dump(R, open(OUT, 'x')); print('mereni hotovo', cas(), 'pocty souhlasi:', R['pocty_souhlasi'], '| sjednoceni', {U: R['M2'][U]['nul'] for U in SJ})

# ---------------------------------------------------------------- modely
def E_gue(x, m=40):
    if x <= 0: return 1.0
    t, w = np.polynomial.legendre.leggauss(m); t = (t + 1) * x / 2; w = w * x / 2
    return float(np.linalg.det(np.eye(m) - np.sqrt(w)[:, None] * np.sinc(t[:, None] - t[None, :]) * np.sqrt(w)[None, :]))
def E_cue(x, N):
    if x <= 0: return 1.0
    if x >= N: return 0.0
    j = np.arange(N); D = j[:, None] - j[None, :]; M = np.where(D == 0, x / N, np.sin(math.pi * D * x / N) / (math.pi * np.where(D == 0, 1, D)))
    return float(np.linalg.det(np.eye(N) - M))
def kose(E, hrany, h=1e-4):
    dE = lambda x: -1.0 if x <= 0 else (E(x + h) - E(x - h)) / (2 * h)
    d = np.array([dE(x) for x in hrany]); return (np.diff(d) / np.diff(hrany))
def vyhodnot(FM, FT, OUT):
    Mz = json.load(open(FM)); Tz = json.load(open(FT)); assert not Mz['zkouska'] and Mz['pocty_souhlasi']; sem = lambda B: np.array(B).std(0, ddof=1) / math.sqrt(len(B))
    R = dict(cas_UTC=cas(), mereni_sha256=hashlib.sha256(open(FM, 'rb').read()).hexdigest(), teorie_sha256=hashlib.sha256(open(FT, 'rb').read()).hexdigest(), M1={}, M2={})
    gue_p = np.array(Tz['GUE']); pg = kose(E_gue, HS); pc = {N: kose(lambda x, N=N: E_cue(x, N), HS) for N in range(2, 7)}
    R['modely'] = dict(norma_GUE_do_3=float(pg.sum() * 0.1), norma_CUE={N: float(pc[N].sum() * 0.1) for N in pc}, p_GUE=pg.tolist())
    S_all = np.zeros(len(HS) - 1); SB_all = np.zeros((BL, len(HS) - 1)); n_all = 0; okB = okG = 0
    for jm, d in Mz['M1'].items():
        q = str(d['q']); bk = np.array(Tz['moduly'][q]['BK']); r = np.array(d['parova']['R2']); e = sem(d['parova']['bloky'])
        cb = float((((r - bk) / e) ** 2).mean()); cg = float((((r - gue_p) / e) ** 2).mean()); a = cb <= 2.0; b = (cg >= 3 * cb) and (cg >= 4.0); okB += a; okG += b
        p = np.array(d['rozestupy']['p']); ep = sem(d['rozestupy']['bloky']); chs = {'GUE': float((((p - pg) / ep) ** 2).mean())}; chs.update({'CUE%d' % N: float((((p - pc[N]) / ep) ** 2).mean()) for N in pc})
        R['M1'][jm] = dict(q=d['q'], nul=d['nul'], chi2_BK=cb, chi2_GUE=cg, BK_v_poradku=bool(a), GUE_horsi=bool(b), rozdil_mereni_minus_GUE_u_0_5=float(r[9:11].mean() - gue_p[9:11].mean()), rozdil_BK_minus_GUE_u_0_5=float(bk[9:11].mean() - gue_p[9:11].mean()),
                           nejvetsi_z_proti_BK=float(np.abs((r - bk) / e).max()), rozestupy_chi2=chs, rozestupy_nejlepsi=min(chs, key=chs.get), N_ef_teorie=Tz['moduly'][q]['N_ef'], rozestupy_max_rozdil_od_GUE=float(np.abs(p - pg).max()))
        S_all += p * d['rozestupy']['n']; SB_all += np.array(d['rozestupy']['bloky']) * d['rozestupy']['n'] / BL; n_all += d['rozestupy']['n']
    pa = S_all / n_all; ea = (SB_all / (n_all / BL)).std(0, ddof=1) / math.sqrt(BL); cha = {'GUE': float((((pa - pg) / ea) ** 2).mean())}; cha.update({'CUE%d' % N: float((((pa - pc[N]) / ea) ** 2).mean()) for N in pc})
    R['M1_dohromady'] = dict(rozestupu=int(n_all), p=pa.tolist(), chyba=ea.tolist(), chi2=cha, nejlepsi=min(cha, key=cha.get), max_rozdil_od_GUE=float(np.abs(pa - pg).max()), kde=float(HS[int(np.abs(pa - pg).argmax())] + 0.05))
    R['M1_verdikt'] = dict(BK_v_poradku=int(okB), GUE_horsi=int(okG), ze=len(Mz['M1']), O=bool(okB >= 14 and okG >= 14))
    pg4 = kose(E_gue, HU); po = -np.diff(np.exp(-HU)) / 0.1; vse = True
    for U, d in Mz['M2'].items():
        f = d['podily']; ES = lambda x, f=f: float(np.prod([E_gue(fi * x) for fi in f])); pS = kose(ES, HU); p = np.array(d['rozestupy']['p']); e = sem(d['rozestupy']['bloky'])
        a = float(np.abs(p - pS).max()); b = float(abs(p[0] - po[0])); c = float(np.abs(p - pg4).max()); ok = (a <= 0.015) and (b >= 0.03) and (c >= 0.2); vse = vse and ok; z = (p - pS) / e
        fd, fh = d['podily_dole'], d['podily_nahore']; pSd = kose(lambda x: float(np.prod([E_gue(fi * x) for fi in fd])), HU); pSh = kose(lambda x: float(np.prod([E_gue(fi * x) for fi in fh])), HU)
        R['M2'][U] = dict(sektoru=len(f), nul=d['nul'], max_rozdil_od_S=a, kde=float(HU[int(np.abs(p - pS).argmax())] + 0.05), prvni_kos=float(p[0]), prvni_kos_S=float(pS[0]), prvni_kos_Poisson=float(po[0]), rozdil_od_Poissona_prvni_kos=b, max_rozdil_od_jednoho_GUE=c,
                          chi2_S=float((z ** 2).mean()), chi2_Poisson=float((((p - po) / e) ** 2).mean()), kosu_nad_3_chyby=int((np.abs(z) > 3).sum()), nejvetsi_z=float(np.abs(z).max()), stred_rozestupu=d['rozestupy']['stred'], nad_4=d['rozestupy']['nad'], S_nad_4=float(1 - pS.sum() * 0.1),
                          citlivost_S_na_vysku=float(np.abs(pSd - pSh).max()), jeden_minus_soucet_f2=float(1 - sum(x * x for x in f)), splneno=bool(ok), p=p.tolist(), p_S=pS.tolist(), chyba=e.tolist())
    R['M2_verdikt'] = dict(O=bool(vse), splneno={U: R['M2'][U]['splneno'] for U in R['M2']})
    json.dump(R, open(OUT, 'x'), ensure_ascii=False, indent=1)
    print('M1', R['M1_verdikt'], '| chi2 BK', [round(R['M1'][j]['chi2_BK'], 2) for j in R['M1']], '| chi2 GUE', [round(R['M1'][j]['chi2_GUE'], 1) for j in R['M1']])
    print('M1 rozestupy dohromady', {k: round(v, 2) for k, v in cha.items()}, 'max rozdil od GUE %.4f u s = %.2f' % (R['M1_dohromady']['max_rozdil_od_GUE'], R['M1_dohromady']['kde']), '| normy', round(R['modely']['norma_GUE_do_3'], 5), {N: round(v, 5) for N, v in R['modely']['norma_CUE'].items()})
    for U in R['M2']: x = R['M2'][U]; print('M2', U, 'sektoru', x['sektoru'], 'max |p - p_S| %.4f u %.2f' % (x['max_rozdil_od_S'], x['kde']), 'prvni kos %.4f (S %.4f, Poisson %.4f)' % (x['prvni_kos'], x['prvni_kos_S'], x['prvni_kos_Poisson']), 'chi2 S %.2f Poisson %.1f' % (x['chi2_S'], x['chi2_Poisson']), 'nejvetsi z %.1f' % x['nejvetsi_z'], x['splneno'])
    print('M2', R['M2_verdikt'])
if __name__ == '__main__':
    if sys.argv[1] == 'mereni': mereni(sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[6]) if len(sys.argv) > 6 and sys.argv[5] == '--zkouska' else None)
    else: vyhodnot(sys.argv[2], sys.argv[3], sys.argv[4])
