"""Studie VII · vyhodnocení doplňku (PRAVIDLA_STUDIE_VII.md, dodatek 2, oddíly 3.6, 3.7, 3.8, 3.10). Zmrazí se otiskem před přenosem dat hr7.
3.6 odv-hr7-teleso  — sjednocení všech charakterů modulu (Dedekindova ζ tělesa Q(ζ_q)): c_k(p) = −w_U·φ(q)·[p^k ≡ 1]·ln p/p^{k/2} (hodnoceno x = p^k ≤ 200, dodatek 3), histogramy p = 2, 3.
3.7 odv-hr7-krizova — křížová kovariance dvou L-funkcí téhož modulu v diagonálním režimu proti prvočíselné stopě s χ_a·conj χ_b; popisně C_ab(δ) na mřížce hřebene.
3.8 odv-hr7-nasyceni — Σ²(ℓ) proti Σ² z dvoubodové funkce BK bez p | q (dodatek 3; hr7_sigma_bk.py), GUE a naivní ζ (BK s q = 1); Berryho diagonální vzorec popisně; stejný modul bez teorie; popisně sjednocení.
3.10 odv-hr7-okno-L  — Hannovo okno a = 3H po 6H, hřeben proti V_BK,q s autokorelací okna, kovariance jackknifem po oknech, Hotelling.
Použití: python3 hr7_eval2.py SLOZKA_HR7 SLOZKA_REF_VI VYSTUP.json [--zkouska T1 T2]
  SLOZKA_REF_VI: nuly χ_D studie VI (formát v1 nebo v2) s pokrytím pásma; nuly ζ z /home/claude/br69/data (zeros_14, zeros_5000, zeros_26000)."""
import gzip, hashlib, json, math, os, sys, time, numpy as np, mpmath as mp
from scipy import stats
from scipy.special import loggamma
t0 = time.time(); PY = os.path.dirname(os.path.abspath(__file__)) + '/'; sys.path.insert(0, PY)
import hr7_berry as HB
VN, VREF, OUT = sys.argv[1:4]; ZK = '--zkouska' in sys.argv
if ZK: i_ = sys.argv.index('--zkouska'); T1, T2Z = float(sys.argv[i_ + 1]), float(sys.argv[i_ + 2])
else: T1, T2Z = 36000.0, None
H = 20.0; A_ = 3 * H; NB = 50; KMAX = 60; NNUL = 200; rng = np.random.default_rng(341072); mp.mp.dps = 20
ZF = ['/home/claude/br69/data/zeros_14_t.npy', '/home/claude/br69/data/zeros_5000_t.npy', '/home/claude/br69/data/zeros_26000_t.npy']
import cypari2
PARI = cypari2.Pari()
def chi_tab(q, chi):
    G = 'znstar(%d,1)' % q; c = '[' + ','.join(str(int(x)) for x in chi) + ']'; T = np.zeros(q, complex)
    for a in range(q):
        if math.gcd(a, q) == 1: T[a] = complex(mp.expjpi(2 * mp.mpf(str(PARI('chareval(%s,%s,%d)' % (G, c, a))))))
    return T, int(PARI('zncharisodd(%s,%s)' % (G, c)))
def lam(n):
    if n < 2: return 0.0
    for p in range(2, int(n ** 0.5) + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return math.log(p) if m == 1 else 0.0
    return math.log(n)
def slozka_D(D): return 'D_%s%d' % ('m' if D < 0 else '', abs(D))
# ---------------------------------------------------------------- načtení (stejné kontroly jako hr7_eval.py)
MN = json.load(open(os.path.join(VN, 'manifest.json'))); assert MN['format'] == 'MAYAN_V1008_FZ_NULY_CHI_v3'
OK = float(MN['okno']); k0 = int(round(T1 / OK)); assert ZK or (abs(float(MN['tmin']) - 36000.0) < 1e-9 and MN.get('sada') == 'hr7')
KD = {}
for jm in MN['L']:
    B = MN['L'][jm]['bloky']; k = k0
    while str(k) in B and os.path.exists(os.path.join(VN, jm, B[str(k)]['soubor'])): k += 1
    KD[jm] = k
K = min(KD.values()); T2 = K * OK
if ZK: K = int(round(T2Z / OK)); T2 = T2Z; assert all(v >= K for v in KD.values())
assert ZK or T2 >= 60000.0
G = {}; INFO = {}
for jm in MN['L']:
    P_ = MN['L'][jm]['popis']; q = int(P_['q']); Tb, odd = chi_tab(q, P_['chi']); B = MN['L'][jm]['bloky']; g = []; sha = 0
    for k in range(k0, K):
        raw = open(os.path.join(VN, jm, B[str(k)]['soubor']), 'rb').read(); sha += hashlib.sha256(raw).hexdigest() == B[str(k)]['sha256']; g.extend(json.loads(gzip.decompress(raw))['g'])
    g = np.array(g); assert sha == K - k0 and np.all(np.diff(g) > 0) and g[0] > T1 and g[-1] <= T2
    G[jm] = g; INFO[jm] = dict(q=q, Tab=Tb, odd=odd)
MR = json.load(open(os.path.join(VREF, 'manifest.json'))); assert MR['format'] in ('MAYAN_V1008_FZ_NULY_L_v1', 'MAYAN_V1008_FZ_NULY_L_v2')
kr0 = int(round(float(MR.get('tmin', 0.0)) / float(MR['okno'])))
for D in (5, -7, 13):
    if abs(D) not in [INFO[j]['q'] for j in list(INFO)]: continue
    B = MR['L'][str(D)]['bloky']; g = []; sha = 0; OKr = float(MR['okno'])
    for k in range(int(round(T1 / OKr)), int(round(T2 / OKr))):
        raw = open(os.path.join(VREF, slozka_D(D), B[str(k)]['soubor']), 'rb').read(); sha += hashlib.sha256(raw).hexdigest() == B[str(k)]['sha256']; g.extend(json.loads(gzip.decompress(raw))['g'])
    g = np.array([v for v in g if T1 < v <= T2]); assert np.all(np.diff(g) > 0)
    q = abs(D); Tb = np.array([float(PARI.kronecker(D, a)) for a in range(q)], complex); jm = 'VI_D%d' % D; G[jm] = g; INFO[jm] = dict(q=q, Tab=Tb, odd=int(D < 0))
z = np.concatenate([np.load(f) for f in ZF]); G['zeta'] = z[(z > T1) & (z <= T2)]; del z; INFO['zeta'] = dict(q=1, Tab=np.ones(1, complex), odd=0)
JM = [j for j in G if j != 'zeta']
out = dict(T1=T1, T2=T2, zkouska=ZK, L=JM, nul={j: int(len(G[j])) for j in G}, skript_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())
print('pásmo %g < γ ≤ %g; L-funkce %s; nul ζ %d | %.0f s' % (T1, T2, {j: len(G[j]) for j in JM}, len(G['zeta']), time.time() - t0), flush=True)
def faze(g, c, Tref): return np.mod((g - Tref) * c + float(mp.frac(mp.mpf(Tref) * mp.mpf(c))), 1.0)
def hist(u): return np.bincount(np.minimum((u * NB).astype(np.int64), NB - 1), minlength=NB).astype(float)
ED = np.linspace(0, 1, NB + 1)
# ---------------------------------------------------------------- 3.6 sjednocení všech charakterů modulu
SJ = {5: ['zeta', 'VI_D5', 'q5_c1', 'q5_c3'], 7: ['zeta', 'VI_D-7', 'q7_c2', 'q7_c4', 'q7_c1', 'q7_c5']}
PRU = {5: [2, 3, 7, 11, 13], 7: [2, 3, 5, 11, 13, 29]}; XMAX = 200      # dodatek 3: hodnotí se jen x = p^k ≤ 200 (chybový člen Landauova vzorce roste s x); větší x popisně
out['teleso'] = {}; dc = []; ph = []; pb = []
for q, mem in SJ.items():
    if not all(m in G for m in mem): continue
    phi = len(mem); NU = sum(len(G[m]) for m in mem); w = (T2 - T1) / (2 * math.pi * NU); res = dict(clenove=mem, N=NU, w=w, c={}, hist={})
    for p in PRU[q]:
        us = [faze(G[m], math.log(p) / (2 * math.pi), T1) for m in mem]
        for k in range(1, 7):
            c = complex(sum(np.sum(np.exp(2j * np.pi * k * u)) for u in us) / NU); one = pow(p, k, q) == 1
            cp = -w * phi * (1.0 if one else 0.0) * math.log(p) / p ** (k / 2); hod = p ** k <= XMAX; res['c']['%d,%d' % (p, k)] = dict(c=[c.real, c.imag], c_pred=cp, p_k_1=one, rozdil=abs(c - cp), hodnoceno=hod)
            if hod: dc.append(abs(c - cp))
        if p in (2, 3):
            O = sum(hist(u) for u in us); nul = np.zeros((NNUL, NB))
            for r in range(NNUL):
                for m in mem:
                    gm = G[m]; gg = gm[0] + np.concatenate([[0.0], np.cumsum(rng.permutation(np.diff(gm)))]); nul[r] += hist(faze(gg, math.log(p) / (2 * math.pi), T1))
            sdn = nul.std(0, ddof=1)
            def Rb(kset):
                R_ = np.zeros(NB)
                for k in kset: R_ -= 2 * w * phi * math.log(p) / p ** (k / 2) * (np.sin(2 * np.pi * k * ED[1:]) - np.sin(2 * np.pi * k * ED[:-1])) / (2 * np.pi * k)
                return R_
            E = NU * (1 / NB + Rb([k for k in range(1, KMAX + 1) if pow(p, k, q) == 1])); Ea = NU * (1 / NB + Rb(range(1, KMAX + 1)))
            err = np.maximum(np.sqrt(E), sdn); erra = np.maximum(np.sqrt(Ea), sdn); c2 = float(np.sum(((O - E) / err) ** 2)); c2a = float(np.sum(((O - Ea) / erra) ** 2))
            r_ = dict(O=O.tolist(), E=E.tolist(), E_bez_interference=Ea.tolist(), err=err.tolist(), chi2=c2, p=float(stats.chi2.sf(c2, NB - 1)), chi2_bez=c2a, p_bez=float(stats.chi2.sf(c2a, NB - 1)))
            res['hist'][str(p)] = r_; ph.append(r_['p'])
            if p == 2: pb.append(r_['p_bez'])
            print('3.6 q %d p %d: histogram χ² %.1f (p %.3g) | bez interference %.1f (p %.2g)' % (q, p, c2, r_['p'], c2a, r_['p_bez']), flush=True)
    mx = max(v['rozdil'] for v in res['c'].values() if v['hodnoceno']); res['max_rozdil'] = mx; res['max_rozdil_x_nad_200'] = max([v['rozdil'] for v in res['c'].values() if not v['hodnoceno']] or [0.0]); out['teleso'][str(q)] = res
    print('3.6 q %d: N_U %d, w_U %.4f | největší |c − c_pred| (x ≤ 200) %.1e | c(2, k) %s | %.0f s' % (q, NU, w, mx, ' '.join('%+.4f' % res['c']['2,%d' % k]['c'][0] for k in range(1, 7)), time.time() - t0), flush=True)
v36 = 'N' if not dc else ('O' if (max(dc) <= 3e-4 and min(ph) >= 1e-3 and len(pb) == 2 and max(pb) < 1e-6) else ('Z' if (max(dc) > 2e-3 or sum(x < 1e-4 for x in ph) >= 2) else 'N'))
# ---------------------------------------------------------------- 3.7 křížová kovariance (okna po H jako 3.5 hr7_eval.py)
KX = [3.0, 4.0, 5.0, 6.0]; TT = T1 + 9 * H; Mw = int((T2 - 9 * H - TT) // H) + 1; Mi = int((T2 - 9 * H - TT) // (6 * H)) + 1; T0s = TT + H * np.arange(Mw)
DEL = np.round(np.arange(-2.5, 2.5 + 1e-9, 0.025), 3)
def okna_S(g, lxs):
    m = np.rint((g - TT) / H).astype(np.int64); S = np.zeros((len(lxs), Mw), complex); N2 = np.zeros(Mw)
    for s in range(-9, 10):
        ix = m + s; ok = (ix >= 0) & (ix < Mw); w = np.where(ok, np.exp(-0.5 * ((g - TT - H * ix) / H) ** 2), 0.0); ix = np.where(ok, ix, 0); N2 += np.bincount(ix, w * w, Mw)
        for j, lx in enumerate(lxs):
            p_ = g * (lx if np.isscalar(lx) else lx[ix]); S[j] += np.bincount(ix, w * np.cos(p_), Mw) + 1j * np.bincount(ix, w * np.sin(p_), Mw)
    return S, N2
SK = {}; SD = {}
for j in JM:
    SK[j] = okna_S(G[j], KX); lxw = np.log(INFO[j]['q'] * T0s / (2 * math.pi)); SD[j] = okna_S(G[j], [lxw + d for d in DEL])
def cov(Sa, Na, Sb, Nb, ix=None):
    if ix is not None: Sa, Na, Sb, Nb = Sa[:, ix], Na[ix], Sb[:, ix], Nb[ix]
    return ((Sa - Sa.mean(1, keepdims=True)) * np.conj(Sb - Sb.mean(1, keepdims=True))).mean(1) / math.sqrt(Na.mean() * Nb.mean())
pairs = [(a, b) for i, a in enumerate(JM) for b in JM[i + 1:] if INFO[a]['q'] == INFO[b]['q']]
pod = np.array_split(np.arange(Mw), 10); out['krizova'] = {}; zz = []; sila = []
for a, b in pairs:
    q = INFO[a]['q']; (Sa, Na), (Sb, Nb) = SK[a], SK[b]; C = cov(Sa, Na, Sb, Nb); Va = cov(Sa, Na, Sa, Na).real; Vb = cov(Sb, Nb, Sb, Nb).real
    Cp = np.array([cov(Sa, Na, Sb, Nb, ix) for ix in pod]); er = np.maximum(np.sqrt(Va * Vb / (2 * Mi)), Cp.real.std(0, ddof=1) / math.sqrt(10)); ei = np.maximum(np.sqrt(Va * Vb / (2 * Mi)), Cp.imag.std(0, ddof=1) / math.sqrt(10))
    rows = []
    for j, kx in enumerate(KX):
        x = math.exp(kx); s = 0j
        for n in range(max(2, int(x * math.exp(-0.7))), int(x * math.exp(0.7)) + 2):
            L_ = lam(n)
            if L_ and math.gcd(n, q) == 1: s += INFO[a]['Tab'][n % q] * np.conj(INFO[b]['Tab'][n % q]) * L_ * L_ / n * math.exp(-H * H * math.log(n / x) ** 2)
        Dp = (H * H / (2 * math.pi)) * s / math.sqrt(Na.mean() * Nb.mean()); zr = float((C[j].real - Dp.real) / er[j]); zi = float((C[j].imag - Dp.imag) / ei[j]); zz += [zr, zi]
        rows.append(dict(k=kx, C=[C[j].real, C[j].imag], D=[Dp.real, Dp.imag], err=[float(er[j]), float(ei[j])], z=[zr, zi], V_a=float(Va[j]), V_b=float(Vb[j]), C_k_chybe=float(abs(C[j]) / math.hypot(er[j], ei[j]))))
    sila.append(rows[0]['C_k_chybe'])
    (Sa, Na), (Sb, Nb) = SD[a], SD[b]; Cd = cov(Sa, Na, Sb, Nb); Vda = cov(Sa, Na, Sa, Na).real; Vdb = cov(Sb, Nb, Sb, Nb).real
    out['krizova']['%s × %s' % (a, b)] = dict(diag=rows, hreben_C=[[float(c.real), float(c.imag)] for c in Cd], hreben_rel_max=float(np.max(np.abs(Cd) / np.sqrt(Vda * Vdb))))
    print('3.7 %-22s z: %s | |C|/chyba u k = 3: %.1f | hřeben: největší |C|/√(VaVb) %.3f' % ('%s × %s' % (a, b), ' '.join('%+.2f' % v for r in rows for v in r['z']), rows[0]['C_k_chybe'], out['krizova']['%s × %s' % (a, b)]['hreben_rel_max']), flush=True)
v37 = 'N' if not zz else ('O' if (max(abs(v) for v in zz) <= 3 and min(sila) >= 5) else ('Z' if sum(abs(v) > 4 for v in zz) >= 2 else 'N'))
# ---------------------------------------------------------------- 3.8 Berryho nasycení (dodatek 3: hlavní předpověď = Σ² z dvoubodové funkce BK bez p | q; Berryho diagonální vzorec popisně)
argv0 = sys.argv; sys.argv = ['x', '/dev/null', str(H), str(T1), str(T2)]; NT = {}
exec(open(PY + 'hr7_teorie.py').read().split("if __name__ == '__main__':")[0], NT); sys.argv = argv0
import hr7_sigma_bk as SB
TABQ = {}
def tabq(q):
    if q not in TABQ: TABQ[q] = NT['tab_q'](q) if q > 1 else NT['QZ']
    return TABQ[q]
LG = np.round(np.logspace(-1, math.log10(200), 30), 6).tolist()
def theta_v(q, odd, t): return 0.5 * t * math.log(q / math.pi) + np.imag(loggamma((0.5 + odd) / 2 + 0.5j * t))
def sig2(x, ell, krok=0.37):
    st = x[0] + krok * np.arange(int(np.ceil((x[-1] - ell - x[0]) / krok))); st = st[st < x[-1] - ell]
    if len(st) < 10: return float('nan')
    c = (np.searchsorted(x, st + ell, 'right') - np.searchsorted(x, st, 'right')).astype(float); return float(np.mean(c * c) - np.mean(c) ** 2)
out['nasyceni'] = dict(ell=LG, L={}); S2 = {}
for j in JM:
    q = INFO[j]['q']; x = theta_v(q, INFO[j]['odd'], G[j]) / math.pi; s2 = np.array([sig2(x, l) for l in LG]); parts = np.array_split(x, 10)
    e = np.array([np.nanstd([sig2(pp, l) for pp in parts], ddof=1) for l in LG]) / math.sqrt(10); S2[j] = (s2, e)
    Hs = [math.log(q * (T1 + (i + 0.5) * (T2 - T1) / 24) / (2 * math.pi)) for i in range(24)]
    Bq = np.array([HB.berry_q(l, Hs, q) for l in LG]); Gu = np.array([HB.gue(l) for l in LG])
    Kq = SB.sigma2_bk(LG, Hs, tabq(q), NT['U'], NT['DU']); Kn = SB.sigma2_bk(LG, Hs, tabq(1), NT['U'], NT['DU'])
    Q = float(np.nanmean(((s2 - Kq) / e) ** 2)); Qg = float(np.nanmean(((s2 - Gu) / e) ** 2)); Qn = float(np.nanmean(((s2 - Kn) / e) ** 2)); Qb = float(np.nanmean(((s2 - Bq) / e) ** 2))
    out['nasyceni']['L'][j] = dict(q=q, sigma2=s2.tolist(), err=e.tolist(), BK_q=Kq.tolist(), BK_naiv=Kn.tolist(), berry_q=Bq.tolist(), gue=Gu.tolist(), Q=Q, Q_GUE=Qg, Q_naiv=Qn, Q_berry_diag=Qb)
    print('3.8 %-10s q %2d: Σ²(ℓ = 1, 10, 100) %s | BK_q %s | Q %.2f, Q_GUE %.1f, Q_naiv %.2f | Berry diag. Q %.2f' % (j, q, ' '.join('%.3f' % s2[min(range(30), key=lambda i: abs(LG[i] - v))] for v in (1, 10, 100)),
          ' '.join('%.3f' % Kq[min(range(30), key=lambda i: abs(LG[i] - v))] for v in (1, 10, 100)), Q, Qg, Qn, Qb), flush=True)
NS = out['nasyceni']['L']; out['nasyceni']['dvojice'] = {}
for a, b in pairs:
    (sa, ea), (sb, eb) = S2[a], S2[b]; out['nasyceni']['dvojice']['%s × %s' % (a, b)] = float(np.nanmean((sa - sb) ** 2 / (ea ** 2 + eb ** 2)))
def hladky(m, t):
    if m == 'zeta': return (t / 2 * np.log(t / (2 * math.pi * math.e)) - math.pi / 8 + 1 / (48 * t)) / math.pi + 1
    return theta_v(INFO[m]['q'], INFO[m]['odd'], t) / math.pi
for q, mem in SJ.items():
    if not all(m in G for m in mem): continue
    gu = np.sort(np.concatenate([G[m] for m in mem])); xu = sum(hladky(m, gu) for m in mem); xs = [hladky(m, G[m]) for m in mem]
    out['nasyceni']['sjednoceni_%d' % q] = dict(sigma2=[sig2(xu, l) for l in LG], soucet_sektoru=[float(np.nansum([sig2(xx, l / len(mem)) for xx in xs])) for l in LG])
HR8 = [j for j in JM if not j.startswith('VI_')]; S6 = [j for j in HR8 if INFO[j]['q'] in (15, 21)]; QA = list(out['nasyceni']['dvojice'].values())
v38 = 'O' if (all(NS[j]['Q'] <= 3 for j in HR8) and all(NS[j]['Q_GUE'] >= 2 * NS[j]['Q'] for j in HR8) and all(NS[j]['Q_naiv'] >= 2 * NS[j]['Q'] for j in S6) and all(x <= 3 for x in QA)) else \
      ('Z' if (sum(NS[j]['Q'] > 10 for j in HR8) >= 2 or sum(NS[j]['Q_naiv'] < NS[j]['Q'] for j in S6) >= 2) else 'N')
# ---------------------------------------------------------------- 3.10 Hannovo okno pro L-funkce
U, DU = NT['U'], NT['DU']; tt = np.linspace(-A_, A_, 24001); dt = tt[1] - tt[0]; wv = np.cos(np.pi * tt / (2 * A_)) ** 2; w0 = float(np.sum(wv * wv) * dt)
UA = np.arange(0, 2 * A_ + 0.05, 0.05); AA = np.array([float(np.sum(wv[np.abs(tt + u) < A_] * np.cos(np.pi * (tt + u)[np.abs(tt + u) < A_] / (2 * A_)) ** 2) * dt) / w0 if u < 2 * A_ else 0.0 for u in UA]); AH = np.interp(U, UA, AA, right=0.0)
def V_A(q):
    qq, h2 = tabq(q); Ts = T1 + 9 * H + (np.arange(24) + 0.5) * (T2 - 18 * H - T1) / 24; Ls = np.log(q * Ts / (2 * math.pi)); acc = np.zeros(len(DEL)); accg = np.zeros(len(DEL))
    for L in Ls:
        C = np.cos(np.outer(L + DEL, U)); Bv = (h2 + 2 * (1 - (np.exp(1j * L * U) * qq).real) / (U * U)) * AH; acc += L * (1 - (C @ Bv) * DU / (math.pi * L))
        Bg = (2 * (1 - np.cos(L * U)) / (U * U)) * AH; accg += L * (1 - (C @ Bg) * DU / (math.pi * L))
    return acc / Ls.sum(), accg / Ls.sum()
Mh = int((T2 - 9 * H - TT) // (6 * H)) + 1; out['okno_L'] = {}; VAc = {}
for j in JM:
    q = INFO[j]['q']; g = G[j]; m = np.rint((g - TT) / (6 * H)).astype(np.int64); d_ = g - TT - 6 * H * m; ok = (m >= 0) & (m < Mh) & (np.abs(d_) < A_); w = np.where(ok, np.cos(np.pi * d_ / (2 * A_)) ** 2, 0.0); ix = np.where(ok, m, 0)
    N2 = np.bincount(ix, w * w, Mh); T0h = TT + 6 * H * np.arange(Mh); lx = np.log(q * T0h / (2 * math.pi)); S = np.zeros((len(DEL), Mh), complex)
    for jd, dl in enumerate(DEL):
        p_ = g * (lx + dl)[ix]; S[jd] = np.bincount(ix, w * np.cos(p_), Mh) + 1j * np.bincount(ix, w * np.sin(p_), Mh)
    V = (np.abs(S - S.mean(1, keepdims=True)) ** 2).mean(1) / N2.mean()
    if q not in VAc: VAc[q] = V_A(q)
    VB, VG = VAc[q]; e = V / math.sqrt(Mh); Q = float(np.mean(((V - VB) / e) ** 2)); Qg = float(np.mean(((V - VG) / e) ** 2)); res = dict(q=q, M=Mh, V=V.tolist(), V_BK=VB.tolist(), V_GUE=VG.tolist(), Q=Q, Q_GUE=Qg)
    if Mh > len(DEL) + 10:
        sS = S.sum(1); sS2 = (np.abs(S) ** 2).sum(1); sN = N2.sum(); J = np.zeros((Mh, len(DEL)))
        for i in range(Mh):
            n_ = Mh - 1; s1 = sS - S[:, i]; s2 = sS2 - np.abs(S[:, i]) ** 2; J[i] = (s2 / n_ - np.abs(s1 / n_) ** 2) / ((sN - N2[i]) / n_)
        Jm = J.mean(0); C = (Mh - 1) / Mh * (J - Jm).T @ (J - Jm); r = V - VB; T2_ = float(r @ np.linalg.solve(C, r)); pdim = len(DEL); F = T2_ * (Mh - pdim) / (pdim * (Mh - 1))
        sd = np.sqrt(np.diag(C)); Rk = C / np.outer(sd, sd); ev = np.linalg.eigvalsh(Rk); res.update(T2=T2_, F=F, p_BK=float(stats.f.sf(F, pdim, Mh - pdim)), N_ef=float(ev.sum() ** 2 / (ev ** 2).sum()))
    out['okno_L'][j] = res; print('3.10 %-10s q %2d: Hann M %d | Q_BK %.2f, Q_GUE %.1f | p_BK %s, N_ef %s' % (j, q, Mh, Q, Qg, ('%.3g' % res['p_BK']) if 'p_BK' in res else '—', ('%.0f' % res['N_ef']) if 'N_ef' in res else '—'), flush=True)
OL = out['okno_L']; HR = [j for j in JM if not j.startswith('VI_')]
v310 = 'N' if not all('p_BK' in OL[j] for j in HR) else ('O' if all(OL[j]['Q'] <= 3 and OL[j]['p_BK'] >= 1e-3 for j in HR) else ('Z' if sum(OL[j]['Q'] > 10 or OL[j]['p_BK'] < 1e-5 for j in HR) >= 2 else 'N'))
out['verdikty'] = {'odv-hr7-teleso': v36, 'odv-hr7-krizova': v37, 'odv-hr7-nasyceni': v38, 'odv-hr7-okno-L': v310}
print('VERDIKTY%s:' % (' (ZKOUŠKA KÓDU, ne výsledek)' if ZK else ''), out['verdikty'], '| %.0f s' % (time.time() - t0))
json.dump(out, open(OUT, 'w'))
