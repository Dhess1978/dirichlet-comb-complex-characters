"""Studie VII · vyhodnocení sady hr7 podle PRAVIDLA_STUDIE_VII.md (oddíly 3.1–3.5). Skript se zmrazí (SHA-256) dřív, než se data hr7 přenesou do cloudu.
Hypotézy: odv-hr7-hreben (hřeben jedné L(s, χ) proti V_BK,q), odv-hr7-nastroje (stejný modul, stejný hřeben — bez teorie), odv-hr7-faze (první Fourierův koeficient
c = −w·χ(p)·ln p/√p a hustota {γ ln p/2π} s fází χ(p)), odv-hr7-presne (explicitní vzorec s Gaussovým oknem, komplexní χ(n)), odv-hr7-diag (rozptyl harmonických
v diagonálním režimu proti součtu s |χ(n)|²).
Měření hřebene je beze změny podle studie VI (fz3_c.py), jen |D| → q; první řád podle fz3_e.py (E1) s komplexním χ(p); explicitní vzorec podle fz3_e.py (E3); diagonální režim podle fz3_e.py (E4).
Použití: python3 hr7_eval.py SLOZKA_HR7 SLOZKA_REF_VI VYSTUP.json [--zkouska T1] [--rychle]
  SLOZKA_REF_VI = nuly χ_D studie VI ve výšce 36 000–80 000 (formát v2, D = 5, −7, 13); '-' = bez reference.
  --zkouska T1: zkouška kódu na generálce (výška od T1; verdikty se nepočítají jako výsledek). --rychle: menší nulové modely a méně ověřených nul (jen zkouška)."""
import gzip, hashlib, json, math, os, sys, time, numpy as np, mpmath as mp
from scipy import stats
t0 = time.time(); PY = os.path.dirname(os.path.abspath(__file__)) + '/'
VN, VREF, OUT = sys.argv[1:4]; ZK = '--zkouska' in sys.argv; RYCHLE = '--rychle' in sys.argv
T1 = float(sys.argv[sys.argv.index('--zkouska') + 1]) if ZK else 36000.0
H = 20.0; PR = [2, 3, 5, 7, 11, 13]; NB = 50; KMAX = 60; NNUL = 40 if RYCHLE else 200; NNUL4 = 20 if RYCHLE else 100; NOV = 8 if RYCHLE else 40; ED = np.linspace(0, 1, NB + 1)
mp.mp.dps = 20; rng = np.random.default_rng(34107)
import cypari2
PARI = cypari2.Pari()
# ---------------------------------------------------------------- charakter: tabulka χ(a), a mod q (PARI chareval), kontroly multiplikativity a parity
def chi_tab(q, chi):
    G = 'znstar(%d,1)' % q; c = '[' + ','.join(str(int(x)) for x in chi) + ']'; T = np.zeros(q, complex)
    for a in range(q):
        if math.gcd(a, q) == 1:
            t = PARI('chareval(%s,%s,%d)' % (G, c, a)); T[a] = complex(mp.expjpi(2 * mp.mpf(str(t))))
    for a in range(q):
        for b in range(q):
            assert abs(T[a * b % q] - T[a] * T[b]) < 1e-12
    odd = int(PARI('zncharisodd(%s,%s)' % (G, c))); assert abs(T[q - 1] - (-1 if odd else 1)) < 1e-12
    return T, odd
def theta(q, odd, T): return 0.5 * T * math.log(q / math.pi) + float(mp.im(mp.loggamma(mp.mpc((0.5 + odd) / 2, T / 2))))
def Lhod(q, Tab, t):
    s = mp.mpc(0.5, t); return float(abs(mp.mpf(q) ** (-s) * mp.fsum(mp.mpc(Tab[a].real, Tab[a].imag) * mp.zeta(s, mp.mpf(a) / q) for a in range(1, q) if abs(Tab[a]) > 0)))
def lam(n):
    if n < 2: return 0.0
    for p in range(2, int(n ** 0.5) + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return math.log(p) if m == 1 else 0.0
    return math.log(n)
# ---------------------------------------------------------------- načtení
def nacti_v3(V, M, jm, k0, k1, OK):
    B = M['L'][jm]['bloky']; g = []; sha = 0; mimo = 0
    for k in range(k0, k1):
        raw = open(os.path.join(V, jm, B[str(k)]['soubor']), 'rb').read(); sha += hashlib.sha256(raw).hexdigest() == B[str(k)]['sha256']
        b = json.loads(gzip.decompress(raw)); assert abs(b['a'] - k * OK) < 1e-9 and abs(b['b'] - (k + 1) * OK) < 1e-9 and b['n'] == len(b['g']), (jm, k)
        mimo += any(not (b['a'] < x <= b['b']) for x in b['g']); g.extend(b['g'])
    return np.array(g, float), sha, mimo
def slozka_D(D): return 'D_%s%d' % ('m' if D < 0 else '', abs(D))
MN = json.load(open(os.path.join(VN, 'manifest.json'))); assert MN['format'] == 'MAYAN_V1008_FZ_NULY_CHI_v3', MN['format']
OK = float(MN['okno']); k0 = int(round(T1 / OK)); assert ZK or (abs(float(MN['tmin']) - 36000.0) < 1e-9 and MN.get('sada') == 'hr7')
JM = list(MN['L'].keys()); KD = {}
for jm in JM:
    B = MN['L'][jm]['bloky']; k = k0
    while str(k) in B and os.path.exists(os.path.join(VN, jm, B[str(k)]['soubor'])): k += 1
    KD[jm] = k
K = min(KD.values()); T2 = K * OK
assert ZK or T2 >= 60000.0, T2
out = dict(slozka=os.path.basename(os.path.normpath(VN)), tmin=T1, T2=T2, okno=OK, divz=MN['divz'], zkouska=ZK, rychle=RYCHLE, H=H, L=JM, kontroly={},
           prostredi={k: v for k, v in MN.get('prostredi', {}).items() if k in ('python', 'system', 'stroj', 'pari', 'skript_sha256')},
           skript_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())
print('pásmo %g < γ ≤ %g (%d oken u všech %d L-funkcí)' % (T1, T2, K - k0, len(JM)), flush=True)
G = {}; INFO = {}; vse = True
for jm in JM:
    P_ = MN['L'][jm]['popis']; q = int(P_['q']); Tab, odd = chi_tab(q, P_['chi']); assert odd == int(P_['lichy']) and int(P_['vodic']) == q
    g, sha, mimo = nacti_v3(VN, MN, jm, k0, K, OK); rost = bool(np.all(np.diff(g) > 0)); oc = (theta(q, odd, T2) - theta(q, odd, T1)) / math.pi
    ok = sha == K - k0 and rost and mimo == 0 and abs(len(g) - oc) <= 6; vse = vse and ok; G[jm] = g; INFO[jm] = dict(q=q, Tab=Tab, odd=odd, chi=P_['chi'])
    out['kontroly'][jm] = dict(q=q, chi=P_['chi'], rad=P_['rad'], lichy=odd, sha_ok=sha, oken=K - k0, rostouci=rost, mimo_okno=mimo, nul=int(len(g)), hladky_pocet=oc, rozdil=float(len(g) - oc), v_poradku=ok)
    print('%-10s q %2d řád %2d: SHA %d/%d, nul %d (hladký počet %.2f, rozdíl %+.2f), rostoucí %s → %s' % (jm, q, P_['rad'], sha, K - k0, len(g), oc, len(g) - oc, rost, 'OK' if ok else 'POZOR'), flush=True)
import multiprocessing as mpr
POOL = mpr.get_context('fork').Pool(2)
for jm in JM:
    g = G[jm]; I_ = INFO[jm]; ix = np.sort(rng.choice(len(g), size=min(NOV, len(g)), replace=False)); hod = POOL.starmap(Lhod, [(I_['q'], I_['Tab'], float(g[i])) for i in ix]); mx = max(hod); ok = mx < 1e-8; vse = vse and ok
    out['kontroly'][jm].update(overeno_nul=len(ix), max_L=mx, overeno_ok=ok)
    print('%-10s %d náhodných nul v mpmath (Hurwitzova ζ): největší |L(½ + iγ, χ)| = %.2e → %s (%.0f s)' % (jm, len(ix), mx, 'OK' if ok else 'POZOR', time.time() - t0), flush=True)
POOL.close(); POOL.join()
out['kontroly_v_poradku'] = vse
if not vse:
    print('kontroly NEPROŠLY — vyhodnocení se neprovádí, napřed dohledat'); json.dump(out, open(OUT, 'w'), indent=0); sys.exit(2)
# ---------------------------------------------------------------- 3.1 hřeben (beze změny podle fz3_c.py studie VI)
argv0 = sys.argv; sys.argv = ['x', '/dev/null', str(H), str(T1), str(T2)]; NT = {}
exec(open(PY + 'hr7_teorie.py').read().split("if __name__ == '__main__':")[0], NT); sys.argv = argv0
DEL = NT['DEL']; assert len(DEL) == 201
def hreben(g, q):
    TT = T1 + 9 * H; Mw = int((T2 - 9 * H - TT) // H) + 1; Mi = int((T2 - 9 * H - TT) // (6 * H)) + 1; m = np.rint((g - TT) / H).astype(np.int64); W = []; I = []
    for s_ in range(-9, 10):
        ix = m + s_; ok = (ix >= 0) & (ix < Mw); W.append(np.where(ok, np.exp(-0.5 * ((g - TT - H * ix) / H) ** 2), 0.0)); I.append(np.where(ok, ix, 0))
    N2 = np.zeros(Mw)
    for w_, ix in zip(W, I): N2 += np.bincount(ix, w_ * w_, Mw)
    T0s = TT + H * np.arange(Mw); lxw = np.log(q * T0s / (2 * math.pi)); S = np.zeros((len(DEL), Mw), complex)
    for j, dl in enumerate(DEL):
        for w_, ix in zip(W, I):
            ph = g * (lxw + dl)[ix]; S[j] += np.bincount(ix, w_ * np.cos(ph), Mw) + 1j * np.bincount(ix, w_ * np.sin(ph), Mw)
    def Vz(S_, N2_): return (np.abs(S_ - S_.mean(1, keepdims=True)) ** 2).mean(1) / N2_.mean()
    V = Vz(S, N2); pod = np.array_split(np.arange(Mw), 10); err = np.maximum(V / math.sqrt(Mi), np.array([Vz(S[:, ix], N2[ix]) for ix in pod]).std(0, ddof=1) / math.sqrt(10))
    return V, err, Mw, Mi
PRED = {}
def pred(q):
    if q not in PRED:
        Lm, VB, VN_, VG = NT['predikce'](q); PRED[q] = (Lm, VB, VN_, VG)
    return PRED[q]
out['hreben'] = {}; HR = {}
for jm in JM:
    q = INFO[jm]['q']; V, err, Mw, Mi = hreben(G[jm], q); Lm, VB, VNv, VG = pred(q); HR[jm] = (V, err)
    Q = float(np.mean(((V - VB) / err) ** 2)); Qn = float(np.mean(((V - VNv) / err) ** 2)); Qg = float(np.mean(((V - VG) / err) ** 2))
    out['hreben'][jm] = dict(q=q, nul=int(len(G[jm])), oken=Mw, M_nezavislych=Mi, L_stred=Lm, p_deli_q=NT['pdel'](q), V=V.tolist(), err=err.tolist(), V_BK_q=VB.tolist(), V_naiv=VNv.tolist(),
                             V_GUE=VG.tolist(), Q=Q, Q_naiv=Qn, Q_GUE=Qg, max_abs_z=float(np.max(np.abs((V - VB) / err))))
    print('3.1 %-10s q %2d: %d oken (M %d) | Q proti BK_q %.2f, naivní ζ %.2f, GUE %.2f | největší |z| %.2f | %.0f s' % (jm, q, Mw, Mi, Q, Qn, Qg, out['hreben'][jm]['max_abs_z'], time.time() - t0), flush=True)
R = out['hreben']; S6 = [j for j in JM if R[j]['q'] in (5, 15, 21)]
v31 = 'O' if (all(R[j]['Q'] <= 3 for j in JM) and all(R[j]['Q_naiv'] >= 2 * R[j]['Q'] for j in S6) and all(R[j]['Q_GUE'] >= 2 * R[j]['Q'] for j in JM)) else \
      ('Z' if (sum(R[j]['Q'] > 10 for j in JM) >= 2 or sum(R[j]['Q_naiv'] < R[j]['Q'] for j in S6) >= 2) else 'N')
# ---------------------------------------------------------------- 3.2 stejný modul, stejný hřeben (bez teorie); reference χ_D ze studie VI
out['reference_VI'] = {}
if VREF != '-':
    MR = json.load(open(os.path.join(VREF, 'manifest.json'))); assert MR['format'] == 'MAYAN_V1008_FZ_NULY_L_v2' and abs(float(MR['tmin']) - T1) < 1e-9
    for D in (5, -7, 13):
        if abs(D) not in [INFO[j]['q'] for j in JM]: continue
        B = MR['L'][str(D)]['bloky']; g = []; sha = 0
        for k in range(k0, K):
            raw = open(os.path.join(VREF, slozka_D(D), B[str(k)]['soubor']), 'rb').read(); sha += hashlib.sha256(raw).hexdigest() == B[str(k)]['sha256']; g.extend(json.loads(gzip.decompress(raw))['g'])
        g = np.array(g, float); assert sha == K - k0 and np.all(np.diff(g) > 0); V, err, Mw, Mi = hreben(g, abs(D)); jm = 'VI_D%d' % D; HR[jm] = (V, err); INFO[jm] = dict(q=abs(D))
        out['reference_VI'][jm] = dict(D=D, nul=int(len(g)), V=V.tolist(), err=err.tolist())
        print('3.2 reference studie VI D = %d: %d nul, %d oken' % (D, len(g), Mw), flush=True)
def Qab(a, b): (Va, ea), (Vb, eb) = HR[a], HR[b]; return float(np.mean((Va - Vb) ** 2 / (ea ** 2 + eb ** 2)))
dv = []; LL = list(HR.keys())
for i, a in enumerate(LL):
    for b in LL[i + 1:]:
        if INFO[a]['q'] == INFO[b]['q']: dv.append((a, b))
kon = [('q5_c1', 'q15_c1-1'), ('q7_c1', 'q21_c1-1'), ('q13_c1', 'q15_c3-1')]; kon = [k for k in kon if k[0] in HR and k[1] in HR]
out['nastroje'] = dict(dvojice={'%s × %s' % d: Qab(*d) for d in dv}, kontrola={'%s × %s' % d: Qab(*d) for d in kon})
for k_, v_ in out['nastroje']['dvojice'].items(): print('3.2 stejný modul %-24s Q_ab %.2f' % (k_, v_), flush=True)
for k_, v_ in out['nastroje']['kontrola'].items(): print('3.2 kontrola síly %-23s Q_ab %.2f' % (k_, v_), flush=True)
QD = list(out['nastroje']['dvojice'].values()); QK = list(out['nastroje']['kontrola'].values())
v32 = 'O' if (QD and all(x <= 3 for x in QD) and len(QK) == 3 and all(x >= 5 for x in QK)) else ('Z' if sum(x > 10 for x in QD) >= 2 else 'N')
# zuby (popisně)
from fractions import Fraction
out['zuby'] = {}
for jm in JM:
    q = INFO[jm]['q']; Tb = INFO[jm]['Tab']; V, err = HR[jm]; VB = np.array(R[jm]['V_BK_q']); rows = []
    for ab in ('1/3', '1/2', '2/3', '1', '3/2', '2', '3', '5'):
        f = Fraction(ab); i = int(np.argmin(np.abs(DEL - math.log(f)))); A_, B_ = f.numerator, f.denominator; cab = Tb[A_ % q] * np.conj(Tb[B_ % q]) if math.gcd(A_ * B_, q) == 1 else 0j
        rows.append(dict(AB=ab, delta=float(DEL[i]), chiA_chiB=[float(cab.real), float(cab.imag)], V=float(V[i]), err=float(err[i]), V_BK_q=float(VB[i]), z=float((V[i] - VB[i]) / err[i])))
    out['zuby'][jm] = rows
# ---------------------------------------------------------------- 3.3 první řád: fáze χ(p)
def faze(g, c, Tref): return np.mod((g - Tref) * c + float(mp.frac(mp.mpf(Tref) * mp.mpf(c))), 1.0)
def hist(u): return np.bincount(np.minimum((u * NB).astype(np.int64), NB - 1), minlength=NB).astype(float)
def Rfaz(w, p, cp):
    """∫ přes koše z −2w·Σ_k (ln p/p^{k/2})·cos(2πkθ − k·arg χ(p)) dθ"""
    R_ = np.zeros(NB); ph = math.atan2(cp.imag, cp.real); r = abs(cp)
    for k in range(1, KMAX + 1):
        a = r ** k * math.log(p) / p ** (k / 2)
        if a: R_ -= 2 * w * a * (np.sin(2 * np.pi * k * ED[1:] - k * ph) - np.sin(2 * np.pi * k * ED[:-1] - k * ph)) / (2 * np.pi * k)
    return R_
out['faze'] = {}; pP = []; pK = []; dc = []
for jm in JM:
    g = G[jm]; q = INFO[jm]['q']; Tb = INFO[jm]['Tab']; N = len(g); w = (T2 - T1) / (2 * math.pi * N); sp = np.diff(g); ps = [p for p in PR if q % p]; nul = {p: np.zeros((NNUL, NB)) for p in ps}
    for r in range(NNUL):
        gg = g[0] + np.concatenate([[0.0], np.cumsum(rng.permutation(sp))])
        for p in ps: nul[p][r] = hist(faze(gg, math.log(p) / (2 * math.pi), T1))
    for p in ps:
        cp = Tb[p % q]; u = faze(g, math.log(p) / (2 * math.pi), T1); O = hist(u); sdn = nul[p].std(0, ddof=1); Eu = np.full(NB, N / NB)
        E = N * (1 / NB + Rfaz(w, p, cp)); Ek = N * (1 / NB + Rfaz(w, p, np.conj(cp)))
        err = np.maximum(np.sqrt(E), sdn); errk = np.maximum(np.sqrt(Ek), sdn); erru = np.maximum(np.sqrt(Eu), sdn)
        c2 = float(np.sum(((O - E) / err) ** 2)); c2k = float(np.sum(((O - Ek) / errk) ** 2)); c2u = float(np.sum(((O - Eu) / erru) ** 2))
        c = complex(np.mean(np.exp(2j * np.pi * u))); cpr = -w * cp * math.log(p) / math.sqrt(p)
        r_ = dict(q=q, p=p, chi_p=[cp.real, cp.imag], N=N, w=w, c=[c.real, c.imag], c_pred=[cpr.real, cpr.imag], rozdil=abs(c - cpr), chi2=c2, p_hodnota=float(stats.chi2.sf(c2, NB - 1)),
                  chi2_sdruzena=c2k, p_sdruzena=float(stats.chi2.sf(c2k, NB - 1)), chi2_rovn=c2u, p_rovn=float(stats.chi2.sf(c2u, NB - 1)), sd_nul_rel=float(np.mean(sdn / np.sqrt(Eu))),
                  O=O.tolist(), E=E.tolist(), E_sdruzena=Ek.tolist(), err=err.tolist())
        out['faze']['%s,%d' % (jm, p)] = r_; pP.append(r_['p_hodnota']); dc.append(r_['rozdil'])
        if abs(cp.imag) >= 0.5 and p <= 7: pK.append(r_['p_sdruzena'])
        print('3.3 %-10s p %2d χ(p) %+.3f%+.3fi: |c − c_pred| %.1e | χ² %6.1f (p %.3g) | sdružená %9.1f (p %.2g) | rovnoměrné %9.1f' % (jm, p, cp.real, cp.imag, r_['rozdil'], c2, r_['p_hodnota'], c2k, r_['p_sdruzena'], c2u), flush=True)
v33 = 'O' if (max(dc) <= 3e-4 and min(pP) >= 1e-3 and pK and max(pK) < 1e-6) else ('Z' if (max(dc) > 2e-3 or sum(x < 1e-4 for x in pP) >= 2) else 'N')
out['faze_souhrn'] = dict(kombinaci=len(pP), max_rozdil_c=max(dc), min_p=min(pP), max_p_sdruzena=max(pK) if pK else None, kombinaci_sdruzena=len(pK))
# ---------------------------------------------------------------- 3.4 explicitní vzorec s Gaussovým oknem (komplexní χ(n))
XS = [('2', 2), ('3', 3), ('4', 4), ('5', 5), ('6', 6), ('7', 7), ('8', 8), ('9', 9), ('10', 10), ('e^sqrt2', None)]; CUT = 12.0
T0S = (40000, 50000, 60000, 70000) if not ZK else tuple(T1 + f * (T2 - T1) for f in (0.3, 0.5, 0.7))
out['presne'] = []; mx3 = 0.0; mxo = 0.0
for jm in JM:
    g = G[jm]; q = INFO[jm]['q']; Tb = INFO[jm]['Tab']
    for T0 in T0S:
        for H3 in (20, 50):
            if T0 + 40 * H3 > T2 or T0 - 40 * H3 < T1: continue
            a = np.searchsorted(g, T0 - 40 * H3); b = np.searchsorted(g, T0 + 40 * H3); d = g[a:b] - T0; wg = np.exp(-0.5 * (d / H3) ** 2); Ng = float(wg.sum())
            for nm, x in XS:
                xm = mp.e ** mp.sqrt(2) if x is None else mp.mpf(x); lx = mp.log(xm); f0 = float(mp.fmod(T0 * lx, 2 * mp.pi)); ph = d * float(lx) + f0
                Z = complex(float(np.sum(wg * np.cos(ph))), float(np.sum(wg * np.sin(ph))))
                xf = float(xm); n1 = max(2, int(math.floor(xf * math.exp(-CUT / H3)))); n2 = int(math.ceil(xf * math.exp(CUT / H3))); P = mp.mpc(0)
                for n in range(n1, n2 + 1):
                    L_ = lam(n)
                    if L_ and math.gcd(n, q) == 1:
                        ch = Tb[n % q]; l = mp.log(mp.mpf(n) / xm); P += mp.mpc(ch.real, ch.imag) * mp.mpf(L_) / mp.sqrt(n) * mp.exp(-(H3 * l) ** 2 / 2) * mp.expj(-T0 * l)
                P = complex(-H3 / mp.sqrt(2 * mp.pi) * P); rz = abs(Z - P) / Ng; mx3 = max(mx3, rz)
                s = H3 * math.sqrt(2 * math.pi) / 2; a2 = np.searchsorted(g, T0 - s); b2 = np.searchsorted(g, T0 + s); ph2 = (g[a2:b2] - T0) * float(lx) + f0
                Sx = complex(float(np.cos(ph2).sum()), float(np.sin(ph2).sum())); chx = Tb[x % q] if (x and math.gcd(x, q) == 1) else 0
                Lnd = -(2 * s) / (2 * math.pi) * (lam(x) * chx if x else 0.0) / math.sqrt(x if x else 1.0); ro = abs(Sx - Lnd) / max(1, b2 - a2); mxo = max(mxo, ro)
                out['presne'].append(dict(L=jm, T0=T0, H=H3, x=nm, Z=[Z.real, Z.imag], P=[P.real, P.imag], Ng=Ng, rozdil=rz, ostre=ro))
v34 = 'O' if mx3 <= 1e-6 else ('Z' if mx3 > 1e-4 else 'N')
out['presne_souhrn'] = dict(bodu=len(out['presne']), max_rozdil=mx3, ostre_max=mxo, ostre_min=min(r['ostre'] for r in out['presne']) if out['presne'] else None)
print('3.4 %d bodů, největší |Z − P|/N_g = %.2e (ostré okno do %.2e) | %.0f s' % (len(out['presne']), mx3, mxo, time.time() - t0), flush=True)
# ---------------------------------------------------------------- 3.5 diagonální režim (podle fz3_e.py E4, okna po H)
KX = [3.0, 4.0, 5.0, 6.0]
def okna(g, Ta, Tb_):
    TT = Ta + 9 * H; Mw = int((Tb_ - 9 * H - TT) // H) + 1; m = np.rint((g - TT) / H).astype(np.int64); W = []; I = []
    for s in range(-9, 10):
        ix = m + s; ok = (ix >= 0) & (ix < Mw); W.append(np.where(ok, np.exp(-0.5 * ((g - TT - H * ix) / H) ** 2), 0.0)); I.append(np.where(ok, ix, 0))
    return TT, Mw, W, I
def Smat(g, lxs, Mw, W, I):
    S = np.zeros((len(lxs), Mw), complex)
    for j, lx in enumerate(lxs):
        for w, ix in zip(W, I):
            ph = g * lx; S[j] += np.bincount(ix, w * np.cos(ph), Mw) + 1j * np.bincount(ix, w * np.sin(ph), Mw)
    return S
def Vz(S, N2): return (np.abs(S - S.mean(1, keepdims=True)) ** 2).mean(1) / N2.mean()
out['diag'] = {}; zs = []
for jm in JM:
    g = G[jm]; q = INFO[jm]['q']; TT, Mw, W, I = okna(g, T1, T2); N2 = np.zeros(Mw)
    for w, ix in zip(W, I): N2 += np.bincount(ix, w * w, Mw)
    S = Smat(g, KX, Mw, W, I); V = Vz(S, N2); Mi = int((T2 - 9 * H - TT) // (6 * H)) + 1; pod = np.array_split(np.arange(Mw), 10)
    err_pod = np.array([Vz(S[:, ix], N2[ix]) for ix in pod]).std(0, ddof=1) / math.sqrt(10); sp = np.diff(g); Vn = []
    for r in range(NNUL4):
        gg = g[0] + np.concatenate([[0.0], np.cumsum(rng.permutation(sp))]); _, Mw2, W2, I2 = okna(gg, T1, T2); N2n = np.zeros(Mw2)
        for w, ix in zip(W2, I2): N2n += np.bincount(ix, w * w, Mw2)
        Vn.append(Vz(Smat(gg, KX, Mw2, W2, I2), N2n))
    Vn = np.array(Vn); err = np.maximum(np.maximum(V / math.sqrt(Mi), err_pod), V * Vn.std(0, ddof=1) / Vn.mean(0)); T0s = TT + H * np.arange(Mw); rows = []
    for j, kx in enumerate(KX):
        x = math.exp(kx); s_chi = 0.0
        for n in range(max(2, int(x * math.exp(-0.7))), int(x * math.exp(0.7)) + 2):
            L_ = lam(n)
            if L_ and math.gcd(n, q) == 1: s_chi += L_ * L_ / n * math.exp(-H * H * math.log(n / x) ** 2)
        Vd = (H * H / (2 * math.pi)) * s_chi / N2.mean(); tau = float(np.mean(kx / np.log(q * T0s / (2 * math.pi)))); z = float((V[j] - Vd) / err[j]); zs.append(z)
        rows.append(dict(k=kx, tau=tau, V=float(V[j]), err=float(err[j]), V_diag=Vd, z=z))
    out['diag'][jm] = rows; print('3.5 %-10s z: %s' % (jm, ', '.join('%+.2f' % r['z'] for r in rows)), flush=True)
v35 = 'O' if max(abs(z) for z in zs) <= 3 else ('Z' if sum(abs(z) > 4 for z in zs) >= 2 else 'N')
out['diag_souhrn'] = dict(cisel=len(zs), max_abs_z=max(abs(z) for z in zs))
out['verdikty'] = {'odv-hr7-hreben': v31, 'odv-hr7-nastroje': v32, 'odv-hr7-faze': v33, 'odv-hr7-presne': v34, 'odv-hr7-diag': v35}
print('VERDIKTY%s:' % (' (ZKOUŠKA KÓDU, ne výsledek)' if ZK else ''), out['verdikty'], '| %.0f s' % (time.time() - t0))
json.dump(out, open(OUT, 'w'))
