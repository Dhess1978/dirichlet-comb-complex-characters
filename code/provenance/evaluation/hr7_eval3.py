"""Studie VII · hudební čtení a slepé testy „použití“ (PRAVIDLA_STUDIE_VII.md, dodatek 5): 3.11 odv-hr7-nastroj, 3.12 odv-hr7-ladeni, 3.13 odv-hr7-teleso2;
popisně prvočíselná mez a negativní kontrola souzvuku. Zmrazí se otiskem dřív, než se spustí jakékoli vyhodnocení dat hr7.
Použití: python3 hr7_eval3.py SLOZKA_HR7 SLOZKA_HR7B|- SLOZKA_REF_VI VYSTUP.json [--zkouska T1 T2]"""
import gzip, hashlib, json, math, os, sys, time, itertools, numpy as np, mpmath as mp
t0 = time.time(); PY = os.path.dirname(os.path.abspath(__file__)) + '/'
VA, VB_, VREF, OUT = sys.argv[1:5]; ZK = '--zkouska' in sys.argv
if ZK: i_ = sys.argv.index('--zkouska'); T1, T2Z = float(sys.argv[i_ + 1]), float(sys.argv[i_ + 2])
else: T1, T2Z = 36000.0, None
H = 20.0; PR = [2, 3, 5, 7, 11, 13]; mp.mp.dps = 20
ZF = ['/home/claude/br69/data/zeros_14_t.npy', '/home/claude/br69/data/zeros_5000_t.npy', '/home/claude/br69/data/zeros_26000_t.npy']
import cypari2
PARI = cypari2.Pari()
def chi_hodnoty(q, chis, n_list):
    G = 'znstar(%d,1)' % q; return {n: (complex(mp.expjpi(2 * mp.mpf(str(PARI('chareval(%s,%s,%d)' % (G, chis, n)))))) if math.gcd(n, q) == 1 else 0j) for n in n_list}
def slozka_D(D): return 'D_%s%d' % ('m' if D < 0 else '', abs(D))
# ---------------------------------------------------------------- načtení (kontroly jako hr7_eval.py; manifest se bere jen pro seznam a popis, bloky podle jmen souborů)
G = {}; INFO = {}; K = None
def nacti_slozku(V):
    global K
    M = json.load(open(os.path.join(V, 'manifest.json'))); assert M['format'] == 'MAYAN_V1008_FZ_NULY_CHI_v3'; OK = float(M['okno']); k0 = int(round(T1 / OK)); res = {}
    for jm, z in M['L'].items():
        k = k0
        while os.path.exists(os.path.join(V, jm, 'okno_%06d.json.gz' % k)): k += 1
        res[jm] = (z['popis'], k)
    return M, OK, k0, res
SL = [(VA, 'hr7')] + ([(VB_, 'hr7b')] if VB_ != '-' else [])
NAC = [(V, s, *nacti_slozku(V)) for V, s in SL]
Kmin = min(k for _, _, _, _, _, res in NAC for (_, k) in res.values()); OK = NAC[0][3]; k0 = NAC[0][4]; T2 = Kmin * OK
if ZK: T2 = T2Z; Kmin = int(round(T2 / OK))
assert ZK or T2 >= 60000.0
for V, s, M, OK_, k0_, res in NAC:
    for jm, (popis, _) in res.items():
        g = []
        for k in range(k0_, Kmin):
            b = json.loads(gzip.decompress(open(os.path.join(V, jm, 'okno_%06d.json.gz' % k), 'rb').read())); assert abs(b['a'] - k * OK_) < 1e-9 and b['n'] == len(b['g']); g.extend(b['g'])
        g = np.array(g); assert np.all(np.diff(g) > 0)
        q = int(popis['q']); chis = '[' + ','.join(str(int(c)) for c in popis['chi']) + ']'; G[jm] = g; INFO[jm] = dict(q=q, chis=chis, sada=s, chi=chi_hodnoty(q, chis, list(range(q))))
MR = json.load(open(os.path.join(VREF, 'manifest.json'))); OKr = float(MR['okno'])
for D in (-3, 5, -7):
    B = MR['L'][str(D)]['bloky']; g = []
    for k in range(int(round(T1 / OKr)), int(round(T2 / OKr))):
        raw = open(os.path.join(VREF, slozka_D(D), B[str(k)]['soubor']), 'rb').read(); assert hashlib.sha256(raw).hexdigest() == B[str(k)]['sha256']; g.extend(json.loads(gzip.decompress(raw))['g'])
    g = np.array([v for v in g if T1 < v <= T2]); q = abs(D); G['VI_D%d' % D] = g; INFO['VI_D%d' % D] = dict(q=q, chis=None, sada='VI', chi={a: complex(float(PARI.kronecker(D, a))) for a in range(q)})
z = np.concatenate([np.load(f) for f in ZF]); G['zeta'] = z[(z > T1) & (z <= T2)]; del z; INFO['zeta'] = dict(q=1, chis=None, sada='zeta', chi={0: 1 + 0j})
HL = [j for j in G if INFO[j]['sada'] in ('hr7', 'hr7b')]
out = dict(T1=T1, T2=T2, zkouska=ZK, L=HL, nul={j: int(len(G[j])) for j in G}, skript_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())
print('pásmo %g < γ ≤ %g; %s | %.0f s' % (T1, T2, {j: len(G[j]) for j in G}, time.time() - t0), flush=True)
def faze(g, c, Tref): return np.mod((g - Tref) * c + float(mp.frac(mp.mpf(Tref) * mp.mpf(c))), 1.0)
def ck(g, p, k=1): u = faze(g, math.log(p) / (2 * math.pi), T1); return complex(np.mean(np.exp(2j * np.pi * k * u)))
# ---------------------------------------------------------------- 3.11 slyšet nástroj
def prvocinitele(n): return {p for p in PR if n % p == 0}
KAND = []
for q in range(3, 31):
    r = q
    for p in PR:
        while r % p == 0: r //= p
    if r != 1 or q % 4 == 2: continue
    Gs = 'znstar(%d,1)' % q; allc = PARI('my(G=%s,v=List()); forvec(x=vector(#G.cyc,i,[0,G.cyc[i]-1]), listput(v,x)); Vec(v)' % Gs)
    for c in allc:
        cs = str(c).replace(' ', '')
        if int(PARI('zncharconductor(%s,%s)' % (Gs, cs))) != q: continue
        KAND.append((q, cs, chi_hodnoty(q, cs, PR)))
out['nastroj'] = {}; spravne = 0
for j in HL:
    g = G[j]; N = len(g); w = (T2 - T1) / (2 * math.pi * N); c = {p: ck(g, p) for p in PR}; amp = {p: w * math.log(p) / math.sqrt(p) for p in PR}
    Dh = {p for p in PR if abs(c[p]) < 0.5 * amp[p]}; chih = {p: -c[p] / amp[p] for p in PR}
    sc = sorted((sum(abs(v[p] - chih[p]) ** 2 for p in PR if p not in Dh), q, cs) for q, cs, v in KAND if prvocinitele(q) == Dh)
    best = sc[0] if sc else (None, None, None); ok = best[1] == INFO[j]['q'] and best[2] == INFO[j]['chis']; spravne += ok
    out['nastroj'][j] = dict(D_odhad=sorted(Dh), chi_odhad={str(p): [chih[p].real, chih[p].imag] for p in PR}, vybrano=[best[1], best[2]], S=best[0], druhy=list(sc[1]) if len(sc) > 1 else None,
                             pravda=[INFO[j]['q'], INFO[j]['chis']], spravne=bool(ok), kandidatu=len(sc))
    print('3.11 %-10s dělitelé %s → vybráno q = %s, χ = %s (S %.2e; druhý S %.2f) | pravda q = %d, χ = %s → %s' % (j, sorted(Dh), best[1], best[2], best[0] if best[0] is not None else float('nan'),
          sc[1][0] if len(sc) > 1 else float('nan'), INFO[j]['q'], INFO[j]['chis'], 'SPRÁVNĚ' if ok else 'ŠPATNĚ'), flush=True)
v311 = 'O' if spravne == len(HL) else ('Z' if len(HL) - spravne >= 2 else 'N')
# ---------------------------------------------------------------- 3.12 slyšet ladění (hřeben jako hr7_eval.py 3.1; předpovědi pro 16 podmnožin S)
argv0 = sys.argv; sys.argv = ['x', '/dev/null']; NM = {}; exec(open(PY + 'hr7_meze_teorie.py').read().split("\nif __name__ == '__main__':")[0], NM); sys.argv = argv0
DEL = NM['DEL']; tab_S, V_sled = NM['tab_S'], NM['V_sled']
def Ls(q): TT = T1 + 9 * H; Tb = T2 - 9 * H; return np.log(q * (TT + (np.arange(24) + 0.5) * (Tb - TT) / 24) / (2 * math.pi))
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
    return V, err
MAL = [2, 3, 5, 7]; SUBS = [list(s) for r in range(5) for s in itertools.combinations(MAL, r)]; PRED = {}
out['ladeni'] = {}; spr = {}; HRB = {}
for j in HL:
    q = INFO[j]['q']; V, err = hreben(G[j], q); HRB[j] = (V, err)
    if q not in PRED: PRED[q] = {str(S): V_sled(DEL, Ls(q), tab_S(S)) for S in SUBS}
    Qs = {k: float(np.mean(((V - v) / err) ** 2)) for k, v in PRED[q].items()}; best = min(Qs, key=Qs.get); pravda = str([p for p in MAL if q % p == 0]); spr[j] = best == pravda
    out['ladeni'][j] = dict(q=q, Q=Qs, vybrano=best, pravda=pravda, spravne=bool(spr[j]))
    print('3.12 %-10s q %2d: vybráno S = %s (Q %.2f; druhé nejlepší %.2f) | pravda %s → %s | %.0f s' % (j, q, best, Qs[best], sorted(Qs.values())[1], pravda, 'SPRÁVNĚ' if spr[j] else 'ŠPATNĚ', time.time() - t0), flush=True)
v312 = 'O' if all(spr.values()) else ('Z' if sum(1 for j in HL if INFO[j]['q'] in (5, 15, 21) and not spr[j]) >= 2 else 'N')
# ---------------------------------------------------------------- 3.13 slyšet těleso
SJ = {5: ['zeta', 'VI_D5', 'q5_c1', 'q5_c3'], 7: ['zeta', 'VI_D-7', 'q7_c2', 'q7_c4', 'q7_c1', 'q7_c5'],
      15: ['zeta', 'VI_D-3', 'VI_D5', 'q5_c1', 'q5_c3', 'q15_c2-1', 'q15_c1-1', 'q15_c3-1'],
      21: ['zeta', 'VI_D-3', 'VI_D-7', 'q7_c2', 'q7_c4', 'q7_c1', 'q7_c5', 'q21_c3-1', 'q21_c1-1', 'q21_c5-1', 'q21_c2-1', 'q21_c4-1']}
PRU = {5: [2, 3, 7, 11, 13], 7: [2, 3, 5, 11, 13, 29], 15: [2, 7, 11, 13, 31], 21: [2, 5, 11, 13, 43]}
FPRED = {5: {2: 4, 3: 4, 7: 4, 11: 1}, 7: {2: 3, 3: 6, 11: 3, 13: 2, 29: 1}, 15: {2: 4, 7: 4, 11: 2, 31: 1}, 21: {2: 6, 13: 2, 43: 1}}
out['teleso2'] = {}; fok = []; dmax = []
for q, mem in SJ.items():
    if not all(m in G for m in mem): continue
    phi = len(mem); NU = sum(len(G[m]) for m in mem); w = (T2 - T1) / (2 * math.pi * NU); C = {}
    for p in PRU[q]:
        us = [faze(G[m], math.log(p) / (2 * math.pi), T1) for m in mem]
        for k in range(1, 7): C[(p, k)] = complex(sum(np.sum(np.exp(2j * np.pi * k * u)) for u in us) / NU)
    nul = [abs(C[(p, k)]) for (p, k) in C if pow(p, k, q) != 1 and p ** k <= 200]; sig = 1.4826 * float(np.median(nul))
    # dodatek 6: f̂ = nejmenší k ≤ 6, při kterém |c_k| dosáhne aspoň poloviny plné síly w_U·φ(q)·ln p/p^{k/2} (σ̂ jen popisně)
    fh = {p: next((k for k in range(1, 7) if abs(C[(p, k)]) >= 0.5 * w * phi * math.log(p) / p ** (k / 2)), None) for p in PRU[q]}
    res = dict(clenove=mem, N=NU, w=w, sigma=sig, f_odhad={str(p): fh[p] for p in PRU[q]}, f_pred={str(p): f for p, f in FPRED[q].items()},
               c={'%d,%d' % pk: dict(c=[v.real, v.imag], c_pred=-w * phi * (1.0 if pow(pk[0], pk[1], q) == 1 else 0.0) * math.log(pk[0]) / pk[0] ** (pk[1] / 2), x=pk[0] ** pk[1]) for pk, v in C.items()})
    for p, f in FPRED[q].items(): fok.append(fh[p] == f)
    if q in (15, 21): dmax.append(max(abs(complex(*v['c']) - v['c_pred']) for v in res['c'].values() if v['x'] <= 200))
    # negativní kontrola: sjednocení bez jednoho sektoru (popisně)
    if q in (5, 7):
        vyn = 'VI_D5' if q == 5 else 'zeta'; mem2 = [m for m in mem if m != vyn]; N2_ = sum(len(G[m]) for m in mem2); w2 = (T2 - T1) / (2 * math.pi * N2_); rows = []
        for p in PRU[q][:3]:
            for k in (1, 2):
                c_ = complex(sum(np.sum(np.exp(2j * np.pi * k * faze(G[m], math.log(p) / (2 * math.pi), T1))) for m in mem2) / N2_)
                pr = -w2 * sum(INFO[m]['chi'][p % INFO[m]['q']] ** k if INFO[m]['q'] > 1 else 1 for m in mem2) * math.log(p) / p ** (k / 2)
                rows.append(dict(p=p, k=k, c=[c_.real, c_.imag], c_pred=[pr.real, pr.imag], c_pred_uplne=0.0))
        res['negativni_kontrola'] = dict(vynechano=vyn, radky=rows)
    out['teleso2'][str(q)] = res
    print('3.13 Q(ζ%d): %d sektorů, N_U %d, σ̂ %.1e | řády f̂: %s | předpověď %s' % (q, phi, NU, sig, {p: fh[p] for p in PRU[q]}, FPRED[q]), flush=True)
v313 = 'N' if len(fok) < 16 else ('O' if (all(fok) and dmax and max(dmax) <= 3e-4) else ('Z' if (sum(not x for x in fok) >= 2 or (dmax and max(dmax) > 2e-3)) else 'N'))
out['teleso2_souhrn'] = dict(radu=len(fok), spravne=int(sum(fok)), max_rozdil_15_21=max(dmax) if dmax else None)
# ---------------------------------------------------------------- popisně: prvočíselná mez
out['mez'] = {}
for q in (5, 15, 21):
    js = [j for j in HL if INFO[j]['q'] == q]
    if not js: continue
    PD = [p for p in PR if q % p == 0]; VP = {str(P): V_sled(DEL, Ls(q), tab_S(PD, P)) for P in (2, 3, 5, 7, 11, 13)}; VP['vse'] = V_sled(DEL, Ls(q), tab_S(PD))
    out['mez'][str(q)] = {j: {P: float(np.mean(((HRB[j][0] - v) / HRB[j][1]) ** 2)) for P, v in VP.items()} for j in js}
    print('mez q %d: %s' % (q, {j: {P: round(v, 2) for P, v in d.items()} for j, d in out['mez'][str(q)].items()}), flush=True)
out['verdikty'] = {'odv-hr7-nastroj': v311, 'odv-hr7-ladeni': v312, 'odv-hr7-teleso2': v313}
print('VERDIKTY%s:' % (' (ZKOUŠKA KÓDU, ne výsledek)' if ZK else ''), out['verdikty'], '| %.0f s' % (time.time() - t0))
json.dump(out, open(OUT, 'w'))
