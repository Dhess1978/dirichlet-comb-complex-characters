"""Studie VII · část R — robustnost hlavní rezonance studie VI na nulách ζ (pásma P8 a P3e10), PRAVIDLA_STUDIE_VII.md oddíl 4.
R1 (odv-hr7-okno): Hannovo okno w(t) = cos²(πt/2a), a = 3H = 60, okna po 6H (nepřekrývají se); předpověď V_BK s autokorelací okna A(u) = (w⋆w)(u)/(w⋆w)(0).
R2 (odv-hr7-podpasma): v τ = 1 (Gaussovo okno H = 20, okna po 6H jako fz2_c_hreben.py studie VI) chyby z K = 5, 10, 20, 40 podpásem, V/√M a blokového jackknife (1 000 bloků).
R3 (odv-hr7-kovar): kovariance 201 bodů jemného průběhu blokovým jackknifem (B = 1 000), Hotellingova oprava, N_ef = (tr R)²/tr(R²), globální p pro největší |z|.
R4 (popisně): kontrolní oblast τ = 0,5 … 2,0 s chybami z R2.
Teorie: funkce studie VI (fz2_c_teorie1.py, část před sítem) beze změny; pro Hannovo okno se jen nahradí váha GW funkcí A(u).
Použití: python3 hr7_r.py --kontrola VYSTUP.json   (jen teorie: Hannova předpověď proti sítu v diagonálním režimu; nuly se nečtou)
         python3 hr7_r.py PASMO VYSTUP.json          (PASMO = P8 nebo P3e10)"""
import json, math, sys, time, numpy as np
from scipy import stats
t0 = time.time(); KONTROLA = sys.argv[1] == '--kontrola'; OUT = sys.argv[2]; nm = None if KONTROLA else sys.argv[1]
H = 20.0; A_ = 3 * H; Z = '/home/claude/br69/data/zeros_%d_t.npy'
PAS = {'P8': [99146000, 101246000, 103346000], 'P3e10': [30603746000, 30605846000, 30607946000]}
DEL = np.round(np.arange(-2.5, 2.5 + 1e-9, 0.025), 3); TAUK = [0.5, 0.7, 0.9, 1.0, 1.1, 1.3, 1.6, 2.0]; NH = 24; B = 1000; rng = np.random.default_rng(34107)
src = open('/home/claude/br87/py/fz2_c_teorie1.py').read().split('SMAX = ')[0]; argv0 = sys.argv; sys.argv = ['x', '20']; ns = {}; exec(src, ns); sys.argv = argv0
U, DU, GW, QU, H2U = ns['U'], ns['DU'], ns['GW'], ns['QU'], ns['H2U']
# autokorelace Hannova okna, normovaná na 1 v nule: (w⋆w)(u) = ∫ w(t) w(t + u) dt
def hann_A(u):
    tt = np.linspace(-A_, A_, 24001); dt = tt[1] - tt[0]; w = np.cos(np.pi * tt / (2 * A_)) ** 2; w0 = float(np.sum(w * w) * dt); r = np.zeros(len(u))
    for i, uu in enumerate(u):
        if uu >= 2 * A_: continue
        t2 = tt + uu; m = np.abs(t2) < A_; r[i] = float(np.sum(w[m] * np.cos(np.pi * t2[m] / (2 * A_)) ** 2) * dt) / w0
    return r
UA = np.arange(0, 2 * A_ + 0.05, 0.05); AA = hann_A(UA); AH = np.interp(U, UA, AA, right=0.0)
def V_T(k, L, q, h2, Wt):
    return 1 - (1 / (math.pi * L)) * float(np.sum(np.cos(k * U) * (h2 + 2 * (1 - (np.exp(1j * L * U) * q).real) / (U * U)) * Wt) * DU)
def pred(KX, Ls, Wt):
    VB = np.array([float(sum(L * V_T(k, L, QU, H2U, Wt) for L in Ls) / Ls.sum()) for k in KX]); VG = np.array([float(sum(L * V_T(k, L, 1.0, 0.0, Wt) for L in Ls) / Ls.sum()) for k in KX])
    return VB, VG
def what_hann(v):
    """|ŵ(v)|² Hannova okna (numericky)"""
    tt = np.linspace(-A_, A_, 24001); dt = tt[1] - tt[0]; w = np.cos(np.pi * tt / (2 * A_)) ** 2
    return np.array([abs(float(np.sum(w * np.cos(vv * tt)) * dt)) ** 2 for vv in np.atleast_1d(v)])
if KONTROLA:
    # diagonální režim: V_BK s A(u) proti V_diag = Σ Λ²/n·|ŵ(ln(n/x))|²/(2π·L·∫w²) sítem; k = 3 … 6 při L pásma P8
    lam = lambda n: next((math.log(p) if (lambda m: [m := m // p for _ in iter(lambda: m % p == 0, False)] and m == 1)(n) else 0.0 for p in range(2, int(n ** 0.5) + 1) if n % p == 0), math.log(n) if n > 1 else 0.0)
    def lam2(n):
        for p in range(2, int(n ** 0.5) + 1):
            if n % p == 0:
                m = n
                while m % p == 0: m //= p
                return math.log(p) if m == 1 else 0.0
        return math.log(n) if n > 1 else 0.0
    tt = np.linspace(-A_, A_, 24001); w2 = float(np.sum(np.cos(np.pi * tt / (2 * A_)) ** 4) * (tt[1] - tt[0])); L = math.log(1.0e8 / (2 * math.pi)); rows = []
    for k in (3.0, 4.0, 5.0, 6.0):
        x = math.exp(k); ns_ = [n for n in range(max(2, int(x * math.exp(-1.0))), int(x * math.exp(1.0)) + 2) if lam2(n)]
        wv = what_hann(np.log(np.array(ns_, float) / x)); s = float(sum(lam2(n) ** 2 / n * wv[i] for i, n in enumerate(ns_)))
        vd = s / (2 * math.pi * L * w2); vb = V_T(k, L, QU, H2U, AH); vg = V_T(k, L, 1.0, 0.0, AH); vbg = V_T(k, L, QU, H2U, GW)
        rows.append(dict(k=k, tau=k / L, V_BK_hann=vb, V_diag_hann=vd, rozdil=vb - vd, V_GUE_hann=vg, V_BK_gauss=vbg)); print('k %.0f τ %.3f: V_BK(Hann) %.5f, síto %.5f, rozdíl %.1e | GUE %.4f | Gauss %.5f' % (k, k / L, vb, vd, vb - vd, vg, vbg), flush=True)
    json.dump(dict(H=H, a=A_, w2=w2, diag=rows), open(OUT, 'w')); sys.exit(0)
# ---------------------------------------------------------------- data
g = np.concatenate([np.load(Z % b) for b in PAS[nm]]); Ta = float(PAS[nm][0]); Tb = Ta + 2100000.0 * 3; Lm = math.log((Ta + Tb) / 2 / (2 * math.pi))
KX = np.concatenate([Lm + DEL, np.array(TAUK) * Lm]); T1 = Ta + 9 * H; M = int((Tb - 9 * H - T1) // (6 * H)) + 1
SG = np.zeros((len(KX), M), complex); N2G = np.zeros(M); SH = np.zeros((len(KX), M), complex); N2H = np.zeros(M)
for a in range(0, len(g), 4000000):
    gg = g[a:a + 4000000]; m = np.rint((gg - T1) / (6 * H)).astype(np.int64); d0 = gg - Ta; WG = []; IG = []
    for s in (-1, 0, 1):
        ix = m + s; ok = (ix >= 0) & (ix < M); w = np.where(ok, np.exp(-0.5 * ((gg - T1 - 6 * H * ix) / H) ** 2), 0.0); ix = np.where(ok, ix, 0); WG.append(w); IG.append(ix); N2G += np.bincount(ix, w * w, M)
    dt_ = gg - T1 - 6 * H * m; okh = (m >= 0) & (m < M) & (np.abs(dt_) < A_); wh = np.where(okh, np.cos(np.pi * dt_ / (2 * A_)) ** 2, 0.0); ih = np.where(okh, m, 0); N2H += np.bincount(ih, wh * wh, M)
    for j, kx in enumerate(KX):
        ph = d0 * kx; c = np.cos(ph); sn = np.sin(ph)
        for w, ix in zip(WG, IG): SG[j] += np.bincount(ix, w * c, M) + 1j * np.bincount(ix, w * sn, M)
        SH[j] += np.bincount(ih, wh * c, M) + 1j * np.bincount(ih, wh * sn, M)
    print('  data %d / %d nul, %.0f s' % (min(a + 4000000, len(g)), len(g), time.time() - t0), flush=True)
del g
def Vz(S, N2): return (np.abs(S - S.mean(1, keepdims=True)) ** 2).mean(1) / N2.mean()
Ts = Ta + (np.arange(NH) + 0.5) * (Tb - Ta) / NH; Ls = np.log(Ts / (2 * math.pi))
VBG, VGG = pred(KX, Ls, GW); VBH, VGH = pred(KX, Ls, AH)
nd = len(DEL); out = dict(pasmo=nm, H=H, a=A_, L=Lm, M=M, delta=DEL.tolist(), tau_kontrola=TAUK, B=B)
# R1 Hannovo okno
VH = Vz(SH, N2H); eH = VH / math.sqrt(M); zb = (VH - VBH) / eH; zg = (VH - VGH) / eH; i0 = nd // 2
pod = np.array_split(np.arange(M), 10); ePH = np.array([Vz(SH[:, ix], N2H[ix]) for ix in pod]).std(0, ddof=1) / math.sqrt(10); eH2 = np.maximum(eH, ePH)
out['R1'] = dict(V=VH[:nd].tolist(), err=eH[:nd].tolist(), err_podpasma=ePH[:nd].tolist(), V_BK=VBH[:nd].tolist(), V_GUE=VGH[:nd].tolist(),
                 z_BK_tau1=float(zb[i0]), z_GUE_tau1=float(zg[i0]), z_BK_tau1_max_chyba=float((VH[i0] - VBH[i0]) / eH2[i0]), z_GUE_tau1_max_chyba=float((VH[i0] - VGH[i0]) / eH2[i0]),
                 Q_BK=float(np.mean(zb[:nd] ** 2)), Q_GUE=float(np.mean(zg[:nd] ** 2)), V_tau1=float(VH[i0]), V_BK_tau1=float(VBH[i0]), V_GUE_tau1=float(VGH[i0]))
print('R1 %s Hann: V(τ=1) %.4f ± %.4f, BK %.4f (z %+.2f), GUE %.4f (z %+.1f) | Q_BK %.2f, Q_GUE %.1f' % (nm, VH[i0], eH[i0], VBH[i0], zb[i0], VGH[i0], zg[i0], out['R1']['Q_BK'], out['R1']['Q_GUE']), flush=True)
# R2 odhady chyby v τ = 1 (Gaussovo okno)
VG_ = Vz(SG, N2G); est = {'V/sqrt(M)': VG_ / math.sqrt(M)}
for Kp in (5, 10, 20, 40):
    pod = np.array_split(np.arange(M), Kp); est['podpasma_%d' % Kp] = np.array([Vz(SG[:, ix], N2G[ix]) for ix in pod]).std(0, ddof=1) / math.sqrt(Kp)
bl = np.array_split(np.arange(M), B); Sm = SG.mean(1); N2m = N2G.mean(); J = np.zeros((B, len(KX)))
sumS = SG.sum(1); sumS2 = (np.abs(SG) ** 2).sum(1); sumN = N2G.sum()
for b, ix in enumerate(bl):
    n_ = M - len(ix); s1 = sumS - SG[:, ix].sum(1); s2 = sumS2 - (np.abs(SG[:, ix]) ** 2).sum(1); nn = sumN - N2G[ix].sum()
    J[b] = (s2 / n_ - np.abs(s1 / n_) ** 2) * n_ / n_ / (nn / n_)
Jm = J.mean(0); CJ = (B - 1) / B * (J - Jm).T @ (J - Jm); est['jackknife_%d' % B] = np.sqrt(np.diag(CJ))
out['R2'] = {}
for k_, e in est.items():
    out['R2'][k_] = dict(err_tau1=float(e[i0]), z_BK_tau1=float((VG_[i0] - VBG[i0]) / e[i0]), z_GUE_tau1=float((VG_[i0] - VGG[i0]) / e[i0]))
    print('R2 %s %-14s chyba %.5f → z_BK %+.2f, z_GUE %+.1f' % (nm, k_, e[i0], out['R2'][k_]['z_BK_tau1'], out['R2'][k_]['z_GUE_tau1']), flush=True)
out['R2_V_tau1'] = float(VG_[i0]); out['R2_V_BK_tau1'] = float(VBG[i0]); out['R2_V_GUE_tau1'] = float(VGG[i0])
# R3 kovariance jemného průběhu
C = CJ[:nd, :nd]; rB = VG_[:nd] - VBG[:nd]; rG = VG_[:nd] - VGG[:nd]; Ci = np.linalg.inv(C); p_ = nd
def hot(r):
    T2_ = float(r @ Ci @ r); F = T2_ * (B - p_) / (p_ * (B - 1)); return T2_, F, float(stats.f.sf(F, p_, B - p_))
T2B, FB, pB = hot(rB); T2G, FG, pG = hot(rG); sd = np.sqrt(np.diag(C)); Rk = C / np.outer(sd, sd); ev = np.linalg.eigvalsh(Rk); Nef = float(ev.sum() ** 2 / (ev ** 2).sum())
zz = rB / sd; S2 = float(np.sum(zz ** 2)); zmax = float(np.max(np.abs(zz)))
Lc = np.linalg.cholesky(Rk + 1e-12 * np.eye(nd)); sim = np.abs(rng.standard_normal((20000, nd)) @ Lc.T).max(1); pglob = float(np.mean(sim >= zmax))
out['R3'] = dict(T2_BK=T2B, F_BK=FB, p_BK=pB, T2_GUE=T2G, F_GUE=FG, p_GUE=pG, N_ef=Nef, sum_z2_BK=S2, p_sum_z2_chi2_Nef=float(stats.chi2.sf(S2 * Nef / nd, Nef)),
                 max_abs_z_BK=zmax, p_globalni_max_z=pglob, err_jackknife=sd.tolist())
print('R3 %s: T² BK %.1f (F %.3f, p %.3g) | T² GUE %.3g (p %.3g) | N_ef %.1f | Σz² %.1f | největší |z| %.2f, globální p %.3f' % (nm, T2B, FB, pB, T2G, pG, Nef, S2, zmax, pglob), flush=True)
# R4 kontrolní oblast (popisně)
out['R4'] = []
for i, tau in enumerate(TAUK):
    j = nd + i; e = max(float(est[k_][j]) for k_ in est)
    out['R4'].append(dict(tau=tau, V=float(VG_[j]), err=e, V_BK=float(VBG[j]), V_GUE=float(VGG[j]), z_BK=float((VG_[j] - VBG[j]) / e), z_GUE=float((VG_[j] - VGG[j]) / e)))
    print('R4 %s τ %.1f: V %.4f ± %.4f | BK %.4f (z %+.2f) | GUE %.4f (z %+.1f)' % (nm, tau, VG_[j], e, VBG[j], out['R4'][-1]['z_BK'], VGG[j], out['R4'][-1]['z_GUE']), flush=True)
json.dump(out, open(OUT, 'w')); print('hotovo %.0f s' % (time.time() - t0))
