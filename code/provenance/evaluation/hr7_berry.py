"""Studie VII · Berryho vzorec pro rozptyl počtu nul Σ²(ℓ) jedné L(s, χ) (diagonální přiblížení; kód studie I / VI G3 beze změny, jen bez prvočísel p | q
a s H = ln(qT/2π)). Funkce: berry_q(ell, Hs, q), gue(ell). Spuštěno samostatně vypíše předpovědi a odhad síly pro sadu hr7."""
import math, numpy as np
from scipy.special import sici
def sieve(M):
    s = np.ones(M + 1, bool); s[:2] = False
    for i in range(2, int(M ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
P0 = 10 ** 6; U0 = math.log(P0); PR = sieve(P0).astype(float)
def powers(pr, Hmax):
    lp = np.log(pr); u, w = [], []
    for k in range(1, 80):
        m = k * lp <= Hmax
        if not m.any(): break
        u.append(k * lp[m]); w.append(1 / (k * k * np.exp(k * lp[m])))
    u = np.concatenate(u); w = np.concatenate(w); o = np.argsort(u, kind='stable'); return u[o], w[o]
_C = {}
def tab(q):
    if q not in _C: _C[q] = powers(PR[np.array([q % int(p) != 0 for p in PR])] if q > 1 else PR, 60.0)
    return _C[q]
def berry1(L, H, U, Wt):
    m = U <= H; a = 2 * math.pi * L / H; s = float(np.sum(Wt[m] * (1 - np.cos(a * U[m]))))
    if H > U0: s += (math.log(H) - math.log(U0)) - (sici(a * H)[1] - sici(a * U0)[1])
    si, ci = sici(2 * math.pi * L); return s / math.pi ** 2 + (1 - math.cos(2 * math.pi * L)) / math.pi ** 2 + L * (1 - 2 / math.pi * si)
def berry_q(L, Hs, q):
    U, Wt = tab(q); return float(np.mean([berry1(L, H, U, Wt) for H in Hs]))
def gue(L): si, ci = sici(2 * math.pi * L); return (math.log(2 * math.pi * L) + 0.5772156649015329 + 1 - math.cos(2 * math.pi * L) - ci) / math.pi ** 2 + L * (1 - 2 / math.pi * si)
if __name__ == '__main__':
    LG = np.round(np.logspace(-1, math.log10(200), 30), 6).tolist(); T1, T2 = 36000.0, 80000.0
    for q in (5, 7, 13, 15, 21):
        Hs = [math.log(q * (T1 + (i + 0.5) * (T2 - T1) / 24) / (2 * math.pi)) for i in range(24)]
        N = (T2 - T1) / (2 * math.pi) * np.mean(Hs); B = np.array([berry_q(L, Hs, q) for L in LG]); Bn = np.array([berry_q(L, Hs, 1) for L in LG]); G = np.array([gue(L) for L in LG])
        e = B * np.sqrt(2 * np.array(LG) / N) + 1e-4; i = [LG.index(v) for v in LG if v >= 10]
        print('q %2d N ≈ %.0f | Σ² při ℓ = 1, 10, 100, 200: %s | naivní ζ %s | GUE %s | síla Q(GUE) %.1f, Q(naivní) %.1f' % (q, N, ' '.join('%.3f' % B[min(range(len(LG)), key=lambda j: abs(LG[j]-v))] for v in (1, 10, 100, 200)),
              ' '.join('%.3f' % Bn[min(range(len(LG)), key=lambda j: abs(LG[j]-v))] for v in (10, 100, 200)), ' '.join('%.3f' % G[min(range(len(LG)), key=lambda j: abs(LG[j]-v))] for v in (10, 100, 200)),
              float(np.mean(((G - B) / e) ** 2)), float(np.mean(((Bn - B) / e) ** 2))))
