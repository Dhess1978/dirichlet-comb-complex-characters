"""Studie VII · teorie hřebene jedné Dirichletovy L-funkce L(s, χ) s primitivním charakterem modulu q (počítá se PŘED pohledem na jakákoli data VII).
Předpověď V_BK,q(δ): dvoubodová funkce Bogomolného a Keatinga (J. Phys. A 46, 2013) — funkce ζ s odstraněnými Eulerovými činiteli prvočísel p | q, L = ln(qT/2π);
na hodnotách χ nezávisí. Kód je beze změny převzatý ze studie VI (fz3_c.py: tab_D, V_sled; fz2_c_teorie1.py: F_derivs, interp), jen seznam p | q se bere z rozkladu q
(ve studii VI z pevného seznamu p ≤ 23, pro q ≤ 23 totéž).
Nulové předpovědi: naivní ζ (bez odstranění činitelů p | q, s L = ln(qT/2π)) a GUE bez aritmetiky.
Použití: python3 hr7_teorie.py VYSTUP.json [H] [T1] [T2]"""
import json, math, sys, time, numpy as np
from scipy.special import exp1
t0 = time.time()
OUT = sys.argv[1]; H = float(sys.argv[2]) if len(sys.argv) > 2 else 20.0; T1 = float(sys.argv[3]) if len(sys.argv) > 3 else 36000.0; T2 = float(sys.argv[4]) if len(sys.argv) > 4 else 80000.0
DEL = np.round(np.arange(-2.5, 2.5 + 1e-9, 0.025), 3)
SRC = '/home/claude/br87/py/fz2_c_teorie1.py'
src = open(SRC).read().split('SMAX = ')[0]; argv0 = sys.argv; sys.argv = ['x', str(H)]; ns = {}; exec(src, ns); sys.argv = argv0
EG, PRB, LPB, U, GW, DU, interp, F_derivs, PB = ns['EG'], ns['PRB'], ns['LPB'], ns['U'], ns['GW'], ns['DU'], ns['interp'], ns['F_derivs'], ns['PB']
def pdel(q):
    r, n, p = [], q, 2
    while p * p <= n:
        if n % p == 0:
            r.append(p)
            while n % p == 0: n //= p
        p += 1
    if n > 1: r.append(n)
    return r
def tab_q(q):
    """q_q, h''_q na mřížce EG (totéž co tab_D ze studie VI, fz3_c.py, s PD = pdel(q))"""
    PD = pdel(q); m = ~np.isin(PRB, PD); P_ = PRB[m]; LP_ = LPB[m]; Q = np.zeros(len(EG), complex); H2 = np.zeros(len(EG))
    for j0 in range(0, len(EG), 2000):
        e = EG[j0:j0 + 2000]; F, F1, F2 = F_derivs(e); lnF2 = F2 / F - (F1 / F) ** 2; A2 = np.abs(F) ** 2
        lo = np.zeros(len(e), complex); ld2 = np.zeros(len(e))
        for k0 in range(0, len(e), 256):
            ee = e[k0:k0 + 256]; v = np.exp(1j * np.outer(ee, LP_)); w = (v - 1) ** 2 / (P_ - 1) ** 2
            lo[k0:k0 + 256] = np.log(1 - w).sum(1); z = v / P_; ld2[k0:k0 + 256] = (2 * (LP_ ** 2 * (z * z / (1 - z) ** 2)).real).sum(1)
        Y = math.log(PB); c1 = 1 - 1j * e; c2 = 1 - 2j * e
        lo += -(exp1(c2 * Y) - 2 * exp1(c1 * Y) + exp1(Y + 0j)); ld2 += 2 * (np.exp(-c2 * Y) * (Y / c2 + 1 / c2 ** 2)).real
        s = 1 + 1j * e; ex = np.ones(len(e)); hd = np.zeros(len(e))
        for p in PD:
            wp = p ** (-s); ex *= np.abs(1 - wp) ** 2 / (1 - 1 / p) ** 2; hd += (-math.log(p) ** 2 * wp / (1 - wp) ** 2).real
        Q[j0:j0 + 2000] = A2 * ex * np.exp(lo); H2[j0:j0 + 2000] = -2 * lnF2.real - 2 * hd + ld2
    return interp(Q, U, True), interp(H2, U, False)
QZ = (ns['QU'], ns['H2U']); QG = (np.ones(len(U), complex), np.zeros(len(U)))
def V_sled(dls, Ls, qh):
    q, h2 = qh; acc = np.zeros(len(dls))
    for L in Ls:
        C = np.cos(np.outer(L + np.asarray(dls), U)); Bv = (h2 + 2 * (1 - (np.exp(1j * L * U) * q).real) / (U * U)) * GW; acc += L * (1 - (C @ Bv) * DU / (math.pi * L))
    return acc / np.sum(Ls)
def predikce(q):
    TT = T1 + 9 * H; Tb = T2 - 9 * H; Ts = TT + (np.arange(24) + 0.5) * (Tb - TT) / 24; Ls = np.log(q * Ts / (2 * math.pi)); qq = tab_q(q)
    return float(np.mean(Ls)), V_sled(DEL, Ls, qq), V_sled(DEL, Ls, QZ), V_sled(DEL, Ls, QG)
if __name__ == '__main__':
    MODULY = [5, 7, 13, 15, 21]
    out = dict(H=H, T1=T1, T2=T2, delta=DEL.tolist(), q={}, zdroj_kodu=SRC)
    M = int((T2 - 9 * H - (T1 + 9 * H)) // (6 * H)) + 1
    for q in MODULY:
        Lm, VB, VN, VG = predikce(q); e = VB / math.sqrt(M)
        out['q'][str(q)] = dict(p_deli_q=pdel(q), L_stred=Lm, V_BK_q=VB.tolist(), V_naiv=VN.tolist(), V_GUE=VG.tolist(), M_nezavislych=M,
                                Q_naiv_ocekavane=float(np.mean(((VN - VB) / e) ** 2)), Q_GUE_ocekavane=float(np.mean(((VG - VB) / e) ** 2)))
        print('q %3d (p | q: %s) L %.2f | síla: Q(naivní ζ) %.1f, Q(GUE) %.1f při M = %d | %.0f s' % (q, pdel(q), Lm, out['q'][str(q)]['Q_naiv_ocekavane'], out['q'][str(q)]['Q_GUE_ocekavane'], M, time.time() - t0), flush=True)
    json.dump(out, open(OUT, 'w'))
