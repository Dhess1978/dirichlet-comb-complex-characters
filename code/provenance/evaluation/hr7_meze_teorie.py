"""Studie VII · dodatek 5 · teorie bez dat: (a) hřeben V_BK,q s Eulerovými činiteli odstraněnými pro množinu S ⊆ {2, 3, 5, 7} (slyšet ladění — síla rozlišení S),
(b) prvočíselná mez: V_BK s Φ omezeným na prvočísla ≤ P (P = 2, 3, 5, 7, 11, 13) — popisně. Kód tab_q ze studie VI/hr7_teorie.py, jen seznam odstraněných p je parametr."""
import json, math, sys, time, itertools, numpy as np
from scipy.special import exp1
OUT = sys.argv[1]; argv0 = sys.argv; sys.argv = ['x', '/dev/null', '20', '36000', '80000']; NT = {}
exec(open('/home/claude/br107/py/hr7_teorie.py').read().split("if __name__ == '__main__':")[0], NT); sys.argv = argv0
EG, PRB, LPB, U, GW, DU, interp, F_derivs, PB, V_sled, DEL = NT['EG'], NT['PRB'], NT['LPB'], NT['U'], NT['GW'], NT['DU'], NT['interp'], NT['F_derivs'], NT['PB'], NT['V_sled'], NT['DEL']
def tab_S(PD, PMAX=None):
    """jako tab_q (studie VI), ale odstraněná prvočísla PD zadaná přímo; PMAX: Φ jen s prvočísly ≤ PMAX (prvočíselná mez; ocas nad PB se pak vynechá)"""
    m = ~np.isin(PRB, PD)
    if PMAX is not None: m &= PRB <= PMAX
    P_ = PRB[m]; LP_ = LPB[m]; Q = np.zeros(len(EG), complex); H2 = np.zeros(len(EG))
    for j0 in range(0, len(EG), 2000):
        e = EG[j0:j0 + 2000]; F, F1, F2 = F_derivs(e); lnF2 = F2 / F - (F1 / F) ** 2; A2 = np.abs(F) ** 2
        lo = np.zeros(len(e), complex); ld2 = np.zeros(len(e))
        for k0 in range(0, len(e), 256):
            ee = e[k0:k0 + 256]; v = np.exp(1j * np.outer(ee, LP_)); w = (v - 1) ** 2 / (P_ - 1) ** 2
            lo[k0:k0 + 256] = np.log(1 - w).sum(1); z = v / P_; ld2[k0:k0 + 256] = (2 * (LP_ ** 2 * (z * z / (1 - z) ** 2)).real).sum(1)
        if PMAX is None:
            Y = math.log(PB); c1 = 1 - 1j * e; c2 = 1 - 2j * e
            lo += -(exp1(c2 * Y) - 2 * exp1(c1 * Y) + exp1(Y + 0j)); ld2 += 2 * (np.exp(-c2 * Y) * (Y / c2 + 1 / c2 ** 2)).real
        s = 1 + 1j * e; ex = np.ones(len(e)); hd = np.zeros(len(e))
        for p in PD:
            wp = p ** (-s); ex *= np.abs(1 - wp) ** 2 / (1 - 1 / p) ** 2; hd += (-math.log(p) ** 2 * wp / (1 - wp) ** 2).real
        Q[j0:j0 + 2000] = A2 * ex * np.exp(lo); H2[j0:j0 + 2000] = -2 * lnF2.real - 2 * hd + ld2
    return interp(Q, U, True), interp(H2, U, False)
T1, T2 = 36000.0, 80000.0; H = 20.0; M = 364
def Ls(q): TT = T1 + 9 * H; Tb = T2 - 9 * H; return np.log(q * (TT + (np.arange(24) + 0.5) * (Tb - TT) / 24) / (2 * math.pi))
if __name__ == '__main__':
    t0 = time.time(); out = dict(S={}, mez={}); MAL = [2, 3, 5, 7]; SUBS = [list(s) for r in range(5) for s in itertools.combinations(MAL, r)]
    for q in (5, 7, 13, 15, 21):
        L_ = Ls(q); big = [p for p in (11, 13) if q % p == 0]; V = {}
        for S in SUBS: V[str(S)] = V_sled(DEL, L_, tab_S(S + big))
        true = str([p for p in MAL if q % p == 0]); e = V[true] / math.sqrt(M)
        Qs = {k: float(np.mean(((v - V[true]) / e) ** 2)) for k, v in V.items()}; second = min(v for k, v in Qs.items() if k != true)
        out['S'][str(q)] = dict(pravda=true, Q_ocekavane={k: v for k, v in Qs.items()}, nejblizsi_jina=second, V={k: v.tolist() for k, v in V.items()})
        print('q %2d: pravá S = %s | nejbližší jiná S má očekávané Q = %.1f | %.0f s' % (q, true, second, time.time() - t0), flush=True)
    for q in (5, 15, 21):
        L_ = Ls(q); PD = [p for p in (2, 3, 5, 7, 11, 13) if q % p == 0]; out['mez'][str(q)] = {}
        for P in (2, 3, 5, 7, 11, 13):
            out['mez'][str(q)][str(P)] = V_sled(DEL, L_, tab_S(PD, P)).tolist()
        print('q %2d: prvočíselné meze spočteny | %.0f s' % (q, time.time() - t0), flush=True)
    json.dump(out, open(OUT, 'w'))
