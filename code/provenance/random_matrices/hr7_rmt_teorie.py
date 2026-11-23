"""Studie VII · dodatek 9 (RMT), predpoved pred merenim: parova korelace R2(d) v kosich 0,05 na (0, 3] pro moduly 5, 7, 13, 15, 21
z dvoubodove funkce Bogomolneho a Keatinga (tabulky q(u), h''(u) z hr7_teorie.py, tj. ze studie VI beze zmeny), prumer pres kos a pres 24 vysek pasma s vahou hustoty;
GUE v tychz kosich; popisne efektivni N z chovani R2 u nuly. Pouziti: python3 hr7_rmt_teorie.py VYSTUP.json"""
import json, math, sys, hashlib, datetime, numpy as np
OUT = sys.argv[1]; T1, T2 = 36000.0, 80000.0; argv0 = sys.argv; sys.argv = ['x', '/tmp/_hr7_rmt_nic.json']; NT = {}
exec(open('/home/claude/br107/py/hr7_teorie.py').read().split("if __name__ == '__main__':")[0], NT); sys.argv = argv0
U, DU, tab_q = NT['U'], NT['DU'], NT['tab_q']
HR = np.round(np.arange(0, 3.0001, 0.05), 6); NK = len(HR) - 1; JEM = 40
def R2_unf(d, L, q, h2):
    rho = L / (2 * math.pi); u = d / rho; i = (u - U[0]) / DU; i0 = np.clip(np.floor(i).astype(int), 0, len(U) - 2); f = i - i0
    qq = q[i0] * (1 - f) + q[i0 + 1] * f; hh = h2[i0] * (1 - f) + h2[i0 + 1] * f
    return 1 + (-hh / (4 * math.pi ** 2) - (1 - (np.exp(1j * L * u) * qq).real) / (2 * math.pi ** 2 * u * u)) / rho ** 2
res = dict(cas_UTC=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S'), T1=T1, T2=T2, hranice=HR.tolist(), moduly={})
dd = np.array([[HR[k] + (j + 0.5) * 0.05 / JEM for j in range(JEM)] for k in range(NK)])
res['GUE'] = (1 - np.sinc(dd) ** 2).mean(1).tolist()
for q in (5, 7, 13, 15, 21):
    tq, th = tab_q(q); Ts = T1 + (np.arange(24) + 0.5) * (T2 - T1) / 24; Ls = np.log(q * Ts / (2 * math.pi)); W = Ls / Ls.sum()
    B = sum(w * R2_unf(dd, L, tq, th).mean(1) for w, L in zip(W, Ls))
    Lm = float((W * Ls).sum()); d0 = np.array([0.02, 0.04, 0.06]); r0 = sum(w * R2_unf(d0, L, tq, th) for w, L in zip(W, Ls)); c = r0 * 3 / (math.pi ** 2 * d0 ** 2)
    res['moduly'][str(q)] = dict(BK=B.tolist(), L_stred=Lm, hustota=Lm / (2 * math.pi), R2_u_nuly=r0.tolist(), pomer_ke_GUE_u_nuly=c.tolist(), N_ef=float(1 / math.sqrt(max(1e-12, 1 - c[1]))) if c[1] < 1 else None,
        rozdil_BK_minus_GUE_max=float(np.abs(B - np.array(res['GUE'])).max()), rozdil_v_d_0_5=float(B[9:11].mean() - np.mean(res['GUE'][9:11])))
    print(q, 'L %.3f' % Lm, 'N_ef', res['moduly'][str(q)]['N_ef'], 'max |BK-GUE| %.4f' % res['moduly'][str(q)]['rozdil_BK_minus_GUE_max'], 'u d=0,5: %+.4f' % res['moduly'][str(q)]['rozdil_v_d_0_5'], 'pomer u nuly', np.round(c, 4))
res['skript_sha256'] = hashlib.sha256(open(__file__, 'rb').read()).hexdigest()
json.dump(res, open(OUT, 'w')); print('ulozeno', OUT)
