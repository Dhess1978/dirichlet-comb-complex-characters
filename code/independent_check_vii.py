"""Independent check of study VII from the results files alone (no zeros needed): recomputes the test statistics from the stored
measured quantities with new code, recomputes the predictions that have a closed form (first-order coefficients, cyclotomic union,
Hotelling p-values, the superposition model of the unions with a different quadrature), applies the pass marks of the rules
(data/rules, translated in the comments) and compares everything with the numbers stored by the evaluation scripts.
usage: python3 -I code/independent_check_vii.py REPO          exit code 0 = all checks passed"""
import sys, os, json, math
import numpy as np
from scipy import stats, special
REPO = sys.argv[1]; R = os.path.join(REPO, 'data/results'); J = lambda n: json.load(open(os.path.join(R, n)))
E1, E2, E3 = J('hr7_eval_hr7.json'), J('hr7_eval2_hr7.json'), J('hr7_eval3_final.json'); RES = []
def chk(name, cond, info=''):
    RES.append((name, bool(cond), info)); print('%-4s %-62s %s' % ('ok' if cond else 'FAIL', name, info))
close = lambda a, b, t=1e-9: abs(a - b) <= t * max(1.0, abs(b))
HR = list(E1['hreben']); MOD = {j: E1['hreben'][j]['q'] for j in HR}; VERD = {}
cplx = lambda v: complex(*v) if isinstance(v, list) else complex(v)
# 3.1 comb: Q = mean z^2 over 201 points against V_BK,q, naive zeta and GUE
H = E1['hreben']; Qs = {}
for j in HR:
    h = H[j]; V, e = np.array(h['V']), np.array(h['err']); Qs[j] = {k: float((((V - np.array(h[a])) / e) ** 2).mean()) for k, a in (('Q', 'V_BK_q'), ('Q_naiv', 'V_naiv'), ('Q_GUE', 'V_GUE'))}
chk('3.1 Q, Q_naiv, Q_GUE recomputed (36 numbers)', all(close(Qs[j][k], H[j][k]) for j in HR for k in Qs[j]))
O = all(Qs[j]['Q'] <= 3 for j in HR) and all(Qs[j]['Q_naiv'] >= 2 * Qs[j]['Q'] for j in HR if MOD[j] in (5, 15, 21)) and all(Qs[j]['Q_GUE'] >= 2 * Qs[j]['Q'] for j in HR)
VERD['odv-hr7-hreben'] = 'O' if O else '?'
# 3.2 same modulus, same comb: Q_ab = mean (V_a - V_b)^2/(e_a^2 + e_b^2)
ref = dict(E1['reference_VI']); allc = {j: (np.array(H[j]['V']), np.array(H[j]['err'])) for j in HR}; allc.update({j: (np.array(v['V']), np.array(v['err'])) for j, v in ref.items()})
qab = lambda a, b: float(((allc[a][0] - allc[b][0]) ** 2 / (allc[a][1] ** 2 + allc[b][1] ** 2)).mean())
NS = E1['nastroje']; val = lambda x: x['Q'] if isinstance(x, dict) else x
pr = {k: qab(*k.split(' × ')) for k in NS['dvojice']}; kc = {k: qab(*k.split(' × ')) for k in NS['kontrola']}
chk('3.2 Q_ab of 18 pairs and 3 controls recomputed', len(pr) == 18 and all(close(pr[k], val(NS['dvojice'][k])) for k in pr) and all(close(kc[k], val(NS['kontrola'][k])) for k in kc))
mod_of = lambda j: abs(ref[j]['D']) if j in ref else MOD[j]
chk('3.2 pairs are of the same modulus, controls of different moduli', all(mod_of(a) == mod_of(b) for a, b in (k.split(' × ') for k in pr)) and all(mod_of(a) != mod_of(b) for a, b in (k.split(' × ') for k in kc)))
VERD['odv-hr7-nastroje'] = 'O' if all(v <= 3 for v in pr.values()) and all(v >= 5 for v in kc.values()) else '?'
# 3.3 first-order phase: c_pred = -w chi(p) ln p / sqrt p; chi^2 of 50 bins against the prediction and the conjugate phase
def alt_bounds(Ob, E, Ea, er):     # error of a bin = max(sqrt(E), sd of the null model); the sd is stored only through err, so the chi^2 of an alternative is bracketed
    sd = np.where(er > np.sqrt(E) * (1 + 1e-12), er, 0.0); lo = float(((Ob - Ea) ** 2 / np.maximum(Ea, er ** 2)).sum()); hi = float(((Ob - Ea) ** 2 / np.maximum(Ea, sd ** 2)).sum()); return lo, hi
F = E1['faze']; mx = 0; okp = True; okc = True; okchi = True
for k, v in F.items():
    p = v['p']; cp = -v['w'] * cplx(v['chi_p']) * math.log(p) / math.sqrt(p); okp &= abs(cp - cplx(v['c_pred'])) < 1e-12; mx = max(mx, abs(cplx(v['c']) - cp))
    Ob, Ex, Es, er = (np.array(v[a]) for a in ('O', 'E', 'E_sdruzena', 'err')); c2 = float((((Ob - Ex) / er) ** 2).sum()); lo, hi = alt_bounds(Ob, Ex, Es, er)
    okchi &= close(c2, v['chi2'], 1e-9) and lo * (1 - 1e-9) <= v['chi2_sdruzena'] <= hi * (1 + 1e-9) and close(float(stats.chi2.sf(c2, 49)), v['p_hodnota'], 1e-6)
    if abs(cplx(v['chi_p']).imag) >= 0.5 and p <= 7: okc &= stats.chi2.sf(lo, 49) < 1e-6
    v['_p'] = float(stats.chi2.sf(c2, 49))
chk('3.3 predicted coefficients recomputed from w and chi(p) (56)', okp and len(F) == 56); chk('3.3 chi^2 and p of 56 histograms recomputed, 56 alternatives bracketed', okchi)
chk('3.3 max |c - c_pred| = %.2e (stored %.2e)' % (mx, E1['faze_souhrn']['max_rozdil_c']), close(mx, E1['faze_souhrn']['max_rozdil_c'], 1e-9))
VERD['odv-hr7-faze'] = 'O' if mx <= 3e-4 and all(v['_p'] >= 1e-3 for v in F.values()) and okc else '?'
# 3.4 explicit formula with a Gaussian window: |Z - P| / N_g
d = [abs(cplx(r['Z']) - cplx(r['P'])) / r['Ng'] for r in E1['presne']]; chk('3.4 |Z - P|/N_g of 960 points recomputed, max %.2e' % max(d), len(d) == 960 and all(close(a, r['rozdil'], 1e-6) for a, r in zip(d, E1['presne'])))
VERD['odv-hr7-presne'] = 'O' if max(d) <= 1e-6 else '?'
# 3.5 diagonal regime: z = (V - V_diag)/err
zd = [(r['V'] - r['V_diag']) / r['err'] for j in HR for r in E1['diag'][j]]; chk('3.5 z of 48 numbers recomputed', len(zd) == 48 and all(close(a, r['z'], 1e-9) for a, r in zip(zd, [r for j in HR for r in E1['diag'][j]])))
VERD['odv-hr7-diag'] = 'O' if all(abs(z) <= 3 for z in zd) else '?'
# 3.6 cyclotomic union: c_k(p) = -w_U phi(q) [p^k = 1 mod q] ln p / p^(k/2), x <= 200
phi = {5: 4, 7: 6, 15: 8, 21: 12}; T = E2['teleso']; ok = True; mxt = 0
for q, v in T.items():
    for key, c in v['c'].items():
        p, k = map(int, key.split(',')); pred = -v['w'] * phi[int(q)] * (pow(p, k, int(q)) == 1) * math.log(p) / p ** (k / 2); ok &= abs(pred - c['c_pred']) < 1e-12
        if c['hodnoceno']: mxt = max(mxt, abs(cplx(c['c']) - pred))
hp = [h for v in T.values() for h in v['hist'].values()]; chk('3.6 predictions of 66 coefficients recomputed (Landau + orthogonality)', ok and sum(len(v['c']) for v in T.values()) == 66)
okh = all(close(float((((np.array(h['O']) - np.array(h['E'])) / np.array(h['err'])) ** 2).sum()), h['chi2']) and (lambda b: b[0] * (1 - 1e-9) <= h['chi2_bez'] <= b[1] * (1 + 1e-9))(alt_bounds(*(np.array(h[a]) for a in ('O', 'E', 'E_bez_interference', 'err')))) for h in hp)
chk('3.6 chi^2 of 4 histograms recomputed, the alternative bracketed', okh)
VERD['odv-hr7-teleso'] = 'O' if mxt <= 3e-4 and all(h['p'] >= 1e-3 for h in hp) and all(v['hist']['2']['p_bez'] < 1e-6 for v in T.values()) else '?'
# 3.7 cross-covariance: z = (C - D)/err (real and imaginary parts), |C|/err at ln x = 3
K = E2['krizova']; zs = []; okz = True
for v in K.values():
    for r in v['diag']:
        z = [(r['C'][i] - r['D'][i]) / r['err'][i] for i in (0, 1)]; okz &= all(close(a, b, 1e-6) for a, b in zip(z, r['z'])); zs += z
chk('3.7 144 components z recomputed', okz and len(zs) == 144)
VERD['odv-hr7-krizova'] = 'O' if all(abs(z) <= 3 for z in zs) and all(r['C_k_chybe'] >= 5 for v in K.values() for r in v['diag'] if r['k'] == 3) else '?'
# 3.8 saturation (after addendum 3: against the full two-point function BK_q)
NS8 = E2['nasyceni']; L8 = NS8['L']; q8 = {j: {k: float((((np.array(v['sigma2']) - np.array(v[a])) / np.array(v['err'])) ** 2).mean()) for k, a in (('Q', 'BK_q'), ('Q_GUE', 'gue'), ('Q_naiv', 'BK_naiv'), ('Q_berry_diag', 'berry_q'))} for j, v in L8.items()}
chk('3.8 Q against BK, GUE, naive and Berry recomputed (60)', all(close(q8[j][k], L8[j][k]) for j in q8 for k in q8[j]))
s8 = lambda a, b: float(((np.array(L8[a]['sigma2']) - np.array(L8[b]['sigma2'])) ** 2 / (np.array(L8[a]['err']) ** 2 + np.array(L8[b]['err']) ** 2)).mean())
d8 = {k: s8(*k.split(' × ')) for k in NS8['dvojice']}; chk('3.8 Q_ab of 18 pairs recomputed', len(d8) == 18 and all(close(d8[k], NS8['dvojice'][k]) for k in d8))
VERD['odv-hr7-nasyceni'] = 'O' if all(q8[j]['Q'] <= 3 for j in HR) and all(q8[j]['Q_GUE'] >= 2 * q8[j]['Q'] for j in HR) and all(q8[j]['Q_naiv'] >= 2 * q8[j]['Q'] for j in HR if MOD[j] in (15, 21)) and all(v <= 3 for v in d8.values()) else '?'
# 3.10 Hann window for the L-functions: Hotelling p from F with B = M = 364
OL = E2['okno_L']; pf = {j: float(stats.f.sf(v['F'], 201, v['M'] - 201)) for j, v in OL.items()}; chk('3.10 Hotelling p of 15 L-functions recomputed from F', all(close(pf[j], OL[j]['p_BK'], 1e-6) for j in OL))
VERD['odv-hr7-okno-L'] = 'O' if all(OL[j]['Q'] <= 3 and pf[j] >= 1e-3 for j in HR) else '?'
# part R: Hann window, six error estimates, full covariance (B = 1000)
PR = {n: J('hr7_r_%s.json' % n) for n in ('P8', 'P3e10')}; okR = True; okw = okp = okk = True
for n, P in PR.items():
    a = P['R1']; V, e = np.array(a['V']), np.array(a['err']); i0 = P['delta'].index(0.0)
    okR &= close(float((((V - np.array(a['V_BK'])) / e) ** 2).mean()), a['Q_BK']) and close(float((((V - np.array(a['V_GUE'])) / e) ** 2).mean()), a['Q_GUE']) and close((V[i0] - a['V_BK'][i0]) / e[i0], a['z_BK_tau1'])
    okw &= abs(a['z_BK_tau1']) <= 2 and a['z_GUE_tau1'] >= 10 and a['Q_BK'] <= 3
    for k, v in P['R2'].items():
        okR &= close((P['R2_V_tau1'] - P['R2_V_BK_tau1']) / v['err_tau1'], v['z_BK_tau1']) and close((P['R2_V_tau1'] - P['R2_V_GUE_tau1']) / v['err_tau1'], v['z_GUE_tau1']); okp &= abs(v['z_BK_tau1']) <= 2 and v['z_GUE_tau1'] >= 10
    r3 = P['R3']; B = P['B']; fb = r3['T2_BK'] * (B - 201) / (201 * (B - 1)); okR &= close(fb, r3['F_BK'], 1e-9) and close(float(stats.f.sf(fb, 201, B - 201)), r3['p_BK'], 1e-6)
    pg = float(stats.f.sf(r3['T2_GUE'] * (B - 201) / (201 * (B - 1)), 201, B - 201)); okk &= r3['p_BK'] >= 0.01 and pg < 1e-10
chk('R Q, z at tau = 1, six estimates, Hotelling F and p recomputed', okR)
VERD.update({'odv-hr7-okno': 'O' if okw else '?', 'odv-hr7-podpasma': 'O' if okp else '?', 'odv-hr7-kovar': 'O' if okk else '?'})
# 3.11-3.13 blind tests
NJ, LD, T2 = E3['nastroj'], E3['ladeni'], E3['teleso2']
chk('3.11 chosen = truth and score below the next candidate (16)', len(NJ) == 16 and all(v['vybrano'] == v['pravda'] and v['S'] < v['druhy'][0] for v in NJ.values()))
pdiv = lambda q: str([p for p in (2, 3, 5, 7) if q % p == 0])
chk('3.12 best subset recomputed = primes dividing q (16)', len(LD) == 16 and all(min(v['Q'], key=v['Q'].get) == v['vybrano'] == pdiv(v['q']) for v in LD.values()))
def order(p, q):
    k = 1
    while pow(p, k, q) != 1: k += 1
    return k
okf = all(T2[q]['f_odhad'][p] == order(int(p), int(q)) for q in T2 for p in T2[q]['f_pred']); n16 = sum(len(T2[q]['f_pred']) for q in T2)
okc2 = True; mx2 = 0
for q, v in T2.items():
    for key, c in v['c'].items():
        p, k = map(int, key.split(',')); pred = -v['w'] * len(v['clenove']) * (pow(p, k, int(q)) == 1) * math.log(p) / p ** (k / 2); okc2 &= abs(pred - c['c_pred']) < 1e-12
        if c['x'] <= 200 and q in ('15', '21'): mx2 = max(mx2, abs(cplx(c['c']) - pred))
chk('3.13 16 orders = multiplicative orders (recomputed)', okf and n16 == 16); chk('3.13 predictions of the union coefficients recomputed, max %.2e' % mx2, okc2 and close(mx2, E3['teleso2_souhrn']['max_rozdil_15_21'], 1e-9))
VERD.update({'odv-hr7-nastroj': 'O' if all(v['vybrano'] == v['pravda'] for v in NJ.values()) else '?', 'odv-hr7-ladeni': 'O' if all(v['vybrano'] == pdiv(v['q']) for v in LD.values()) else '?',
             'odv-hr7-teleso2': 'O' if okf and mx2 <= 3e-4 else '?'})
# addendum 9: M1 pair correlation chi^2 per bin from the stored R2 and its ten blocks; M2 with an independent superposition model
M, V1, TH = J('hr7_rmt_mereni.json'), J('hr7_rmt_vysledek.json'), J('hr7_rmt_teorie.json'); sem = lambda B: np.array(B).std(0, ddof=1) / math.sqrt(len(B)); g = np.array(TH['GUE']); ok1 = True; nb = ng = 0
for j, d in M['M1'].items():
    r = np.array(d['parova']['R2']); e = sem(d['parova']['bloky']); bk = np.array(TH['moduly'][str(d['q'])]['BK']); cb = float((((r - bk) / e) ** 2).mean()); cg = float((((r - g) / e) ** 2).mean())
    ok1 &= close(cb, V1['M1'][j]['chi2_BK'], 1e-9) and close(cg, V1['M1'][j]['chi2_GUE'], 1e-9); nb += cb <= 2.0; ng += (cg >= 3 * cb and cg >= 4)
chk('M1 chi^2 per bin against [8] and GUE recomputed (32)', ok1 and len(M['M1']) == 16); VERD['odv-hr7-rmt-par'] = 'O' if nb >= 14 and ng >= 14 else '?'
X, W = special.roots_legendre(64)                                      # 64 nodes (the evaluation used numpy's 40)
def Egue(s):
    if s <= 0: return 1.0
    t = (X + 1) * s / 2; w = W * s / 2; K = np.sinc(t[:, None] - t[None, :])
    return float(np.linalg.det(np.eye(len(t)) - np.sqrt(w)[:, None] * K * np.sqrt(w)[None, :]))
def dE(E, s, h=2e-3):                                                   # five-point derivative
    return -1.0 if s <= 0 else (E(s - 2 * h) - 8 * E(s - h) + 8 * E(s + h) - E(s + 2 * h)) / (12 * h)
HU = np.round(np.arange(41) * 0.1, 10); ok2 = True; dmax = 0; okm2 = True
for U, d in M['M2'].items():
    f = d['podily']; ES = lambda s: float(np.prod([Egue(fi * s) for fi in f])); de = np.array([dE(ES, s) for s in HU]); pS = np.diff(de) / 0.1
    st = V1['M2'][U]; dmax = max(dmax, float(np.abs(pS - np.array(st['p_S'])).max())); p = np.array(st['p']); e = np.array(st['chyba'])
    ok2 &= close(float(np.abs(p - np.array(st['p_S'])).max()), st['max_rozdil_od_S']) and close(float((((p - np.array(st['p_S'])) / e) ** 2).mean()), st['chi2_S'])
    okm2 &= float(np.abs(p - pS).max()) <= 0.015 and (-np.expm1(-0.1) / 0.1 - p[0]) >= 0.03 and st['max_rozdil_od_jednoho_GUE'] >= 0.2
chk('M2 superposition model recomputed with 64 nodes, largest difference %.1e' % dmax, dmax < 1e-4); chk('M2 largest |p - p_S| and chi^2 recomputed', ok2)
VERD['odv-hr7-rmt-orchestr'] = 'O' if okm2 else '?'
# addendum 10: amplitude A along (BK - GUE), size V = A * rms, pooled over the four L-functions, ratio
O3, V3 = J('hr7_rmt_vyska_obr.json'), J('hr7_rmt_vyska_vysledek.json'); T0, T1 = J('hr7_rmt_teorie_vyska_036.json'), J('hr7_rmt_teorie_vyska_400.json'); g0 = np.array(T0['GUE'])
def amp(par, bk, gue):
    r = np.array(par['R2']); B = np.array(par['bloky']); e = sem(B)
    e = np.where(e > 0, e, np.sqrt(np.maximum(bk, 1e-4) / (par['n'] / len(B) * 0.05) / len(B))); x = bk - gue; w = x / e ** 2; jm = float((x * w).sum())
    A = float((r - gue) @ w / jm); Ab = np.array([float((np.array(b) - gue) @ w / jm) for b in B]); return A, float(Ab.std(ddof=1) / math.sqrt(len(B))), float(np.sqrt((x ** 2).mean()))
okA = True; acc = {'nahore': [], 'dole': []}
for j, x in V3['L'].items():
    q = str(x['q']); A1, s1, r1 = amp(O3['L'][j]['parova'], np.array(T1['moduly'][q]['BK']), g0); A0, s0, r0 = amp(M['M1'][j]['parova'], np.array(T0['moduly'][q]['BK']), g0)
    okA &= close(A1, x['A_nahore'], 1e-8) and close(s1, x['chyba_A_nahore'], 1e-6) and close(A0, x['A_dole'], 1e-8) and close(r1, x['rms_teorie_nahore'], 1e-8)
    acc['nahore'].append((A1, s1, r1)); acc['dole'].append((A0, s0, r0))
def pool(L):
    w = np.array([1 / s ** 2 for _, s, _ in L]); A = float(sum(a * wi for (a, _, _), wi in zip(L, w)) / w.sum()); sA = float(1 / math.sqrt(w.sum()))
    Vv = np.array([a * r for a, _, r in L]); sv = np.array([s * r for _, s, r in L]); wv = 1 / sv ** 2; return A, sA, float((Vv * wv).sum() / wv.sum()), float(1 / math.sqrt(wv.sum()))
An, sAn, Vn, sVn = pool(acc['nahore']); Ad, sAd, Vd, sVd = pool(acc['dole']); SD = V3['sdruzeno']; ratio = Vn / Vd; sr = ratio * math.hypot(sVn / Vn, sVd / Vd)
chk('M3 amplitudes A of 8 measurements recomputed from the blocks', okA)
chk('M3 pooled A %.4f +- %.4f, V %.5f -> %.5f, ratio %.3f +- %.3f' % (An, sAn, Vd, Vn, ratio, sr), close(An, SD['A_nahore'], 1e-3) and close(Vn, SD['V_nahore'], 1e-3) and close(Vd, SD['V_dole'], 1e-3) and abs(ratio - SD['pomer_V']) < 2e-3)
okv = sum(1 for x in V3['L'].values() if x['chi2_BK_nahore'] <= 2.0) >= 3 and abs(An - 1) <= 3 * sAn and (Vd - Vn) >= 3 * math.hypot(sVn, sVd) and abs(ratio - SD['pomer_teorie']) <= 3 * sr
VERD['odv-hr7-rmt-vyska'] = 'O' if okv else '?'
# all verdicts
st = dict(**E1['verdikty'], **E2['verdikty'], **J('hr7_r_verdikty.json')['verdikty'], **E3['verdikty'], **{'odv-hr7-rmt-par': 'O' if V1['M1_verdikt']['O'] else 'Z', 'odv-hr7-rmt-orchestr': 'O' if V1['M2_verdikt']['O'] else 'Z', 'odv-hr7-rmt-vyska': SD['verdikt']})
chk('18 verdicts recomputed from the pass marks = stored verdicts (all O)', len(VERD) == 18 and VERD == st and all(v == 'O' for v in VERD.values()), {k: v for k, v in VERD.items() if v != st.get(k)})
nf = sum(1 for _, c, _ in RES if not c); print('\n%d checks, %d failed' % (len(RES), nf)); sys.exit(1 if nf else 0)
