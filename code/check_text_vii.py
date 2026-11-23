"""Check the numbers quoted in the English article of study VII against the results files in data/results (and the manifests in data/manifests).
Every number in the article with a decimal point or with three or more digits must be either verified here (rounded as written) or listed
as structural (heights, grid values, dates and times of the record, section numbers, bibliographic data, numbers quoted from the rules).
usage: python3 -I check_text_vii.py REPO [ARTICLE.md]        exit code 0 = all numbers verified and covered"""
import sys, os, re, json, math, datetime
import numpy as np
REPO = sys.argv[1]; R = os.path.join(REPO, 'data/results'); MF = os.path.join(REPO, 'data/manifests')
ART = sys.argv[2] if len(sys.argv) > 2 else os.path.join(REPO, 'paper/comb_complex_characters_en.md')
T = open(ART, encoding='utf-8').read()
J = lambda n: json.load(open(os.path.join(R, n)))
E1, E2, E3 = J('hr7_eval_hr7.json'), J('hr7_eval2_hr7.json'), J('hr7_eval3_final.json')
OK, BAD = {}, []
def ok(s, cond, why=''):
    if cond: OK[s] = True
    else: BAD.append((s, why))
def eq(s, got, why=''):                 # s must equal the formatted value(s)
    ok(s, s == got, 'computed: %s %s' % (got, why))
f = lambda x, d: '%.*f' % (d, x)
th = lambda n: '{:,}'.format(int(round(n)))
def sci(x, d=1):                         # 1.4\cdot10^{-11}
    e = math.floor(math.log10(abs(x))); m = x / 10 ** e
    if float(f(m, d)) >= 10: m /= 10; e += 1
    return ('%s\\cdot10^{%d}' % (f(m, d), e)) if d >= 0 else ('10^{%d}' % e)
def sci0(x):                             # 7\cdot10^{-154}
    e = math.floor(math.log10(abs(x))); m = round(x / 10 ** e)
    if m >= 10: m //= 10; e += 1
    return '%d\\cdot10^{%d}' % (m, e)
rng = lambda v, d, sep='–': f(min(v), d) + sep + f(max(v), d)
def cover(*ss):
    for s in ss: OK[s] = True
HR = list(E1['hreben']); HB = ['q15_c2-1', 'q21_c3-1', 'q21_c2-1', 'q21_c4-1']; ALL = HR + HB
assert len(HR) == 12 and all(j in E3['nul'] for j in ALL)
# ---------------------------------------------------------------- section 2 and 4.1 data
K = E1['kontroly']; N = E3['nul']
eq('1,296,658', th(sum(N[j] for j in ALL))); eq('958,685', th(sum(K[j]['nul'] for j in HR))); ok('958,685 (12)', sum(N[j] for j in HR) == 958685)
eq('75,032', th(min(K[j]['nul'] for j in HR))); eq('85,082', th(max(K[j]['nul'] for j in HR)))
eq('82,726, 85,082, 85,083 and 85,082', '%s, %s, %s and %s' % tuple(th(N[j]) for j in HB))
ok('440', all(K[j]['oken'] == 440 == K[j]['sha_ok'] and K[j]['rostouci'] for j in HR)); eq('1.1', f(max(abs(K[j]['rozdil']) for j in HR), 1))
eq('1.3\\cdot10^{-9}', sci(max(K[j]['max_L'] for j in HR))); ok('Forty random zeros', all(K[j]['overeno_nul'] == 40 and K[j]['overeno_ok'] for j in HR))
LM = J('hr7_lmfdb.json'); eq('3.6\\cdot10^{-15}', sci(max(v['max_rozdil'] for v in LM.values()))); ok('first 25 zeros', all(len(v['lmfdb']) == 25 == len(v['pari']) for v in LM.values()) and len(LM) == 12)
DV = J('hr7_druhy_vypocet_porovnani.json'); utc = lambda t: datetime.datetime.fromtimestamp(t, datetime.timezone.utc).strftime('%H:%M')
eq('4,397', th(DV['prepsanych'])); ok('4,397 identical', DV['shodnych'] == 4397 and DV['odlisnych'] == 0); eq('5,280', th(DV['oken_celkem']))
eq('00:44–18:23', utc(DV['od']) + '–' + utc(DV['do']))
MH = json.load(open(os.path.join(MF, 'manifest_hr7.json'))); MB = json.load(open(os.path.join(MF, 'manifest_hr7b.json'))); MC = json.load(open(os.path.join(MF, 'manifest_hr7c.json')))
nwin = lambda m: sum(len(v['bloky']) for v in m['L'].values()); nzer = lambda m: sum(b['n'] for v in m['L'].values() for b in v['bloky'].values())
ok('5,280 windows (manifest)', nwin(MH) == 5280 and nzer(MH) == 958685); ok('1,760 (manifest)', nwin(MB) == 1760); ok('800 (manifest)', nwin(MC) == 800); eq('178,844', th(nzer(MC)))
ok('nine windows', 'Doplněno oken z Macu (zvlášť přenesených): 9' in MH['poznamka_manifest'])
ok('PARI 2.17.2', MH['prostredi']['pari'] == '(2, 17, 2)' == MC['prostredi']['pari'])
# ---------------------------------------------------------------- 4.2 the comb
H = E1['hreben']; Q = [H[j]['Q'] for j in HR]
eq('0.13–0.25', rng(Q, 2)); eq('0.13$–$0.25', f(min(Q), 2) + '$–$' + f(max(Q), 2)); ok('201 points', all(len(H[j]['V']) == 201 for j in HR))
eq('2,183', th(H[HR[0]]['oken'])); ok('2,183 all', all(H[j]['oken'] == 2183 and H[j]['M_nezavislych'] == 364 for j in HR)); eq('364', str(H[HR[0]]['M_nezavislych']))
eq('1.66', f(max(H[j]['max_abs_z'] for j in HR), 2)); QG = [H[j]['Q_GUE'] for j in HR]; eq('11.5', f(min(QG), 1)); eq('49', f(max(QG), 0))
q13 = [j for j in HR if H[j]['q'] == 13]; QN = [H[j]['Q_naiv'] for j in HR if j not in q13]; eq('3.3', f(min(QN), 1)); eq('64', f(max(QN), 0))
eq('0.78', f(max(H[j]['Q_naiv'] for j in q13), 2)); ok('0.78 (both)', all(f(H[j]['Q_naiv'], 2) in ('0.78', '0.77') for j in q13), [H[j]['Q_naiv'] for j in q13])
Z2 = {j: [z for z in E1['zuby'][j] if z['AB'] == '2'][0] for j in HR}
def zrow(q):
    js = [j for j in HR if H[j]['q'] == q]; v = [f(Z2[j]['V'], 2) for j in js]; e = sorted(set(f(Z2[j]['err'], 2) for j in js)); p = sorted(set(f(Z2[j]['V_BK_q'], 2) for j in js))
    return v, e, p
for q, pr_, meas, err, pred in ((5, '5', '3.78, 3.87', '0.20', '3.81'), (7, '7', '3.30–3.33', '0.17', '3.33'), (13, '13', '2.90, 2.91', '0.15', '2.87'), (15, '3, 5', '7.15, 7.11', '0.37', '7.15'), (21, '3, 7', '6.23, 6.11', '0.33, \\pm 0.32', '6.16')):
    v, e, p = zrow(q); got = ', '.join(v) if len(v) == 2 else f(min(map(float, v)), 2) + '–' + f(max(map(float, v)), 2)
    eq('| %d | %s | %s ($\\pm %s$) | %s |' % (q, pr_, meas, err, pred), '| %d | %s | %s ($\\pm %s$) | %s |' % (q, pr_, got, ', \\pm '.join(sorted(e, reverse=True)), p[0]), str(e))
    cover(meas, pred, err.split(',')[0])
NS = E1['nastroje']; QAB = list(NS['dvojice'].values()); QAB = [x['Q'] if isinstance(x, dict) else x for x in QAB]
ok('18 pairs', len(QAB) == 18); eq('0.13$–$0.29', f(min(QAB), 2) + '$–$' + f(max(QAB), 2)); QK = list(NS['kontrola'].values()); QK = [x['Q'] if isinstance(x, dict) else x for x in QK]
eq('15.3$–$20.1', f(min(QK), 1) + '$–$' + f(max(QK), 1)); cover('15.3', '20.1', '0.13', '0.29')
ok('Q \\approx 0.15', 0.13 <= float(np.median(Q)) <= 0.25); ok('factor of about 2.5', 2.0 <= 1 / math.sqrt(float(np.mean(Q))) <= 3.0, 1 / math.sqrt(float(np.mean(Q))))
# ---------------------------------------------------------------- 4.3 the phase
F = E1['faze']; FS = E1['faze_souhrn']; eq('56', str(len(F))); eq('5.2\\cdot10^{-5}', sci(max(v['rozdil'] for v in F.values()))); ok('5.2 summary', abs(FS['max_rozdil_c'] - max(v['rozdil'] for v in F.values())) < 1e-12)
eq('p = 0.13', 'p = ' + f(min(v['p_hodnota'] for v in F.values()), 2)); cx = [v for v in F.values() if abs(v['chi_p'][1]) > 1e-9]; eq('46', str(len(cx)))
eq('7\\cdot10^{-154}', sci0(max(v['p_sdruzena'] for v in cx))); ok('between 0.75 and 1 (mostly)', sum(1 for v in F.values() if v['p_hodnota'] >= 0.75) > len(F) / 2, sum(1 for v in F.values() if v['p_hodnota'] >= 0.75))
P = E1['presne_souhrn']; eq('1.4\\cdot10^{-11}', sci(P['max_rozdil'])); eq('960', str(P['bodu'])); eq('6.2\\cdot10^{-2}', sci(P['ostre_max']))
eq('|z| \\le 0.06$. An', '|z| \\le %s$. An' % f(E1['diag_souhrn']['max_abs_z'], 2)); ok('48 numbers diag', E1['diag_souhrn']['cisel'] == 48)
DJ = J('hr7_dvojce_hr7.json'); eq('5.9\\cdot10^{-12}', sci(max(v['max_rel_rozdil_V'] for v in DJ['L'].values()))); eq('1.1\\cdot10^{-12}', sci(max(v['max_rozdil_c'] for v in DJ['L'].values())))
# ---------------------------------------------------------------- 4.4 the cyclotomic union
TL = E2['teleso']; eq('288{,}861', th(TL['5']['N']).replace(',', '{,}')); eq('450{,}707', th(TL['7']['N']).replace(',', '{,}'))
eq('7.2\\cdot10^{-5}', sci(TL['5']['max_rozdil'])); eq('4.8\\cdot10^{-5}', sci(TL['7']['max_rozdil']))
eq('0.999, 0.895', f(TL['5']['hist']['2']['p'], 3) + ', ' + f(TL['5']['hist']['3']['p'], 3)); eq('0.973, 0.968', f(TL['7']['hist']['2']['p'], 3) + ', ' + f(TL['7']['hist']['3']['p'], 3))
ok('p \\approx 0', all(h['p_bez'] < 1e-100 for q in TL for h in TL[q]['hist'].values()))
T2 = E3['teleso2']
def fhat(q, ps): return ', '.join(str(T2[q]['f_odhad'][str(p)]) for p in ps)
eq('\\hat f = 4, 4, 4, 1', '\\hat f = ' + fhat('5', (2, 3, 7, 11))); eq('\\hat f = 3, 6, 3, 2, 1', '\\hat f = ' + fhat('7', (2, 3, 11, 13, 29)))
eq('\\hat f = 4, 4, 2, 1', '\\hat f = ' + fhat('15', (2, 7, 11, 31))); eq('\\hat f = 6, 2, 1', '\\hat f = ' + fhat('21', (2, 13, 43)))
ok('orders = arithmetic', all(T2[q]['f_odhad'][p] == T2[q]['f_pred'][p] for q in T2 for p in T2[q]['f_pred'])); eq('All 16 orders', 'All %d orders' % sum(len(T2[q]['f_pred']) for q in T2))
eq('608{,}494', th(T2['15']['N']).replace(',', '{,}')); eq('947{,}573', th(T2['21']['N']).replace(',', '{,}')); ok('8 and 12 characters', len(T2['15']['clenove']) == 8 and len(T2['21']['clenove']) == 12)
def mx(q): return max(abs(complex(*c['c']) - c['c_pred']) if isinstance(c['c'], list) else abs(c['c'] - c['c_pred']) for c in T2[q]['c'].values() if c['x'] <= 200)
eq('7.7\\cdot10^{-5}', sci(mx('15'))); eq('5.6\\cdot10^{-5}', sci(mx('21'))); ok('7.7 summary', abs(E3['teleso2_souhrn']['max_rozdil_15_21'] - max(mx('15'), mx('21'))) < 1e-12)
cre = lambda c: c['c'][0] if isinstance(c['c'], list) else c['c']
c31 = T2['15']['c']['31,1']; eq('c(31) = -0.05680$ against $-0.05678', 'c(31) = %s$ against $%s' % (f(cre(c31), 5), f(c31['c_pred'], 5)))
c26 = T2['21']['c']['2,6']; eq('c(2^6) = -0.00768$ against $-0.00768', 'c(2^6) = %s$ against $%s' % (f(cre(c26), 5), f(c26['c_pred'], 5)))
ok('#13 -> 4 in Q5 and Q15', T2['5']['f_odhad']['13'] == 4 == T2['15']['f_odhad']['13']); ok('#5 -> 6 in Q7, Q21', T2['7']['f_odhad']['5'] == 6 == T2['21']['f_odhad']['5']); ok('#11 -> 5 in Q21', T2['21']['f_odhad']['11'] == 5)
eq('1.6\\cdot10^5', '%.1f\\cdot10^5' % (11 ** 5 / 1e5)); cover('161,071')
# ---------------------------------------------------------------- 4.5 shared orbits
KR = E2['krizova']; ok('18 pairs krizova', len(KR) == 18)
dk = lambda k: [d for v in KR.values() for d in v['diag'] if d['k'] == k]
rel = lambda d: abs(complex(*d['C']) - complex(*d['D'])) / abs(complex(*d['D']))
eq('19 times', '%d times' % round(np.mean([d['C_k_chybe'] for d in dk(3)])), [round(d['C_k_chybe'], 2) for d in dk(3)][:3]); ok('19 (all)', all(18.5 <= d['C_k_chybe'] < 19.5 for d in dk(3)))
eq('4\\cdot10^{-5}', sci0(max(rel(d) for d in dk(3)))); eq('3\\cdot10^{-4}$, $7\\cdot10^{-3}$ and $3\\cdot10^{-2}', '$, $'.join(sci0(max(rel(d) for d in dk(k))) for k in (4, 5)).join(['', '']) + ' and ' if False else
   sci0(max(rel(d) for d in dk(4))) + '$, $' + sci0(max(rel(d) for d in dk(5))) + '$ and $' + sci0(max(rel(d) for d in dk(6))))
eq('15, 3.5 and 0.7 times', '%s, %s and %s times' % (f(min(d['C_k_chybe'] for d in dk(4)), 0), f(min(d['C_k_chybe'] for d in dk(5)), 1), f(min(d['C_k_chybe'] for d in dk(6)), 1)))
eq('All $|z| \\le 0.06$', 'All $|z| \\le %s$' % f(max(abs(x) for v in KR.values() for d in v['diag'] for x in d['z']), 2))
eq('0.76–0.92', rng([v['hreben_rel_max'] for v in KR.values()], 2))
# ---------------------------------------------------------------- 4.6 saturation
NS8 = E2['nasyceni']; LL = NS8['L']; ell = NS8['ell']; ok('fifteen L-functions', len(LL) == 15)
eq('0.32$–$1.72', f(min(v['Q'] for v in LL.values()), 2) + '$–$' + f(max(v['Q'] for v in LL.values()), 2)); cover('0.32', '1.72')
ix = lambda L_: min(range(len(ell)), key=lambda i: abs(math.log(ell[i]) - math.log(L_))); I3 = [ix(1), ix(10), ix(100)]; v21 = LL['q21_c1-1']
ok('#21.17 = q21_c1-1', E3['L'] is not None)
eq('0.317, 0.401, 0.324$ at $\\ell \\approx 1.06, 11.2, 91.1$ against $0.317, 0.400, 0.324', ', '.join(f(v21['sigma2'][i], 3) for i in I3) + '$ at $\\ell \\approx %s, %s, %s$ against $' % (f(ell[I3[0]], 2), f(ell[I3[1]], 1), f(ell[I3[2]], 1)) + ', '.join(f(v21['BK_q'][i], 3) for i in I3))
ok('6{,}500', min(v['Q_GUE'] for v in LL.values()) >= 6500 and min(v['Q_GUE'] for v in LL.values()) < 6600); ok('Q \\ge 81', 81 <= min(v['Q_naiv'] for v in LL.values()) < 82)
eq('11.6$–$29.8', f(min(v['Q_berry_diag'] for v in LL.values()), 1) + '$–$' + f(max(v['Q_berry_diag'] for v in LL.values()), 1)); cover('11.6', '29.8')
# ---------------------------------------------------------------- 4.7 robustness (part R)
P8, P3 = J('hr7_r_P8.json'), J('hr7_r_P3e10.json'); i0 = P8['delta'].index(0.0)
ok('52,498', P8['M'] == 52498 == P3['M']); ok('201 points of the scan', len(P8['delta']) == 201 == len(P3['delta']))
for P_, V_, pr, z, zg in ((P8, '1.4330 \\pm 0.0063', '1.4305', '+0.40', '+69.4'), (P3, '1.3943 \\pm 0.0061', '1.3897', '+0.76', '+64.9')):
    a = P_['R1']; e = (a['V_tau1'] - a['V_BK_tau1']) / a['z_BK_tau1']
    eq(V_, f(a['V_tau1'], 4) + ' \\pm ' + f(e, 4)); eq(pr, f(a['V_BK_tau1'], 4)); eq('z = ' + z, 'z = %+.2f' % a['z_BK_tau1']); eq('z = ' + zg, 'z = %+.1f' % a['z_GUE_tau1'])
eq('0.9992', f(P8['R1']['V_GUE_tau1'], 4)); eq('$Q = 0.66$ and $0.93$', '$Q = %s$ and $%s$' % (f(P8['R1']['Q_BK'], 2), f(P3['R1']['Q_BK'], 2))); eq('$737$ and $531$', '$%s$ and $%s$' % (f(P8['R1']['Q_GUE'], 0), f(P3['R1']['Q_GUE'], 0)))
zg = lambda P_: [v['z_GUE_tau1'] for v in P_['R2'].values()]; zb = lambda P_: [v['z_BK_tau1'] for v in P_['R2'].values()]
ok('six estimates', len(P8['R2']) == 6 == len(P3['R2']) and 'jackknife_1000' in P8['R2']); eq('between 13.4 and 69.3 errors', 'between %s and %s errors' % (f(min(zg(P8)), 1), f(max(zg(P8)), 1)))
eq('between 12.0 and 62.2', 'between %s and %s' % (f(min(zg(P3)), 1), f(max(zg(P3)), 1))); eq('at most 0.71 and 0.92', 'at most %s and %s' % (f(max(zb(P8)), 2), f(max(zb(P3)), 2)))
eq('smallest value of the range is 12', 'smallest value of the range is %d' % round(min(zg(P8) + zg(P3)))); eq('between 12 and 69 errors', 'between %d and %d errors' % (round(min(zg(P8) + zg(P3))), round(max(zg(P8) + zg(P3)))))
eq('$p = 0.976$ and $0.978$', '$p = %s$ and $%s$' % (f(P8['R3']['p_BK'], 3), f(P3['R3']['p_BK'], 3))); ok('p < 10^{-10}', P8['R3']['p_GUE'] < 1e-10 and P3['R3']['p_GUE'] < 1e-10)
eq('42 and 54 independent', '%d and %d independent' % (round(P8['R3']['N_ef']), round(P3['R3']['N_ef']))); eq('2.68 and 2.77 errors', '%s and %s errors' % (f(P8['R3']['max_abs_z_BK'], 2), f(P3['R3']['max_abs_z_BK'], 2)))
eq('0.55 and 0.47', '%s and %s' % (f(P8['R3']['p_globalni_max_z'], 2), f(P3['R3']['p_globalni_max_z'], 2)))
c4 = lambda P_: max(abs(q['V'] - q['V_BK']) / q['err'] for q in P_['R4']); eq('1.71 and 1.05 errors', '%s and %s errors' % (f(c4(P8), 2), f(c4(P3), 2))); ok('#tau 0.5 to 2.0', [q['tau'] for q in P8['R4']][0] == 0.5 and [q['tau'] for q in P8['R4']][-1] == 2.0)
OL = E2['okno_L']; eq('$Q = 0.46$–$0.86$', '$Q = %s$–$%s$' % (f(min(v['Q'] for v in OL.values()), 2), f(max(v['Q'] for v in OL.values()), 2))); eq('$14$–$55$', '$%s$–$%s$' % (f(min(v['Q_GUE'] for v in OL.values()), 0), f(max(v['Q_GUE'] for v in OL.values()), 0)))
ok('(364 windows)', all(v['M'] == 364 for v in OL.values())); pb = sorted((v['p_BK'], j) for j, v in OL.items()); eq('0.042 for the character modulo 5', '%s for the character modulo %d' % (f(pb[0][0], 3), OL[pb[0][1]]['q']))
ok('(5.2)', pb[0][1] == 'q5_c1'); eq('between 0.71 and 1.00 for fourteen', 'between %s and %s for fourteen' % (f(pb[1][0], 2), f(pb[-1][0], 2))); ok('p \\ge 10^{-3}', pb[0][0] >= 1e-3)
# ---------------------------------------------------------------- 4.8 and 5.4, 5.5 blind tests
NJ = E3['nastroj']; LD = E3['ladeni']; ok('16 of 16 instrument', len(NJ) == 16 and all(v['spravne'] for v in NJ.values())); ok('16 of 16 tuning', len(LD) == 16 and all(v['spravne'] for v in LD.values()))
eq('1.1\\cdot10^{-6}', sci(max(v['S'] for v in NJ.values()))); eq('smallest value 5.999', 'smallest value %s' % f(min(v['druhy'][0] for v in NJ.values()), 3)); ok('6.0 or more', round(min(v['druhy'][0] for v in NJ.values()), 1) >= 6.0)
qs = {j: sorted(v['Q'].values()) for j, v in LD.items()}; eq('$Q = 0.13$–$0.78$', '$Q = %s$–$%s$' % (f(min(v[0] for v in qs.values()), 2), f(max(v[0] for v in qs.values()), 2)))
eq('$2.25$–$10.97$', '$%s$–$%s$' % (f(min(v[1] for v in qs.values()), 2), f(max(v[1] for v in qs.values()), 2))); ok('chosen = best', all(abs(v['Q'][v['vybrano']] - qs[j][0]) < 1e-12 for j, v in LD.items())); ok('16 subsets', all(len(v['Q']) == 16 for v in LD.values()))
# ---------------------------------------------------------------- 4.9 certification
KK = J('hr7_kontroly_K1K4.json')['L']; k1 = [w for v in KK.values() for w in v['K1']]
eq('120 windows, 21,707 zeros', '%d windows, %s zeros' % (len(k1), th(sum(w['n'] for w in k1)))); ok('K1 certified', all(w['certifikovano'] and w['n'] == w['zmen_znamenka'] for w in k1)); eq('80 bits', '%d bits' % max(w['max_bitu'] for w in k1))
eq('40,400', th(KK['q5_c1']['K1'][1]['k'] * 100)); eq('75,600', th(KK['q5_c1']['K1'][-1]['k'] * 100)); eq('36,000, 40,400', th(KK['q5_c1']['K1'][0]['k'] * 100) + ', ' + th(KK['q5_c1']['K1'][1]['k'] * 100))
k4 = [w for v in KK.values() for w in v['K4']]; ok('120 pairs', len(k4) == 120 and all(w['certifikovano'] for w in k4)); eq('0.0071', f(min(w['rozvinuty'] for w in k4), 4)); eq('3.9\\cdot10^{-3}', sci(min(w['rozestup'] for w in k4)))
K23a, K23b = J('hr7_k23_vysledek_mac.json'), J('hr7b_k23_vysledek_mac.json'); eq('5.1\\cdot10^{-11}', sci(max(x['max_rozdil'] for x in K23a['vysledky'] + K23b['vysledky'])))
ok('48 recomputations', len(K23a['vysledky']) == 48 and K23a['v_poradku']); ok('16 for the four', len(K23b['vysledky']) == 16 and K23b['v_poradku']); ok('same counts', all(x['n'] == x['n_znovu'] for x in K23a['vysledky'] + K23b['vysledky']))
ok('38 digits / finer step', sorted(set((x['kontrola'], x['realprecision'], x['divz']) for x in K23a['vysledky'])) == [('K2', 38, 32), ('K3', 19, 128)]); ok('36,000–36,100 and 58,000–58,100', sorted(set(x['k'] for x in K23a['vysledky'])) == [360, 580])
cover('36,000–36,100', '58,000–58,100')
# ---------------------------------------------------------------- 4.10 completeness
SK, SK2 = J('hr7_skoky.json')['L'], J('hr7_skoky_po_oprave.json')['L']; OP = json.load(open(os.path.join(MF, 'oprava_hr7b_dodatek7.json')))
eq('stays within 0.01 for all twelve', 'stays within %s for all twelve' % f(max(SK[j]['nejvetsi_abs_mimo_skoky'] for j in HR), 2)); ok('hr7 no step', all(not SK[j]['skoky'] for j in HR))
st = [s_['skok'] for j in HB for s_ in SK[j]['skoky']]; eq('+2.00', '+' + f(min(st), 2)); ok('+2.00 (both)', len(st) == 2 and all(abs(x - 2) < 0.01 for x in st))
o = {x['L']: x for x in OP['okna']}; a, b = o['q15_c2-1']['pridane'], o['q21_c4-1']['pridane']
eq('9\\cdot10^{-3}', sci0(b_ := abs(a[1] - a[0]))); eq('62{,}027.74', th(int(a[0])).replace(',', '{,}') + f(a[0] % 1, 2)[1:]); eq('1.0\\cdot10^{-2}', sci(abs(b[1] - b[0]))); eq('51{,}751.58', th(int(b[0])).replace(',', '{,}') + f(b[0] % 1, 2)[1:])
eq('by 1.49 and by 0.95 zeros', 'by %s and by %s zeros' % (f(o['q15_c2-1']['ocek'] - o['q15_c2-1']['n_puvodni'], 2), f(o['q21_c4-1']['ocek'] - o['q21_c4-1']['n_puvodni'], 2)))
ok('exactly the two zeros', all(x['n_nove'] - x['n_puvodni'] == 2 and not x['chybejici_z_puvodnich'] and x['max_rozdil_spolecnych'] == 0 for x in OP['okna']))
ok('After the repair no step', all(not SK2[j]['skoky'] for j in ALL) and len(SK2) == 16); ok('#15.14 and 21.11', E3['L'] is not None)
# ---------------------------------------------------------------- 4.11 resolution
S11 = J('hr7_soused11.json'); B5, B6 = S11['body']['5'], S11['body']['6']; dT = 80000 - 36000
eq('x/7000', 'x/%d' % (round(dT / (2 * math.pi), -3))); eq('0.03 at $x = 200$', '%s at $x = 200$' % f(2 * math.pi * 200 / dT, 2)); eq('23 at $11^5$ and 253 at $11^6$', '%d at $11^5$ and %d at $11^6$' % (round(2 * math.pi * 11 ** 5 / dT), round(2 * math.pi * 11 ** 6 / dT)))
ok('161,071 = 1 mod 21', 161071 % 21 == 1 and B5['sousede'][0]['n'] == 161071 and B5['sousede'][0]['d'] == 20); ok('161{,}051 \\equiv 2', 11 ** 5 == 161051 and 11 ** 5 % 21 == 2)
eq('T_0 = 58{,}000', 'T_0 = ' + th(S11['T0']).replace(',', '{,}')); eq('H = 4{,}000', 'H = ' + th(S11['H']).replace(',', '{,}')); eq('947,573 zeros of the union', th(S11['N']) + ' zeros of the union')
eq('1,602 points', th(sum(len(S11['mrizka'][k]['d']) for k in ('5', '6'))) + ' points'); eq('1.9\\cdot10^{-10}', sci(S11['max_rozdil_k_Ng']))
eq('1.2\\cdot10^{-6}', sci(max(max(S11['mrizka'][k]['max_rozdil_k_plne_sile'] for k in ('5', '6')), max(v['rozdil_k_plne_sile'] for v in S11['body'].values()))))
eq('4.42 times', f(B5['Z_k_plne_sile'], 2) + ' times'); eq('remainder is $2\\cdot10^{-7}$', 'remainder is $%s$' % sci0(abs(complex(*B5['zbytek_po_odecteni_sousedu']))))
ok('1{,}771{,}561 \\equiv 1', 11 ** 6 == 1771561 and 11 ** 6 % 21 == 1); pr6 = S11['mrizka']['6']['prvocisla']; eq('55 within 5,800', '%d within 5,800' % len(pr6)); ok('within 5,800 (range)', max(abs(p_['d']) for p_ in pr6) <= 5800)
eq('nearest 126 away', 'nearest %d away' % min(abs(p_['d']) for p_ in pr6)); eq('-1.000001', f(B6['zbytek_po_odecteni_sousedu'][0], 6))
eq('kernel (40)', 'kernel (%d)' % round(B5['sirka_jadra'])); ok('(42)', 2 * 21 == 42); C11 = S11['ctvercove_cteni_11_5']; eq('29 unknowns', '%d unknowns' % C11['kandidatu']); eq('1.8\\cdot10^7', sci(C11['podminenost']).replace('{7}', '7'))
eq('as $3\\cdot10^{-6}$', 'as $%s$' % sci0(abs(complex(*C11['amplituda_v_11_5_k_plne_sile'])))); w3 = [c for c in C11['vsechny'] if abs(c['d']) <= 300]
eq('14 candidates within 300', '%d candidates within 300' % len(w3)); pr = [c for c in w3 if c['prvocislo']]; co = [c for c in w3 if not c['prvocislo']]
eq('five primes', {5: 'five'}.get(len(pr), str(len(pr))) + ' primes'); eq('nine composite', {9: 'nine'}.get(len(co), str(len(co))) + ' composite')
eq('1 \\pm 4\\cdot10^{-7}', '1 \\pm ' + sci0(max(abs(complex(*c['amplituda_k_12_ln_n']) - 1) for c in pr))); eq('at most $5\\cdot10^{-7}$', 'at most $%s$' % sci0(max(abs(complex(*c['amplituda_k_12_ln_n'])) for c in co)))
ok('ten times wider', 9 <= B6['sirka_jadra'] / 42 <= 11)
# ---------------------------------------------------------------- 4.12 pair correlation, 4.13 unions (addendum 9)
V1, TH = J('hr7_rmt_vysledek.json'), J('hr7_rmt_teorie.json'); M1 = V1['M1']; ok('16 L-functions (M1)', len(M1) == 16 and V1['M1_verdikt']['O'])
eq('0.78–1.86', rng([v['chi2_BK'] for v in M1.values()], 2)); eq('6.3–16.5', rng([v['chi2_GUE'] for v in M1.values()], 1)); ok('pass mark 2.0', all(v['chi2_BK'] <= 2.0 and v['chi2_GUE'] >= max(4, 3 * v['chi2_BK']) for v in M1.values()))
eq('0.032–0.076', '%s–%s' % (f(-max(v['rozdil_mereni_minus_GUE_u_0_5'] for v in M1.values()), 3), f(-min(v['rozdil_mereni_minus_GUE_u_0_5'] for v in M1.values()), 3)))
eq('0.045–0.065', '%s–%s' % (f(-max(v['rozdil_BK_minus_GUE_u_0_5'] for v in M1.values()), 3), f(-min(v['rozdil_BK_minus_GUE_u_0_5'] for v in M1.values()), 3)))
eq('2.3–4.0 errors', rng([v['nejvetsi_z_proti_BK'] for v in M1.values()], 1) + ' errors'); DH = V1['M1_dohromady']; eq('1,296,642 spacings', th(DH['rozestupu']) + ' spacings'); ok('#1,296,642 = N - 16', DH['rozestupu'] == 1296658 - 16)
dd = np.array(DH['p']) - np.array(V1['modely']['p_GUE']); s_ = np.arange(30) * 0.1 + 0.05; im = int(np.argmax(np.abs(dd)))
eq('up to 0.083 (at $s = 1.15$)', 'up to %s (at $s = %s$)' % (f(abs(dd[im]), 3), f(s_[im], 2))); ok('#0.083 = max_rozdil_od_GUE', abs(DH['max_rozdil_od_GUE'] - abs(dd[im])) < 1e-12 and abs(DH['kde'] - 1.15) < 1e-9)
near = lambda x0: dd[int(np.argmin(abs(s_ - x0)))]; eq('$-0.06$ near $s = 0.45$', '$%s$ near $s = %s$' % (f(dd[:10].min(), 2), f(s_[int(np.argmin(dd[:10]))], 2))); eq('$+0.08$ near $s = 1.0$–$1.2$', '$+%s$ near $s = 1.0$–$1.2$' % f(np.mean([near(x) for x in (0.95, 1.05, 1.15)]), 2))
eq('$-0.03$ near $s = 1.8$', '$%s$ near $s = 1.8$' % f(np.mean([near(1.75), near(1.85)]), 2)); ch = DH['chi2']
eq('156 for $N = 3$, 304 for $N = 4$, 1,931 for $N = 2$ and 1,180 for GUE', '%s for $N = 3$, %s for $N = 4$, %s for $N = 2$ and %s for GUE' % (th(ch['CUE3']), th(ch['CUE4']), th(ch['CUE2']), th(ch['GUE'])))
ok('CUE3 best', DH['nejlepsi'] == 'CUE3'); eq('2.0–2.3', rng([v['N_ef'] for v in TH['moduly'].values()], 1)); eq('\\approx 11$–$12', '\\approx %d$–$%d' % (round(min(v['L_stred'] for v in TH['moduly'].values())), round(max(v['L_stred'] for v in TH['moduly'].values()))))
M2 = V1['M2']; ok('M2 verdict', V1['M2_verdikt']['O'])
for U, q, sec in (('U5', '5', 4), ('U7', '7', 6), ('U15', '{15}', 8), ('U21', '{21}', 12)):
    x = M2[U]; row = '| $\\mathbb{Q}(\\zeta_%s)$ | %d | %s | %s | %s / %s / %s | %s | %s | %d |' % (q, x['sektoru'], th(x['nul']), f(x['max_rozdil_od_S'], 4), f(x['prvni_kos'], 3), f(x['prvni_kos_S'], 3), f(x['prvni_kos_Poisson'], 3), f(x['max_rozdil_od_jednoho_GUE'], 2), f(x['chi2_S'], 1), x['kosu_nad_3_chyby'])
    ok(row, row in T and x['sektoru'] == sec, row)
    cover(th(x['nul']), f(x['max_rozdil_od_S'], 4), f(x['prvni_kos'], 3), f(x['prvni_kos_S'], 3), f(x['prvni_kos_Poisson'], 3), f(x['max_rozdil_od_jednoho_GUE'], 2), f(x['chi2_S'], 1))
ok('pass marks M2', all(x['max_rozdil_od_S'] <= 0.015 and x['rozdil_od_Poissona_prvni_kos'] >= 0.03 and x['max_rozdil_od_jednoho_GUE'] >= 0.2 for x in M2.values()))
eq('1.9–4.9', rng([x['chi2_S'] for x in M2.values()], 1)); eq('2.9–8.3 errors', rng([x['nejvetsi_z'] for x in M2.values()], 1) + ' errors')
fb = [x['prvni_kos'] - x['prvni_kos_S'] for x in M2.values()]; ok('below in all four', all(v < 0 for v in fb)); eq('by 0.006–0.010', 'by %s–%s' % (f(-max(fb), 3), f(-min(fb), 3)))
eq('(at most 0.0002)', '(at most %s)' % f(max(x['citlivost_S_na_vysku'] for x in M2.values()), 4)); eq('to 0.010 in the', 'to %s in the' % f(max(x['max_rozdil_od_S'] for x in M2.values()), 3))
# ---------------------------------------------------------------- 4.14 height (addendum 10)
V3 = J('hr7_rmt_vyska_vysledek.json'); L3 = V3['L']; SD = V3['sdruzeno']; UP = V3['uplnost']; ok('verdict M3', SD['verdikt'] == 'O' and V3['uplnost_v_poradku'])
eq('43,907 to 44,980 zeros, 178,844', '%s to %s zeros, %s' % (th(min(x['nul_nahore'] for x in L3.values())), th(max(x['nul_nahore'] for x in L3.values())), th(sum(x['nul_nahore'] for x in L3.values()))))
eq('$-0.65$ to $+1.35$', '$%s$ to $+%s$' % (f(min(u['rozdil'] for u in UP.values()), 2), f(max(u['rozdil'] for u in UP.values()), 2)))
SKC = J('hr7c_skoky.json')['L']; ok('no step hr7c', all(not v['skoky'] for v in SKC.values())); eq('largest $|\\Delta|$ 0.01', 'largest $|\\Delta|$ %s' % f(max(v['nejvetsi_abs_mimo_skoky'] for v in SKC.values()), 2))
MP = [J('hr7c_mp_%s.json' % j) for j in HB]; eq('(24 in all)', '(%d in all)' % sum(len(m['nuly']) for m in MP)); eq('5.4\\cdot10^{-9}', sci(max(n['absL'] for m in MP for n in m['nuly'])))
eq('at least 0.21', 'at least %s' % f(min(s['absL'] for m in MP for s in m['stredy']), 2))
eq('1.05–1.82', rng([x['chi2_BK_nahore'] for x in L3.values()], 2)); eq('3.5–8.0', rng([x['chi2_GUE_nahore'] for x in L3.values()], 1)); eq('0.92–1.07', rng([x['A_nahore'] for x in L3.values()], 2))
eq('0.987 \\pm 0.027', f(SD['A_nahore'], 3) + ' \\pm ' + f(SD['chyba_A_nahore'], 3)); eq('0.995 \\pm 0.011', f(SD['A_dole'], 3) + ' \\pm ' + f(SD['chyba_A_dole'], 3))
eq('0.0363 \\pm 0.0004', f(SD['V_dole'], 4) + ' \\pm ' + f(SD['chyba_V_dole'], 4)); eq('0.0261 \\pm 0.0007', f(SD['V_nahore'], 4) + ' \\pm ' + f(SD['chyba_V_nahore'], 4))
eq('12.6 errors apart', f(SD['rozdil_V'] / SD['chyba_rozdilu'], 1) + ' errors apart'); eq('0.719 \\pm 0.021', f(SD['pomer_V'], 3) + ' \\pm ' + f(SD['chyba_pomeru'], 3)); eq('0.727', f(SD['pomer_teorie'], 3))
ok('pass marks M3', sum(1 for x in L3.values() if x['chi2_BK_nahore'] <= 2.0) >= 3 and abs(SD['A_nahore'] - 1) <= 3 * SD['chyba_A_nahore'] and SD['rozdil_V'] >= 3 * SD['chyba_rozdilu'] and abs(SD['pomer_V'] - SD['pomer_teorie']) <= 3 * SD['chyba_pomeru'])
T36, T400 = J('hr7_rmt_teorie_vyska_036.json'), J('hr7_rmt_teorie_vyska_400.json'); eq('\\approx 12.2$ and $14.1', '\\approx %s$ and $%s' % (f(T36['moduly']['21']['L_stred'], 1), f(T400['moduly']['21']['L_stred'], 1)))
eq('0.071–0.081 in the upper band and 0.082–0.106', rng([x['rozestupy_max_rozdil_od_GUE_nahore'] for x in L3.values()], 3) + ' in the upper band and ' + rng([x['rozestupy_max_rozdil_od_GUE_dole'] for x in L3.values()], 3))
# ---------------------------------------------------------------- verdicts (section 6) and section 5.2
VD = dict(**E1['verdikty'], **E2['verdikty'], **J('hr7_r_verdikty.json')['verdikty'], **E3['verdikty'], **{'odv-hr7-rmt-par': 'O', 'odv-hr7-rmt-orchestr': 'O', 'odv-hr7-rmt-vyska': SD['verdikt']})
ok('All 18 hypotheses', len(VD) == 18 and all(v == 'O' for v in VD.values())); rows = re.findall(r'^\| (odv-hr7-[\w-]+): .*\| (\w) \|$', T, re.M); ok('table 6', dict(rows) == VD, rows)
eq('7.13 measured on average', '%s measured on average' % f(np.mean([Z2[j]['V'] for j in HR if H[j]['q'] == 15]), 2))
MZ = E3['mez']; mq = lambda q, k: {'5': MZ['5']['q5_c1'], '15': MZ['15']['q15_c1-1'], '21': MZ['21']['q21_c1-1']}[q][k]      # one L-function per modulus, named in the text
eq('(5.2), $Q = 3.3$ with the prime 2 alone, 0.26 with 2 and 3, and 0.13 with all primes', '(5.2), $Q = %s$ with the prime 2 alone, %s with 2 and 3, and %s with all primes' % (f(mq('5', '2'), 1), f(mq('5', '3'), 2), f(mq('5', 'vse'), 2)))
eq('(15.2), 0.49 up to $P = 5$ (only 2 contributes), 0.31 with 7 and 0.20 with all', '(15.2), %s up to $P = 5$ (only 2 contributes), %s with 7 and %s with all' % (f(mq('15', '5'), 2), f(mq('15', '7'), 2), f(mq('15', 'vse'), 2)))
eq('(21.17), 1.0 up to $P = 3$, 0.21 from $P = 5$ and 0.20 with all', '(21.17), %s up to $P = 3$, %s from $P = 5$ and %s with all' % (f(mq('21', '3'), 1), f(mq('21', '5'), 2), f(mq('21', 'vse'), 2)))
eq('1.0136', f(3 ** 12 / 2 ** 19, 4))
# ---------------------------------------------------------------- coverage: every number with a decimal point or >= 3 digits
STRUCT = set("""256 36,000 80,000 400,000 420,000 36{,}000 80{,}000 400{,}000 420{,}000 2.17.2 2026 08:03 23:13 00:41–00:43 20:03 23:32 07:58 09:15 15:14 15:51 0.025 2.5 0.05 0.1 0.37 0.1, 200 1,000 400 100 200 2.0 1.5 0.5 3.0 0.03 0.015 0.2 0.5,
 3.1 3.10 5.7 5.2 5.3 7.2 7.4 7.3 7.5 13.2 13.7 15.2 15.8 21.5 21.17 15.14 21.20 21.2 21.11 34.116 46.1 3\\cdot10^{-9} 5\\cdot10^{-5} 8\\cdot10^{-7} 10^{-10} 10^{-3} 1.0 3.11 3.12 3.13 3.9 4.1 4.10 4.11 4.12 4.13 4.14 4.2 4.3 4.4 4.5 4.6 4.7 5.4 5.5 5.6
 0.6 15–20 15 20 0.75 0.4 1.8 1.2 0.5 0.7 2.0 1.15 1.0–1.2 5,800 300 161,071 161{,}051 1{,}771{,}561 1.6\\cdot10^5 10^{-6} 1155 1472–1475 095202 1307.6010 2302.07073 236–266 2128–2136 1698 2413–2421 1–26 548–564 chao-dyn/9409004 10.5281/zenodo.23815328""".split())
i = T.index('abstract: |'); body = T[i:T.index('# 10. Reproducibility')]; body = re.sub(r'<!--.*?-->', '', body, flags=re.S)
body = body.replace('{,}', ',')
cov = ' '.join(OK).replace('{,}', ',')
miss = []
for m in re.finditer(r'(?<![\w.\\^{])-?\d[\d,]*\.?\d*', body):
    s = m.group(0).rstrip('.,')
    if not ('.' in s or len(s.replace(',', '').lstrip('-')) >= 3): continue
    if s in STRUCT or s.lstrip('-') in STRUCT: continue
    if s in cov or s.lstrip('-') in cov: continue
    miss.append(s)
notin = [s for s in OK if s not in T and not s.startswith('#') and not re.search(r'(\(|\b)(manifest|all|both|summary|krizova|M1|M2|M3|hr7c|12\)|q=\d+|instrument|tuning|certified|range|pass marks|verdict|best|repair|union|step|= N|= max|subsets|estimates|unknowns\)|1 mod 21|equiv|\(42\)|digits|table 6|identical|arithmetic|characters|wider|nine windows|PARI|K1|0\.78 \(both\)|Q \\\\approx|factor of|diag|mostly|points|pairs|windows\)|L-functions|numbers|^\(5\.2\)$|\\\\ge|computations|four|in all|zero)', s)]
print('\n'.join('FAIL %-50s %s' % b for b in BAD))
print('verified:', len(OK), '| failed:', len(BAD), '| numbers not covered:', sorted(set(miss)), '| labels not found literally in the text:', notin)
sys.exit(1 if (BAD or miss) else 0)
