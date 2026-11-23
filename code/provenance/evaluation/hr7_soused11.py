"""Studie VII · dodatek 8: sousední prvočísla a mez rozlišení u 11^5 a 11^6 ve sjednocení U21 (popisně, bez verdiktu). Jen čte data.
Z(x) = Σ_{γ ∈ U21} exp(−(γ − T0)²/2H²)·x^{iγ} proti P(x) = −(H/√(2π))·Σ_n a(n)Λ(n)n^{−1/2}·exp(−H² ln²(n/x)/2)·exp(−iT0 ln(n/x)).
Použití: python3 hr7_soused11.py SLOZKA_HR7 SLOZKA_HR7B SLOZKA_REF_VI VYSTUP.json"""
import gzip, hashlib, json, math, os, sys, time, numpy as np
t0 = time.time(); VA, VB, VREF, OUT = sys.argv[1:5]; T1, T2, T0, H, Q = 36000.0, 80000.0, 58000.0, 4000.0, 21
ZF = ['/home/claude/br69/data/zeros_14_t.npy', '/home/claude/br69/data/zeros_5000_t.npy', '/home/claude/br69/data/zeros_26000_t.npy']
CLEN = dict(hr7=['q7_c2', 'q7_c4', 'q7_c1', 'q7_c5', 'q21_c1-1', 'q21_c5-1'], hr7b=['q21_c3-1', 'q21_c2-1', 'q21_c4-1'])
G = {}
for V, js in ((VA, CLEN['hr7']), (VB, CLEN['hr7b'])):
    M = json.load(open(os.path.join(V, 'manifest.json')))
    for jm in js:
        g = []
        for k in range(360, 800):
            b = M['L'][jm]['bloky'][str(k)]; raw = open(os.path.join(V, jm, b['soubor']), 'rb').read(); assert hashlib.sha256(raw).hexdigest() == b['sha256']; g.extend(json.loads(gzip.decompress(raw))['g'])
        G[jm] = np.array(g)
MR = json.load(open(os.path.join(VREF, 'manifest.json'))); OKr = float(MR['okno'])
for D in (-3, -7):
    B = MR['L'][str(D)]['bloky']; g = []
    for k in range(int(round(T1 / OKr)), int(round(T2 / OKr))):
        raw = open(os.path.join(VREF, 'D_m%d' % abs(D), B[str(k)]['soubor']), 'rb').read(); assert hashlib.sha256(raw).hexdigest() == B[str(k)]['sha256']; g.extend(json.loads(gzip.decompress(raw))['g'])
    G['VI_D%d' % D] = np.array([v for v in g if T1 < v <= T2])
z = np.concatenate([np.load(f) for f in ZF]); G['zeta'] = z[(z > T1) & (z <= T2)]; del z
assert len(G) == 12; U = np.sort(np.concatenate(list(G.values()))); NU = len(U); wg = np.exp(-0.5 * ((U - T0) / H) ** 2); Ng = float(wg.sum())
print('U21: 12 sektorů, %d nul; N_g = %.1f | %.0f s' % (NU, Ng, time.time() - t0), flush=True)
# prvočísla a mocniny do 1,8e6: Λ(n) a koeficient a(n) sjednocení
NMAX = 1_790_000; sito = np.ones(NMAX + 1, bool); sito[:2] = False
for i in range(2, int(NMAX ** 0.5) + 1):
    if sito[i]: sito[i * i::i] = False
PRV = np.nonzero(sito)[0]; LAM = {}
for p in PRV[PRV <= int(NMAX ** 0.5) + 1]:
    m = int(p) * int(p)
    while m <= NMAX: LAM[m] = math.log(int(p)); m *= int(p)
def a_lam(n):
    """a(n)·Λ(n) pro sjednocení U21 (0, když n není mocnina prvočísla nebo se sektory vyruší)"""
    L = math.log(n) if sito[n] else LAM.get(n, 0.0)
    if not L: return 0.0
    if n % 3 == 0: return 6.0 * L if pow(n, 1, 7) == 1 else 0.0
    if n % 7 == 0: return 2.0 * L
    return 12.0 * L if n % Q == 1 else 0.0
def Zm(xs): return np.array([complex(np.sum(wg * np.exp(1j * U * math.log(x)))) for x in xs])
def cleny(x, vyn=()):
    """seznam (n, příspěvek k P(x)) pro n do 8 šířek jádra; vyn = vynechaná n"""
    W = int(8 * x / H) + 2; out = []
    for n in range(max(2, x - W), x + W + 1):
        al = a_lam(n)
        if al and n not in vyn:
            u = math.log1p((n - x) / x); out.append((n, -(H / math.sqrt(2 * math.pi)) * al / math.sqrt(n) * math.exp(-0.5 * H * H * u * u) * complex(math.cos(T0 * u), -math.sin(T0 * u))))
    return out
def P(x, vyn=()): return sum(c for _, c in cleny(x, vyn))
A = lambda x: (H / math.sqrt(2 * math.pi)) * 12 * math.log(11) / math.sqrt(x)
out = dict(T0=T0, H=H, N=NU, Ng=Ng, skript_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest(), body={}, mrizka={})
dmax = 0.0
for k in range(1, 7):
    x = 11 ** k; Z = Zm([x])[0]; Pp = P(x); Pb = P(x, vyn=(x,)); R = Z - Pb; dmax = max(dmax, abs(Z - Pp) / Ng)
    sous = sorted(cleny(x, vyn=(x,)), key=lambda t: -abs(t[1]))[:6]
    out['body'][str(k)] = dict(x=x, x_mod_21=x % Q, sirka_jadra=x / H, Z=[Z.real, Z.imag], P=[Pp.real, Pp.imag], Z_k_plne_sile=abs(Z) / A(x), P_k_plne_sile=abs(Pp) / A(x), rozdil_k_Ng=abs(Z - Pp) / Ng, rozdil_k_plne_sile=abs(Z - Pp) / A(x),
                               zbytek_po_odecteni_sousedu=[R.real / A(x), R.imag / A(x)], sousede=[dict(n=n, d=n - x, podil_plne_sily=abs(c) / A(x)) for n, c in sous])
    print('11^%d = %-8d (mod 21 = %2d, šířka jádra %.3g): |Z|/A %.4f, |P|/A %.4f, |Z − P|/A %.1e | zbytek po odečtení sousedů %+.4f%+.4fi (čeká se %s) | sousedé: %s | %.0f s'
          % (k, x, x % Q, x / H, abs(Z) / A(x), abs(Pp) / A(x), abs(Z - Pp) / A(x), R.real / A(x), R.imag / A(x), '−1' if x % Q == 1 else '0', ', '.join('%d (d %+d): %.2f' % (n, n - x, abs(c) / A(x)) for n, c in sous[:4]), time.time() - t0), flush=True)
for k, ds in ((5, range(-400, 401)), (6, range(-4000, 4001, 10))):
    x0 = 11 ** k; xs = [x0 + d for d in ds]; Z = Zm(xs); Pp = np.array([P(x) for x in xs]); Pb = np.array([P(x, vyn=(x0,)) for x in xs]); dmax = max(dmax, float(np.max(np.abs(Z - Pp))) / Ng)
    pr = [n for n in range(xs[0] - int(4 * x0 / H), xs[-1] + int(4 * x0 / H) + 1) if a_lam(n) and n != x0]
    out['mrizka'][str(k)] = dict(x0=x0, d=list(ds), Z=[[v.real, v.imag] for v in Z], P=[[v.real, v.imag] for v in Pp], P_bez_x0=[[v.real, v.imag] for v in Pb], A=A(x0), prvocisla=[dict(n=n, d=n - x0) for n in pr], max_rozdil_k_plne_sile=float(np.max(np.abs(Z - Pp))) / A(x0))
    print('mřížka kolem 11^%d: %d bodů, největší |Z − P|/A %.1e, největší |Z|/A %.2f; prvočísel ≡ 1 (mod 21) v dosahu %d | %.0f s' % (k, len(xs), np.max(np.abs(Z - Pp)) / A(x0), np.max(np.abs(Z)) / A(x0), len(pr), time.time() - t0), flush=True)
    if k == 5:
        # (4) čtení „z nul a znalosti q“: neznámé amplitudy na lichých n ≡ 1 (mod 21) + zkušební člen v 11^5
        kand = [n for n in range(x0 - 400 - int(5 * x0 / H), x0 + 400 + int(5 * x0 / H) + 1) if n % Q == 1 and n % 2 == 1] + [x0]
        K = np.array([[-(H / math.sqrt(2 * math.pi)) / math.sqrt(n) * math.exp(-0.5 * (H * math.log1p((n - x) / x)) ** 2) * complex(math.cos(T0 * math.log1p((n - x) / x)), -math.sin(T0 * math.log1p((n - x) / x))) for n in kand] for x in xs])
        sv = np.linalg.svd(K, compute_uv=False); amp, res, rk, _ = np.linalg.lstsq(K, Z, rcond=1e-10)
        jed = 12 * math.log(11)     # amplituda zubu jedenáctky plné síly
        top = sorted(((abs(v), n, v) for n, v in zip(kand[:-1], amp[:-1])), reverse=True)[:6]
        out['ctvercove_cteni_11_5'] = dict(kandidatu=len(kand), podminenost=float(sv[0] / sv[-1]), amplituda_v_11_5_k_plne_sile=[float(amp[-1].real / jed), float(amp[-1].imag / jed)],
                                           nejvetsi=[dict(n=int(n), d=int(n - x0), amplituda_k_12_ln_n=[float(v.real / (12 * math.log(n))), float(v.imag / (12 * math.log(n)))], prvocislo=bool(sito[n])) for _, n, v in top],
                                           zbytek_rms_k_plne_sile=float(np.sqrt(np.mean(np.abs(K @ amp - Z) ** 2)) / A(x0)),
                                           vsechny=[dict(n=int(n), d=int(n - x0), amplituda_k_12_ln_n=[float(v.real / (12 * math.log(n))), float(v.imag / (12 * math.log(n)))], prvocislo=bool(sito[n])) for n, v in zip(kand[:-1], amp[:-1])])
        vn = [r for r in out['ctvercove_cteni_11_5']['vsechny'] if abs(r['d']) <= 300]
        print('  uvnitř |d| ≤ 300: prvočísel %d, největší |amplituda − 1| %.1e; složených %d, největší |amplituda| %.1e' % (sum(r['prvocislo'] for r in vn), max(abs(complex(*r['amplituda_k_12_ln_n']) - 1) for r in vn if r['prvocislo']),
              sum(not r['prvocislo'] for r in vn), max(abs(complex(*r['amplituda_k_12_ln_n'])) for r in vn if not r['prvocislo'])), flush=True)
        print('čtení z nul a znalosti q u 11^5: %d kandidátů, podmíněnost %.1e, amplituda v 11^5 = %+.4f%+.4fi plné síly (čeká se 0); největší amplitudy: %s' % (len(kand), sv[0] / sv[-1], amp[-1].real / jed, amp[-1].imag / jed,
              ', '.join('%d (d %+d, %s): %.3f' % (n, n - x0, 'prvočíslo' if sito[n] else 'složené', v.real / (12 * math.log(n))) for _, n, v in top)), flush=True)
out['max_rozdil_k_Ng'] = dmax; json.dump(out, open(OUT, 'w')); print('největší |Z − P|/N_g = %.1e | %.0f s' % (dmax, time.time() - t0))
