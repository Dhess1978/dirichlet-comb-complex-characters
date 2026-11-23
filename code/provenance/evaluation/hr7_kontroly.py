"""Studie VII · kontroly kvality dat (PRAVIDLA_STUDIE_VII.md, dodatek 4): K1 certifikace vzorku oken intervalovou aritmetikou (Arb přes python-flint),
K2 přepočet vybraných oken v PARI s realprecision 38, K3 přepočet s divz 128, K4 certifikace 10 nejtěsnějších dvojic každé L-funkce.
Použití: python3 hr7_kontroly.py SLOZKA VYSTUP.json [--k1 360,404,...] [--k23 360,580] [--tmin T1]
Výchozí okna: K1 = 360, 404, …, 756 (pod T₂), K2 a K3 = 360 a 580."""
import gzip, hashlib, json, math, os, sys, time, numpy as np
from flint import acb, arb, ctx, dirichlet_char
from scipy.special import loggamma
import cypari2
t0 = time.time(); V, OUT = sys.argv[1:3]
def arg(nm, dflt):
    return [int(x) for x in sys.argv[sys.argv.index(nm) + 1].split(',')] if nm in sys.argv else dflt
K1 = arg('--k1', [360, 404, 448, 492, 536, 580, 624, 668, 712, 756]); K23 = arg('--k23', [360, 580])
T1 = float(sys.argv[sys.argv.index('--tmin') + 1]) if '--tmin' in sys.argv else 36000.0
P = cypari2.Pari(); P.allocatemem(2 * 10 ** 9, silent=True)
M = json.load(open(os.path.join(V, 'manifest.json'))); OK = float(M['okno'])
def blok(jm, k):
    b = M['L'][jm]['bloky'].get(str(k))
    if not b: return None
    raw = open(os.path.join(V, jm, b['soubor']), 'rb').read(); assert hashlib.sha256(raw).hexdigest() == b['sha256']; return json.loads(gzip.decompress(raw))
def conrey(q, chi):
    G = 'znstar(%d,1)' % q; c = '[' + ','.join(str(int(x)) for x in chi) + ']'
    m = [m for m in range(1, q) if math.gcd(m, q) == 1 and str(P('znconreychar(%s,%d)' % (G, m))).replace(' ', '') == c]
    assert len(m) == 1; return m[0]
def znak(ch, t):
    for prec in (80, 160, 320):
        ctx.prec = prec; v = ch.hardy_z(acb(t)).real
        if v > 0: return 1, prec
        if v < 0: return -1, prec
    return 0, prec
def certifikuj(ch, body):
    zn = [znak(ch, t) for t in body]; s = [z for z, _ in zn]
    ok = all(x != 0 for x in s) and all(s[i] * s[i + 1] < 0 for i in range(len(s) - 1))
    return ok, sum(1 for i in range(len(s) - 1) if s[i] * s[i + 1] < 0), max(p for _, p in zn)
out = dict(slozka=os.path.basename(os.path.normpath(V)), K1=K1, K23=K23, L={}, skript_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())
for jm, z in M['L'].items():
    q = int(z['popis']['q']); chi = z['popis']['chi']; m = conrey(q, chi); ch = dirichlet_char(q, m); r = dict(q=q, chi=chi, conrey='%d.%d' % (q, m), K1=[], K2=[], K3=[])
    # K1
    for k in K1:
        b = blok(jm, k)
        if b is None: continue
        g = b['g']; body = [b['a']] + [(g[i] + g[i + 1]) / 2 for i in range(len(g) - 1)] + [b['b']]; t1 = time.time(); ok, ns, pm = certifikuj(ch, body)
        r['K1'].append(dict(k=k, n=len(g), zmen_znamenka=ns, certifikovano=ok, max_bitu=pm, sekund=round(time.time() - t1, 1)))
        print('K1 %-10s (%s) okno %d: %d nul, certifikovaných změn znaménka %d → %s (nejvýš %d bitů, %.0f s)' % (jm, r['conrey'], k, len(g), ns, 'OK' if ok else 'POZOR', pm, time.time() - t1), flush=True)
    # K2, K3
    for k in K23:
        b = blok(jm, k)
        if b is None: continue
        g = np.array(b['g'])
        for nm, rp, dz in (('K2', 38, 32), ('K3', 19, 128)):
            P.default('realprecision', rp); t1 = time.time()
            z2 = sorted(float(v) for v in P('lfunzeros(lfuncreate([znstar(%d,1),[%s]]), [%s, %s], %d)' % (q, ','.join(map(str, chi)), repr(float(b['a'])), repr(float(b['b'])), dz)))
            z2 = np.array([v for v in z2 if b['a'] < v <= b['b']]); d = float(np.max(np.abs(z2 - g))) if len(z2) == len(g) else None
            r[nm].append(dict(k=k, n=len(g), n_znovu=len(z2), max_rozdil=d, v_poradku=(d is not None and d <= 1e-9), sekund=round(time.time() - t1, 1)))
            print('%s %-10s okno %d: %d / %d nul, největší |Δγ| %s → %s (%.0f s)' % (nm, jm, k, len(g), len(z2), ('%.1e' % d) if d is not None else '—', 'OK' if (d is not None and d <= 1e-9) else 'POZOR', time.time() - t1), flush=True)
        P.default('realprecision', 19)
    # K4: 10 nejtěsnějších dvojic v pásmu (rozvinuto θ_χ/π)
    g = []; k = int(round(T1 / OK))
    while True:
        b = blok(jm, k)
        if b is None: break
        g.extend(b['g']); k += 1
    g = np.array(g); odd = int(z['popis']['lichy']); x = (0.5 * g * math.log(q / math.pi) + np.imag(loggamma((0.5 + odd) / 2 + 0.5j * g))) / math.pi; dx = np.diff(x); ix = np.argsort(dx)[:10]; K4 = []
    for i in sorted(ix):
        lo = g[i - 1] if i > 0 else g[i] - 0.5; hi = g[i + 2] if i + 2 < len(g) else g[i + 1] + 0.5
        body = [(lo + g[i]) / 2, (g[i] + g[i + 1]) / 2, (g[i + 1] + hi) / 2]; ok, ns, pm = certifikuj(ch, body)
        K4.append(dict(gamma=[float(g[i]), float(g[i + 1])], rozestup=float(g[i + 1] - g[i]), rozvinuty=float(dx[i]), certifikovano=ok, max_bitu=pm))
    r['K4'] = K4; r['K4_nejmensi_rozvinuty'] = float(dx.min()); r['nul_v_pasmu'] = int(len(g))
    print('K4 %-10s: %d nul, nejmenší rozvinutý rozestup %.4f, certifikováno %d / %d dvojic | %.0f s' % (jm, len(g), dx.min(), sum(d['certifikovano'] for d in K4), len(K4), time.time() - t0), flush=True)
    out['L'][jm] = r
    json.dump(out, open(OUT, 'w'))
vse = all(all(w['certifikovano'] for w in r['K1']) and all(w['v_poradku'] for w in r['K2'] + r['K3']) and all(d['certifikovano'] for d in r['K4']) for r in out['L'].values())
out['kontroly_v_poradku'] = vse; json.dump(out, open(OUT, 'w')); print('KONTROLY K1–K4:', 'v pořádku' if vse else 'POZOR — dohledat', '| %.0f s' % (time.time() - t0))
