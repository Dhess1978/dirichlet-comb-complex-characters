#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MAYAN ALFA V1008 · studie VII · kontroly K2 a K3 (PRAVIDLA_STUDIE_VII.md, dodatek 4) na Macu — jen PARI, bez Arb.
K2: vybraná okna znovu s realprecision 38; K3: táž okna s divz 128. Porovná počty a polohy nul s uloženými okny (|Δγ| ≤ 1e-9).
Logika je táž jako v hr7_kontroly.py (SHA-256 95ebce69…), jen bez části K1/K4, která potřebuje knihovnu Arb.
Použití: python3 hr7_k23.py SLOZKA VYSTUP.json [--okna 360,580] [--jadra 9]
Nic jiného než VYSTUP.json nezapisuje, nic neposílá po síti."""
import gzip, hashlib, json, os, sys, time
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
def uloha(args):
    V, jm, b_info, q, chi, k, rp, dz = args
    import cypari2, numpy as np
    P = cypari2.Pari(); P.allocatemem(2 * 2 ** 30, silent=True); P.default('realprecision', rp)
    raw = open(os.path.join(V, jm, b_info['soubor']), 'rb').read(); b = json.loads(gzip.decompress(raw)); g = np.array(b['g']); t = time.time()
    z = sorted(float(v) for v in P('lfunzeros(lfuncreate([znstar(%d,1),[%s]]), [%s, %s], %d)' % (q, ','.join(map(str, chi)), repr(float(b['a'])), repr(float(b['b'])), dz)))
    z = np.array([v for v in z if b['a'] < v <= b['b']]); d = float(np.max(np.abs(z - g))) if len(z) == len(g) else None
    return dict(L=jm, k=k, kontrola='K2' if rp == 38 else 'K3', realprecision=rp, divz=dz, n=len(g), n_znovu=len(z), max_rozdil=d, v_poradku=(d is not None and d <= 1e-9),
                sha_bloku=hashlib.sha256(raw).hexdigest(), sekund=round(time.time() - t, 1))
def main():
    V, OUT = sys.argv[1:3]; ok = [int(x) for x in sys.argv[sys.argv.index('--okna') + 1].split(',')] if '--okna' in sys.argv else [360, 580]
    jad = int(sys.argv[sys.argv.index('--jadra') + 1]) if '--jadra' in sys.argv else 9
    M = json.load(open(os.path.join(V, 'manifest.json'))); U = []
    for jm, z in M['L'].items():
        q = int(z['popis']['q']); chi = z['popis']['chi']
        for k in ok:
            f = 'okno_%06d.json.gz' % k
            if not os.path.exists(os.path.join(V, jm, f)): continue
            for rp, dz in ((38, 32), (19, 128)): U.append((V, jm, dict(soubor=f), q, chi, k, rp, dz))
    print('Úloh %d (12 L-funkcí × %d okna × K2, K3), procesů %d' % (len(U), len(ok), jad), flush=True); t0 = time.time(); R = []
    with ProcessPoolExecutor(max_workers=jad, mp_context=mp.get_context('spawn')) as ex:
        for r in ex.map(uloha, U):
            R.append(r); print('  %s %-10s okno %d: %d / %d nul, největší |Δγ| %s → %s (%.0f s)' % (r['kontrola'], r['L'], r['k'], r['n'], r['n_znovu'], ('%.1e' % r['max_rozdil']) if r['max_rozdil'] is not None else '—', 'OK' if r['v_poradku'] else 'POZOR', r['sekund']), flush=True)
    vse = all(r['v_poradku'] for r in R); json.dump(dict(vysledky=R, v_poradku=vse, skript_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()), open(OUT, 'w'), indent=1)
    print('K2 a K3:', 'v pořádku' if vse else 'POZOR', '| %.2f h' % ((time.time() - t0) / 3600))
if __name__ == '__main__':
    main()
