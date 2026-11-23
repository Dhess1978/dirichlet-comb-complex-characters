#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MAYAN ALFA V1008 · studie VII · doplnění 9 oken sady hr7 na Macu (okna rozpracovaná v okamžiku Ctrl+C 20. 11. v 10:27 UTC).
Volá funkci pocitej ze zmrazeného skriptu fz_nuly_chi_v3.py (ze stejné složky, SHA-256 musí začínat cf90cceb) beze změny: PARI, realprecision 19, divz 32,
hlídání úplnosti. Zapisuje jen chybějící soubory okno_XXXXXX.json.gz do výstupní složky (existující soubory nepřepisuje, manifest nemění).
Použití: python3 doplnek9_mac.py SLOZKA_HR7"""
import hashlib, importlib.util, json, os, sys, time
TU = os.path.dirname(os.path.abspath(__file__)); FZ = os.path.join(TU, 'fz_nuly_chi_v3.py')
sha = hashlib.sha256(open(FZ, 'rb').read()).hexdigest(); assert sha.startswith('cf90cceb'), 'nesouhlasí otisk fz_nuly_chi_v3.py: ' + sha
spec = importlib.util.spec_from_file_location('fz', FZ); fz = importlib.util.module_from_spec(spec); spec.loader.exec_module(fz)
OKNA = [('q13_c1', 432), ('q13_c11', 432), ('q15_c1-1', 432), ('q21_c1-1', 431), ('q21_c5-1', 431), ('q7_c1', 432), ('q7_c2', 432), ('q7_c4', 432), ('q7_c5', 432)]
def main():
    V = sys.argv[1]; fz._init(32, 2.0); L = {fz.jmeno(q, c): (q, c) for q, c in fz.SADY['hr7']}; out = []; t0 = time.time()
    for jm, k in OKNA:
        cesta = os.path.join(V, jm, 'okno_%06d.json.gz' % k)
        if os.path.exists(cesta): print('  %s okno %d už existuje — přeskočeno' % (jm, k), flush=True); continue
        q, chi = L[jm]; P = fz.popis(q, chi); r = fz.pocitej(q, chi, P['lichy'], k, k * 100.0, (k + 1) * 100.0, cesta); out.append(r)
        print('  %s okno %d: %d nul (hladké očekávání %.2f), krok divz %d, %.0f s' % (jm, k, r['n'], r['ocek'], r['divz'], r['sekund']), flush=True)
    json.dump(dict(okna=out, prostredi=fz.prostredi(), skript_sha256=sha), open(os.path.join(TU, 'doplnek9_mac_vysledek.json'), 'w'), indent=1)
    print('Hotovo: doplněno %d oken za %.1f min.' % (len(out), (time.time() - t0) / 60))
if __name__ == '__main__':
    main()
