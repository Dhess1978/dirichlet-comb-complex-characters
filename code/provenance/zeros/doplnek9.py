"""Doplnění 9 oken sady hr7, která přerušil Ctrl+C na Macu 20. 11. v 10:27 UTC (úlohy rozpracované v tu chvíli): týž zmrazený skript fz_nuly_chi_v3.py
(funkce pocitej beze změny, PARI 2.17.2, realprecision 19, divz 32), jen zavolaný na zadaná okna. Použití: python3 doplnek9.py VYSTUP_SLOZKA JMENO:K [JMENO:K ...]"""
import sys, os, json, importlib.util
spec = importlib.util.spec_from_file_location('fz', '/home/claude/br107/mac/fz_nuly_chi_v3.py'); fz = importlib.util.module_from_spec(spec); spec.loader.exec_module(fz)
V = sys.argv[1]; fz._init(32, 2.0); L = {fz.jmeno(q, c): (q, c) for q, c in fz.SADY['hr7']}; out = []
for arg in sys.argv[2:]:
    jm, k = arg.split(':'); k = int(k); q, chi = L[jm]; P = fz.popis(q, chi); os.makedirs(os.path.join(V, jm), exist_ok=True)
    r = fz.pocitej(q, chi, P['lichy'], k, k * 100.0, (k + 1) * 100.0, os.path.join(V, jm, 'okno_%06d.json.gz' % k)); out.append(r); print(r, flush=True)
json.dump(out, open(os.path.join(V, 'doplnek_%s.json' % os.getpid()), 'w'), indent=1)
