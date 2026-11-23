"""Studie VII · kontrola úplnosti nul (Turingova typu): rozvinuté polohy x_n = θ_χ(γ_n)/π; chybějící nula = trvalý skok x_n − n o +1.
Vypíše největší rozvinuté mezery a polohu skoků (klouzavý medián x_n − n přes 2 000 nul). Jen čte data.
Použití: python3 hr7_uplnost.py SLOZKA [SLOZKA ...] VYSTUP.json"""
import json, os, gzip, math, sys, numpy as np
from scipy.special import loggamma
from scipy.ndimage import median_filter
out = {}
for V in sys.argv[1:-1]:
    M = json.load(open(os.path.join(V, 'manifest.json')))
    for jm, z in M['L'].items():
        g = []; okna = []
        for k in sorted(int(x) for x in z['bloky']):
            b = json.loads(gzip.decompress(open(os.path.join(V, jm, 'okno_%06d.json.gz' % k), 'rb').read())); g.extend(b['g']); okna.extend([k] * len(b['g']))
        g = np.array(g); q = int(z['popis']['q']); odd = int(z['popis']['lichy'])
        x = (0.5 * g * math.log(q / math.pi) + np.imag(loggamma((0.5 + odd) / 2 + 0.5j * g))) / math.pi
        r = x - np.arange(len(g)); r -= np.median(r); m = median_filter(r, size=2001, mode='nearest')
        dx = np.diff(x); ix = np.argsort(dx)[::-1][:5]
        skoky = [int(i) for i in np.where(np.abs(np.diff(np.round(m[::500]))) >= 1)[0] * 500]
        out[jm] = dict(n=int(len(g)), median_min=float(m.min()), median_max=float(m.max()), nejvetsi_mezery=[dict(gamma=[float(g[i]), float(g[i + 1])], rozvinuta=float(dx[i]), okno=int(okna[i])) for i in ix])
        print('%-10s n=%d | klouzavý medián x−n: %+.2f … %+.2f | největší rozvinuté mezery: %s' % (jm, len(g), m.min(), m.max(), ', '.join('%.2f (γ %.3f, okno %d)' % (dx[i], g[i], okna[i]) for i in ix[:3])), flush=True)
json.dump(out, open(sys.argv[-1], 'w'), indent=1)
