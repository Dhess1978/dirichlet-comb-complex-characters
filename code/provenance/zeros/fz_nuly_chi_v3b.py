#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""MAYAN ALFA V1008 · studie VII — vysoké nuly Dirichletových L-funkcí L(s, χ) s primitivními (i komplexními) charaktery, po oknech výšky.

Verze 3 skriptu ze studie VI (fz_nuly_L_v2.py). Jediné změny: L-funkce je dána modulem q a charakterem χ v souřadnicích znstar(q, 1) z PARI
(místo Kroneckerova symbolu χ_D), hladký počet nul používá paritu χ a do manifestu se zapisují hodnoty χ(p). Výpočet nul je týž: PARI/GP lfunzeros
(změny znaménka Hardyho funkce) s krokem = střední rozestup / divz, hlídání úplnosti proti (θ_χ(b) − θ_χ(a))/π, při schodku dvou a více nul jemnější krok (divz × 4, nejvýš dvakrát).
Komplexní charakter χ a jeho sdružený χ̄ jsou dvě různé L-funkce; nuly L(s, χ) pod reálnou osou jsou nuly L(s, χ̄) nad ní (kontroluje --test).

Sada hr7 (studie VII): q = 5 řád 4; q = 7 řád 3 a řád 6; q = 13 řád 12; q = 15 řád 4; q = 21 řád 6 — vždy χ i χ̄, celkem 12 L-funkcí.
Sada gen (generálka, mimo sadu hr7): q = 11 řád 5, χ i χ̄.
Sada hr7b (verze 3b, 20. 11. večer; doplnění celých těles Q(ζ15) a Q(ζ21)): reálné χ₋₁₅ a χ₂₁, komplexní pár řádu 6 mod 21 ([2,1] a [4,1]).

Použití:
  python3 fz_nuly_chi_v3.py --test                                                kontrola instalace (asi minuta)
  python3 fz_nuly_chi_v3.py --odhad --sada hr7 --tmin 36000 --tmax 80000 --jadra 9 odhad času (asi 1–3 minuty)
  python3 fz_nuly_chi_v3.py --sada hr7 --tmin 36000 --tmax 80000 --jadra 9 --vystup SLOZKA   hlavní běh (lze přerušit; hotová okna se příště přeskočí)
  python3 fz_nuly_chi_v3.py --souhrn --vystup SLOZKA                               stav složky

Výstup: SLOZKA/manifest.json a SLOZKA/<jméno L-funkce>/okno_XXXXXX.json.gz s obsahem {"q", "chi", "a", "b", "divz", "n", "ocek", "g": [nuly]}.
Potřebuje: Python 3 a cypari2 (PARI/GP). Nic jiného než výstupní složku nemění, nic neposílá po síti."""
import argparse, gzip, hashlib, json, math, os, platform, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp

FMT = 'MAYAN_V1008_FZ_NULY_CHI_v3'      # verze 3b mění jen sady (hr7b), formát výstupu je týž
# (q, χ v souřadnicích znstar(q, 1)); jméno složky q<q>_c<složky χ>
SADY = {'hr7': [(5, [1]), (5, [3]), (7, [2]), (7, [4]), (7, [1]), (7, [5]), (13, [1]), (13, [11]), (15, [1, 1]), (15, [3, 1]), (21, [1, 1]), (21, [5, 1])],
        'gen': [(11, [2]), (11, [8])],
        'hr7b': [(15, [2, 1]), (21, [3, 1]), (21, [2, 1]), (21, [4, 1])]}
TMAX_VYCH, OKNO_VYCH, DIVZ_VYCH = 36000.0, 100.0, 32
PRV = (2, 3, 5, 7, 11, 13)
# reference pro --test (PARI 2.17, realprecision 38, počítáno v cloudu 20. 11. 2026; nezávisle ověřeno v mpmath přes Hurwitzovu ζ: |L| < 4e-28)
TEST_PRVNI = {(5, (1,)): [6.18357819545085, 8.45722917442323], (7, (2,)): [4.35640162473628, 8.78555471449908]}
TEST_OKNO = (5, [1], 1000.0, 1005.0, [1000.7050248872, 1001.9392850865, 1002.4349854037, 1003.2889195433, 1004.0906652162])
_P = None; DIVZ = DIVZ_VYCH; ZAS = 2.0

def jmeno(q, chi): return 'q%d_c%s' % (q, '-'.join(str(int(c)) for c in chi))
def chis(chi): return '[' + ','.join(str(int(c)) for c in chi) + ']'

def _init(divz=None, zas=None):
    global _P, DIVZ, ZAS
    if divz is not None: DIVZ = divz
    if zas is not None: ZAS = zas
    import cypari2
    _P = cypari2.Pari()
    try: _P.allocatemem(64 * 2**20, int(ZAS * 2**30), silent=True)
    except TypeError: _P.allocatemem(512 * 2**20)
    try: _P.default('nbthreads', 1)
    except Exception: pass
    _P.default('realprecision', 19)

def popis(q, chi):
    """řád, parita, sdružený charakter, primitivita a χ(p) = exp(2πi·t) pro p = 2, 3, 5, 7, 11, 13 (t jako zlomek; '-' pro p | q)"""
    G = 'znstar(%d,1)' % q; c = chis(chi)
    t = [str(_P('chareval(%s,%s,%d)' % (G, c, p))) if q % p else '-' for p in PRV]
    return dict(q=q, chi=[int(x) for x in chi], rad=int(_P('charorder(%s,%s)' % (G, c))), lichy=int(_P('zncharisodd(%s,%s)' % (G, c))),
                vodic=int(_P('zncharconductor(%s,%s)' % (G, c))), konj=[int(x) for x in _P('charconj(%s,%s)' % (G, c))], chi_p_t=t)

def theta(q, lichy, T):
    """θ_χ(T) = T/2 · ln(q/π) + Im lnΓ((½ + κ)/2 + iT/2), κ = 1 pro lichý χ; počet nul v (a, b] je (θ(b) − θ(a))/π + (S(b) − S(a))"""
    if T <= 0: return 0.0
    return 0.5 * T * math.log(q / math.pi) + float(_P.imag(_P.lngamma(_P('(1/2+%d)/2' % lichy) + _P('I') * T / 2)))

def nuly_okna(q, chi, a, b, divz):
    for pokus in range(3):
        try:
            z = _P('lfunzeros(lfuncreate([znstar(%d,1),%s]), [%s, %s], %d)' % (q, chis(chi), repr(float(a)), repr(float(b)), divz))
            return sorted(float(v) for v in z)
        except Exception as e:      # přetečení zásobníku PARI → větší zásobník a znovu
            if 'stack' not in str(e) or pokus == 2: raise RuntimeError('q = %d, χ = %s, okno [%g, %g]: %s (zkuste větší --zasobnik nebo méně --jadra)' % (q, chis(chi), a, b, str(e)[:200]))
            try: _P.allocatemem(64 * 2**20, int(ZAS * 2**(31 + pokus)), silent=True)
            except TypeError: _P.allocatemem(int(ZAS * 2**(31 + pokus)))

def okno(q, chi, lichy, a, b):
    """nuly v (a, b] s hlídáním úplnosti: při schodku dvou a více nul proti hladkému očekávání jemnější krok"""
    ocek = (theta(q, lichy, b) - theta(q, lichy, a)) / math.pi; divz = DIVZ; g = nuly_okna(q, chi, a, b, divz); pokusu = 1
    while len(g) < ocek - 1.5 and pokusu < 3:
        divz *= 4; g2 = nuly_okna(q, chi, a, b, divz); pokusu += 1
        if len(g2) >= len(g): g = g2
        else: divz //= 4; break
    g = [v for v in g if a < v <= b] if a > 0 else [v for v in g if 1e-9 < v <= b]
    return g, ocek, divz, pokusu

def pocitej(q, chi, lichy, k, a, b, cesta):
    if _P is None: _init()
    t0 = time.time(); g, ocek, divz, pokusu = okno(q, chi, lichy, a, b)
    data = (json.dumps(dict(q=q, chi=list(chi), a=a, b=b, divz=divz, n=len(g), ocek=round(ocek, 4), g=[float('%.15g' % v) for v in g]), separators=(',', ':')) + '\n').encode('utf-8'); tmp = cesta + '.tmp'
    with gzip.GzipFile(tmp, 'wb', mtime=0) as f: f.write(data)
    os.replace(tmp, cesta)
    return dict(j=jmeno(q, chi), k=k, soubor=os.path.basename(cesta), sha256=hashlib.sha256(open(cesta, 'rb').read()).hexdigest(), n=len(g), ocek=round(ocek, 4), divz=divz, pokusu=pokusu, sekund=round(time.time() - t0, 2))

def manifest_cti(V):
    p = os.path.join(V, 'manifest.json')
    return json.load(open(p, encoding='utf-8')) if os.path.exists(p) else None

def manifest_pis(V, M):
    p = os.path.join(V, 'manifest.json'); tmp = p + '.tmp'
    json.dump(M, open(tmp, 'w', encoding='utf-8'), ensure_ascii=False, indent=1); os.replace(tmp, p)

def prostredi():
    import cypari2
    if _P is None: _init()
    return dict(python=platform.python_version(), system=platform.system(), stroj=platform.machine(), pari=str(_P.version()), cypari2=getattr(cypari2, '__version__', '?'),
                skript_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest())

def test():
    _init(); chyba = 0.0
    for (q, chi), ref in TEST_PRVNI.items():
        g = nuly_okna(q, list(chi), 0.0, 12.0, DIVZ)[:2]; d = max(abs(x - y) for x, y in zip(g, ref)) if len(g) == 2 else 1.0; chyba = max(chyba, d)
        print('  q = %d, χ = %s: první nuly %s (reference %s) %s' % (q, list(chi), [round(v, 10) for v in g], ref, 'OK' if d < 1e-9 else 'ROZDÍL'))
    if TEST_OKNO:
        q, chi, a, b, ref = TEST_OKNO; P = popis(q, chi); g, ocek, divz, pokusu = okno(q, chi, P['lichy'], a, b); d = max(abs(x - y) for x, y in zip(g, ref)) if len(g) == len(ref) else 1.0; chyba = max(chyba, d)
        print('  q = %d, χ = %s, okno [%g, %g]: %d nul (hladké očekávání %.2f) %s' % (q, chi, a, b, len(g), ocek, 'OK' if d < 1e-6 else 'ROZDÍL'))
    # nuly L(s, χ) pod osou = nuly L(s, χ̄) nad osou
    q, chi = 5, [1]; P = popis(q, chi); dole = sorted(-v for v in nuly_okna(q, chi, -30.0, -0.5, DIVZ)); nahore = nuly_okna(q, P['konj'], 0.5, 30.0, DIVZ)
    d = max(abs(x - y) for x, y in zip(dole, nahore)) if len(dole) == len(nahore) and dole else 1.0; chyba = max(chyba, d)
    print('  q = 5: nuly χ = [1] pod osou (%d) proti nulám χ̄ = %s nad osou (%d): největší rozdíl %.1e %s' % (len(dole), P['konj'], len(nahore), d, 'OK' if d < 1e-9 else 'ROZDÍL'))
    for nm, S in SADY.items():
        for q, chi in S:
            P = popis(q, chi); assert P['vodic'] == q and (P['rad'] > 2 or nm == 'hr7b'), (q, chi, P)
        print('  sada %s:' % nm, {jmeno(q, chi): 'řád %d, %s, χ(p) pro p = %s: t = %s' % (popis(q, chi)['rad'], 'lichý' if popis(q, chi)['lichy'] else 'sudý', PRV, popis(q, chi)['chi_p_t']) for q, chi in S})
    print('TEST', 'v pořádku' if chyba < 1e-6 else 'NESOUHLASÍ — pošlete prosím tento výpis', '(největší rozdíl %.1e)' % chyba)
    return chyba < 1e-6

def odhad(L, tmax, sirka, jadra, tmin=0.0):
    _init(); cs = []; qs = sorted(set(q for q, _ in L))
    print('Odhad rychlosti (okno šířky 20 v polovině výšky u nejmenšího a největšího q; čas na nulu roste zhruba s odmocninou z q; skutečný čas se může lišit o desítky procent):')
    for q in dict.fromkeys([qs[0], qs[-1]]):
        chi = [c for qq, c in L if qq == q][0]; P = popis(q, chi); a = 0.5 * (tmin + tmax); t0 = time.time(); g, ocek, divz, pokusu = okno(q, chi, P['lichy'], a, a + 20.0); dt = time.time() - t0
        s1 = dt / max(1, len(g)); cs.append(s1 / math.sqrt(q)); print('  q = %d, okno [%g, %g]: %d nul za %.1f s (%.3f s na nulu)' % (q, a, a + 20.0, len(g), dt, s1), flush=True)
    c = sum(cs) / len(cs); nul = {jmeno(q, chi): (theta(q, popis(q, chi)['lichy'], tmax) - theta(q, popis(q, chi)['lichy'], tmin)) / math.pi for q, chi in L}
    sek = sum(nul[jmeno(q, chi)] * c * math.sqrt(q) for q, chi in L); h = sek / 3600 / max(1, jadra)
    print('  celkem asi %d nul u %d L-funkcí ve výšce %g až %g · asi %.1f h na %d jádrech (%.0f procesorových hodin)' % (round(sum(nul.values())), len(L), tmin, tmax, h, jadra, sek / 3600))
    print('  (okna šířky 20 mají větší režii než okna šířky %g v hlavním běhu; běh lze kdykoli přerušit — okna se počítají odspodu u všech L-funkcí zároveň)' % sirka)

def souhrn(V):
    M = manifest_cti(V)
    if not M: print('Ve složce', V, 'zatím není manifest.json.'); return
    print('Složka %s · formát %s · výška %g až %g, okna %g, divz %d' % (V, M['format'], M.get('tmin', 0.0), M['tmax'], M['okno'], M['divz']))
    for j, z in M['L'].items():
        B = z['bloky']; ok = 0
        for k, b in B.items():
            p = os.path.join(V, j, b['soubor'])
            if os.path.exists(p) and hashlib.sha256(open(p, 'rb').read()).hexdigest() == b['sha256']: ok += 1
        n = sum(b['n'] for b in B.values()); e = sum(b['ocek'] for b in B.values()); jem = sum(1 for b in B.values() if b['divz'] != M['divz'])
        print('  %-12s oken %d z %d (SHA-256 v pořádku %d), nul %d (hladké očekávání %.1f, rozdíl %+.1f), okna s jemnějším krokem %d, %.2f h procesoru' % (j, len(B), z['oken'], ok, n, e, n - e, jem, sum(b['sekund'] for b in B.values()) / 3600))

def main():
    ap = argparse.ArgumentParser(description='MAYAN ALFA V1008 · studie VII — vysoké nuly Dirichletových L-funkcí s primitivními charaktery (PARI/GP)')
    ap.add_argument('--sada', default='', help='sada L-funkcí: ' + ', '.join(SADY))
    ap.add_argument('--tmin', type=float, default=0.0, help='výška, od které se nuly počítají (výchozí 0; musí být násobkem --okno)')
    ap.add_argument('--tmax', type=float, default=TMAX_VYCH, help='výška, do které se nuly počítají (výchozí 36 000)'); ap.add_argument('--okno', type=float, default=OKNO_VYCH, help='šířka okna = jednoho bloku (výchozí 100)')
    ap.add_argument('--divz', type=int, default=DIVZ_VYCH, help='krok hledání lfunzeros = střední rozestup / divz (výchozí 32)'); ap.add_argument('--zasobnik', type=float, default=2, help='největší zásobník PARI na proces v GB')
    ap.add_argument('--jadra', type=int, default=max(1, (os.cpu_count() or 2) - 1)); ap.add_argument('--vystup', default=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hr7_vystup'))
    ap.add_argument('--test', action='store_true'); ap.add_argument('--odhad', action='store_true'); ap.add_argument('--souhrn', action='store_true')
    a = ap.parse_args()
    if a.test: _init(a.divz, a.zasobnik); sys.exit(0 if test() else 1)
    if a.souhrn: souhrn(a.vystup); return
    if not a.sada: ap.print_help(); return
    L = SADY[a.sada]
    if a.odhad: _init(a.divz, a.zasobnik); odhad(L, a.tmax, a.okno, a.jadra, a.tmin); return
    k0 = int(round(a.tmin / a.okno))
    if abs(k0 * a.okno - a.tmin) > 1e-9 or a.tmin >= a.tmax: sys.exit('--tmin musí být násobkem --okno a menší než --tmax')
    _init(a.divz, a.zasobnik); PO = {jmeno(q, chi): popis(q, chi) for q, chi in L}
    for j, P in PO.items(): assert P['vodic'] == P['q'], '%s není primitivní' % j
    V = a.vystup; os.makedirs(V, exist_ok=True); E = prostredi()
    M = manifest_cti(V) or dict(format=FMT, prostredi=E, sada=a.sada, tmin=a.tmin, tmax=a.tmax, okno=a.okno, divz=a.divz, L={})
    if M['format'] != FMT or M.get('tmin') != a.tmin or M['tmax'] != a.tmax or M['okno'] != a.okno or M['divz'] != a.divz: sys.exit('Složka %s patří jinému nastavení (tmin, tmax, okno nebo divz) — použijte jinou složku --vystup.' % V)
    M['prostredi_posledni'] = E; K = int(math.ceil(a.tmax / a.okno)); ulohy = []
    for q, chi in L:
        j = jmeno(q, chi); z = M['L'].setdefault(j, dict(popis=PO[j], oken=K - k0, bloky={})); os.makedirs(os.path.join(V, j), exist_ok=True)
        for k in range(k0, K):
            cesta = os.path.join(V, j, 'okno_%06d.json.gz' % k); hot = z['bloky'].get(str(k))
            if hot and os.path.exists(cesta) and hashlib.sha256(open(cesta, 'rb').read()).hexdigest() == hot['sha256']: continue
            ulohy.append((q, list(chi), PO[j]['lichy'], k, k * a.okno, min((k + 1) * a.okno, a.tmax), cesta))
    ulohy.sort(key=lambda u: (u[3], u[0]))      # napřed nízká okna všech L-funkcí: při přerušení jsou všechny hotové do podobné výšky
    manifest_pis(V, M); print('Oken celkem %d (výška %g až %g), zbývá spočítat %d; procesů %d; výstup %s' % ((K - k0) * len(L), a.tmin, a.tmax, len(ulohy), a.jadra, V), flush=True)
    if not ulohy: souhrn(V); return
    t0 = time.time(); hotovo = 0
    try:
        with ProcessPoolExecutor(max_workers=a.jadra, mp_context=mp.get_context('spawn'), initializer=_init, initargs=(a.divz, a.zasobnik)) as ex:
            fut = [ex.submit(pocitej, *u) for u in ulohy]
            for f in as_completed(fut):
                r = f.result(); hotovo += 1; M['L'][r['j']]['bloky'][str(r['k'])] = {x: r[x] for x in ('soubor', 'sha256', 'n', 'ocek', 'divz', 'pokusu', 'sekund')}
                if hotovo % 20 == 0 or hotovo == len(ulohy):
                    manifest_pis(V, M); el = time.time() - t0
                    print('  %d / %d oken · %.2f h · zbývá asi %.2f h' % (hotovo, len(ulohy), el / 3600, el / hotovo * (len(ulohy) - hotovo) / 3600), flush=True)
    except KeyboardInterrupt:
        print('\nPřerušeno — hotová okna jsou uložená, příště se přeskočí.'); manifest_pis(V, M); sys.exit(1)
    manifest_pis(V, M); print('Hotovo za %.2f h.' % ((time.time() - t0) / 3600)); souhrn(V)

if __name__ == '__main__':
    main()
