"""Studie VII · oprava dvou oken sady hr7b (PRAVIDLA_STUDIE_VII.md, dodatek 7).
Kontrola úplnosti Turingova typu našla v okně 620 L-funkce q15_c2-1 a v okně 517 L-funkce q21_c4-1 po jedné chybějící těsné dvojici nul
(lfunzeros s krokem střední rozestup / 32 přeskočil dvě změny znaménka v jednom kroku; schodek v okně byl menší než práh 1,5 nuly).
Skript přepočte každé z obou oken celé s jemnějším krokem (divz 128) funkcemi skriptu fz_nuly_chi_v3b.py (týž formát okna)
a úsek kolem chybějící dvojice ještě jednou s divz 512. Původní složka HR7B_nuly_2026_11_20 se nemění; opravená okna jdou do nové složky.
Použití (ze složky MAYAN ALFA V1008, s aktivním prostředím .venv):
    python3 "TEST/STUDIE VII - Mac běh/hr7b_oprava.py"
Výstup: DATABASE/KS_NULY/HR7B_oprava_2026_11_21/<L-funkce>/okno_XXXXXX.json.gz a DATABASE/KS_NULY/HR7B_oprava_2026_11_21/oprava.json
Volba --zkouska: místo celých oken jen krátké úseky kolem dvojic (zkouška kódu, výstup do složky zadané za --vystup)."""
import argparse, gzip, hashlib, json, os, sys, time
TU = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, TU)
import fz_nuly_chi_v3b as F

OPR = [((15, [2, 1]), 620, (62027.0, 62028.7)), ((21, [4, 1]), 517, (51750.6, 51752.5))]
V0_VYCH, V_VYCH = 'DATABASE/KS_NULY/HR7B_nuly_2026_11_20', 'DATABASE/KS_NULY/HR7B_oprava_2026_11_21'
DIVZ_OKNO, DIVZ_USEK = 128, 512

def nacti(cesta): return json.loads(gzip.decompress(open(cesta, 'rb').read()))

def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--puvodni', default=V0_VYCH); ap.add_argument('--vystup', default=V_VYCH); ap.add_argument('--zkouska', action='store_true')
    a_ = ap.parse_args(); V0, V = a_.puvodni, a_.vystup
    if os.path.exists(os.path.join(V, 'oprava.json')): sys.exit('Složka %s už obsahuje oprava.json — oprava je hotová, nic se nepřepisuje.' % V)
    F._init(DIVZ_OKNO, 2.0); t0 = time.time()
    out = dict(skript_sha256=hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest(),
               fz_skript_sha256=hashlib.sha256(open(os.path.join(TU, 'fz_nuly_chi_v3b.py'), 'rb').read()).hexdigest(),
               zkouska=a_.zkouska, divz_okno=DIVZ_OKNO, divz_usek=DIVZ_USEK, okna=[])
    print('Oprava dvou oken hr7b (dodatek 7)%s · začátek %s' % (' — ZKOUŠKA KÓDU' if a_.zkouska else '', time.strftime('%H:%M:%S')), flush=True)
    for (q, chi), k, (u, v) in OPR:
        jm = F.jmeno(q, chi); P = F.popis(q, chi); a, b = (u - 2.0, v + 2.0) if a_.zkouska else (k * 100.0, k * 100.0 + 100.0)
        os.makedirs(os.path.join(V, jm), exist_ok=True); cesta = os.path.join(V, jm, 'okno_%06d.json.gz' % k); t1 = time.time()
        print('  %s okno %d [%g, %g]: počítám s divz %d …' % (jm, k, a, b, DIVZ_OKNO), flush=True)
        r = F.pocitej(q, chi, P['lichy'], k, a, b, cesta)
        nove = nacti(cesta)['g']; stare = [x for x in nacti(os.path.join(V0, jm, 'okno_%06d.json.gz' % k))['g'] if a < x <= b]
        pridane = [x for x in nove if min((abs(x - y) for y in stare), default=1.0) > 1e-8]
        chybi = [y for y in stare if min((abs(x - y) for x in nove), default=1.0) > 1e-8]
        par = [(x, min(stare, key=lambda y: abs(x - y))) for x in nove if x not in pridane]
        usek = [x for x in F.nuly_okna(q, chi, u, v, DIVZ_USEK) if u < x <= v]; nove_usek = [x for x in nove if u < x <= v]
        shoda_usek = len(usek) == len(nove_usek) and all(abs(x - y) < 1e-8 for x, y in zip(usek, nove_usek))
        w = dict(L=jm, q=q, chi=chi, k=k, a=a, b=b, n_puvodni=len(stare), n_nove=len(nove), ocek=r['ocek'], divz=r['divz'], pokusu=r['pokusu'],
                 pridane=pridane, chybejici_z_puvodnich=chybi, max_rozdil_spolecnych=max((abs(x - y) for x, y in par), default=0.0),
                 usek=[u, v], usek_divz512=usek, usek_shoda=shoda_usek, soubor=os.path.relpath(cesta, V), sha256=r['sha256'], sekund=round(time.time() - t1, 1))
        out['okna'].append(w)
        print('  %s okno %d: původně %d nul, nově %d (hladké očekávání %.2f); přidané %s; z původních chybí %d; společné nuly se liší nejvýš o %.1e; úsek s divz 512 %s | %.0f s'
              % (jm, k, len(stare), len(nove), r['ocek'], ['%.6f' % x for x in pridane], len(chybi), w['max_rozdil_spolecnych'], 'souhlasí' if shoda_usek else 'NESOUHLASÍ', time.time() - t1), flush=True)
    ok = all(len(w['pridane']) == 2 and not w['chybejici_z_puvodnich'] and w['usek_shoda'] and w['max_rozdil_spolecnych'] < 1e-8 for w in out['okna'])
    out['v_poradku'] = ok; out['sekund'] = round(time.time() - t0, 1)
    json.dump(out, open(os.path.join(V, 'oprava.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    print('OPRAVA %s · %s · %.0f s' % ('V POŘÁDKU (v každém okně přibyly právě 2 nuly, nic nezmizelo)' if ok else 'POZOR — výsledek pošlete do cloudu ke kontrole', time.strftime('%H:%M:%S'), time.time() - t0), flush=True)

if __name__ == '__main__':
    main()
