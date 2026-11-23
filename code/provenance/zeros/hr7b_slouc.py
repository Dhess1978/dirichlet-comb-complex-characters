"""Studie VII · sloučení opravených oken hr7b (dodatek 7): kopie původní složky + okna ze složky opravy; manifest se upraví podle souborů (SHA-256, n, ocek, divz, zdroj).
Použití: python3 hr7b_slouc.py PUVODNI OPRAVA VYSTUP 'popis zdroje'"""
import json, os, sys, shutil, gzip, hashlib
V0, VO, V, ZDROJ = sys.argv[1:5]
assert not os.path.exists(V); shutil.copytree(V0, V)
M = json.load(open(os.path.join(V, 'manifest.json'))); O = json.load(open(os.path.join(VO, 'oprava.json'))); assert O['v_poradku'] and not O['zkouska']
M['oprava_dodatek_7'] = dict(oprava_skript_sha256=O['skript_sha256'], okna=[])
for w in O['okna']:
    src = os.path.join(VO, w['soubor']); raw = open(src, 'rb').read(); assert hashlib.sha256(raw).hexdigest() == w['sha256']
    b = json.loads(gzip.decompress(raw)); jm = w['L']; k = str(w['k']); dst = os.path.join(V, jm, 'okno_%06d.json.gz' % w['k'])
    stary = dict(M['L'][jm]['bloky'][k]); shutil.copyfile(src, dst)
    M['L'][jm]['bloky'][k].update(sha256=w['sha256'], n=b['n'], ocek=b['ocek'], divz=b['divz'], zdroj=ZDROJ)
    M['oprava_dodatek_7']['okna'].append(dict(L=jm, k=w['k'], puvodni_sha256=stary['sha256'], puvodni_n=stary['n'], nove_n=b['n'], pridane=w['pridane']))
    print(jm, 'okno', k, ':', stary['n'], '→', b['n'], 'nul')
json.dump(M, open(os.path.join(V, 'manifest.json'), 'w'), ensure_ascii=False, indent=1); print('sloučeno do', V)
