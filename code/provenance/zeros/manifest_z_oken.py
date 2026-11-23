"""Sestaví manifest sady hr7 v cloudové kopii ze souborů oken (původní manifest z Macu se po Ctrl+C v 10:27 UTC aktualizoval jen do ~880 oken).
Původní manifest se uloží jako manifest_mac.json; každý blok dostane SHA-256 souboru, n, ocek, divz ze souboru a zdroj (Mac / cloud – doplnění 9 oken).
Použití: python3 manifest_z_oken.py SLOZKA [SLOZKA_DOPLNKU]"""
import gzip, hashlib, json, os, shutil, sys
V = sys.argv[1]; DOP = sys.argv[2] if len(sys.argv) > 2 else None
M = json.load(open(os.path.join(V, 'manifest.json'))); assert M['format'] == 'MAYAN_V1008_FZ_NULY_CHI_v3'
if not os.path.exists(os.path.join(V, 'manifest_mac.json')): shutil.copy2(os.path.join(V, 'manifest.json'), os.path.join(V, 'manifest_mac.json'))
puv = {j: set(z['bloky']) for j, z in M['L'].items()}; dopl = 0; nesoulad = 0
for jm, z in M['L'].items():
    if DOP and os.path.isdir(os.path.join(DOP, jm)):
        for f in os.listdir(os.path.join(DOP, jm)):
            if f.endswith('.json.gz') and not os.path.exists(os.path.join(V, jm, f)): shutil.copy2(os.path.join(DOP, jm, f), os.path.join(V, jm, f)); dopl += 1
    B = {}
    for f in sorted(os.listdir(os.path.join(V, jm))):
        if not f.endswith('.json.gz'): continue
        raw = open(os.path.join(V, jm, f), 'rb').read(); b = json.loads(gzip.decompress(raw)); k = int(round(b['a'] / float(M['okno'])))
        assert f == 'okno_%06d.json.gz' % k and b['q'] == z['popis']['q'] and list(b['chi']) == list(z['popis']['chi']) and b['n'] == len(b['g'])
        sha = hashlib.sha256(raw).hexdigest(); old = z['bloky'].get(str(k))
        if old and old['sha256'] != sha: nesoulad += 1
        B[str(k)] = dict(soubor=f, sha256=sha, n=b['n'], ocek=b['ocek'], divz=b['divz'], pokusu=(old or {}).get('pokusu'), sekund=(old or {}).get('sekund'),
                         zdroj='Mac' if (DOP is None or not os.path.exists(os.path.join(DOP, jm, f))) else 'Mac, doplněno 21. 11. 00:41–00:43 UTC týmž příkazem (přeneseno zvlášť)')
    z['bloky'] = B
M['poznamka_manifest'] = 'Sestaveno v cloudu ze souborů oken (manifest_z_oken.py); původní manifest z Macu je manifest_mac.json. Doplněno oken z Macu (zvlášť přenesených): %d. Nesoulad SHA s původním manifestem: %d.' % (dopl, nesoulad)
json.dump(M, open(os.path.join(V, 'manifest.json'), 'w'), ensure_ascii=False, indent=1)
print(M['poznamka_manifest'], {j: len(z['bloky']) for j, z in M['L'].items()})
