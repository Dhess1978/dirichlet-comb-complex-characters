"""Figures 1-9 of study VII (numbered in the order of their first mention in the article) in English (_en) and Czech (_cz), PNG 300 dpi and PDF, drawn from the results files.
The plotting scripts in code/figures/ are the working scripts of the study (Czech labels; the _en versions differ only in the translated
labels and two legend positions, see preloz_grafy.py). This driver runs each of them unchanged in a separate process, points their input
paths to data/results, and saves every figure at 300 dpi as PNG and PDF under the published names.
usage: python3 -I code/figures_vii.py REPO OUTDIR"""
import sys, os, subprocess, tempfile, shutil
REPO, OUT = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]); os.makedirs(OUT, exist_ok=True)
F, RES, RMT = os.path.join(REPO, 'code/figures'), os.path.join(REPO, 'data/results'), os.path.join(REPO, 'code/provenance/random_matrices')
BOOT = r'''
import sys, os, runpy, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
src, ren, argv = sys.argv[1], dict(a.split('=', 1) for a in sys.argv[2].split('|') if a), sys.argv[3:]
_orig = plt.savefig
def _save(fn, *a, **k):
    base = os.path.basename(str(fn)); tgt = ren.get(base)
    if tgt is None: raise SystemExit('unexpected output ' + base)
    k['dpi'] = 300; _orig(tgt + '.png', *a, **k); k.pop('dpi'); _orig(tgt + '.pdf', *a, **k); print('  ->', os.path.basename(tgt))
plt.savefig = _save
code = open(src, encoding='utf-8').read().replace("'/home/claude/br107/data/hr7_teorie.json'", repr(os.path.join(os.environ['RES'], 'hr7_teorie.json')))
code = code.replace("f'data/hr7_r_{nm}.json'", "os.path.join(os.environ['RES'], f'hr7_r_{nm}.json')").replace('import json, sys,', 'import json, sys, os,')
sys.argv = [src] + argv; exec(compile(code, src, 'exec'), {'__name__': '__main__'})
'''
def run(script, ren, *args):
    with tempfile.TemporaryDirectory() as tmp:
        r = subprocess.run([sys.executable, '-I', '-c', BOOT, os.path.join(F, script), '|'.join('%s=%s' % (k, os.path.join(OUT, v)) for k, v in ren.items())] + list(args),
                           cwd=tmp, env=dict(os.environ, RES=RES, MPLBACKEND='Agg'), capture_output=True, text=True)
        print(script, ' '.join(args[:1]), r.stdout.strip()); assert r.returncode == 0, r.stderr[-2000:]
for lg, sfx in (('cz', ''), ('en', '_en')):
    run('graf_hr7%s.py' % sfx, {'x_hreben.png': '01_comb_%s' % lg, 'x_faze.png': '02_phase_%s' % lg}, os.path.join(RES, 'hr7_eval_hr7.json'), 'x')
    run('graf_hr7_hudba%s.py' % sfx, {'x_ladeni.png': '05_tuning_%s' % lg, 'x_teleso.png': '03_body_%s' % lg}, os.path.join(RES, 'hr7_eval3_final.json'), 'x')
    run('graf_hr7_soused11%s.py' % sfx, {'x.png': '06_resolution_%s' % lg}, os.path.join(RES, 'hr7_soused11.json'), 'x.png')
    for nm, ab in (('P8', 'a'), ('P3e10', 'b')):
        run('graf_r%s.py' % sfx, {'hr7_R_%s%s.png' % (nm, sfx): '04%s_robustness_%s_%s' % (ab, nm, lg)}, nm)
    L2 = 'cs' if lg == 'cz' else 'en'
    run('graf_hr7_rmt.py', {'x_parova_%s.png' % L2: '07_pair_correlation_%s' % lg, 'x_orchestr_%s.png' % L2: '08_unions_%s' % lg},
        *[os.path.join(RES, n) for n in ('hr7_rmt_mereni.json', 'hr7_rmt_vysledek.json', 'hr7_rmt_teorie.json')], RMT, 'x', L2)
    run('graf_hr7_rmt_vyska.py', {'x_vyska_%s.png' % L2: '09_height_%s' % lg},
        *[os.path.join(RES, n) for n in ('hr7_rmt_mereni.json', 'hr7_rmt_vyska_obr.json', 'hr7_rmt_teorie_vyska_036.json', 'hr7_rmt_teorie_vyska_400.json', 'hr7_rmt_vyska_vysledek.json')], 'x', L2)
print('ok', len(os.listdir(OUT)), 'files in', OUT)
