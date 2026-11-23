"""Compare the numbers of the English and the Czech article of study VII section by section.
Numbers are normalised (Czech decimal comma -> point, thousands separators removed); the reference list is compared by its numbers too.
Differences that are only words in one language (for example "19 times" against "devatenáctkrát") are listed in KNOWN and accepted.
usage: python3 -I check_numbers_en_cs.py EN.md CS.md        exit code 0 = no other difference"""
import re, sys
from collections import Counter
KNOWN = {('4.5', '19', 'EN'), ('4.11', '11', 'EN')}        # "19 times" = "devatenáctkrát"; "a tooth of 11" = "zubu jedenáctky" (figure 6)
def sections(t):
    t = re.sub(r'<!--.*?-->', '', t, flags=re.S); t = re.sub(r'\ndate: "[^"]*"', '', t); parts = re.split(r'\n(?=#{1,2} )', t); out = {}
    for c in parts:
        m = re.match(r'#{1,2} (\d+(?:\.\d+)?)\.? ', c)
        key = m.group(1) if m else ('header' if c.startswith('---') else ('references' if re.match(r'# (References|Literatura)', c) else c[:20]))
        out[key] = re.sub(r'\]\([^)]*\.png\)', ']', c)                # figure file names are not numbers of the text
    return out
def numbers(t, czech):
    t = t.replace('{,}', ',' if czech else '').replace('\\,', '')
    if czech:
        t = re.sub(r'(?<=\d)[  ](?=\d{3}(?!\d))', '', t); t = re.sub(r'(?<=\d),(?=\d)', '.', t)
    else:
        t = re.sub(r'(?<=\d),(?=\d{3}(?!\d))', '', t)
    return Counter(x.rstrip('.') for x in re.findall(r'\d+(?:\.\d+)?', t))
E = sections(open(sys.argv[1], encoding='utf-8').read()); C = sections(open(sys.argv[2], encoding='utf-8').read()); bad = 0
for k in E:
    if k not in C: print('section missing in the Czech text:', k); bad += 1; continue
    a, b = numbers(E[k], False), numbers(C[k], True); onlyE, onlyC = a - b, b - a
    rest = {x for x in onlyE if (k, x, 'EN') not in KNOWN} | {x for x in onlyC if (k, x, 'CS') not in KNOWN}
    if onlyE or onlyC: print('section', k, '| only EN:', dict(onlyE), '| only CS:', dict(onlyC), '| accepted' if not rest else '| DIFFERENCE')
    bad += bool(rest)
for k in C:
    if k not in E: print('extra section in the Czech text:', k); bad += 1
print('sections EN %d, CS %d, with an unexplained difference %d' % (len(E), len(C), bad)); sys.exit(1 if bad else 0)
