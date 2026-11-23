# PRAVIDLA · studie VII, dodatek 10: M3 — slábne odchylka od náhodných matic s výškou tak, jak říká teorie?

Zapsáno 22. 11. 2026 v 19:56 UTC (`date -u`), **před výpočtem nul vyššího pásma** (nuly ještě neexistují); vlákno STUDIE. David 22. 11.: „Chceš provést měření na MAC a využít třeba pár hodin“ a „připrav mi příkaz pro terminal“.

## Proč
Měření M1 (dodatek 9, 19:54 UTC) ukázalo u všech 16 L-funkcí ve výšce 36 000–80 000 odchylku párové korelace od GUE, kterou popisuje dvoubodová funkce Bogomolného a Keatinga (χ² na koš 0,78–1,86 proti 6,3–16,5 u GUE). Teorie říká, že odchylka s výškou slábne zhruba jako 1/ln²(qT/2π). To se dá změřit jen na nulách z jiné výšky.

## Data (nový výpočet na Macu, asi 3 hodiny)
- Sada hr7b (4 L-funkce: q = 15, χ = [2, 1]; q = 21, χ = [3, 1], [2, 1], [4, 1]), pásmo **400 000 < γ ≤ 420 000**, okna 100, divz 32, skript `fz_nuly_chi_v3b.py` (SHA-256 f350356b…0d8b3080, beze změny), výstup `DATABASE/KS_NULY/HR7C_nuly_2026_11_22`. Čeká se asi 179 000 nul.
- Před vyhodnocením: úplnost (počty proti hladkému počtu po oknech, kontrola skoků jako v dodatku 7), druhá knihovna na vzorku, případné chybějící dvojice doplnit postupem dodatku 7. Vyhodnocovací skript dostane otisk dřív, než nuly otevřu.

## Měření a předpověď (bez volného parametru)
- Párová korelace jako v M1 (60 košů po 0,05 do d = 3, chyba z 10 bloků), rozvinutí w = (γ/2π)(ln(qγ/2π) − 1).
- Teorie pro obě pásma: `hr7_rmt_teorie_vyska.py` (8b2ba85f…) → `hr7_rmt_teorie_vyska_036.json` (1099f3d7…), `hr7_rmt_teorie_vyska_400.json` (653bbd26…). Střední kvadratická odchylka BK − GUE přes 60 košů: q = 15 dole 0,0404, nahoře 0,0280, poměr **0,693**; q = 21 dole 0,0353, nahoře 0,0257, poměr **0,727** — přesná čísla jsou v souborech a berou se z nich. (Oprava v 19:57 UTC, před výpočtem nul: v prvním znění z 19:56 UTC stálo omylem 0,78 u obou; to je poměr průměrné odchylky v prohlubni d = 0,3–0,7 (0,774 a 0,779), ne střední kvadratické odchylky přes všech 60 košů. Kritéria se nemění.)
- Pro každou L-funkci a pásmo: A = amplituda měřené odchylky (R₂ − GUE) ve tvaru teorie (BK − GUE) téhož pásma, vážený součet přes koše s chybami z bloků; velikost odchylky V = A × střední kvadratická odchylka teorie. Sdruženo přes 4 L-funkce váženým průměrem.

## Kritéria (předem)
- **odv-hr7-rmt-vyska** = O, když platí všechno: (a) nahoře χ² na koš proti BK nejvýš 2,0 u aspoň 3 ze 4 L-funkcí; (b) sdružená amplituda nahoře |A − 1| ≤ 3 chyby; (c) sdružená velikost nahoře je menší než dole (tytéž 4 L-funkce, data M1) o aspoň 3 chyby rozdílu; (d) poměr velikostí nahoře / dole se od poměru teorie liší nejvýš o 3 chyby.
- (a) nebo (b) nesplněno → Z (teorie v nové výšce nesedí); (c) nesplněno při splněném (a), (b) → N (slábnutí pod rozlišením, uvést čísla); (d) nesplněno → Z.
- Popisně: rozestupy sousedních nul nahoře proti GUE a CUE(N), největší rozdíl od GUE proti hodnotě z M1.

## Co se netvrdí
Opravy konečné výšky pro ζ jsou známé (Bogomolny, Keating 1996; Bogomolny, Bohigas, Leboeuf, Monastra 2006), pro Dirichletovy L-funkce je dal Bogomolny a Keating 2013 bez numeriky. Tady jde o měření na komplexních charakterech ve dvou výškách. Kritéria se po datech nemění.

---
## Záznam běhu a výsledek (zapsáno 22. 11. 2026 v 23:39 UTC)
- Otisky před otevřením nul horního pásma: vyhodnocení `hr7_rmt_vyska.py` (a48f0512…, 20:00:04 UTC, `SHA256_PRED_BEHEM.txt`); pomocné skripty `hr7_rmt_vyska_obr.py` (9cbe9c0f…, data pro obrázek, bez verdiktu), `graf_hr7_rmt_vyska.py` (6250cf1e…), `hr7_skoky.py` (a121fe52…, beze změny z dodatku 7) ve 20:26:08 UTC (`SHA256_PRED_BEHEM_2.txt`); `hr7c_vzorek.py` (77e5d947…) a `hr7c_mpmath.py` (630665…) ve 20:27:56 UTC (`SHA256_PRED_BEHEM_3.txt`).
- Nuly: `DATABASE/KS_NULY/HR7C_nuly_2026_11_22` (skript `fz_nuly_chi_v3b.py` f350356b…, beze změny), 800 oken, první 20:03:43, poslední a manifest 23:32:52 UTC; SHA-256 všech oken v pořádku (výpis skriptu na Macu). Nul: q15_c2-1 43 907, q21_c3-1 44 979, q21_c2-1 44 978, q21_c4-1 44 980 (celkem 178 844); okna s jemnějším krokem 1, 1, 0, 2.
- Úplnost (23:35:09–23:35:30 UTC, `hr7c_skoky.json` 157f82a1…): žádný skok, největší |Δ| 0,01 u všech čtyř. Kontrola ve vyhodnocení: rozdíl od hladkého počtu −0,62; −0,65; +0,35; +1,35, rozpětí kumulativního rozdílu 2,0–2,4 (meze 2,5 a 4) → v pořádku.
- Postupová poznámka: `hr7_skoky.py` s `python3 -I` na Macu nenašel scipy (je jen v uživatelském site-packages); spuštěn s `python3 -E` z domovské složky (bez proměnných prostředí Pythonu, bez cesty k aktuální složce dat). Skript ani data se nezměnily.
- Vyhodnocení (zmrazený `hr7_rmt_vyska.py`, 23:35:39–23:35:41 UTC, `hr7_rmt_vyska_vysledek.json` 9a392d96…, log `beh_vyska.log`):

| L-funkce | nul nahoře | A nahoře | A dole | χ² na koš proti BK nahoře | proti GUE nahoře | rozestupy: největší rozdíl od GUE nahoře / dole |
|---|---|---|---|---|---|---|
| q15_c2-1 (15.14) | 43 907 | 1,069 ± 0,042 | 0,983 ± 0,021 | 1,82 | 8,0 | 0,080 / 0,106 |
| q21_c2-1 (21.2) | 44 978 | 0,917 ± 0,056 | 0,975 ± 0,024 | 1,07 | 4,0 | 0,071 / 0,082 |
| q21_c3-1 (21.20) | 44 979 | 0,938 ± 0,058 | 1,025 ± 0,017 | 1,05 | 3,5 | 0,081 / 0,090 |
| q21_c4-1 (21.11) | 44 980 | 0,950 ± 0,067 | 0,962 ± 0,030 | 1,24 | 5,5 | 0,076 / 0,093 |

Sdruženo: A nahoře **0,987 ± 0,027** (dole 0,995 ± 0,011); velikost odchylky V dole **0,03632 ± 0,00040**, nahoře **0,02612 ± 0,00071** (teorie 0,0365 a 0,0265); rozdíl 0,0102 ± 0,0008 (12,6 chyby); poměr **0,719 ± 0,021**, teorie **0,727** (rozdíl 0,37 chyby).

**odv-hr7-rmt-vyska: O.** (a) χ² proti BK nahoře 1,05–1,82, 4 ze 4 pod 2,0; (b) |A − 1| = 0,48 chyby; (c) pokles o 12,6 chyby; (d) poměr do 0,4 chyby od teorie. Odchylka párové korelace od GUE slábne s výškou tak, jak říká dvoubodová funkce Bogomolného a Keatinga, bez volného parametru (zhruba jako 1/L², L = ln(qT/2π) z 12,2 na 14,1 pro q = 21), a tvar zůstává.
- Odchylka od postupu (popisná část): pravidla uvádějí popisně i srovnání rozestupů nahoře s CUE(N); zmrazený vyhodnocovací skript počítá jen největší rozdíl od GUE (tabulka výše: nahoře 0,071–0,081, dole 0,082–0,106 — rozestupy jsou nahoře blíž GUE). Srovnání s CUE(N) nahoře jsem nedopočítával; verdiktu se netýká.
- Data pro obrázek (`hr7_rmt_vyska_obr.py`, 23:35:51 UTC, `hr7_rmt_vyska_obr.json` c8f02515…): R₂ a rozestupy shodné s vyhodnocením u všech čtyř („shoda vse: True“). Graf `grafy/hr7_rmt_vyska_cs.png`, `_en.png`.
- Kontrola druhou knihovnou: cloud, `hr7c_mpmath.py` (630665…), vzorek `hr7c_vzorek.json` (11c3059e…, 6 nul na L-funkci, semínko 34110, vybráno na Macu ve 23:34 UTC), 23:36–00:03 UTC: |L(½ + iγ, χ)| v nulách nejvýš 5,4·10⁻⁹ (24 nul), ve středech k sousední nule nejméně 0,21 → v pořádku (výsledky `hr7c_mp_*.json`: 6b79af40…, 2417d650…, 279cbd36…, 28a5226b…).

**Studie VII: 18 hypotéz, 18 × O.**
