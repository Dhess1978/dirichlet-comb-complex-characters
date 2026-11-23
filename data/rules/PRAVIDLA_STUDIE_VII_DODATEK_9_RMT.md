# PRAVIDLA · studie VII, dodatek 9: navázání na RMT (část M) — dvě měření na nulách, které už máme

Zapsáno 22. 11. 2026 v 19:51 UTC (`date -u`), **před výpočtem**; vlákno STUDIE. David 22. 11.: „rád bych, aby tam bylo navázání i na RMT“ → volba „text + dvě měření“. Data: nuly 16 L-funkcí sad hr7 a hr7b (po opravě podle dodatku 7), 36 000 < γ ≤ 80 000, k tomu ζ a reálné charaktery modulů 3, 5, 7 ze studie VI pro sjednocení. Tatáž data už prošla testy 3.1–3.13 (hřeben, nasycení, sjednocení); rozdělení rozestupů a párová korelace na krátké vzdálenosti se na nich dosud neměřily. Počítá se na Macu (skript `hr7_rmt.py`), předpověď teorie v cloudu předem (`hr7_rmt_teorie.py` → `hr7_rmt_teorie.json`, otisky do záznamu před měřením).

## M1 — jedna L-funkce proti náhodným maticím (párová korelace na krátké vzdálenosti)
- Rozvinutí w = (γ/2π)·(ln(qγ/2π) − 1). Párová korelace R₂(d) pro 0 < d ≤ 3 v 60 koších po 0,05 (páry s levým koncem aspoň 3 před koncem pásma); chyba z 10 bloků.
- Srovnání bez volného parametru: (i) GUE, 1 − (sin πd/πd)²; (ii) dvoubodová funkce Bogomolného a Keatinga pro modul q (tabulky studie VI beze změny, bez Eulerových činitelů p | q), průměr přes koš a přes 24 výšek pásma s vahou hustoty.
- **odv-hr7-rmt-par:** χ² na koš proti (ii) nejvýš 2,0 u aspoň 14 ze 16 L-funkcí **a** GUE zřetelně horší: χ²_GUE ≥ 3·χ²_BK a χ²_GUE ≥ 4 u aspoň 14 ze 16. Splněno → O (na krátké vzdálenosti už v této výšce nestačí náhodné matice, stačí aritmetická oprava); nesplněno → Z s čísly.
- Popisně (bez verdiktu): rozdělení rozestupů sousedních nul (30 košů po 0,1) proti GUE (Fredholmův determinant sinového jádra) a proti CUE(N) pro N = 2 až 6; efektivní N z teorie (z chování R₂ u nuly) vedle N, které rozestupům sedí nejlépe; totéž pro všech 16 L-funkcí dohromady.

## M2 — celý orchestr proti složení nezávislých spekter
- Sjednocení nul všech charakterů modulu (tj. nuly Dedekindovy funkce kruhového tělesa): U₅ (4 sektory), U₇ (6), U₁₅ (8), U₂₁ (12). Rozvinutí součtem hladkých počtů sektorů. Rozdělení rozestupů sousedních nul, 40 košů po 0,1 na [0, 4].
- Model S bez volného parametru: nezávislá spektra GUE s podíly hustot f_i: pravděpodobnost prázdného úseku E_S(s) = Π E_GUE(f_i·s), p_S = E_S″. Srovnávací modely: Poisson e^{−s}, jedno spektrum GUE.
- **odv-hr7-rmt-orchestr:** u všech čtyř sjednocení největší |p − p_S| přes koše nejvýš 0,015 **a** první koš se od Poissona liší aspoň o 0,03 (odpuzování uvnitř sektorů je vidět) **a** jedno GUE je vedle aspoň o 0,2 v některém koši. Splněno → O; nesplněno → Z s čísly.
- Popisně: zbytek p − p_S s chybou z 10 bloků a χ²; kde zbytek přesahuje 3 chyby, uvést (čeká se malý vliv konečné výšky a společných prvočísel sektorů, test 3.7).

## Co se netvrdí
- Že nuly jedné L-funkce vysoko na přímce mají lokálně statistiku GUE, je známá domněnka (Montgomery, Rudnick–Sarnak) s numerikou pro ζ (Odlyzko); opravy konečné výšky pro ζ: Bogomolny, Bohigas, Leboeuf, Monastra 2006. Tady jde o měření na komplexních charakterech a na sjednocení, ne o nový zákon. Třída symetrie rodiny (unitární pro komplexní charaktery) se týká nízkých nul přes rodinu a tímto měřením se netestuje.
- Kritéria se po datech nemění; odchylky od postupu dodatkem s časem. Zákon běhu, ρ(i), konstanty, perioda 260, krok 18,45°, zrcadlo 300 a ζ beze změny.

---
## Záznam běhu a výsledek (zapsáno 22. 11. 2026 v 19:58 UTC)
- Teorie předem v cloudu 19:51 UTC: `hr7_rmt_teorie.py` (a1dc4d4e…) → `hr7_rmt_teorie.json` (27dc883d…). Skript měření `hr7_rmt.py` (220ce175…): zkouška funkce v cloudu jen na 20 oknech (36 000–38 000; čteny jen počty a střední rozestup), otisky uloženy na Macu v 19:54:12 UTC před během.
- Běh na Macu: měření 19:54:12–19:54:19 UTC (počty nul všech 20 sektorů stejné jako ve vyhodnocení 3.11–3.13; 2 opravená okna hr7b použita), vyhodnocení 19:54:23–19:54:29 UTC. Složka `TEST/STUDIE VII - Mac běh/dodatek9_rmt_2026_11_22/` (`hr7_rmt_mereni.json` 76a4d47f…, `hr7_rmt_vysledek.json` 679c4807…).

**M1 · odv-hr7-rmt-par: O.** χ² na koš proti Bogomolnému a Keatingovi 0,78–1,86 (16 ze 16 pod 2,0); proti GUE 6,3–16,5 (16 ze 16 splňuje obě podmínky). Měření − GUE u d = 0,5: −0,032 až −0,076 (teorie −0,045 až −0,066); největší odchylka koše od teorie 2,3–4,0 chyby.

Popisně, rozestupy sousedních nul (1 296 642 rozestupů 16 L-funkcí dohromady): největší rozdíl od GUE 0,083 u s = 1,15; méně malých rozestupů (−0,05 u s ≈ 0,4), víc kolem středu (+0,08 u s ≈ 1,0–1,2), méně velkých (−0,03 u s ≈ 1,8) — spektrum je tužší než GUE. χ² na koš: GUE 1 180, CUE(2) 1 931, CUE(3) 156, CUE(4) 304, CUE(5) 513, CUE(6) 672: nejblíž je CUE(3), ale žádná celá velikost matice v mezích chyby nesedí. Efektivní N z teorie (z chování R₂ u nuly) je 2,0–2,3. Popis rozestupů jednou maticí konečné velikosti tedy v této výšce nestačí; dvoubodová teorie ano.
Odchylka od postupu: χ² rozestupů po jednotlivých L-funkcích obsahuje koše s nulovou chybou bloků (dělení nulou), proto ho neuvádím; souhrn přes 16 L-funkcí tím není dotčen a verdikt se ho netýká.

**M2 · odv-hr7-rmt-orchestr: O.**
| těleso | sektorů | nul | největší \|p − p_S\| | první koš: měření / S / Poisson | největší rozdíl od jednoho GUE | χ² na koš proti S | košů nad 3 chyby |
|---|---|---|---|---|---|---|---|
| U₅ | 4 | 288 861 | 0,0096 | 0,721 / 0,731 / 0,952 | 0,71 | 3,4 | 3 |
| U₇ | 6 | 450 707 | 0,0075 | 0,798 / 0,805 / 0,952 | 0,79 | 4,9 | 3 |
| U₁₅ | 8 | 608 494 | 0,0068 | 0,837 / 0,842 / 0,952 | 0,83 | 3,2 | 2 |
| U₂₁ | 12 | 947 573 | 0,0073 | 0,872 / 0,879 / 0,952 | 0,86 | 1,9 | 0 |

Všechna tři kritéria splněna u všech čtyř. Popisně: zbytek proti S není šum (χ² na koš 1,9–4,9; největší odchylka 2,9–8,3 chyby) a v prvním koši je u všech čtyř záporný (−0,006 až −0,010): nejmenších rozestupů je o něco méně, než dává složení nezávislých spekter. Znaménkem to odpovídá vazbě mezi L-funkcemi téhož modulu přes společná prvočísla (test 3.7); velikost jsem nepočítal a netvrdím ji. Citlivost modelu S na výšku (podíly dole a nahoře v pásmu) je nejvýš 0,0002.

*Oprava zápisu (23. 11. 2026, 00:49 UTC):* v M1 stojí „teorie −0,045 až −0,066“; správně zaokrouhleno **−0,045 až −0,065** (krajní hodnota −0,065498 v `hr7_rmt_vysledek.json`). Nalezeno při stavbě listu laboratoře 34.116, kde list tuto hodnotu přepočítává živě. Verdikt se nemění. Stejná oprava je v článku (oddíl 4.12).
