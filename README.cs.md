# Různé nástroje, stejný akord: hřeben rezonancí a fáze prvního řádu u Dirichletových L-funkcí s komplexními charaktery, sjednocení kruhového tělesa a robustnost hlavní rezonance studie VI (studie VII)

**David Hess (MAYAN ALFA V1008 LAB v.34.116) · 2026** · [English version](README.md)

Sedmá studie řady, po studiích I–VI ([doi:10.5281/zenodo.23815328](https://doi.org/10.5281/zenodo.23815328), [.23928017](https://doi.org/10.5281/zenodo.23928017), [.23991165](https://doi.org/10.5281/zenodo.23991165), [.24094318](https://doi.org/10.5281/zenodo.24094318), [.24256378](https://doi.org/10.5281/zenodo.24256378), [.24381946](https://doi.org/10.5281/zenodo.24381946)). Studie VI změřila hřeben rezonancí v rozptylu součtů s oknem Σ g x^{iγ} přes nuly ζ a Dirichletových L-funkcí s reálnými charaktery. Tato studie měření rozšiřuje na komplexní charaktery a ptá se, co závisí na charakteru a co jen na jeho modulu. Každé pravidlo bylo zapsáno dřív, než se nuly spočítaly, nebo — u pozdějších doplňků — dřív, než se veličina vyhodnotila.

| výsledek | hodnota |
|---|---|
| data · nuly 16 Dirichletových L-funkcí modulo 5, 7, 13, 15, 21 mezi výškami 36 000 a 80 000, spočtené pro tuto studii v PARI/GP | 1 296 658 nul (dvanáct L-funkcí s komplexními charaktery, 958 685 nul; čtyři, které doplňují rodiny modulo 15 a 21); 178 844 nul čtyř z nich mezi 400 000 a 420 000 |
| hřeben jedné L-funkce proti dvoubodové funkci bez Eulerových činitelů p \| q | průměrné z² přes 201 bodů 0,13–0,25 u dvanácti; GUE bez aritmetiky 11,5–49; 18 dvojic L-funkcí téhož modulu má stejný hřeben (0,13–0,29) |
| první řád | první Fourierův koeficient {γ ln p/2π} je −w χ(p) ln p/√p s přesností 5,2·10⁻⁵ v 56 kombinacích; sdružená fáze vyloučena (p ≤ 7·10⁻¹⁵⁴) |
| sjednocení všech charakterů modulu (nuly Dedekindovy funkce zeta tělesa Q(ζ_q), q = 5, 7, 15, 21, až 947 573 nul) | prvočíslo zní jen v harmonických, které jsou násobky jeho řádu modulo q; koeficienty na 7,7·10⁻⁵ |
| křížová kovariance dvou L-funkcí téhož modulu | rovná se prvočíselné stopě s relativní fází χ_a χ̄_b (19 chyb u ln x = 3, relativní rozdíl 4·10⁻⁵) |
| rozptyl počtu nul | nasycuje se na úrovni dané modulem (Q = 0,32–1,72 proti dvoubodové funkci) |
| robustnost hlavní rezonance studie VI (nuly ζ u t ≈ 10⁸ a 3·10¹⁰) | platí pro Hannovo okno, šest odhadů chyby i plnou kovarianci; GUE 12 až 69 chyb daleko podle odhadu chyby |
| tři slepé testy | modul a charakter, prvočísla dělící modul a řády 16 prvočísel ve čtyřech kruhových tělesech určeny ze samotných nul, 16 ze 16 v každém testu |
| náhodné matice | párová korelace každé L-funkce do tří středních rozestupů: χ² na koš 0,78–1,86 proti dvoubodové funkci Bogomolného a Keatinga, 6,3–16,5 proti samotnému GUE; rozestupy sjednocení = složení nezávislých spekter GUE na 0,010; mezi dvěma pásmy klesne odchylka od GUE na 0,719 ± 0,021 své hodnoty, předpověď 0,727 |
| úplnost | počítací test Turingova typu našel dvě dvojice těsných nul, které vyhledávací mřížka přeskočila; opraveno před vyhodnocením |
| verdikty 18 hypotéz zapsaných předem | všechny PODPOŘENO |

Tvrzení o prvním řádu, o sjednocení a o křížové kovarianci jsou věty nebo jejich přímé důsledky (Landauův vzorec, explicitní vzorec, ortogonalita charakterů); dvoubodová funkce za hřebenem, rozptylem počtu nul a párovou korelací je domněnka. Nové jsou data pro komplexní charaktery s nulovými modely a chybami, nuly, které v době zápisu pravidel neexistovaly, a slepé testy. Shoda je numerická podpora, ne důkaz; s Riemannovou hypotézou nic z toho nesouvisí.

Jak byla použita předem zapsaná pravidla, laboratoř a asistent s umělou inteligencí, popisuje oddíl 3.10 článku.

## Obsah

| složka | co obsahuje |
|---|---|
| `paper/` | článek `comb_complex_characters_en` (PDF, DOCX, Markdown, LaTeX `main.tex`); česká verze v `paper/cs/` |
| `figures/` | obrázky 1–9 anglicky (`_en`) a česky (`_cz`), PNG 300 dpi a PDF |
| `data/results/` | výsledkové soubory vyhodnocovacích skriptů beze změny (každé číslo článku; který soubor patří ke kterému oddílu, říká `data/results/README.md`) |
| `data/manifests/` | manifesty nul (SHA-256 každého okna, prostředí, otisk skriptu): 5 280 oken sady hr7, 1 760 sady hr7b se dvěma opravenými okny (`oprava_hr7b_dodatek7.json`), 800 sady hr7c |
| `data/rules/` | pravidla zapsaná před výpočtem s deseti datovanými dodatky |
| `code/` | `independent_check_vii.py` s logem, `check_text_vii.py`, `check_numbers_en_cs.py`, `figures_vii.py`; `provenance/` obsahuje skripty tak, jak běžely (viz `code/provenance/README.md`) |

## Jak výsledky zopakovat

```
pip install -r requirements.txt
python3 -I code/independent_check_vii.py .
python3 -I code/check_text_vii.py .
python3 -I code/check_numbers_en_cs.py paper/comb_complex_characters_en.md paper/cs/comb_complex_characters_cs.md
python3 -I code/figures_vii.py . figures
```

Nuly spočetl PARI/GP 2.17.2 (cypari2) skripty `code/provenance/zeros/fz_nuly_chi_v3.py` (sada hr7) a `fz_nuly_chi_v3b.py` (sady hr7b a hr7c) na stolním počítači; samotné nuly nejsou přiloženy, manifesty se SHA-256 každého okna ano, takže přepočet lze porovnat okno po okně. Nuly ζ jsou dostupné od A. Odlyzka a z LMFDB (nuly D. Platta).

Text, obrázky a data: CC BY 4.0; kód: MIT. Shrnutí licencí česky: `LICENCE.cs.md`.
