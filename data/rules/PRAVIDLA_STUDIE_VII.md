# PRAVIDLA · studie VII — „různé nástroje, stejný akord“: hřeben rezonancí L-funkcí s komplexními charaktery a robustnost hlavní rezonance studie VI

Zapsáno 20. 11. 2026 v 07:58 UTC (8:58 našeho času), vlákno STUDIE. **Před výpočtem nul sady hr7.** Vyhodnocovací skripty se zmrazí otiskem SHA-256 dřív, než se data z Macu přenesou do cloudu; otisky se zapíšou dodatkem k těmto pravidlům.

## 0. Zadání
- David 20. 11. 2026: „pokračuj na studii VII“. Z nabídnutých směrů zvolil **hřeben rezonancí pro L-funkce s komplexními charaktery a většími vodiči**; jako součást VII i **tabulku robustnosti hlavní rezonance studie VI** (náměty k 46 σ, `claude/POZNAMKA_NAMETY_46_SIGMA_19_11.md`).
- David 20. 11.: „rád bych využil na 1 den MAC“ → výpočet nul sady hr7 na Macu (odhad 16–17 h na 9 jádrech, přesné číslo dá `--odhad` na Macu).
- **Oprava očekávání (sděleno Davidovi před psaním pravidel):** v nabídce stálo, že fáze χ(A)·χ̄(B) rozhodují, které zuby zní. Teorie (Bogomolny a Keating 2013) i data studie VI to nepodporují: dvoubodová funkce jedné L(s, χ) se od funkce ζ liší jen konečnými součiny přes prvočísla p | q. Studie VII proto testuje, že **hřeben určuje modul, ne charakter**.

## 1. Co je známo
- **Bogomolny, Keating, „Two-point correlation function for Dirichlet L-functions“, J. Phys. A 46 (2013) 095202, arXiv:1307.6010:** dvoubodová funkce nul jedné L(s, χ) s primitivním χ (heuristika Hardyho a Littlewoodova typu pro dvojice prvočísel v aritmetických posloupnostech); od ζ se liší „konečnými součiny přes prvočísla, která dělí modul“. Číselné srovnání pro L-funkce v článku není (ověřeno čtením 20. 11.).
- **Studie VI, část II (DOI 10.5281/zenodo.24381946), hypotéza odv-fz3-hreben PODPOŘENO:** hřeben osmi L-funkcí s reálnými χ_D (D = −3, −4, 5, −7, 8, 13, 17, −23) ve výšce 36 000–80 000 souhlasí s V_BK,D. Popisně (přečteno 20. 11. ze zveřejněných výsledků): zuby s χ(A·B) = −1 nejsou obrácené (např. D = −4, poměr 3: 3,81 ± 0,20 proti 3,79).
- **Landauův vzorec pro L(s, χ)** (odvození z −L′/L(s, χ) = Σ χ(n)Λ(n)n^{−s} obdélníkovou integrací): Σ_{0<γ≤T} x^{ρ} = −(T/2π)·χ(x)Λ(x) + O(log T). Odtud první Fourierův koeficient rozdělení {γ·ln p/2π}: **c = −w·χ(p)·ln p/√p**, w = (T₂ − T₁)/(2πN). Pro komplexní χ(p) = e^{iφ} se prohlubeň hustoty posune do θ = φ/2π („melodie nese fázi“). Pro reálná χ je to odv-fz3-w studie VI (PODPOŘENO, největší rozdíl 2,1·10⁻⁵).
- Rešerše, zda někdo měřil hřeben nebo první řád u komplexních charakterů numericky, **zatím neproběhla do konce**; tvrzení o novosti se do článku napíše až po ní.

## 2. Data
### 2.1 Sada hr7 (nová, nikým neviděná)
Skript `fz_nuly_chi_v3.py` (SHA-256 `cf90cceb5a96f25372ae42703b22a653a4532d4156c6a2f44b635effdcd3be43`; verze 3 skriptu studie VI, PARI/GP lfunzeros, krok = střední rozestup/32, okna 100, hlídání úplnosti proti (θ_χ(b) − θ_χ(a))/π). Výška **36 000 < γ ≤ T₂**, T₂ = výška hotová u všech dvanácti (cíl 80 000; vyhodnotí se, co bude hotové, nejméně 60 000).

| L-funkce | q | χ v znstar(q,1) | řád | parita | χ(p) = e^{2πit}: t pro p = 2, 3, 5, 7, 11, 13 |
|---|---|---|---|---|---|
| q5_c1 | 5 | [1] | 4 | lichý | 1/4, 3/4, –, 1/4, 0, 3/4 |
| q5_c3 | 5 | [3] = χ̄ | 4 | lichý | 3/4, 1/4, –, 3/4, 0, 1/4 |
| q7_c2 | 7 | [2] | 3 | sudý | 2/3, 1/3, 2/3, –, 1/3, 0 |
| q7_c4 | 7 | [4] = χ̄ | 3 | sudý | 1/3, 2/3, 1/3, –, 2/3, 0 |
| q7_c1 | 7 | [1] | 6 | lichý | 1/3, 1/6, 5/6, –, 2/3, 1/2 |
| q7_c5 | 7 | [5] = χ̄ | 6 | lichý | 2/3, 5/6, 1/6, –, 1/3, 1/2 |
| q13_c1 | 13 | [1] | 12 | lichý | 1/12, 1/3, 3/4, 11/12, 7/12, – |
| q13_c11 | 13 | [11] = χ̄ | 12 | lichý | 11/12, 2/3, 1/4, 1/12, 5/12, – |
| q15_c1-1 | 15 | [1,1] | 4 | sudý | 3/4, –, –, 1/4, 1/2, 3/4 |
| q15_c3-1 | 15 | [3,1] = χ̄ | 4 | sudý | 1/4, –, –, 3/4, 1/2, 1/4 |
| q21_c1-1 | 21 | [1,1] | 6 | sudý | 5/6, –, 1/3, –, 1/6, 1/2 |
| q21_c5-1 | 21 | [5,1] = χ̄ | 6 | sudý | 1/6, –, 2/3, –, 5/6, 1/2 |

Všechny charaktery jsou primitivní (vodič = q, kontroluje skript). Nuly L(s, χ) pod reálnou osou jsou nuly L(s, χ̄) nad ní (`--test`: shoda 0,0 na 11 nulách q = 5); proto „obě strany výšky“ = dvojice χ, χ̄.

### 2.2 Reálné reference ze studie VI (už viděné)
Nuly χ₅, χ₋₇, χ₁₃ ve výšce 36 000–80 000 (`FZ_L_nuly_2026_11_16`, studie VI). Slouží jen jako reference v testu bez teorie (3.2); jejich hřeben se přepočte tímtéž kódem a musí bit po bitu souhlasit se zveřejněnými čísly studie VI (kontrola kódu).

### 2.3 Generálka (mimo sadu hr7, do verdiktů nevstupuje)
q = 11, řád 5, χ = [2] a χ̄ = [8], 0 < γ ≤ 8 000, počítáno v cloudu týmž skriptem 20. 11. Slouží ke zkoušce vyhodnocovacích skriptů (včetně konvence χ(p) proti χ̄(p) v prvním řádu). Kdyby generálka ukázala chybu v odvození nebo kódu, oprava se zapíše dodatkem **před** přenosem dat hr7.

### 2.4 Co už bylo spočteno z pásma hr7 (zveřejňuje se)
Jen kvůli času a funkčnosti, bez jakékoli statistiky (vypsány jen počty nul a časy): q5_c1 v [36 000, 36 040] a [58 000, 58 020], q13_c1 v [36 000, 36 040], q15_c1-1 a q21_c1-1 v [58 000, 58 020] (`--odhad`); mimo pásmo q5_c1 v [1 000, 1 200], [80 000, 80 040], q7 a q13 v [1 000, 1 200] a q = 61 v [36 000, 36 040]. Tyto nuly se nikam neukládaly.

## 3. Hypotézy sady hr7 a meze verdiktů
Verdikty: **O** = PODPOŘENO, **Z** = ZAMÍTNUTO, **N** = NEROZHODNUTO.

### 3.1 odv-hr7-hreben — hřeben jedné L(s, χ) podle modulu
Měření beze změny podle studie VI (`fz3_c.py`), jen |D| → q: Gaussova okna H = 20, středy T₁ + 9H … T₂ − 9H po H, ln x = ln(qT₀/2π) + δ, δ = −2,5 … 2,5 po 0,025 (201 bodů), V(δ) = ⟨|S − ⟨S⟩|²⟩/⟨N₂⟩, chyba = větší z V/√M (M nezávislých oken po 6H) a směrodatné odchylky deseti podpásem/√10. Q = průměr z² přes 201 bodů.
Předpovědi bez volného parametru (`hr7_teorie.py`, SHA-256 `ed7c3612…`, výstup `b3a620d6…`): V_BK,q (ζ bez Eulerových činitelů p | q; pro q = 5, 7, 13 bit po bitu totožné s V_BK,D studie VI), naivní ζ, GUE bez aritmetiky. Očekávaná síla při M = 364: Q(naivní ζ) 10,9 / 3,3 / 0,6 / 62,0 / 50,2 a Q(GUE) 24,9 / 18,3 / 11,9 / 48,7 / 41,8 pro q = 5, 7, 13, 15, 21.
- **O**, jestliže Q ≤ 3 u všech 12 L-funkcí, Q_naiv ≥ 2Q u všech šesti L-funkcí s q ∈ {5, 15, 21} a Q_GUE ≥ 2Q u všech 12.
- **Z**, jestliže Q > 10 aspoň u dvou L-funkcí, nebo Q_naiv < Q aspoň u dvou ze šesti L-funkcí s q ∈ {5, 15, 21}.
- jinak **N**.

### 3.2 odv-hr7-nastroje — stejný modul, stejný hřeben (bez teorie)
Pro dvojici L-funkcí a, b téhož modulu: Q_ab = průměr přes 201 bodů (V_a − V_b)²/(e_a² + e_b²). Dvojice: q = 5 {q5_c1, q5_c3, χ₅}: 3; q = 7 {q7_c2, q7_c4, q7_c1, q7_c5, χ₋₇}: 10; q = 13 {q13_c1, q13_c11, χ₁₃}: 3; q = 15: 1; q = 21: 1 — celkem 18 dvojic.
Kontrola síly (různé moduly, očekávané Q_ab z teorie): q5_c1 × q15_c1-1 (15,8), q7_c1 × q21_c1-1 (14,7), q13_c1 × q15_c3-1 (20,0).
- **O**, jestliže Q_ab ≤ 3 u všech 18 dvojic téhož modulu a Q_ab ≥ 5 u všech tří kontrolních dvojic.
- **Z**, jestliže Q_ab > 10 aspoň u dvou dvojic téhož modulu.
- jinak **N** (mimo jiné, když kontrola síly neprojde: pak test neumí rozlišit a verdikt nelze vynést).
Popisně: tabulka zubů A/B ∈ {1/3, 1/2, 2/3, 1, 3/2, 2, 3, 5} s hodnotami χ(A)·χ̄(B) a V − V_BK,q u každé L-funkce.

### 3.3 odv-hr7-faze — první řád nese fázi charakteru
Pro každou L-funkci a p ∈ {2, 3, 5, 7, 11, 13}, p ∤ q (celkem 56 kombinací: q = 5 a 13 po 5 prvočíslech, q = 7 po 5, q = 15 a 21 po 4): u = {γ·ln p/2π} (fáze počítaná jako ve studii VI, `faze` s Tref = T₁), c = ⟨e^{2πiu}⟩ (komplexní), předpověď c_pred = −w·χ(p)·ln p/√p, w = (T₂ − T₁)/(2πN). Histogram 50 košů proti hustotě f(θ) = 1 − 2w·Σ_{k=1}^{60} (ln p/p^{k/2})·Re(χ(p)^k·e^{−2πikθ}); chyba koše = větší z √E a sd nulového modelu (200 promíchání rozestupů), jako ve studii VI. Alternativy: (i) sdružená fáze χ̄(p) (prohlubeň zrcadlená), (ii) rovnoměrné rozdělení.
- **O**, jestliže max |c − c_pred| ≤ 3·10⁻⁴ přes všech 56 kombinací, p (χ², 49 st. v.) proti předpovědi ≥ 10⁻³ u všech, a p proti alternativě (i) < 10⁻⁶ u všech kombinací s |Im χ(p)| ≥ 0,5 a p ≤ 7.
- **Z**, jestliže max |c − c_pred| > 2·10⁻³, nebo p < 10⁻⁴ proti předpovědi aspoň u dvou kombinací.
- jinak **N**.

### 3.4 odv-hr7-presne — explicitní vzorec s Gaussovým oknem (kontrola kvality dat)
Jako odv-fz3-presne studie VI, s komplexním χ(n): Z = Σ_γ e^{−(γ−T₀)²/2H²}·x^{iγ}, P = −(H/√(2π))·Σ_n χ(n)Λ(n)n^{−1/2}·e^{−H² ln²(n/x)/2}·e^{−iT₀ ln(n/x)}; x ∈ {2, 3, 4, 5, 6, 7, 8, 9, 10, e^{√2}}, T₀ ∈ {40 000, 50 000, 60 000, 70 000} (je-li T₀ + 40H ≤ T₂), H ∈ {20, 50}.
- **O**, jestliže max |Z − P|/N_g ≤ 10⁻⁶; **Z**, jestliže > 10⁻⁴; jinak **N**. Popisně ostré okno (počet nul v intervalu) jako ve studii VI.

### 3.5 odv-hr7-diag — rozptyl harmonických v diagonálním režimu
Jako odv-fz3-diag studie VI: V v ln x = k ∈ {3, 4, 5, 6} (τ ≤ 0,6) proti diagonálnímu součtu s |χ(n)|² (hlavní charakter mod q), okna po H, chyba = největší z V/√M, rozptylu podpásem a nulového modelu (100 promíchání); 48 čísel.
- **O**, jestliže |z| ≤ 3 u všech 48; **Z**, jestliže |z| > 4 aspoň u dvou; jinak **N**.

## 4. Část R — robustnost hlavní rezonance studie VI (nuly ζ, pásma P8 a P3e10)
Data už byla viděna (studie VI); nová jsou jen čísla níže. Je to slabší než slepý test a článek to tak řekne. Měření jako `fz2_c_hreben.py` studie VI (H = 20, δ = −2,5 … 2,5 po 0,025, okna po 6H, M = 52 498), pokud není řečeno jinak.

### 4.1 odv-hr7-okno — jiný tvar okna
Okno Hannovo w(t) = cos²(πt/2a) pro |t| < a, a = 3H = 60 (okna po 6H se pak nepřekrývají). Předpověď V_BK s autokorelací okna A(u) = (w⋆w)(u)/(w⋆w)(0) místo e^{−u²/4H²}; před použitím se ověří v diagonálním režimu (τ ≤ 0,6) proti sítu (bez dat), jako ve studii VI.
- **O**, jestliže v obou pásmech |z_BK(τ = 1)| ≤ 2, z_GUE(τ = 1) ≥ 10 a Q_BK ≤ 3 na 201 bodech.
- **Z**, jestliže v některém pásmu |z_BK(τ = 1)| > 4 nebo Q_BK > 10.
- jinak **N**.

### 4.2 odv-hr7-podpasma — jiné odhady chyby
V τ = 1 (Gaussovo okno, H = 20, zveřejněné V): chyba z K ∈ {5, 10, 20, 40} podpásem (sd/√K), V/√M a blokový jackknife (1 000 bloků souvislých oken). Pro každý odhad z_BK a z_GUE.
- **O**, jestliže |z_BK| ≤ 2 a z_GUE ≥ 10 pro všech šest odhadů v obou pásmech.
- **Z**, jestliže |z_BK| > 4 nebo z_GUE < 3 pro některý odhad.
- jinak **N**. Popisně: rozpětí „46,1“ přes odhady.

### 4.3 odv-hr7-kovar — celý průběh s kovariancí a efektivní počet bodů
Kovarianční matice C vektoru V(δ) (201 bodů) blokovým jackknifem s B = 1 000 bloky; T² = rᵀC⁻¹r pro r = V − V_BK a r = V − V_GUE; p-hodnota Hotellingovou opravou F = T²(B − 201)/(201(B − 1)) ~ F(201, B − 201). Efektivní počet nezávislých bodů N_ef = (tr R)²/tr(R²) z korelační matice R.
- **O**, jestliže p_BK ≥ 0,01 a p_GUE < 10⁻¹⁰ v obou pásmech.
- **Z**, jestliže p_BK < 10⁻⁴ v některém pásmu.
- jinak **N**. Popisně: N_ef, Σz² proti χ²(N_ef) a globální p-hodnota pro největší |z| přes mřížku (efekt hledání jinde).

### 4.4 Kontrolní oblast (popisně)
τ = 0,5; 0,7; 0,9; 1,3; 1,6; 2,0 mimo rezonance: z_BK a z_GUE s chybami z 4.2, výslovně jako kontrola (čísla tabulky 3 studie VI).

## 5. Postup
1. Mac: `--test`, `--odhad`, hlavní běh sady hr7 (příkazy v `NAVOD_HR7.md`); výsledek se během běhu neprohlíží.
2. Cloud mezitím: generálka q = 11 (2.3), vyhodnocovací skripty `hr7_eval.py` (části 3.1–3.5) a `hr7_r.py` (část 4), nezávislé dvojče výpočtu V(δ) a c (jiný kód), zmrazení otiskem a dodatek s otisky — **před** přenosem dat hr7.
3. Přenos dat (tar na Macu → cloud), kontroly úplnosti a 40 náhodných nul každé L-funkce v mpmath (|L(½ + iγ, χ)| < 10⁻⁸), pak vyhodnocení zmrazenými skripty; průběžně grafy Davidovi.
4. Laboratoř (verze = poslední v TEST + 1, PRAVIDLA_SPOLECNA_LABORATOR), článek EN/CS, balíček, vydání.

## 6. Co by znamenalo neúspěch
- Q > 3 u komplexních L-funkcí při souhlasu reálných (studie VI) by znamenalo, že dvoubodová funkce Bogomolného a Keatinga pro komplexní charaktery neplatí v této podobě, nebo chybu v datech; rozhodne 3.4 (kvalita dat) a 3.2 (bez teorie).
- Q_ab > 3 u dvojice téhož modulu by znamenalo, že hřeben závisí na charakteru — proti teorii.
- Selhání 3.3 při úspěchu 3.4 by ukázalo chybu v odvození fáze prvního řádu.
Každý verdikt se zveřejní, jak padne.

## Dodatek 1 — zmrazení vyhodnocovacích skriptů (20. 11. 2026, 09:15 UTC = 10:15 našeho času)
Zapsáno dřív, než se jakákoli data sady hr7 přenesla do cloudu. Z běhu na Macu (start 08:03 UTC) se do této chvíle četly jen počty souborů oken (průběh), žádný obsah.

**Zmrazené skripty (SHA-256):**
- `hr7_eval.py` `a222f2bfc7655edd22301596242b83bc6665fc291e9b213a5f0188a7c4e0c8c8` — oddíly 3.1–3.5;
- `hr7_dvojce.py` `2d19c65638c06a4d45a5d2a817e92b511de2243b44234ca86f8cd08e1cb774ab` — nezávislé dvojče V(δ) a c;
- `hr7_teorie.py` `ed7c361287c777d76f203d754fea4fed401aa359542667f332089d4a7f2c263d` (výstup `b3a620d6…`);
- `graf_hr7.py` `5e2ce96762006ed2cb8e622f2d6265a22b33921599bc77d7781cde62c78f58df` — jen grafy, do verdiktů nevstupuje.

**Zkouška kódu na generálce q = 11** (nuly 1 000–3 000, `--zkouska 1000 --rychle`; mimo hodnocení): kontroly v pořádku (40 → 8 nul v mpmath, |L| ≤ 5·10⁻¹¹); explicitní vzorec s χ(n) shoda 8,8·10⁻¹³ (konvence χ(n), ne χ̄(n)); první řád sedí s fází χ(p) (p 0,58–0,99), sdružená fáze zamítnuta (p 10⁻⁸–10⁻²²); diagonální režim |z| ≤ 0,31; dvojče shodné (V 4·10⁻¹³, c 3·10⁻¹⁴). Konvence χ(x) v Landauově vzorci souhlasí i s literaturou (Banks, arXiv:2302.07073, hlavní člen −(T₂ − T₁)/(2π)·Λ(x)χ(x)). Pravidla oddílu 3 se nemění.

**Část R** byla vyhodnocena skriptem `hr7_r.py` (`a41766ae…`, napsán po pravidlech, na nulách ζ už viděných ve studii VI) v 08:10–08:27 UTC: odv-hr7-okno O, odv-hr7-podpasma O, odv-hr7-kovar O (`claude/POZNAMKA_STUDIE_VII_CAST_R_20_11.md`).

**Generálka** se dopočítává do výšky 8 000 a vyhodnotí se týmž zmrazeným skriptem popisně (do verdiktů nevstupuje).

## Dodatek 2 — doplnění podle Davida (20. 11. 2026, 15:14 UTC = 16:14 našeho času; před přenosem dat hr7 do cloudu)
David 20. 11.: „Máš ještě návrhy, co by bylo dobré do této studie doplnit pro srovnání s pomocí již existujících dat a teorií“ a „I něco, co nás povede k hamiltoniánu“; na návrhy 1–4 (a test nezávislosti sektorů) „Souhlasím s návrhy“. Rozdělení mezer podle Bogomolného a spol. (2006) se odkládá do studie VIII (konstantu pro L-funkce bez činitelů p | q je třeba nejdřív odvodit). Hamiltonián se tu nehledá; měří se podmínky, které musí splnit každý kandidát (prvočísla jako periodické dráhy s délkou ln p, explicitní vzorec jako stopový vzorec, charakter jako fáze na dráze, L-funkce jako symetrické sektory).

Data: nuly sady hr7 (pásmo 3.1), nuly χ₅ a χ₋₇ studie VI v témže pásmu, nuly ζ v témže pásmu (laboratoř, soubor `zeros_26000_t.npy`, 26 000–236 000, studie I). Všechny výpočty v cloudu, Mac nic dalšího nepočítá.

### 3.6 odv-hr7-teleso — celý orchestr: sjednocení všech charakterů modulu (Dedekindova ζ tělesa Q(ζ_q))
U₅ = nuly {ζ, χ₅, q5_c1, q5_c3}, U₇ = nuly {ζ, χ₋₇, q7_c2, q7_c4, q7_c1, q7_c5} v pásmu T₁ < γ ≤ T₂ (to jsou netriviální nuly ζ_{Q(ζ_q)}; hlavní charakter mod q má tytéž nuly jako ζ). Pro p ∈ {2, 3, 7, 11, 13} (q = 5) a p ∈ {2, 3, 5, 11, 13, 29} (q = 7) a k = 1 … 6 (66 kombinací): c_k(p) = ⟨e^{2πik·γ ln p/2π}⟩_U. Předpověď (Landauův vzorec každé L-funkce + ortogonalita charakterů): **c_k(p) = −w_U·φ(q)·[p^k ≡ 1 (mod q)]·ln p/p^{k/2}**, w_U = (T₂ − T₁)/(2πN_U). Harmonické s p^k ≢ 1 (mod q) mají vymizet: pro q = 7 a p = 2 zbývají jen k = 3, 6 (tři prohlubně v {γ ln 2/2π}), pro q = 5 a p = 2 jen k = 4 (čtyři prohlubně), p = 11 ≡ 1 (mod 5) a p = 29 ≡ 1 (mod 7) mají jednu prohlubeň plné síly φ(q)×.
Histogramy (50 košů) pro p = 2 a 3 v obou sjednoceních proti hustotě f(θ) = 1 − 2w_U·φ(q)·Σ_{k ≤ 60, p^k ≡ 1} (ln p/p^{k/2})·cos(2πkθ); chyba koše = větší z √E a sd nulového modelu (200×: v každém sektoru zvlášť promíchané rozestupy, pak sjednocení). Alternativa „bez interference“ (každý sektor jako ζ, fáze charakterů ignorovány): f(θ) = 1 − 2w_U·φ(q)·Σ_k (ln p/p^{k/2})·cos(2πkθ).
- **O**, jestliže max |c − c_pred| ≤ 3·10⁻⁴ přes všech 66 kombinací, p ≥ 10⁻³ u všech čtyř histogramů a alternativa „bez interference“ má p < 10⁻⁶ u p = 2 v obou sjednoceních.
- **Z**, jestliže max |c − c_pred| > 2·10⁻³, nebo p < 10⁻⁴ aspoň u dvou histogramů. Jinak **N**.

### 3.7 odv-hr7-krizova — sektory sdílejí dráhy: křížová kovariance dvou L-funkcí téhož modulu
Pro 18 dvojic (a, b) téhož modulu z 3.2 a ln x = k ∈ {3, 4, 5, 6} (diagonální režim, okna po H jako 3.5): C_ab(k) = ⟨(S_a − ⟨S_a⟩)·conj(S_b − ⟨S_b⟩)⟩/√(⟨N₂,a⟩⟨N₂,b⟩) (komplexní). Předpověď z explicitního vzorce (diagonála, n = m): **D_ab(k) = (H²/2π)·Σ_n χ_a(n)·conj χ_b(n)·Λ(n)²/n·e^{−H² ln²(n/x)}/√(⟨N₂,a⟩⟨N₂,b⟩)**. Teoreticky (bez dat) |D|/V = 1,00 pro k = 3 (jediné prvočíslo 19 v okně), medián 0,79 / 0,44 / 0,17 pro k = 4 / 5 / 6 — sektory nejsou nezávislé, korelují přes společné dráhy s relativní fází χ_a(p)·conj χ_b(p). Chyba každé složky (Re, Im) = větší z √(V_a V_b/(2M)) a sd deseti podpásem/√10.
- **O**, jestliže všech 144 složek má |z| ≤ 3 a u k = 3 je |C|/chyba ≥ 5 u všech 18 dvojic (nulový model „nezávislé sektory“, C = 0, zamítnut).
- **Z**, jestliže |z| > 4 aspoň u dvou složek. Jinak **N**.
Popisně: C_ab(δ) na mřížce hřebene kolem x = qT/2π (201 bodů) — pro τ ≈ 1 tu předpověď nemáme (mimodiagonální členy s Jacobiho součty a kořenovými čísly; námět pro studii VIII).

### 3.8 odv-hr7-nasyceni — Berryho nasycení podle modulu
Σ²(ℓ) každé L-funkce: nuly rozvinuté hladkým počtem θ_χ(γ)/π, počty v intervalech délky ℓ se začátky po 0,37 (jako studie VI, G3), ℓ = 30 hodnot logaritmicky 0,1 … 200, chyba = sd deseti podpásem/√10. Předpověď: Berryho vzorec (diagonální přiblížení, kód studie I/VI beze změny, `hr7_berry.py` `berry_q`) bez prvočísel p | q, H = ln(qT/2π) průměrně přes 24 výšek pásma (pro q = 1 bit po bitu totožný se studií VI). Nulové předpovědi: GUE (bez nasycení) a naivní ζ (bez odstranění p | q). Hrubý odhad síly (bez dat): Q(GUE) 435–1 041, Q(naivní ζ) 5,1 / 14,6 / 23,1 / 113 / 153 pro q = 13 / 7 / 5 / 21 / 15.
- **O**, jestliže Q_Berry ≤ 3 u všech 12 L-funkcí, Q_GUE ≥ 2Q u všech 12, Q_naiv ≥ 2Q u čtyř L-funkcí s q ∈ {15, 21} a Q_ab ≤ 3 u všech 18 dvojic téhož modulu.
- **Z**, jestliže Q_Berry > 10 aspoň u dvou L-funkcí, nebo Q_naiv < Q aspoň u dvou ze čtyř s q ∈ {15, 21}. Jinak **N**.
Popisně: sjednocení U₅ a U₇ — Σ² sjednocení proti součtu Σ² sektorů (při nezávislých sektorech by se sčítaly).

### 3.9 Kontrola proti LMFDB (ne hypotéza)
První nuly L(s, χ) všech 12 charakterů (Conreyho značení podle PARI znconreylog) z veřejné databáze LMFDB proti nulám z `fz_nuly_chi_v3.py` v cloudu. Shoda na 10⁻⁶ u všech uvedených nul = v pořádku; jinak se vyhodnocení zastaví a rozdíl se dohledá (včetně přiřazení χ a χ̄).

### 3.10 odv-hr7-okno-L — robustnost i pro L-funkce
Hannovo okno a = 3H po 6H (nepřekrývající se okna) jako v části R: hřeben každé L-funkce proti V_BK,q s autokorelací Hannova okna; kovariance 201 bodů jackknifem po jednotlivých oknech (B = M = 364 > 201) s Hotellingovou opravou.
- **O**, jestliže Q_BK ≤ 3 a p_BK ≥ 10⁻³ u všech 12 L-funkcí. **Z**, jestliže Q_BK > 10 nebo p_BK < 10⁻⁵ aspoň u dvou. Jinak **N**. Popisně N_ef.

### Postup
Kód (`hr7_eval2.py`) se vyzkouší na generálkách: q = 11 (χ, χ̄; dvojice pro 3.7) a nová generálka 2 — q = 5, χ = [1] a χ̄ = [3], 0 < γ ≤ 8 000 (cloud, kopie skriptu s doplněnou sadou `gen5`) spolu s nulami ζ a χ₅ studie VI ve výšce 2 000–8 000 (sjednocení U₅ pro 3.6). Pak se zmrazí otiskem (dodatek 3) — stále před přenosem dat hr7.

## Dodatek 3 — dvě úpravy po generálkách a zmrazení doplňku (20. 11. 2026, 15:51 UTC = 16:51 našeho času; před přenosem dat hr7)
Obě úpravy vycházejí jen z generálek (q = 11 a q = 5 nízko) a z nul ζ už viděných ve studiích I a VI; data hr7 nikdo neviděl.

**(a) 3.8 — hlavní předpověď je Σ² z dvoubodové funkce Bogomolného a Keatinga bez p | q, ne Berryho diagonální vzorec.** Na generálce q = 11 (2 000–8 000) dal Berryho vzorec Q = 3,26 a 8,81 (odchylka u ℓ ≈ 0,3–0,4, konečná výška — tentýž jev jako ve studii VI, kde Berryho vzorec měl zbytek až 52 chyb a BK ne), dvoubodová funkce BK Q = 0,32 a 0,37. Výpočet `hr7_sigma_bk.py` (Σ² = ℓ + 2∫₀^λ (λ − u)(R₂ − ρ̄²) du s tabulkami q(u), h''(u) z hr7_teorie.py) souhlasí se zveřejněnými čísly BK studie VI (pásmo P8) na ≤ 5·10⁻⁵. Naivní alternativa = táž dvoubodová funkce s q = 1 (bez odstranění p | q); Berryho vzorec se uvádí popisně. Meze verdiktu 3.8 beze změny.

**(b) 3.6 — hodnotí se jen x = p^k ≤ 200.** Chybový člen Landauova vzorce roste s x: na nulách ζ v témže pásmu 36 000–80 000 (už viděných) je |c − c_pred| ≤ 10⁻⁴ pro x ≤ 125, 2·10⁻⁴–6·10⁻⁴ pro x = 169–3 125 a 10⁻³–6·10⁻³ pro x ≥ 16 807. Hodnocené kombinace: q = 5: p = 2 (k = 1…6), 3 (k = 1…4), 7 (k = 1, 2), 11 (k = 1, 2), 13 (k = 1, 2); q = 7: p = 2 (k = 1…6), 3 (k = 1…4), 5 (k = 1…3), 11 (k = 1, 2), 13 (k = 1, 2), 29 (k = 1). Nenulové předpovědi mezi nimi: 2⁴ = 16, 3⁴ = 81, 11, 11² (q = 5) a 2³ = 8, 2⁶ = 64, 13² = 169, 29 (q = 7). Ostatní (x > 200) popisně. Meze verdiktu 3.6 beze změny.

**Generálky doplňku (popisně, do verdiktů nevstupují):**
- q = 11 (2 000–8 000): 3.7 všech 8 složek |z| ≤ 0,02, |C|/chyba u k = 3 = 6,9 (M = 48); 3.8 BK Q 0,32 / 0,37, Q_GUE ≥ 3 000, Q_naiv 13,6–19,7; 3.10 Q_BK 1,21 / 0,85 (Hotelling nelze, M = 48 < 201).
- q = 5 (2 000–4 000; nuly do 8 000 se dopočítávají): 3.6 sjednocení {ζ, χ₅, q5_c1, q5_c3}, N_U = 9 366: histogramy p = 2 a 3 p = 0,99 a 0,34, „bez interference“ p = 10⁻²⁴ a 10⁻³⁵; c pro 2⁴ = −0,0233 (předpověď −0,0236), pro 11² = −0,0290 (−0,0296); největší |c − c_pred| 1,3·10⁻³ u x = 2 — při N_U ≈ 9 000 odpovídá okrajovým členům řádu několika jednotek, pro data hr7 (N_U ≈ 3–5·10⁵) o řád menší. 3.7 |z| ≤ 0,31. 3.8 BK Q 1,45–1,76.
- 3.9 LMFDB: všech 12 L-funkcí sady hr7, prvních 25 nul každé (kladné nuly χ a záporné nuly = kladné nuly χ̄), největší rozdíl 3,6·10⁻¹⁵ (`data/hr7_lmfdb.json`, SHA-256 `f67f3176…`) — **v pořádku**.

**Zmrazené skripty doplňku (SHA-256):**
- `hr7_eval2.py` `0858c049a1c387f1eef4e4b9a4a420d8ca35089eb45f260f62d0f9128c5fb7f1` (oddíly 3.6, 3.7, 3.8, 3.10);
- `hr7_sigma_bk.py` `688b2813415662ff3410067cf9e86939a23cb4d9d58c7b1d581bbde49c24d93d`;
- `hr7_berry.py` `aca8f54f7a2ab8cf2a49444c5d26b734d3b8ca5b7f14b248fb09845acf5fd356`.
Skripty dodatku 1 (`hr7_eval.py`, `hr7_dvojce.py`, `hr7_teorie.py`) beze změny.

## Dodatek 4 — kontroly kvality dat (20. 11. 2026, 17:47 UTC = 18:47 našeho času; před přenosem dat hr7)
David 20. 11. poslal náměty jiného AI k „ultrapřesnému měření“ (data, ne zadání); po prosetí souhlasil s kontrolami níže („Souhlas“). Verdikty oddílů 3 a 4 se nemění; kontroly rozhodují jen o tom, zda jsou data v pořádku. Selže-li některá, vyhodnocení se zveřejní až po dohledání příčiny a nález se uvede v článku.

**K1 Certifikace vzorku intervalovou aritmetikou (Arb, python-flint 0.9.0).** Pro každou z 12 L-funkcí okna k ∈ {360, 404, 448, 492, 536, 580, 624, 668, 712, 756} (ta, která leží pod T₂): Hardyho funkce Z_χ(t) (`dirichlet_char(q, Conrey).hardy_z`, přesnost 80 bitů, při neurčeném znaménku 160 a 320) v bodech a, středech mezi sousedními nulami z PARI a b. Okno je certifikováno, jestliže znaménka všech sousedních bodů jsou určena a střídají se — pak v každém intervalu leží lichý počet nul, tedy aspoň jedna, a každá nula z PARI je izolována. Úplnost (že v intervalu nejsou tři nuly místo jedné) se dál opírá o hladký počet; Turingovu metodu pro L(s, χ) tu nemáme a článek to řekne.
**K2 Přesnost.** Táž okna znovu v PARI s realprecision 38 (místo 19): stejné počty a největší |Δγ| ≤ 10⁻⁹.
**K3 Krok.** Táž okna s divz 128 (místo 32): stejné počty a největší |Δγ| ≤ 10⁻⁹.
**K4 Těsné dvojice.** Pro každou L-funkci 10 nejmenších rozvinutých rozestupů v pásmu; obě nuly každé dvojice certifikovány jako v K1 (body vlevo, mezi a vpravo se střídajícím znaménkem). Popisně: nejmenší rozestup proti očekávání GUE.
**K5 Poloviny výšky (replikace, popisně).** Zmrazené skripty `hr7_eval.py` a `hr7_eval2.py` (beze změny, v režimu `--zkouska`) zvlášť na 36 000–58 000 a 58 000–T₂; uvedou se Q a odchylky obou polovin vedle celku. Čeká se stejný obraz s chybami asi √2 většími.
Skript kontrol `hr7_kontroly.py` se vyzkouší na generálce q = 11 a zmrazí otiskem před přenosem dat (dodatek 5 nebo záznam ve STAV).
*Doplněk k dodatku 4 (18:11 UTC):* skript kontrol vyzkoušen na generálce q = 11 (K1: 4 okna, 130–152 nul, vše certifikováno při 80 bitech; K2 a K3: |Δγ| 5·10⁻¹², stejné počty; K4: 10 / 10 dvojic certifikováno, nejmenší rozvinutý rozestup 0,085) a zmrazen: `hr7_kontroly.py` SHA-256 `95ebce692b226fa2507e0e67bedff5d9ed8e1ab6af45a1a8a5bbe8b6b19dcb01`. Protože K2 a K3 trvají v cloudu 4–7 min na okno už ve výšce 4 000, poběží K2 a K3 (okna 360 a 580 všech 12 L-funkcí) na Macu po hlavním běhu týmž skriptem; K1 a K4 v cloudu.

## Dodatek 5 — hudební čtení a slepé testy „použití“; doplnění celých těles Q(ζ₁₅) a Q(ζ₂₁) (20. 11. 2026, 21:37 UTC = 22:37 našeho času; před přenosem a jakýmkoli vyhodnocením dat hr7)
David 20. 11.: „Rád bych tato témata víc ve VII rozvinul, hlavně jak to použít, a udělal i podkapitoly“; na otázku, co ještě spočítat, zvolil doplnění celých těles Q(ζ₁₅) i Q(ζ₂₁). Hudební čtení (alikvotní řada, prvočíselná mez, pythagorejské koma, souzvuk celého orchestru) dostane v článku samostatný oddíl s podkapitolami; tvrdá čísla dávají tři slepé testy „použití“: z nul samotných poznat, co jsme dosud jen předpovídali.

### Nová data hr7b (Mac, po hlavním běhu)
Skript `fz_nuly_chi_v3b.py` (SHA-256 `f350356b4e1533734a20433021216a3f5973cb0cf82dcdd8d138e9ed0d8b3080`; od verze 3 se liší jen sadou hr7b a podmínkou testu, formát výstupu týž), sada hr7b: q15_c2-1 (reálný χ₋₁₅, Conrey 15.14), q21_c3-1 (reálný χ₂₁, Conrey 21.20), q21_c2-1 a q21_c4-1 (komplexní pár řádu 6, Conrey 21.2 a 21.11). Pásmo 36 000 < γ ≤ 80 000 (vyhodnotí se výška hotová u všech; nejméně 60 000). Spolu s nulami ζ, χ₋₃, χ₅, χ₋₇ (studie VI), sadou hr7 a hr7b jsou pak úplné všechny charaktery mod 15 (8 = φ(15)) i mod 21 (12 = φ(21)).

### 3.11 odv-hr7-nastroj — slyšet nástroj: který charakter hraje
Pro každou L-funkci sad hr7 a hr7b (16) z prvního řádu bez znalosti štítku: (1) dělitele modulu D̂ = {p ∈ {2, 3, 5, 7, 11, 13} : |c_p| < 0,5·w·ln p/√p}; (2) kandidáti = všechny primitivní charaktery modulů q′ ≤ 30 (q′ ≢ 2 mod 4, všechna prvočísla q′ ≤ 13), jejichž množina prvočinitelů v {2, …, 13} je D̂; (3) χ̂(p) = −c_p·√p/(w·ln p), skóre S = Σ_{p ∉ D̂} |χ′(p) − χ̂(p)|², vybere se nejmenší. Bez dat: pravý charakter má S = 0, nejbližší jiný kandidát S ≥ 6,0 (q = 5 proti charakterům mod 25: 6,75; q = 7: 12; q = 13: 8; q = 15: 6; q = 21: 8–9).
- **O**, jestliže je správně určen (q i charakter) u všech 16; **Z**, jestliže špatně aspoň u dvou; jinak **N**.

### 3.12 odv-hr7-ladeni — slyšet ladění: která malá prvočísla dělí modul
Pro každou L-funkci (16) z hřebene samotného (měření 3.1): Ŝ = argmin přes 16 podmnožin S ⊆ {2, 3, 5, 7} průměru z² dat proti V_BK s odstraněnými Eulerovými činiteli právě v S (kód `hr7_meze_teorie.py`, `tab_S`; činitele 11 a 13 se v rozsahu δ neprojeví). Bez dat: nejbližší jiná podmnožina má očekávané Q = 7,3 (q = 5), 3,3 (q = 7), 3,6 (q = 13), 9,4 (q = 15), 10,5 (q = 21).
- **O**, jestliže Ŝ = {p ≤ 7 : p | q} u všech 16; **Z**, jestliže je Ŝ špatně aspoň u dvou L-funkcí s q ∈ {5, 15, 21} (rezerva ≥ 7); jinak **N**.

### 3.13 odv-hr7-teleso2 — slyšet těleso: jak se prvočísla rozkládají v Q(ζ_q)
Sjednocení U₅, U₇ (jako 3.6), U₁₅ = {ζ, χ₋₃, χ₅, q5_c1, q5_c3, χ₋₁₅, q15_c1-1, q15_c3-1}, U₂₁ = {ζ, χ₋₃, χ₋₇, q7_c2, q7_c4, q7_c1, q7_c5, χ₂₁, q21_c1-1, q21_c5-1, q21_c2-1, q21_c4-1}. Pro každé sjednocení c_k(p) pro k = 1 … 6; šum σ̂_U = 1,4826·medián |c_k(p)| přes kombinace s p^k ≢ 1 (mod q) a p^k ≤ 200; f̂(p) = nejmenší k ≤ 6 s |c_k(p)| ≥ 10·σ̂_U. Předpověď (řád p modulo q): q = 5: f(2) = 4, f(3) = 4, f(7) = 4, f(11) = 1; q = 7: f(2) = 3, f(3) = 6, f(11) = 3, f(13) = 2, f(29) = 1; q = 15: f(2) = 4, f(7) = 4, f(11) = 2, f(31) = 1; q = 21: f(2) = 6, f(13) = 2, f(43) = 1 (16 čísel). Navíc pro U₁₅ a U₂₁ táž kontrola c_k(p) = −w_U·φ(q)·[p^k ≡ 1]·ln p/p^{k/2} pro p^k ≤ 200 jako v 3.6.
- **O**, jestliže všech 16 řádů f̂(p) souhlasí a pro U₁₅, U₂₁ je max |c − c_pred| (x ≤ 200) ≤ 3·10⁻⁴; **Z**, jestliže nesouhlasí aspoň dva řády nebo max |c − c_pred| > 2·10⁻³; jinak **N**.

### Popisně (bez verdiktu)
- **Prvočíselná mez:** Q dat proti V_BK s Φ omezeným na prvočísla ≤ P (P = 2, 3, 5, 7, 11, 13, všechna) pro q = 5, 15, 21 („kolik prvočísel potřebuje ladění nul“).
- **Negativní kontrola souzvuku:** sjednocení bez jednoho sektoru (U₇ bez ζ; U₅ bez χ₅) — harmonické k = 1, 2 se pak nevyruší; data proti součtu přes přítomné sektory.
- **Komma:** text (iracionalita ln 3/ln 2, racionální a iracionální násobky ze studie VI), bez nového měření.

### Kontroly K2 a K3 na Macu
Skript `hr7_k23.py` (SHA-256 `bc1e6e4437d93ccde8c2b5ee2accc13992503cfa40ee5b38d38750268f520a41`) — táž logika K2/K3 jako `hr7_kontroly.py`, jen bez části Arb; vyzkoušen na generálce q = 11 (|Δγ| 5·10⁻¹²). Okna 360 a 580 všech 12 L-funkcí hr7 (a 4 hr7b).

### Postup
Vyhodnocovací skript `hr7_eval3.py` (3.11–3.13 a popisné části) se napíše, vyzkouší na generálkách a zmrazí otiskem **dřív, než se spustí jakékoli vyhodnocení dat hr7** (ani 3.1 se nespustí dřív, aby výsledek hřebene neovlivnil test 3.12).

### K poznámkám jiných AI, které David 20. 11. večer poslal (data, ne zadání)
Náměty k režimu vysokých L-funkcí (přechod od jednotlivých intervalů k obálce, téměř-rezonance, kontrapunkt Eulerových činitelů, krátké proti dlouhým škálám) se zčásti kryjí s tím, co studie VII měří (nasycení, křížová kovariance, prvočíselná mez); nové hypotézy z nich se do studie VII nepřidávají. Tvrzení o čerstvých výsledcích a ohlášených důkazech (2026) se ve studii nepoužívají; spekulativní důkazy RH se podle dřívějšího rozhodnutí Davida nesledují ani necitují.

## Dodatek 6 — úprava odhadu řádu f̂(p) a zmrazení hudebních testů (20. 11. 2026, 21:49 UTC = 22:49 našeho času; před jakýmkoli vyhodnocením dat hr7)
Na generálce 2 (q = 5, 2 000–4 000, N_U = 9 366) dal odhad f̂ podle 10·σ̂ z dodatku 5 falešné výsledky na vysokých mocninách: f̂(7) nenalezen a f̂(13) = 5 (u x = 13⁵ je chybový člen Landauova vzorce velký, viz dodatek 3). **Nově:** f̂(p) = nejmenší k ≤ 6, při kterém |c_k(p)| dosáhne aspoň poloviny plné síly w_U·φ(q)·ln p/p^{k/2}; σ̂ se uvádí jen popisně. Na téže generálce pak f̂ = 4, 4, 4, 1 pro p = 2, 3, 7, 11 (předpověď 4, 4, 4, 1). Meze verdiktu 3.13 i seznam 16 řádů beze změny.
Zkouška celého skriptu na generálce 2: 3.11 obě L-funkce určeny správně (S 7·10⁻⁴ a 1·10⁻³ proti nejbližšímu jinému kandidátovi 6,7); 3.12 obě S = [5] správně (Q 0,24 a 0,17 proti 0,81 a 0,74).
**Zmrazeno (SHA-256):** `hr7_eval3.py` `f57c30fe9878af8de98010cfc420505fa3470d2fb5091a342d87c085cb7b84ea`; `hr7_meze_teorie.py` `816690289fddcab44196377d3284d0768310d73bba9897770d5cf732d6381604` (výstup teorie bez dat `a709ec5e…`).
Pořadí vyhodnocení po přenosu: kontroly úplnosti (hr7_eval.py část 1), pak `hr7_eval3.py`, `hr7_eval.py`, `hr7_eval2.py`, `hr7_dvojce.py`, kontroly K1 a K4; hr7b až po druhém běhu na Macu.

## Dodatek 7 — úplnost nul sady hr7b a oprava dvou oken (21. 11. 2026, 09:52 UTC = 10:52 našeho času; před jakýmkoli vyhodnocením dat hr7b)
Sada hr7b doběhla na Macu v 09:36 UTC (1 760 oken, SHA-256 všech oken v pořádku, skript `fz_nuly_chi_v3b.py` `f350356b…`). Před vyhodnocením jsem přidal kontrolu úplnosti Turingova typu, která v pravidlech dosud nebyla: r_n = θ_χ(γ_n)/π − n; chybějící nula se projeví trvalým skokem r o +1. Skript `hr7_skoky.py` (SHA-256 `a121fe52f1f2b2b7e2a377d4cdd4da69fbfde81e7583c9b8933575f61eac7e92`; skok = průměr r přes 400 následujících minus 400 předchozích nul, hlásí se |skok| > 0,6). Výsledek (`data/hr7_skoky.json`):
- **hr7 (12 L-funkcí): bez skoku** (největší |Δ| 0,01) — hlavní data jsou úplná; dosavadní verdikty platí.
- **hr7b: q15_c2-1 skok +2,00 v okně 620 (γ ≈ 62 028), q21_c4-1 skok +2,00 v okně 517 (γ ≈ 51 752); q21_c3-1 a q21_c2-1 bez skoku.** V obou oknech chybí těsná dvojice nul: lfunzeros s krokem střední rozestup/32 přeskočil dvě změny znaménka v jednom kroku a schodek okna (−1,49 u okna 620, −0,95 u okna 517; upřesněno 21. 11. ve 13:55 UTC) nepřekročil práh skriptu 1,5 nuly. Kontrola v cloudu na krátkém úseku obě dvojice našla (62 027,7402 a 62 027,7491; 51 751,5759 a 51 751,5860; rozestupy 0,009 a 0,010).
- **Co už bylo z hr7b vidět (poctivě):** zmrazený `hr7_eval3.py` s hr7b jsem spustil v 09:38 UTC souběžně s kontrolou úplnosti a po nálezu skoků ho zastavil. Do té doby stihl vypsat jen oddíl 3.11, i pro čtyři L-funkce hr7b — na neúplných datech: všechny čtyři určeny správně (S 4·10⁻⁷ až 1,3·10⁻⁶; `log/hr7_eval3_hr7b.log`). Nejdřív jsem log četl jen zčásti a mylně zapsal, že o hr7b nic nevypsal; opraveno v 09:55 UTC. Oddíly 3.12 a 3.13 pro hr7b se nespočítaly. Dvě chybějící nuly z 85 000 mění koeficienty prvního řádu řádově o 10⁻⁵ relativně, kdežto rezerva testu 3.11 je S ≈ 10⁻⁶ proti ≥ 6; konečný verdikt se přesto spočítá až z opravených dat.

**Postup opravy (data na Macu, podle Davidova pravidla):** skript `hr7b_oprava.py` (SHA-256 `76d6de4b96ee9d93cdade96f87112fb3c848a09caf42bb2d6321c20718e930c6`) přepočte obě okna celá funkcemi `fz_nuly_chi_v3b.py` s divz 128 (týž formát okna) a úsek kolem dvojice ještě s divz 512; výstup do nové složky `DATABASE/KS_NULY/HR7B_oprava_2026_11_21` (původní složka se nemění). Zkouška kódu v cloudu (`--zkouska`, úseky kolem dvojic): v každém okně přibyly právě 2 nuly, žádná původní nezmizela, společné nuly totožné. Cloud přepočte obě celá okna nezávisle jako kontrolu.
**Přijetí opravy:** v každém okně přibudou právě 2 nuly (dvojice výše), žádná původní nezmizí, úsek s divz 512 souhlasí; po sloučení (manifest sestavený ze souborů oken, zdroj „Mac, oprava 21. 11. (dodatek 7)“) `hr7_skoky.py` bez skoku u všech 16 L-funkcí. Teprve pak se spustí zmrazený `hr7_eval3.py` s hr7b. Meze verdiktů 3.11–3.13 beze změny.

*Doplněk k dodatku 7 (21. 11. 2026, 13:21 UTC = 14:21 našeho času) — oprava přijata.* Mac spočítal obě okna skriptem `hr7b_oprava.py` (`76d6de4b…`) v 13:03–13:18 UTC (`DATABASE/KS_NULY/HR7B_oprava_2026_11_21`, `oprava.json` SHA-256 `ebf11a59…`): q15_c2-1 okno 620 188 → 190 nul (+62 027,7401508586 a 62 027,7491002939), q21_c4-1 okno 517 191 → 193 nul (+51 751,5759459649 a 51 751,5860134948); žádná původní nula nezmizela, společné nuly totožné (rozdíl 0), úseky s divz 512 souhlasí. Obě okna jsou bit po bitu shodná s nezávislým přepočtem v cloudu (SHA-256 `60547477…`, `beba8c68…`). Po sloučení (`hr7b_slouc.py`, složka `HR7B_slouceno`, zdroj „Mac, oprava 21. 11. (dodatek 7)“) dává `hr7_skoky.py` u všech 16 L-funkcí „bez skoku“ (největší |Δ| 0,01); počty nul hr7b: 82 726, 85 082, 85 083, 85 082. Všechny podmínky přijetí splněny → spouští se zmrazený `hr7_eval3.py` (`f57c30fe…`) se sadami hr7 a hr7b (sloučenou).

## Dodatek 8 — proč popisné čtení řádu 11 v Q(ζ₂₁) vyšlo špatně: sousední prvočísla a mez rozlišení (21. 11. 2026, 16:08 UTC = 17:08 našeho času; popisně, bez verdiktu, po vyhodnocení)
David 21. 11.: „řád prvočísla 11 v Q(ζ₂₁) — můžeme tomu nějak pomoct nějakým výpočtem“. Je to analýza **po** datech: hodnoty |c₅(11)| = 0,80 a |c₆(11)| = 5,68 plné síly jsou známé (dodatek 7, konečné vyhodnocení). Verdikty se nemění, nic se nepřidává k 16 řádům testu 3.13.

**Domněnka (z aritmetiky, bez dat):** ve sjednocení U₂₁ přispívá každá mocnina prvočísla n ≡ 1 (mod 21) vahou φ(21)·Λ(n)/√n. Součet přes nuly v pásmu šířky ΔT = 44 000 rozliší dvě čísla až od vzdálenosti asi 2πx/ΔT (23 u x = 11⁵, 253 u x = 11⁶). Vedle 11⁵ = 161 051 (≡ 2) leží prvočíslo 161 071 ≡ 1 (mod 21) ve vzdálenosti 20 (odhad ozvěny při ostrém okně 0,73 plné síly) a 160 903 ve vzdálenosti 148 (0,24); v okolí ±700 kolem 11⁶ = 1 771 561 je šest prvočísel ≡ 1 (mod 21), nejbližší 1 771 687 (odhad 3,8× plné síly). Chybné čtení f̂(11) = 5 je tedy ozvěna sousedů, ne šum.

**Výpočet (skript `hr7_soused11.py`, jen čte data):** místo ostrého okna Gaussovo okno přes celé pásmo, g(γ) = exp(−(γ − T₀)²/2H²), T₀ = 58 000, H = 4 000 (okraje pásma ve vzdálenosti 5,5 H; zanedbaná váha 4·10⁻⁸). Z(x) = Σ_{γ ∈ U₂₁} g(γ)·x^{iγ}, N_g = Σ g. Předpověď z explicitního vzorce bez volného parametru (jako oddíl 3.4): P(x) = −(H/√(2π))·Σ_n a(n)·Λ(n)·n^{−1/2}·e^{−H² ln²(n/x)/2}·e^{−iT₀ ln(n/x)}, a(n) = součet hodnot 12 primitivních charakterů sjednocení (pro n nesoudělné s 21: 12·[n ≡ 1 mod 21]; pro n = 3^j: 6·[3^j ≡ 1 mod 7]; pro n = 7^j: 2). Body: x = 11^k (k = 1 … 6), mřížka 11⁵ + d (d = −400 … 400 po 1) a 11⁶ + d (d = −4 000 … 4 000 po 10). Jednotka „plná síla“ A(x) = (H/√(2π))·12·ln 11/√x.
Uvede se: (1) největší |Z − P|/N_g; (2) |Z|/A v 11⁵ a 11⁶ a rozklad předpovědi na jednotlivé sousedy; (3) zbytek po odečtení všech členů kromě n = 11^k: u 11⁵ se čeká 0, u 11⁶ se čeká −1 (zub jedenáctky plné síly) — to využívá znalost prvočísel, není to čtení „jen z nul“; (4) zkusmo čtení „z nul a znalosti q“: nejmenší čtverce na mřížce kolem 11⁵ s neznámými amplitudami na lichých číslech ≡ 1 (mod 21) a zkušebním členem přímo v 11⁵ (čeká se amplituda 0); uvede se i číslo podmíněnosti — vyjde-li nestabilně, řekne se to; (5) pravidlo rozlišení pro článek (šířka jádra x/H proti rozteči 42 lichých kandidátů ≡ 1 mod 21).
Čtení řádu 6 jen z nul by vyžadovalo pásmo výšek asi desetkrát širší (rozlišit 11⁶ od kandidátů ve vzdálenosti 42 žádá ΔT řádu 2πx/42·2 ≈ 5·10⁵); to se v této studii nedělá.

*Výsledek dodatku 8 (21. 11. 2026, 16:13 UTC; `hr7_soused11.py` SHA-256 `02392bd4…`, `data/hr7_soused11.json` `a25e0808…`, graf `hr7_soused11.png`).* Domněnka potvrzena. (1) Shoda měření s předpovědí na 1 608 bodech: |Z − P|/N_g ≤ 1,9·10⁻¹⁰, tj. ≤ 1,2·10⁻⁶ plné síly zubu jedenáctky. (2) V 11⁵: |Z| = 4,416 plné síly, celé od prvočísla 161 071 (4,42); v 11⁶: |Z| = 6,22, nejbližší sousedé 1 771 687 (5,76), 1 771 183 (4,17), 1 771 981 (3,83). (3) Zbytek po odečtení sousedů: v 11⁵ 2·10⁻⁷ (čekáno 0), v 11⁶ −1,000001 − 0,0000003i (čekáno −1). (4) Čtení „z nul a znalosti q“ u 11⁵ (29 neznámých, podmíněnost 1,8·10⁷): zkušební člen v 11⁵ má amplitudu 3·10⁻⁶ plné síly; mezi 14 lichými kandidáty ≡ 1 (mod 21) do vzdálenosti 300 vyšlo 5 prvočísel s amplitudou 1 ± 4·10⁻⁷ a 9 složených čísel s amplitudou ≤ 5·10⁻⁷ — ticho v páté harmonické tedy z nul a modulu přečíst lze. U 11⁶ (jádro ±443 proti rozteči 42) takové čtení v tomto pásmu možné není; řád 6 „jen z nul“ zůstává nepřečtený. (5) Pravidlo rozlišení: dvě čísla x a x + d rozliší nuly z pásma šířky ΔT, když d ≳ 2πx/ΔT ≈ x/7 000 (0,03 u x = 200; 23 u 11⁵; 253 u 11⁶). Popisné čtení „11 → 5“ z konečného vyhodnocení zůstává v záznamu jako chybné; tato analýza ho vysvětluje, nepřepisuje.
