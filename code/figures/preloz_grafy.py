"""Studie VII: anglické verze grafů 1–6 — kopie skriptů graf_*.py s přeloženými popisky (přesné náhrady řetězců v kódu; výpočty a data beze změny). Výstup: graf_*_en.py"""
import re
T = {
 'graf_hr7.py': [("'GUE bez prvočísel'", "'GUE without primes'"), ("'naivní ζ (bez odstranění p | q)'", "'naive ζ (p | q not removed)'"), ("'předpověď s prvočísly, modul %d'", "'prediction with primes, modulus %d'"),
   ("'%s (reálný, studie VI)'", "'%s (real, study VI)'"), ("'q = %d · Q proti předpovědi: %s'", "'q = %d · Q against the prediction: %s'"), ("'předpověď: arg c = arg(−χ(p))'", "'prediction: arg c = arg(−χ(p))'"),
   ("'sdružená fáze: arg c = −arg(−χ(p))'", "'conjugate phase: arg c = −arg(−χ(p))'"), ("'předpovězená fáze prvního koeficientu (°)'", "'predicted phase of the first coefficient (°)'"), ("'naměřená fáze (°)'", "'measured phase (°)'"),
   ("'melodie nese fázi: %d kombinací L-funkce × prvočíslo\\nnejvětší |c − c_pred| = %.1e'", "'the melody carries the phase: %d combinations L-function × prime\\nlargest |c − c_pred| = %.1e'"),
   ("label='data')", "label='data')"), ("'předpověď s fází χ(p)'", "'prediction with the phase χ(p)'"), ("'sdružená fáze χ̄(p)'", "'conjugate phase χ̄(p)'"),
   ("'%s, p = %d, χ(p) = %+.3f%+.3fi: χ² p = %.2g, sdružená p = %.1g'", "'%s, p = %d, χ(p) = %+.3f%+.3fi: χ² p = %.2g, conjugate p = %.1g'"), ("b.set_ylabel('hustota')", "b.set_ylabel('density')")],
 'graf_hr7_hudba.py': [("'%s: vybráno %s (Q %.2f)'", "'%s: chosen %s (Q %.2f)'"), ("'q = %d · skutečně chybí %s'", "'q = %d · actually missing %s'"), ("'žádné z kandidátů 2, 3, 5, 7 (13 mezi kandidáty není)'", "'none of the candidates 2, 3, 5, 7 (13 is not a candidate)'"),
   ("'kandidátní „ladění“: prvočísla vyňatá z hřebene (S)'", "'candidate tuning: primes removed from the comb (S)'"), ("'3.12 slyšet ladění: z nul samotných vybráno správné S u %d / %d L-funkcí (zelený pruh = pravda)'", "'3.12 hearing the tuning: the correct S chosen from the zeros alone for %d / %d L-functions (green band = truth)'"),
   ("'harmonická k (x = p^k)'", "'harmonic k (x = p^k)'"), ("'\\npopisně mimo test: '", "'\\ndescriptive, outside the test: '"), ("'správně'", "'correct'"), ("'chybně, řád je %d; x = %s^%s > 200'", "'wrong, the order is %d; x = %s^%s > 200'"),
   ("'Q(ζ%s): %d sektorů, N = %d nul\\nřád z dat: %s · předpověď %s%s'", "'Q(ζ%s): %d sectors, N = %d zeros\\norder from the data: %s · prediction %s%s'"), ("'splněna'", "'met'"), ("'NESPLNĚNA'", "'NOT MET'"),
   ("'p^k ≡ 1 (mod q): zub předpovězen (−1)'", "'p^k ≡ 1 (mod q): tooth predicted (−1)'"), ("'p^k ≢ 1: předpověď 0 (sektory se vyruší)'", "'p^k ≢ 1: prediction 0 (the sectors cancel)'"), ("'prázdný kroužek: x = p^k > 200 (mimo spolehlivý rozsah)'", "'open circle: x = p^k > 200 (outside the reliable range)'"),
   ("'3.13 slyšet těleso: prvočíslo p „zazní“ v celém orchestru až v harmonické rovné jeho řádu mod q — %d / %d řádů správně\\n(úsečky ±1 σ; body mimo osu přitlačeny k okraji)'", "'3.13 hearing the body: a prime p sounds in the whole orchestra only in the harmonic equal to its order mod q — %d / %d orders correct\\n(bars ±1 σ; points off the axis pressed to the edge)'")],
 'graf_hr7_soused11.py': [("'předpověď z prvočísel (explicitní vzorec, bez volného parametru)'", "'prediction from the primes (explicit formula, no free parameter)'"), ("'měření z nul (12 L-funkcí, 947 573 nul)'", "'measurement from the zeros (12 L-functions, 947,573 zeros)'"),
   ("'měření po odečtení sousedních prvočísel'", "'measurement minus the neighbouring primes'"), ("'11⁵ = 161 051\\nsama mlčí (zbytek %.0e)'", "'11⁵ = 161,051\\nsilent itself (residual %.0e)'"), ("'soused 161 071\\n(o 20 dál, ≡ 1 mod 21)'", "'neighbour 161,071\\n(20 further, ≡ 1 mod 21)'"),
   ("'kolem 11⁵: jedenáctka tu neleží v tělese (11⁵ ≡ 2 mod 21); naměřený „zub“ je ozvěna sousedního prvočísla · šířka rozlišení ±%.0f'", "'around 11⁵: 11 does not sound here (11⁵ ≡ 2 mod 21); the measured tooth is the echo of the neighbouring prime · resolution width ±%.0f'"),
   ("'zub jedenáctky po odečtení sousedů:\\n%.4f plné síly (čeká se 1)'", "'tooth of 11 after subtracting the neighbours:\\n%.4f of the full strength (expected 1)'"),
   ("'kolem 11⁶ = 1 771 561 (≡ 1 mod 21): zub jedenáctky je schovaný pod %d sousedními prvočísly (čárky dole) · šířka rozlišení ±%.0f'", "'around 11⁶ = 1,771,561 (≡ 1 mod 21): the tooth of 11 is hidden under %d neighbouring primes (ticks below) · resolution width ±%.0f'"),
   ("'|Z(x)| / plná síla jedenáctky'", "'|Z(x)| / full strength of 11'"), ("'Proč čtení řádu 11 v Q(ζ₂₁) selhalo: sousední prvočísla a mez rozlišení (shoda měření s předpovědí %.0e plné síly)'", "'Why reading the order of 11 in Q(ζ₂₁) failed: neighbouring primes and the limit of resolution (measurement = prediction to %.0e of the full strength)'")],
 'graf_r.py': [("'GUE bez prvočísel'", "'GUE without primes'"), ("'předpověď s prvočísly (Hannovo okno)'", "'prediction with primes (Hann window)'"), ("'data, Hannovo okno a = 3H'", "'data, Hann window a = 3H'"),
   (": hřeben s jiným tvarem okna — z_BK(τ=1) ", ": comb with another window shape — z_BK(τ=1) "), ("'z proti GUE'", "'z against GUE'"), ("'z proti předpovědi s prvočísly'", "'z against the prediction'"), ("b.legend(fontsize=8)", "b.legend(fontsize=8, loc='center right', bbox_to_anchor=(1, 0.66), framealpha=1)"), ("'podpásma '", "'sub-bands '"),
   ("τ = 1: šest odhadů chyby; R3: N_ef", "τ = 1: six error estimates; R3: N_ef"), ("f'out/hr7_R_{nm}.png'", "f'out/hr7_R_{nm}_en.png'")]}
for f, L in T.items():
    s = open(f, encoding='utf-8').read()
    for a, b in L:
        n = s.count(a); assert n >= 1, (f, a[:60]); s = s.replace(a, b)
    rest = [m for m in re.findall(r"'[^'\n]*'", s.split('\n', 2)[2]) if re.search(r'[ěščřžýáíéůúťďňó]', m) and 'Použití' not in m]
    open(f.replace('.py', '_en.py'), 'w', encoding='utf-8').write(s); print(f, 'zbyva cesky:', rest[:5])
