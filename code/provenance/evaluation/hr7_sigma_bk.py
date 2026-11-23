"""Studie VII · Σ²(ℓ) z dvoubodové funkce Bogomolného a Keatinga bez Eulerových činitelů p | q (tabulky q(u), h''(u) z hr7_teorie.py, tj. ze studie VI beze změny):
Σ²(ℓ) = ℓ + 2∫_0^{ℓ/ρ̄} (ℓ/ρ̄ − u)·[R₂(u) − ρ̄²] du·ρ̄/ρ̄ … v nerozvinutých jednotkách: Σ² = ℓ + 2∫_0^λ (λ − u)(R₂ − ρ̄²) du, λ = ℓ/ρ̄, ρ̄ = L/2π,
R₂(u) − ρ̄² = −h''(u)/(4π²) − (1 − Re[e^{iLu} q(u)])/(2π²u²). Průměr přes výšky pásma (váha L jako ve studii VI)."""
import math, numpy as np
def sigma2_bk(ells, Ls, qh, U, DU):
    q, h2 = qh; out = np.zeros(len(ells)); W = 0.0
    for L in Ls:
        rho = L / (2 * math.pi); Y = -h2 / (4 * math.pi ** 2) - (1 - (np.exp(1j * L * U) * q).real) / (2 * math.pi ** 2 * U * U)
        cY = np.cumsum(Y) * DU; cUY = np.cumsum(U * Y) * DU
        for i, l in enumerate(ells):
            lam_ = l / rho; n = min(len(U), int(lam_ / DU)); out[i] += L * (l + 2 * (lam_ * cY[n - 1] - cUY[n - 1]) * rho / rho) if n > 0 else L * l
        W += L
    return out / W
