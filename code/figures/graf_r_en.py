import json, sys, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
nm = sys.argv[1]; R = json.load(open(f'data/hr7_r_{nm}.json')); d = np.array(R['delta']); r1 = R['R1']
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw=dict(width_ratios=[2.2, 1]))
a = ax[0]; a.plot(d, r1['V_GUE'], '--', color='#d95926', lw=1.5, label='GUE without primes'); a.plot(d, r1['V_BK'], color='#199e70', lw=2, label='prediction with primes (Hann window)')
a.errorbar(d, r1['V'], yerr=r1['err'], fmt='o', ms=2.5, color='#3987e5', elinewidth=0.6, label='data, Hann window a = 3H')
a.set_xlabel('δ = ln x − ln(T/2π)'); a.set_ylabel('V(x)'); a.set_title(f'{nm}: comb with another window shape — z_BK(τ=1) {r1["z_BK_tau1"]:+.2f}, z_GUE {r1["z_GUE_tau1"]:+.1f}, Q_BK {r1["Q_BK"]:.2f}', fontsize=10); a.legend(fontsize=8); a.grid(alpha=.3)
b = ax[1]; ks = list(R['R2'].keys()); zg = [R['R2'][k]['z_GUE_tau1'] for k in ks]; zb = [R['R2'][k]['z_BK_tau1'] for k in ks]; y = np.arange(len(ks))
b.barh(y + 0.2, zg, 0.4, color='#d95926', label='z against GUE'); b.barh(y - 0.2, zb, 0.4, color='#199e70', label='z against the prediction')
b.set_yticks(y); b.set_yticklabels([k.replace('podpasma_', 'sub-bands ').replace('jackknife_', 'jackknife ') for k in ks], fontsize=8); b.axvline(10, color='k', lw=.6, ls=':'); b.axvline(2, color='k', lw=.6, ls=':')
b.set_title(f'τ = 1: six error estimates; R3: N_ef {R["R3"]["N_ef"]:.0f}, p_BK {R["R3"]["p_BK"]:.2f}', fontsize=10); b.legend(fontsize=8, loc='center right', bbox_to_anchor=(1, 0.66), framealpha=1); b.grid(alpha=.3, axis='x')
plt.tight_layout(); plt.savefig(f'out/hr7_R_{nm}_en.png', dpi=130); print('ok')
