"""Draws the validation charts into docs/. python3 bench/plots.py"""
import json, glob, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 140})
INK, ORANGE, GREY, RED = '#1d2433', '#e8622c', '#9aa3b2', '#c4382b'

# 1. inner solver vs Clarabel
d = json.load(open('bench/out/inner_cmp.json'))
both = [r for r in d if r['js']['feasible'] and r['CLARABEL']['status'] == 'optimal']
gap = np.array([(r['js']['fuel'] - r['CLARABEL']['fuel']) / r['CLARABEL']['fuel'] * 100 for r in both])
fig, ax = plt.subplots(1, 2, figsize=(9, 3.2))
ax[0].hist(gap, bins=30, color=ORANGE); ax[0].set_xlabel('fuel difference vs Clarabel (%)'); ax[0].set_ylabel('problems'); ax[0].set_title(f'Same problem, {len(both)} cases', loc='left')
ms = [r['js']['ms'] for r in d]; cl = [r['CLARABEL']['solve_ms'] for r in d]; ec = [r['ECOS']['solve_ms'] for r in d]
ax[1].boxplot([ms, cl, ec], tick_labels=['in-house ADMM\n(JavaScript)', 'Clarabel\n(Rust)', 'ECOS\n(C)'], showfliers=False)
ax[1].set_yscale('log'); ax[1].set_ylabel('solve time (ms)'); ax[1].set_title('Solve time per problem', loc='left')
fig.tight_layout(); fig.savefig('docs/solver-vs-reference.png'); plt.close(fig)

# 2. full planner vs G-FOLD reference
f = [r for r in json.load(open('bench/out/full_cmp.json')) if r['ref_fuel'] and r['js_fuel']]
fig, ax = plt.subplots(figsize=(9, 3.2))
gap2 = [(r['js_fuel'] - r['ref_fuel']) / r['ref_fuel'] * 100 for r in f]; names = ['Mars' if 'Mars' in r['name'] else r['name'].replace('Booster ', 'B') for r in f]
ax.bar(range(len(f)), gap2, color=[ORANGE if g >= 0 else GREY for g in gap2]); ax.axhline(0, color=INK, lw=0.8)
ax.set_xticks(range(len(f))); ax.set_xticklabels(names, fontsize=8); ax.set_ylabel('extra propellant vs G-FOLD (%)')
ax.set_title('Free-final-time plans against the G-FOLD reference (negative = this planner used less)', loc='left', fontsize=10)
fig.tight_layout(); fig.savefig('docs/planner-vs-gfold.png'); plt.close(fig)

# 3. Monte Carlo
mc = sum((json.load(open(p)) for p in sorted(glob.glob('bench/out/mc_*.json'))), [])
if mc:
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    ok = [r for r in mc if r['ok']]; bad = [r for r in mc if not r['ok']]
    ax[0].hist([r['miss'] for r in ok], bins=np.arange(0, 26, 1), color=ORANGE, label=f'on the pad ({sum(r["miss"] <= 25 for r in ok)})')
    off = [r['miss'] for r in bad]
    ax[0].set_xlabel('distance from pad centre at touchdown (m)'); ax[0].set_ylabel('flights')
    ax[0].set_title(f'Landed flights; {len(bad)} of {len(mc)} failed (beyond 25 m or too hard)', loc='left', fontsize=10)
    ax[1].hist([r['vs'] for r in mc], bins=np.arange(0, max(4, max(r['vs'] for r in mc)) + 0.25, 0.25), color=ORANGE)
    ax[1].axvline(3, color=RED, ls='--', lw=1); ax[1].set_xlabel('touchdown speed (m/s)'); ax[1].set_ylabel('flights'); ax[1].set_title('Touchdown speed (hard landing above 3 m/s)', loc='left', fontsize=10)
    fig.tight_layout(); fig.savefig('docs/monte-carlo.png'); plt.close(fig)
print('charts written')
