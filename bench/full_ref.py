"""Reference: the exact lossless-convexification G-FOLD problem (log-mass variable,
Acikmese & Ploen 2007 / Blackmore et al. 2010), solved with CVXPY + Clarabel and a
bounded scalar search over flight time. Compares fuel and flight time with the
in-house planner.  python3 bench/full_ref.py"""
import json, subprocess, time, numpy as np, cvxpy as cp
from scipy.optimize import minimize_scalar
G0 = 9.80665

def gfold(c, tf, N=30):
    dt = tf / N; g = np.array([0, -c['g'], 0]); alpha = 1 / (c['isp'] * G0)
    r, v, u, z, s = cp.Variable((N + 1, 3)), cp.Variable((N + 1, 3)), cp.Variable((N, 3)), cp.Variable(N + 1), cp.Variable(N)
    cons = [r[0] == c['r0'], v[0] == c['v0'], z[0] == np.log(c['mass']), r[N] == 0, v[N] == 0, z[N] >= np.log(c['dry'])]
    tg, tp = np.tan(np.radians(c['glideDeg'])), np.tan(np.radians(c['pointDeg']))
    for k in range(N):
        cons += [r[k + 1] == r[k] + v[k] * dt + 0.5 * (u[k] + g) * dt ** 2, v[k + 1] == v[k] + (u[k] + g) * dt, z[k + 1] == z[k] - alpha * s[k] * dt,
                 cp.norm(u[k]) <= s[k], cp.norm(cp.hstack([u[k, 0], u[k, 2]])) <= tp * u[k, 1]]
        z0 = np.log(max(c['dry'], c['mass'] - alpha * c['Tmax'] * k * dt)); mu1, mu2 = c['Tmin'] * np.exp(-z0), c['Tmax'] * np.exp(-z0); dz = z[k] - z0
        cons += [s[k] >= mu1 * (1 - dz + cp.square(dz) / 2), s[k] <= mu2 * (1 - dz), z[k] >= z0 if c['mass'] - alpha * c['Tmax'] * k * dt > c['dry'] else z[k] >= np.log(c['dry'])]
        if k > 0: cons.append(cp.norm(cp.hstack([r[k, 0], r[k, 2]])) <= r[k, 1] / tg)
    prob = cp.Problem(cp.Maximize(z[N]), cons)
    try: prob.solve(solver='CLARABEL')
    except Exception: return None
    if prob.status not in ('optimal', 'optimal_inaccurate'): return None
    return c['mass'] - np.exp(z.value[N])

def search(c):
    est = 2 * c['r0'][1] / max(-c['v0'][1], 20); lo, hi = max(3, est * 0.5), est * 2.2
    grid = np.linspace(lo, hi, 25); vals = [(t, gfold(c, t)) for t in grid]; ok = [(t, f) for t, f in vals if f is not None]
    if not ok: return None, None
    t0 = min(ok, key=lambda x: x[1])[0]; step = grid[1] - grid[0]
    res = minimize_scalar(lambda t: gfold(c, t) or 1e9, bounds=(max(lo, t0 - step), min(hi, t0 + step)), method='bounded', options={'xatol': 0.02})
    return res.x, res.fun

cases = [
    # Blackmore, Acikmese & Scharf 2010, Mars lander, case 1 (altitude axis is y here). No pointing limit in the paper: 89 deg used.
    dict(name='Mars lander (Blackmore et al. 2010, case 1)', r0=[500, 1500, 2000], v0=[0, -75, 100], mass=1905, dry=1505, Tmax=13260, Tmin=4972, isp=225, g=3.7114, glideDeg=4, pointDeg=89),
    dict(name='Mars lander, 30 deg glide slope', r0=[500, 1500, 2000], v0=[0, -75, 100], mass=1905, dry=1505, Tmax=13260, Tmin=4972, isp=225, g=3.7114, glideDeg=30, pointDeg=89),
]
rng = np.random.default_rng(7)
for i in range(20):
    cases.append(dict(name=f'Booster {i + 1}', r0=[float(rng.uniform(-500, 500)), float(rng.uniform(1200, 2500)), float(rng.uniform(-300, 300))],
                      v0=[float(rng.uniform(-30, 30)), float(rng.uniform(-100, -70)), float(rng.uniform(-15, 15))], mass=30000, dry=24000,
                      Tmax=450e3 * 0.86, Tmin=180e3 * 1.08, isp=282, g=9.81, glideDeg=20, pointDeg=30))
json.dump(cases, open('bench/out/cases.json', 'w'))
subprocess.run(['node', 'bench/full_js.js', 'bench/out/cases.json', 'bench/out/full_js.json'], check=True)
js = {r['name']: r for r in json.load(open('bench/out/full_js.json'))}
rows = []
for c in cases:
    t = time.perf_counter(); tf, fuel = search(c); wall = time.perf_counter() - t
    j = js[c['name']]
    rows.append(dict(name=c['name'], ref_tf=tf, ref_fuel=fuel, js_tf=j.get('tf'), js_fuel=j.get('fuelKg'), js_feasible=j['feasible'], js_thrust_violation=j.get('thrustBoundViolation'), js_ms=j['ms'], ref_s=wall))
    r = rows[-1]
    print(f"{c['name'][:44]:44s} ref tf {tf if tf is None else round(tf, 2)} fuel {fuel if fuel is None else round(fuel, 1)} | ours tf {r['js_tf'] and round(r['js_tf'], 2)} fuel {r['js_fuel'] and round(r['js_fuel'], 1)} viol {r['js_thrust_violation']}")
json.dump(rows, open('bench/out/full_cmp.json', 'w'), indent=1)
