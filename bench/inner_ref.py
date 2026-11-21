"""Solves the same fixed-time landing problems with CVXPY (Clarabel and ECOS)
and compares them with the in-house ADMM solver.
python3 bench/inner_ref.py bench/out/inner.json bench/out/inner_cmp.json"""
import json, sys, time
import numpy as np, cvxpy as cp

def build(p):
    N, tf = p['N'], p['tf']; dt = tf / N
    r0, v0, g = np.array(p['r0']), np.array(p['v0']), np.array([0, -1.0, 0])
    a, s = cp.Variable((N, 3)), cp.Variable(N)
    def pos(k):
        tk = k * dt
        e = r0 + v0 * tk + 0.5 * g * tk * tk
        return e + sum(dt * dt * (k - j - 0.5) * a[j] for j in range(k)) if k else e
    def vel(k):
        return v0 + g * k * dt + (sum(dt * a[j] for j in range(k)) if k else 0)
    c = [pos(N) == 0, vel(N) == 0]
    for k in range(1, N):
        r = pos(k); c.append(cp.norm(cp.hstack([r[0], r[2]])) <= r[1] / p['tanGlide'])
    for k in range(N):
        c += [cp.norm(cp.hstack([a[k, 0], a[k, 2]])) <= p['tanPoint'] * a[k, 1], cp.norm(a[k]) <= s[k], s[k] >= p['amin'][k], s[k] <= p['amax'][k]]
    return cp.Problem(cp.Minimize(dt * cp.sum(s)), c), a, s

def check(p, x):
    """Worst constraint violation of a candidate solution x (in-house layout)."""
    N, dt = p['N'], p['tf'] / p['N']; x = np.array(x); a = x[:3 * N].reshape(N, 3); s = x[3 * N:]
    r, v, g, worst = np.array(p['r0'], float), np.array(p['v0'], float), np.array([0, -1.0, 0]), 0
    for k in range(N):
        worst = max(worst, np.linalg.norm(a[k]) - s[k], p['amin'][k] - s[k], s[k] - p['amax'][k], np.hypot(a[k, 0], a[k, 2]) - p['tanPoint'] * a[k, 1])
        r = r + v * dt + 0.5 * (a[k] + g) * dt * dt; v = v + (a[k] + g) * dt
        if k < N - 1: worst = max(worst, np.hypot(r[0], r[2]) - r[1] / p['tanGlide'])
    return max(worst, np.linalg.norm(r), np.linalg.norm(v))

data = json.load(open(sys.argv[1])); out = []
for i, d in enumerate(data):
    p = d['p']; row = {'i': i, 'js': {k: d['js'][k] for k in ('feasible', 'fuel', 'iters', 'ms')}}
    row['js']['violation'] = float(check(p, d['js']['x']))
    for solver in ('CLARABEL', 'ECOS'):
        prob, a, s = build(p)
        t = time.perf_counter(); prob.solve(solver=solver); wall = (time.perf_counter() - t) * 1e3
        st = prob.solver_stats
        row[solver] = {'status': prob.status, 'fuel': float(prob.value) if prob.value is not None and np.isfinite(prob.value) else None,
                       'solve_ms': (st.solve_time or 0) * 1e3, 'wall_ms': wall}
    out.append(row)
json.dump(out, open(sys.argv[2], 'w'), indent=1)
both = [r for r in out if r['js']['feasible'] and r['CLARABEL']['status'] == 'optimal']
gap = np.array([(r['js']['fuel'] - r['CLARABEL']['fuel']) / r['CLARABEL']['fuel'] * 100 for r in both])
print('reference feasible:', sum(r['CLARABEL']['status'] == 'optimal' for r in out), 'of', len(out))
print('in-house feasible :', sum(r['js']['feasible'] for r in out))
print('disagree on feasibility:', [(r['i'], r['js']['feasible'], r['CLARABEL']['status']) for r in out if r['js']['feasible'] != (r['CLARABEL']['status'] == 'optimal')])
print('fuel gap vs Clarabel %%: median %.4f  max |gap| %.4f' % (np.median(gap), np.abs(gap).max()))
print('in-house worst violation: %.2e' % max(r['js']['violation'] for r in both))
print('ECOS vs Clarabel max diff %%: %.5f' % max(abs(r['ECOS']['fuel'] - r['CLARABEL']['fuel']) / r['CLARABEL']['fuel'] * 100 for r in both if r['ECOS']['fuel']))
for k in ('CLARABEL', 'ECOS'): print(k, 'median solve ms %.1f' % np.median([r[k]['solve_ms'] for r in out]))
print('in-house median ms %.1f' % np.median([r['js']['ms'] for r in out]))
