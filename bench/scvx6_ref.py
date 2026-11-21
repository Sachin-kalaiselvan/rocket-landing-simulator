"""Independent checks of the C++ 6-DOF solver (cpp/build/rocket6dof):
1. a separate SCvx implementation in Python: CVXPY with the Clarabel solver (not ECOS), scipy's
   adaptive RK45 for the discretisation (not fixed-step RK4);
2. the C++ controls flown through the nonlinear model with scipy's DOP853 at 1e-10 tolerance;
3. a lower bound: the 3-DOF (point-mass) fuel-optimal landing with the same thrust, glide-slope and
   pointing limits, solved exactly by lossless convexification. Attitude dynamics can only cost fuel,
   so the 6-DOF answer must use at least this much.
python3 bench/scvx6_ref.py"""
import json, subprocess, time, numpy as np, cvxpy as cp
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar

P = dict(m_wet=2.0, m_dry=1.0, alpha=0.05, Tmin=0.3, Tmax=5.0, gimbal=np.radians(20), tilt=np.radians(75), glide=np.radians(20), wmax=np.radians(60),
         g=np.array([-1.0, 0, 0]), rT=np.array([-0.25, 0, 0]), J=np.diag([0.02] * 3))
x0 = np.array([2, 4, 4, 0, -1, -2, 0, 1, 0, 0, 0, 0, 0, 0.0]); xf = np.array([1, 0, 0, 0, -0.1, 0, 0, 1, 0, 0, 0, 0, 0, 0.0])
NX, NU, K = 14, 3, 30

def dcm(q):
    w, x, y, z = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)], [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)], [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)]])

def f(x, u):
    m, q, w = x[0], x[7:11], x[11:14]
    O = np.array([[0, -w[0], -w[1], -w[2]], [w[0], 0, w[2], -w[1]], [w[1], -w[2], 0, w[0]], [w[2], w[1], -w[0], 0]])
    return np.concatenate([[-P['alpha'] * np.linalg.norm(u)], x[4:7], dcm(q) @ u / m + P['g'], 0.5 * O @ q, np.linalg.solve(P['J'], np.cross(P['rT'], u) - np.cross(w, P['J'] @ w))])

def jac(x, u, h=1e-7):
    A = np.column_stack([(f(x + h * e, u) - f(x - h * e, u)) / (2 * h) for e in np.eye(NX)])
    B = np.column_stack([(f(x, u + h * e) - f(x, u - h * e)) / (2 * h) for e in np.eye(NU)])
    return A, B

def discretise(X, U, sig):
    out = []
    for k in range(K - 1):
        def rhs(t, y):
            s = t * (K - 1); x = y[:NX]; Phi = y[NX:NX + NX * NX].reshape(NX, NX)
            u = (1 - s) * U[k] + s * U[k + 1]; A, B = jac(x, u); fx = f(x, u); Pi = np.linalg.inv(Phi)
            return np.concatenate([sig * fx, (sig * A @ Phi).ravel(), (Pi @ (sig * B) * (1 - s)).ravel(), (Pi @ (sig * B) * s).ravel(), Pi @ fx, Pi @ (-sig * A @ x - sig * B @ u)])
        y0 = np.concatenate([X[k], np.eye(NX).ravel(), np.zeros(2 * NX * NU + 2 * NX)])
        y = solve_ivp(rhs, (0, 1 / (K - 1)), y0, method='RK45', rtol=1e-8, atol=1e-10).y[:, -1]
        Phi = y[NX:NX + NX * NX].reshape(NX, NX); o = NX + NX * NX
        Bm = Phi @ y[o:o + NX * NU].reshape(NX, NU); o += NX * NU; Bp = Phi @ y[o:o + NX * NU].reshape(NX, NU); o += NX * NU
        out.append((Phi, Bm, Bp, Phi @ y[o:o + NX], Phi @ y[o + NX:o + 2 * NX], np.linalg.norm(y[:NX] - X[k + 1])))
    return out

def scvx(tf_guess=5.0, iters=25):
    a = np.linspace(0, 1, K)[:, None]; X = (1 - a) * x0 + a * xf; X[:, 0] = (1 - a[:, 0]) * 2 + a[:, 0]; X[:, 7:] = [1, 0, 0, 0, 0, 0, 0]
    U = np.column_stack([X[:, 0], np.zeros(K), np.zeros(K)]); sig = tf_guess; log = []
    for it in range(1, iters + 1):
        t0 = time.perf_counter(); D = discretise(X, U, sig)
        Xv, Uv, nu, s, tr, ds = cp.Variable((K, NX)), cp.Variable((K, NU)), cp.Variable((K - 1, NX)), cp.Variable(), cp.Variable(K), cp.Variable()
        miss = cp.Variable()
        # same end conditions as the C++ solver: altitude, velocity, upright, no pitch/yaw rate; roll free;
        # horizontal miss penalised (w = 50) rather than forced to zero
        c = [Xv[0] == x0, Xv[K - 1, 1] == xf[1], Xv[K - 1, 4:7] == xf[4:7], Xv[K - 1, 9:11] == 0, Xv[K - 1, 12:14] == 0,
             cp.norm(Xv[K - 1, 2:4] - xf[2:4]) <= miss, cp.abs(s - sig) <= ds, s >= 0.1]
        for k in range(K - 1):
            A, Bm, Bp, S, z, _ = D[k]; c.append(Xv[k + 1] == A @ Xv[k] + Bm @ Uv[k] + Bp @ Uv[k + 1] + S * s + z + nu[k])
        uh = U / np.linalg.norm(U, axis=1, keepdims=True)
        for k in range(K):
            c += [cp.norm(Uv[k]) <= Uv[k, 0] / np.cos(P['gimbal']), cp.norm(Uv[k]) <= P['Tmax'], uh[k] @ Uv[k] >= P['Tmin'], cp.norm(cp.hstack([Xv[k] - X[k], Uv[k] - U[k]])) <= tr[k]]
            if k == 0: continue
            c += [Xv[k, 0] >= P['m_dry'], cp.norm(Xv[k, 2:4]) <= Xv[k, 1] / np.tan(P['glide']), cp.norm(Xv[k, 9:11]) <= np.sqrt((1 - np.cos(P['tilt'])) / 2),
                  cp.norm(Xv[k, 11:14]) <= P['wmax']]
        prob = cp.Problem(cp.Minimize(-Xv[K - 1, 0] + 50 * miss + 1e5 * cp.sum(cp.abs(nu)) + 0.1 * cp.sum(tr) + 0.1 * ds), c); prob.solve(solver='CLARABEL')
        X, U, sig = Xv.value, Uv.value, float(s.value)
        log.append(dict(it=it, defect=max(d[5] for d in D), trust=float(tr.value.max()), nu=float(np.abs(nu.value).sum()), m_final=float(X[-1, 0]), s=time.perf_counter() - t0))
        if tr.value.max() < 1e-3 and np.abs(nu.value).sum() < 1e-6: break
    return X, U, sig, log

def fly(U, sig, x=x0):
    """Nonlinear model, first-order-hold controls, DOP853 at tight tolerance."""
    xs = [x]; dt = sig / (K - 1)
    for k in range(K - 1):
        r = solve_ivp(lambda t, y: f(y, (1 - t / dt) * U[k] + t / dt * U[k + 1]), (0, dt), xs[-1], method='DOP853', rtol=1e-10, atol=1e-12)
        xs.append(r.y[:, -1])
    return np.array(xs)

def three_dof(tf, N=60):
    """Lossless-convexified point-mass problem (log mass), same limits; returns fuel used."""
    dt = tf / N; a = P['alpha']; g = P['g']
    r, v, u, z, s = cp.Variable((N + 1, 3)), cp.Variable((N + 1, 3)), cp.Variable((N, 3)), cp.Variable(N + 1), cp.Variable(N)
    c = [r[0] == x0[1:4], v[0] == x0[4:7], z[0] == np.log(2), r[N] == 0, v[N] == xf[4:7], z[N] >= 0]
    for k in range(N):
        c += [r[k + 1] == r[k] + v[k] * dt + 0.5 * (u[k] + g) * dt ** 2, v[k + 1] == v[k] + (u[k] + g) * dt, z[k + 1] == z[k] - a * s[k] * dt,
              cp.norm(u[k]) <= s[k], cp.norm(u[k, 1:]) <= np.tan(P['tilt']) * u[k, 0]]
        z0 = np.log(max(1.0, 2 - a * P['Tmax'] * k * dt)); mu1, mu2 = P['Tmin'] * np.exp(-z0), P['Tmax'] * np.exp(-z0); dz = z[k] - z0
        c += [s[k] >= mu1 * (1 - dz + cp.square(dz) / 2), s[k] <= mu2 * (1 - dz)]
        if k: c.append(cp.norm(r[k, 1:]) <= r[k, 0] / np.tan(P['glide']))
    p = cp.Problem(cp.Maximize(z[N]), c); p.solve(solver='CLARABEL')
    return 2 - np.exp(z.value[N]) if p.status == 'optimal' else None

if __name__ == '__main__':
    t = time.perf_counter(); cpp = json.loads(subprocess.run(['cpp/build/rocket6dof'], capture_output=True, text=True).stdout); t_cpp_wall = time.perf_counter() - t
    t = time.perf_counter(); X, U, sig, log = scvx(); t_py = time.perf_counter() - t
    Xc, Uc = np.array(cpp['X']), np.array(cpp['U'])
    flown = fly(Uc, cpp['tf'])
    res = minimize_scalar(lambda tf: three_dof(tf) or 1e9, bounds=(2.0, 8.0), method='bounded', options={'xatol': 1e-3})
    out = dict(cpp=dict(tf=cpp['tf'], fuel=2 - cpp['m_final'], iterations=cpp['iterations'], solve_ms=cpp['solve_ms']),
               python=dict(tf=sig, fuel=2 - X[-1, 0], iterations=len(log), solve_s=t_py),
               agreement=dict(fuel_diff_pct=abs((2 - cpp['m_final']) - (2 - X[-1, 0])) / (2 - X[-1, 0]) * 100, tf_diff=abs(cpp['tf'] - sig),
                              max_state_diff=float(np.abs(Xc - X).max())),
               nonlinear_check=dict(final_position_error=float(np.linalg.norm(flown[-1, 1:4] - xf[1:4])), final_velocity_error=float(np.linalg.norm(flown[-1, 4:7] - xf[4:7])),
                                    max_node_error=float(np.abs(flown - Xc).max())),
               three_dof_bound=dict(tf=float(res.x), fuel=float(res.fun)))
    out['attitude_cost_pct'] = ((2 - cpp['m_final']) / res.fun - 1) * 100
    print(json.dumps(out, indent=1))
    json.dump(dict(out, cpp_traj=dict(X=cpp['X'], U=cpp['U'], tf=cpp['tf']), py_traj=dict(X=X.tolist(), U=U.tolist(), tf=sig), py_log=log, cpp_log=cpp['log']), open('bench/out/scvx6.json', 'w'))
