// Powered-descent guidance: lossless-convexified fuel-optimal landing (G-FOLD style),
// solved with an operator-splitting (ADMM) conic solver written from scratch.
// Units inside the solver are scaled so gravity = 1.
function cholesky(A, n) {
  const L = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      let s = A[i * n + j];
      for (let k = 0; k < j; k++) s -= L[i * n + k] * L[j * n + k];
      if (i === j) { if (s <= 0) s = 1e-12; L[i * n + i] = Math.sqrt(s); }
      else L[i * n + j] = s / L[j * n + j];
    }
  }
  return L;
}
function cholSolve(L, n, b, out) {
  for (let i = 0; i < n; i++) { let s = b[i]; for (let k = 0; k < i; k++) s -= L[i * n + k] * out[k]; out[i] = s / L[i * n + i]; }
  for (let i = n - 1; i >= 0; i--) { let s = out[i]; for (let k = i + 1; k < n; k++) s -= L[k * n + i] * out[k]; out[i] = s / L[i * n + i]; }
  return out;
}
// projection onto circular cone {x : |x_perp| <= tan(alpha) * (x . axis)}; axis = last component index ax
function projCone(s, off, dim, ax, tanA) {
  let y = s[off + ax], n2 = 0;
  for (let i = 0; i < dim; i++) if (i !== ax) n2 += s[off + i] * s[off + i];
  const n = Math.sqrt(n2);
  if (n <= tanA * y) return;
  if (tanA * n <= -y || (n < 1e-14 && y <= 0)) { for (let i = 0; i < dim; i++) s[off + i] = 0; return; }
  const ca = 1 / Math.sqrt(1 + tanA * tanA), sa = tanA * ca;
  const t = n * sa + y * ca;
  for (let i = 0; i < dim; i++) if (i !== ax) s[off + i] = t * sa * s[off + i] / n;
  s[off + ax] = t * ca;
}

// p: {r0:[3], v0:[3], tf, N, amin:[N], amax:[N], tanGlide, tanPoint, warm}
// all in scaled units (g = 1 pointing -y). Returns plan + diagnostics.
function solveLanding(p) {
  const N = p.N, dt = p.tf / N, n = 4 * N;
  const g = [0, -1, 0];
  // rows: blocks {type, rows:[{coef:Map idx->val, h}], ...}
  const rows = []; // each row: {idx:Int32Array, val:Float64Array, h}
  const blocks = [];
  const addRow = (coefs, h) => { rows.push({ idx: Int32Array.from(coefs.map(c => c[0])), val: Float64Array.from(coefs.map(c => c[1])), h }); return rows.length - 1; };
  const posRow = (k, d) => { // r_k component d as affine in a
    const tk = k * dt, coefs = [];
    for (let j = 0; j < k; j++) coefs.push([3 * j + d, dt * dt * (k - j - 0.5)]);
    return addRow(coefs, p.r0[d] + p.v0[d] * tk + 0.5 * g[d] * tk * tk);
  };
  const velRow = (k, d) => {
    const tk = k * dt, coefs = [];
    for (let j = 0; j < k; j++) coefs.push([3 * j + d, dt]);
    return addRow(coefs, p.v0[d] + g[d] * tk);
  };
  // terminal equalities
  blocks.push({ type: 'zero', start: rows.length, dim: 3, rho: 400 }); for (let d = 0; d < 3; d++) posRow(N, d);
  blocks.push({ type: 'zero', start: rows.length, dim: 3, rho: 400 }); for (let d = 0; d < 3; d++) velRow(N, d);
  // glide slope cone on r_k (axis y): |(x,z)| <= y / tanGlide
  for (let k = 1; k < N; k++) { blocks.push({ type: 'cone', start: rows.length, dim: 3, ax: 1, tanA: 1 / p.tanGlide, rho: 1 }); for (let d = 0; d < 3; d++) posRow(k, d); }
  for (let k = 0; k < N; k++) {
    // thrust pointing cone on a_k
    blocks.push({ type: 'cone', start: rows.length, dim: 3, ax: 1, tanA: p.tanPoint, rho: 1 });
    for (let d = 0; d < 3; d++) addRow([[3 * k + d, 1]], 0);
    // |a_k| <= sigma_k  (4-D second-order cone, axis = sigma)
    blocks.push({ type: 'cone', start: rows.length, dim: 4, ax: 3, tanA: 1, rho: 1 });
    for (let d = 0; d < 3; d++) addRow([[3 * k + d, 1]], 0);
    addRow([[3 * N + k, 1]], 0);
    // amin <= sigma_k <= amax
    blocks.push({ type: 'box', start: rows.length, dim: 1, lo: p.amin[k], hi: p.amax[k], rho: 1 });
    addRow([[3 * N + k, 1]], 0);
  }
  const m = rows.length;
  let rho0 = p.rho || 0.1;
  const sig = 1e-6, alpha = 1.6, eps = 1e-7;
  const R = new Float64Array(m);
  const K = new Float64Array(n * n);
  let L;
  const factor = () => {
    for (const b of blocks) for (let i = 0; i < b.dim; i++) R[b.start + i] = rho0 * b.rho;
    K.fill(0);
    for (let i = 0; i < n; i++) K[i * n + i] = eps + sig;
    for (let r = 0; r < m; r++) {
      const { idx, val } = rows[r], w = R[r];
      for (let a = 0; a < idx.length; a++) for (let b = 0; b < idx.length; b++) K[idx[a] * n + idx[b]] += w * val[a] * val[b];
    }
    L = cholesky(K, n);
  };
  factor();
  const c = new Float64Array(n); for (let k = 0; k < N; k++) c[3 * N + k] = dt;
  let cInf = 0; for (let i = 0; i < n; i++) cInf = Math.max(cInf, Math.abs(c[i]));
  let x = new Float64Array(n), z = new Float64Array(m), y = new Float64Array(m);
  if (p.warm && p.warm.x && p.warm.x.length === n) { x.set(p.warm.x); }
  const Fx = (xv, out) => { for (let r = 0; r < m; r++) { const { idx, val, h } = rows[r]; let s = h; for (let a = 0; a < idx.length; a++) s += val[a] * xv[idx[a]]; out[r] = s; } return out; };
  const FTv = (vv, out) => { out.fill(0); for (let r = 0; r < m; r++) { const { idx, val } = rows[r]; for (let a = 0; a < idx.length; a++) out[idx[a]] += val[a] * vv[r]; } return out; };
  const project = (s) => {
    for (const b of blocks) {
      if (b.type === 'zero') for (let i = 0; i < b.dim; i++) s[b.start + i] = 0;
      else if (b.type === 'box') s[b.start] = Math.min(b.hi, Math.max(b.lo, s[b.start]));
      else projCone(s, b.start, b.dim, b.ax, b.tanA);
    }
  };
  const st = new Float64Array(m); Fx(x, st); z.set(st); project(z);
  const rhs = new Float64Array(n), sv = new Float64Array(m), zprev = new Float64Array(m), tmp = new Float64Array(m);
  const dz = new Float64Array(m), ftv = new Float64Array(n), fty = new Float64Array(n);
  const hist = []; let it = 0, rp = Infinity, rd = Infinity, refactors = 0;
  const maxIt = p.maxIt || 6000;
  for (it = 1; it <= maxIt; it++) {
    for (let i = 0; i < n; i++) rhs[i] = sig * x[i] - c[i];
    for (let r = 0; r < m; r++) {
      const { idx, val, h } = rows[r]; const w = R[r] * (z[r] - h) - y[r];
      for (let a = 0; a < idx.length; a++) rhs[idx[a]] += val[a] * w;
    }
    cholSolve(L, n, rhs, x);
    Fx(x, sv);
    zprev.set(z);
    for (let r = 0; r < m; r++) tmp[r] = alpha * sv[r] + (1 - alpha) * zprev[r] + y[r] / R[r];
    project(tmp);
    z.set(tmp);
    for (let r = 0; r < m; r++) y[r] += R[r] * (alpha * sv[r] + (1 - alpha) * zprev[r] - z[r]);
    if (it % 10 === 0) {
      rp = 0; let sInf = 0, zInf = 0;
      for (let r = 0; r < m; r++) { rp = Math.max(rp, Math.abs(sv[r] - z[r])); sInf = Math.max(sInf, Math.abs(sv[r])); zInf = Math.max(zInf, Math.abs(z[r])); }
      for (let r = 0; r < m; r++) dz[r] = R[r] * (z[r] - zprev[r]);
      FTv(dz, ftv); rd = 0; for (let i = 0; i < n; i++) rd = Math.max(rd, Math.abs(ftv[i]));
      if (it % 20 === 0) hist.push([it, rp, rd]);
      if (rp < 2e-5 && rd < 2e-5) break;
      if (it % 50 === 0) {
        FTv(y, fty); let yInf = 0; for (let i = 0; i < n; i++) yInf = Math.max(yInf, Math.abs(fty[i]));
        const rpn = rp / Math.max(1e-9, sInf, zInf), rdn = rd / Math.max(1e-9, yInf, cInf);
        const ratio = Math.sqrt(rpn / Math.max(1e-12, rdn));
        if (ratio > 4 || ratio < 0.25) { rho0 = Math.min(1e5, Math.max(1e-5, rho0 * ratio)); factor(); refactors++; }
      }
    }
  }
  // unpack
  const a = [], sigma = [], r = [], v = [];
  for (let k = 0; k < N; k++) { a.push([x[3 * k], x[3 * k + 1], x[3 * k + 2]]); sigma.push(x[3 * N + k]); }
  let rk = p.r0.slice(), vk = p.v0.slice();
  r.push(rk.slice()); v.push(vk.slice());
  for (let k = 0; k < N; k++) {
    for (let d = 0; d < 3; d++) { rk[d] += vk[d] * dt + 0.5 * (a[k][d] + g[d]) * dt * dt; vk[d] += (a[k][d] + g[d]) * dt; }
    r.push(rk.slice()); v.push(vk.slice());
  }
  let fuel = 0; for (let k = 0; k < N; k++) fuel += sigma[k] * dt;
  const term = Math.hypot(...r[N]) + Math.hypot(...v[N]);
  const feasible = rp < 5e-4 && term < 5e-3;
  return { tf: p.tf, N, dt, a, sigma, r, v, fuel, iters: Math.min(it, maxIt), primal: rp, dual: rd, hist, feasible, x: Array.from(x), terminalError: term, refactors, rho: rho0 };
}

// Scaled wrapper with free-final-time golden search.
// Real units: r0,v0 m, m/s; vehicle {mass, Tmax, Tmin, isp}; limits {glideDeg, pointDeg}
function planLanding(req) {
  const L = req.L || Math.max(400, Math.hypot(...req.r0));
  const G = req.g || 9.81, T = Math.sqrt(L / G), V = L / T, A = L / (T * T);
  const N = req.N || 30;
  const r0 = req.r0.map(q => q / L), v0 = req.v0.map(q => q / V);
  const tanGlide = Math.tan(req.glideDeg * Math.PI / 180), tanPoint = Math.tan(req.pointDeg * Math.PI / 180);
  const g0 = 9.80665;
  const bounds = (tf, prev) => {
    // mass profile from previous solution (successive approximation), else linear guess
    const amin = [], amax = []; let mk = req.mass; const dt = tf / N * T;
    for (let k = 0; k < N; k++) {
      amin.push(req.Tmin / mk / A); amax.push(req.Tmax / mk / A);
      const acc = prev ? prev.sigma[k] * A : 0.6 * req.Tmax / req.mass;
      mk = Math.max(req.dry, mk - dt * mk * acc / (req.isp * g0));
    }
    return { amin, amax };
  };
  const solveAt = (tfReal, warm) => {
    const tf = tfReal / T;
    let b = bounds(tf, null), sol = null;
    for (let pass = 0; pass < 2; pass++) {
      sol = solveLanding({ r0, v0, tf, N, amin: b.amin, amax: b.amax, tanGlide, tanPoint, warm: warm || (sol && { x: sol.x }), maxIt: req.maxIt || 2500 });
      if (!sol.feasible) break;
      b = bounds(tf, sol);
    }
    return sol;
  };
  const evals = [];
  let best = null;
  const cost = (tfReal) => {
    const s = solveAt(tfReal, best && best.sol && { x: best.sol.x });
    s._tf = tfReal;
    const score = s.feasible ? s.fuel : 1e3 + s.primal * 1e3 + s.terminalError * 10;
    evals.push({ tf: tfReal, fuel: s.feasible ? req.mass * (1 - Math.exp(-s.fuel * V / (req.isp * g0))) : null, feasible: s.feasible, iters: s.iters });
    if (req.onEval) req.onEval(evals[evals.length - 1]);
    if (s.feasible && (!best || s.fuel < best.sol.fuel)) best = { tf: tfReal, sol: s };
    return score;
  };
  let tfFixed = req.tfFixed;
  if (tfFixed) { cost(tfFixed); }
  else {
    // rough time-to-land estimate brackets the search
    const h = Math.max(1, req.r0[1]), vy = -req.v0[1];
    let est = 2 * h / Math.max(vy, 20);
    est = Math.max(est, Math.sqrt(Math.hypot(req.r0[0], req.r0[2])) * 1.5, 4);
    let lo = req.tfLo || Math.max(1.5, est * 0.5), hi = req.tfHi || est * 2.0;
    // coarse scan then golden-section refine
    const scan = [];
    for (let i = 0; i <= 6; i++) { const t = lo + (hi - lo) * i / 6; scan.push([t, cost(t)]); }
    let bi = 0; for (let i = 1; i < scan.length; i++) if (scan[i][1] < scan[bi][1]) bi = i;
    let a = scan[Math.max(0, bi - 1)][0], b = scan[Math.min(6, bi + 1)][0];
    const gr = 0.6180339887;
    let c1 = b - gr * (b - a), c2 = a + gr * (b - a), f1 = cost(c1), f2 = cost(c2);
    for (let i = 0; i < 9; i++) {
      if (f1 < f2) { b = c2; c2 = c1; f2 = f1; c1 = b - gr * (b - a); f1 = cost(c1); }
      else { a = c1; c1 = c2; f1 = f2; c2 = a + gr * (b - a); f2 = cost(c2); }
    }
  }
  if (!best) {
    const last = evals[evals.length - 1];
    return { feasible: false, evals };
  }
  const s = best.sol;
  const plan = {
    feasible: true, tf: best.tf, N, dt: s.dt * T,
    r: s.r.map(q => q.map(e => e * L)), v: s.v.map(q => q.map(e => e * V)),
    a: s.a.map(q => q.map(e => e * A)), sigma: s.sigma.map(e => e * A),
    iters: s.iters, primal: s.primal, hist: s.hist, evals,
    fuelKg: req.mass * (1 - Math.exp(-s.fuel * V / (req.isp * g0))), dv: s.fuel * V, warm: s.x, L,
  };
  return plan;
}
if (typeof module !== 'undefined') module.exports = { solveLanding, planLanding };
