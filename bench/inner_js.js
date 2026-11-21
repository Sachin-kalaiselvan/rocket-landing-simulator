// Generates random fixed-final-time landing problems in the solver's scaled units,
// solves each with the in-house ADMM solver and writes problems + results to JSON.
// node bench/inner_js.js 100 bench/out/inner.json
const fs = require('fs');
const { solveLanding } = require('../src/solver.js');
const n = +(process.argv[2] || 100), out = process.argv[3] || 'bench/out/inner.json';
let seed = 12345; const R = () => (seed = seed * 16807 % 2147483647) / 2147483647;
const res = [];
for (let i = 0; i < n; i++) {
  // physical scenario -> scaled problem exactly like planLanding does (linear mass guess)
  const g0 = 9.80665, G = 9.81, mass = 30000, dry = 24000, Tmax = 450e3 * 0.86, Tmin = 180e3 * 1.08, isp = 282;
  const r0p = [(R() * 2 - 1) * 500, 1200 + R() * 1300, (R() * 2 - 1) * 300], v0p = [(R() * 2 - 1) * 30, -(70 + R() * 30), (R() * 2 - 1) * 15];
  const L = Math.max(400, Math.hypot(...r0p)), T = Math.sqrt(L / G), V = L / T, A = L / (T * T), N = 30;
  const tfReal = 2 * r0p[1] / -v0p[1] * (0.95 + R() * 0.3), tf = tfReal / T, dt = tf / N * T;
  const amin = [], amax = []; let mk = mass;
  for (let k = 0; k < N; k++) { amin.push(Tmin / mk / A); amax.push(Tmax / mk / A); mk = Math.max(dry, mk - dt * 0.6 * Tmax / (isp * g0)); }
  const p = { r0: r0p.map(q => q / L), v0: v0p.map(q => q / V), tf, N, amin, amax, tanGlide: Math.tan(20 * Math.PI / 180), tanPoint: Math.tan(30 * Math.PI / 180), maxIt: +(process.argv[4] || 2500) };
  const t0 = process.hrtime.bigint(); const s = solveLanding(p); const ms = Number(process.hrtime.bigint() - t0) / 1e6;
  res.push({ p, js: { feasible: s.feasible, fuel: s.fuel, x: s.x, iters: s.iters, primal: s.primal, terminalError: s.terminalError, ms } });
}
fs.mkdirSync(require('path').dirname(out), { recursive: true });
fs.writeFileSync(out, JSON.stringify(res));
console.log(`${n} problems, ${res.filter(r => r.js.feasible).length} feasible by the in-house solver, median ${res.map(r => r.js.ms).sort((a, b) => a - b)[n >> 1].toFixed(1)} ms`);
