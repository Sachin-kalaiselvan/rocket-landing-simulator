// Runs the in-house planner (free final time, successive mass approximation) on the
// cases in bench/cases.json and writes fuel, flight time and thrust-bound checks.
const fs = require('fs');
const { planLanding } = require('../src/solver.js');
const cases = JSON.parse(fs.readFileSync(process.argv[2], 'utf8')), out = [];
for (const c of cases) {
  const t0 = process.hrtime.bigint();
  const p = planLanding({ r0: c.r0, v0: c.v0, mass: c.mass, dry: c.dry, Tmax: c.Tmax, Tmin: c.Tmin, isp: c.isp, g: c.g, glideDeg: c.glideDeg, pointDeg: c.pointDeg, N: c.N || 30 });
  const ms = Number(process.hrtime.bigint() - t0) / 1e6;
  if (!p.feasible) { out.push({ name: c.name, feasible: false, ms }); continue; }
  // fly the plan's thrust-acceleration profile with true mass flow and check the engine limits
  let m = c.mass, worst = 0; const g0 = 9.80665;
  for (let k = 0; k < p.N; k++) {
    const T = m * p.sigma[k]; worst = Math.max(worst, (T - c.Tmax) / c.Tmax, (c.Tmin - T) / c.Tmin);
    m *= Math.exp(-p.sigma[k] * p.dt / (c.isp * g0));
  }
  out.push({ name: c.name, feasible: true, tf: p.tf, fuelKg: c.mass - m, thrustBoundViolation: Math.max(0, worst), ms, evals: p.evals.length });
}
fs.writeFileSync(process.argv[3], JSON.stringify(out, null, 1));
