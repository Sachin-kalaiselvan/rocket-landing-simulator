// Closed-loop Monte Carlo: full flights (120 Hz physics, re-planning every 1.5 s) with
// dispersed start state, wind, gusts and engine faults. Same dispersions as the app's
// Reliability view.  node bench/monte_carlo.js <first> <count> <out.json>
const fs = require('fs');
const { planLanding } = require('../src/solver.js');
const { flyHeadless } = require('../src/sim.js');
const first = +process.argv[2], count = +process.argv[3], out = process.argv[4];
const base = { world: 'earth', alt: 1500, x: 350, z: -120, vx: -25, vy: -90, vz: 10, glide: 20, point: 30 };
const res = [];
for (let i = first; i < first + count; i++) {
  let st = 1000 + i * 7919; const rnd = () => (st = st * 16807 % 2147483647) / 2147483647, pm = k => 1 + k * (2 * rnd() - 1);
  const alt = base.alt * (1 + 0.15 * rnd()), faulty = rnd() < 0.2;
  const cfg = { ...base, seed: i + 1, alt, x: base.x * pm(0.4), z: base.z * pm(0.4), vx: base.vx * pm(0.3), vz: base.vz * pm(0.3), vy: base.vy * (1.05 - 0.15 * rnd()),
    wind: 12 * rnd(), gust: 3 * rnd(), windDir: 360 * rnd(), faultAlt: faulty ? alt * (0.35 + 0.4 * rnd()) : 0, faultLevel: faulty ? 0.85 + 0.1 * rnd() : 1 };
  const t0 = Date.now(), r = flyHeadless(cfg, planLanding);
  res.push({ i, cfg, x: r.x, z: r.z, ok: r.ok, text: r.text, vs: r.vs, hs: r.hs, miss: r.miss, fuel: r.fuel, t: r.t, replans: r.replans, planned: r.planned, fault: !!r.fault, wall_s: (Date.now() - t0) / 1000 });
  fs.writeFileSync(out, JSON.stringify(res));
}
