// Flies each preset with no graphics and checks the result: node test/landings.js
const { planLanding } = require('../src/solver.js');
const { flyHeadless } = require('../src/sim.js');
const base = { gust: 1, windDir: 30, faultAlt: 0, faultLevel: 1, glide: 20, point: 30, seed: 3 };
const cases = [
  ['Routine return', true, { ...base, world: 'earth', alt: 1500, x: 350, z: -120, vx: -25, vy: -90, vz: 10, wind: 3 }],
  ['Engine fault', true, { ...base, world: 'earth', alt: 2400, x: 300, z: -120, vx: -25, vy: -85, vz: 10, wind: 4, faultAlt: 2000, faultLevel: 0.85 }],
  ['Mars lander', true, { ...base, world: 'mars', wind: 8, gust: 2, vehicle: { mass: 1900, dry: 1500, Tmax: 16000, Tmin: 4000, isp: 225, area: 12 }, alt: 1500, x: 300, z: 100, vx: -20, vy: -75, vz: 0 }],
  ['Lunar lander', true, { ...base, world: 'moon', wind: 0, gust: 0, vehicle: { mass: 15000, dry: 7000, Tmax: 45000, Tmin: 4700, isp: 311, area: 40 }, alt: 2000, x: 800, z: -200, vx: -30, vy: -40, vz: 5, glide: 10, point: 45 }],
  ['Too late to stop', false, { ...base, world: 'earth', alt: 450, x: 60, z: 0, vx: -5, vy: -110, vz: 0, wind: 2 }],
];
let bad = 0;
for (const [name, expect, cfg] of cases) {
  const r = flyHeadless(cfg, planLanding), pass = r.ok === expect; if (!pass) bad++;
  console.log(`${pass ? 'PASS' : 'FAIL'}  ${name.padEnd(18)} ${r.text}, ${r.vs.toFixed(2)} m/s, ${r.miss.toFixed(1)} m from centre`);
}
process.exit(bad ? 1 : 0);
