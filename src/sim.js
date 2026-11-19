// Flight dynamics + guidance tracking. Shared by the page and the headless test.
const G0 = 9.80665, GRAV = 9.81;
const VEHICLE = { mass: 30000, dry: 24000, Tmax: 450e3, Tmin: 180e3, isp: 282, area: 3.7 * 40 };
const WORLDS = { earth: { g: 9.81, rho: 1.225 }, mars: { g: 3.71, rho: 0.020 }, moon: { g: 1.62, rho: 0 } };
function makeSim(cfg) {
  return {
    veh: Object.assign({}, cfg.vehicle || VEHICLE), g: (WORLDS[cfg.world] || WORLDS.earth).g, rho: (WORLDS[cfg.world] || WORLDS.earth).rho,
    phase: 'idle', t: 0, r: [cfg.x, cfg.alt, cfg.z], v: [cfg.vx, cfg.vy, cfg.vz], m: (cfg.vehicle || VEHICLE).mass,
    TmaxEff: (cfg.vehicle || VEHICLE).Tmax, T: [0, 0, 0], plan: null, planT0: 0, tEnd: 0, lastReplan: -9, replans: 0,
    pending: false, mode: 'guided', faultDone: false, gust: [0, 0], trace: [], events: [], solveMs: 0, solveIters: 0,
    result: null, rng: cfg.seed || 1,
  };
}
function rand(sim) { sim.rng = (sim.rng * 16807) % 2147483647; return sim.rng / 2147483647; }
function gauss(sim) { return Math.sqrt(-2 * Math.log(rand(sim) + 1e-12)) * Math.cos(2 * Math.PI * rand(sim)); }
function guidanceRequest(sim, cfg, extra) {
  return Object.assign({
    r0: sim.r.slice(), v0: sim.v.slice(), mass: sim.m, dry: sim.veh.dry, Tmax: sim.TmaxEff * (extra && extra.tight ? 0.985 : sim.faultDone ? 0.95 : 0.86), Tmin: sim.veh.Tmin * (extra && extra.tight ? 1.02 : 1.08),
    isp: sim.veh.isp, g: sim.g, glideDeg: cfg.glide, pointDeg: Math.max(5, cfg.point - 4), N: 30,
  }, extra || {}, { tight: undefined });
}
// full search: with the usual thrust reserve first, then using almost all the thrust there is
function searchRequests(sim, cfg) { return [guidanceRequest(sim, cfg), guidanceRequest(sim, cfg, { tight: true })]; }
function searchPlan(sim, cfg, solve) { let p; for (const r of searchRequests(sim, cfg)) { p = solve(r); if (p.feasible) return p; } return p; }
function acceptPlan(sim, plan, tReq) {
  sim.plan = plan; sim.planT0 = tReq; sim.tEnd = tReq + plan.tf; sim.solveMs = plan.ms || 0; sim.solveIters = plan.iters;
}
// commanded thrust-acceleration (m/s^2) at current time
function commandAccel(sim) {
  const p = sim.plan, r = sim.r, v = sim.v;
  const tau = sim.t - sim.planT0;
  const terminal = !p || sim.mode === 'emergency' || tau > p.tf - 0.6 || r[1] < 2.0;
  if (terminal) {
    if (sim.mode === 'emergency') {
      // no safe solution: brake as hard as needed to come down slowly, kill sideways speed
      const vyT = -Math.max(1, Math.min(6, r[1] * 0.08));
      return [-0.4 * v[0], sim.g + 2.0 * (vyT - v[1]), -0.4 * v[2]];
    }
    const vyT = -Math.max(0.7, Math.min(4, r[1] * 0.5));
    return [-0.5 * r[0] - 2.0 * v[0], sim.g + 3.0 * (vyT - v[1]), -0.5 * r[2] - 2.0 * v[2]];
  }
  const k = Math.min(p.N - 1, Math.max(0, Math.floor(tau / p.dt)));
  const s = Math.max(0, tau - k * p.dt), a = p.a[k];
  const out = [0, 0, 0], g = [0, -sim.g, 0];
  for (let d = 0; d < 3; d++) {
    const rr = p.r[k][d] + p.v[k][d] * s + 0.5 * (a[d] + g[d]) * s * s;
    const vr = p.v[k][d] + (a[d] + g[d]) * s;
    out[d] = a[d] + 0.30 * (rr - r[d]) + 1.1 * (vr - v[d]);
  }
  return out;
}
function limitThrust(sim, cfg, acc) {
  if (sim.m <= sim.veh.dry + 1) return [0, 0, 0];
  let T = acc.map(q => q * sim.m);
  if (T[1] < 1) T[1] = 1;
  const h = Math.hypot(T[0], T[2]), maxTan = Math.tan(cfg.point * Math.PI / 180);
  if (h > maxTan * T[1]) { const s = maxTan * T[1] / h; T[0] *= s; T[2] *= s; }
  const mag = Math.hypot(...T), lo = sim.veh.Tmin, hi = sim.TmaxEff;
  const target = Math.min(hi, Math.max(lo, mag));
  return T.map(q => q * target / mag);
}
function step(sim, cfg, dt) {
  if (sim.phase !== 'flight') return;
  // wind: steady + Ornstein-Uhlenbeck gusts
  const tauG = 3;
  for (let i = 0; i < 2; i++) sim.gust[i] += (-sim.gust[i] / tauG) * dt + cfg.gust * Math.sqrt(2 * dt / tauG) * gauss(sim);
  const dir = (cfg.windDir || 0) * Math.PI / 180;
  const w = [cfg.wind * Math.cos(dir) + sim.gust[0], cfg.wind * Math.sin(dir) + sim.gust[1]];
  const rel = [w[0] - sim.v[0], w[1] - sim.v[2]], relM = Math.hypot(...rel);
  const kd = 0.5 * sim.rho * 0.8 * (sim.veh.area || 100) / sim.m;
  sim.wind = w;
  // engine fault
  if (cfg.faultAlt > 0 && !sim.faultDone && sim.r[1] < cfg.faultAlt) {
    sim.faultDone = true; sim.TmaxEff = sim.veh.Tmax * cfg.faultLevel; sim.needSearch = true;
    sim.events.push({ t: sim.t, kind: 'fault', text: `Engine thrust dropped to ${Math.round(cfg.faultLevel * 100)}%` });
  }
  const acc = commandAccel(sim);
  sim.T = limitThrust(sim, cfg, acc);
  const Tm = Math.hypot(...sim.T);
  const a = [sim.T[0] / sim.m + kd * relM * rel[0], sim.T[1] / sim.m - sim.g, sim.T[2] / sim.m + kd * relM * rel[1]];
  for (let d = 0; d < 3; d++) { sim.v[d] += a[d] * dt; sim.r[d] += sim.v[d] * dt; }
  sim.m = Math.max(sim.veh.dry, sim.m - Tm / (sim.veh.isp * G0) * dt);
  sim.t += dt;
  if (sim.r[1] <= 0) touchdown(sim);
}
function touchdown(sim) {
  sim.r[1] = 0;
  const vs = -sim.v[1], hs = Math.hypot(sim.v[0], sim.v[2]), miss = Math.hypot(sim.r[0], sim.r[2]);
  let result, ok = false;
  if (vs > 8 || hs > 4) result = 'Vehicle lost on impact';
  else if (miss > 25) result = 'Landed off the pad';
  else if (vs > 3 || hs > 1.5) result = 'Hard landing';
  else { result = 'Landed on the pad'; ok = true; }
  sim.result = { text: result, ok, vs, hs, miss, fuel: sim.veh.mass - sim.m, t: sim.t };
  sim.phase = ok ? 'landed' : 'crashed';
  sim.T = [0, 0, 0]; sim.v = [0, 0, 0];
}
// Re-plan policy: keep the landing time, else slip it later a little, else full free-time search.
function replanSequence(sim, cfg, solve) {
  const rem = sim.tEnd - sim.t;
  for (const slip of [0, 1.5, 3.5]) {
    const q = solve(guidanceRequest(sim, cfg, { tfFixed: rem + slip }));
    if (q.feasible) return q;
  }
  return null;
}
// One complete flight with no graphics: used by the batch reliability test (runs in a worker).
function flyHeadless(cfg, solve) {
  const sim = makeSim(cfg); const p0 = searchPlan(sim, cfg, solve);
  if (!p0.feasible) sim.mode = 'emergency'; else acceptPlan(sim, p0, 0);
  sim.phase = 'flight';
  while (sim.phase === 'flight' && sim.t < 600) {
    if (sim.needSearch) { sim.needSearch = false; const q = searchPlan(sim, cfg, solve); if (q.feasible) { acceptPlan(sim, q, sim.t); sim.replans++; } else sim.mode = 'emergency'; }
    else if (sim.mode === 'guided' && sim.t - sim.lastReplan >= 1.5 && sim.tEnd - sim.t > 3) {
      sim.lastReplan = sim.t; const q = replanSequence(sim, cfg, solve);
      if (q) { acceptPlan(sim, q, sim.t); sim.replans++; } else if (!sim.lastSearch || sim.t - sim.lastSearch > 6) { sim.lastSearch = sim.t; const q2 = searchPlan(sim, cfg, solve); if (q2.feasible) { acceptPlan(sim, q2, sim.t); sim.replans++; } }
    }
    step(sim, cfg, 1 / 60);
  }
  const r = sim.result || { text: 'Did not land', ok: false, vs: 0, hs: 0, miss: 0, fuel: sim.veh.mass - sim.m, t: sim.t };
  return { ok: r.ok, text: r.text, vs: r.vs, hs: r.hs, miss: r.miss, x: sim.r[0], z: sim.r[2], fuel: r.fuel, t: r.t, replans: sim.replans, planned: !!p0.feasible, fault: sim.faultDone, emergency: sim.mode === 'emergency', wind: cfg.wind };
}
if (typeof module !== 'undefined') module.exports = { searchPlan, searchRequests, flyHeadless, VEHICLE, WORLDS, makeSim, step, guidanceRequest, acceptPlan, commandAccel, replanSequence };
