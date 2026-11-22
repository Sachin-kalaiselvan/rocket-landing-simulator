# rocket-landing-simulator

fuel-optimal powered descent guidance, solved live in your browser tab. a booster (or a hopper, mars lander, lunar lander) comes in high and fast, the planner works out the thrust profile that puts it on the pad using the least propellant, then re-solves about once a second while wind, gusts, shoves and engine faults try to knock it off.

## what's in here

- **convex guidance, G-FOLD style.** the min-throttle limit isn't convex, so it gets relaxed with a slack variable (lossless convexification). relaxed problem is convex and the optimum still sits on the real limits. dynamics w/ gravity, thrust bounds, engine tilt, glide slope all in there
- **my own conic solver.** zero optimisation libraries. ADMM over 120 thrust variables, cached factorisation, adaptive step. flight time found by golden-section search
- **closed loop flight.** 120 Hz point mass w/ drag, ornstein-uhlenbeck gusts, thrust faults. re-plans keep the touchdown time if they can, slip it if they can't. if there's no safe landing left it tells you and brakes as hard as it can
- **bring your own vehicle.** mass, thrust range, isp, area. land on earth, mars or the moon
- **reliability runs.** fly a scenario 25 to 250 times with random start, wind and faults across parallel web workers. success rate, touchdown scatter, touchdown speed histogram, why each failure happened, and you can replay any failure in 3D
- **csv export** for flight data, latest plan, reliability results

## 6-DOF in c++

`cpp/` has full rigid-body landing guidance: attitude, body rates, gimbal limits. successive convexification (szmuk & açıkmeşe 2018) on top of ECOS.

- matches a separate python implementation to 0.0001% of propellant
- solves in 0.29s
- closed loop with model errors + re-planning: lands upright on target 91/100 runs

method, checks and the failures are in [docs/SIXDOF.md](docs/SIXDOF.md).

## does it actually work

every layer got checked against something that shares zero code with it. full writeup in [docs/VALIDATION.md](docs/VALIDATION.md).

- **solver vs clarabel + ECOS**, 200 problems: fuel agrees within 0.0025%. never calls an infeasible problem feasible
- **planner vs G-FOLD** (açıkmeşe & ploen 2007, blackmore et al. 2010) on the published mars lander case + 19 boosters: within 0.14% of optimal propellant. mars case lands in 79.2s on 362.5 kg vs 79.4s / 362.5 kg for the reference. this check also caught a mass model bug that was burning up to 3.6% extra fuel, fixed now
- **400 monte carlo flights** w/ wind, gusts, faults: 96.3% on the pad, and every single no-fault flight made it. faults leaving ≥87% thrust: all survived. deeper than that eats past the planner's 14% thrust reserve

## run it

open `index.html`. that's it.

```bash
python3 build.py                   # rebuild index.html from src/
node test/landings.js              # fly the presets headless + check them
python3 test/reference_check.py    # solver vs clarabel & ECOS (pip install cvxpy clarabel ecos)
```

`src/solver.js` = optimiser, `src/sim.js` = flight model + guidance loop, `src/index.src.html` = the UI (three.js).

## what's simplified

point mass with ideal attitude control, one engine, flat ground, basic drag. mass loss handled by successive approximation inside each solve.

MIT
