# Rocket Landing Simulator

Fuel-optimal powered-descent guidance, solved live in the browser. A booster (or a hopper, Mars lander or lunar lander) comes in high and fast; the planner finds the thrust profile that lands it on the pad using the least propellant, then re-solves about once a second as wind, gusts, shoves and engine faults push it off course.

## What is in it

- **Convex guidance (lossless convexification, G-FOLD style).** The non-convex minimum-throttle limit is relaxed with a slack variable; the relaxed problem is convex and its optimum sits on the real limits. Constraints: exact dynamics with gravity, thrust bounds, engine tilt, glide slope.
- **A conic solver written from scratch.** No optimisation library. ADMM over 120 thrust variables with cached factorisation and adaptive step size. Flight time is found by golden-section search.
- **Closed-loop flight.** 120 Hz point-mass model with drag, Ornstein–Uhlenbeck gusts and thrust faults. Re-plans keep the touchdown time if they can and slip it otherwise. If no safe landing exists, it says so and brakes as hard as it can.
- **Your own vehicle.** Enter mass, thrust range, Isp and area; land on Earth, Mars or the Moon.
- **Reliability testing.** Flies a scenario 25 to 250 times with random start state, wind and engine faults across parallel Web Workers. Reports success rate, touchdown scatter, touchdown-speed histogram and failure causes, with replay of any failure in 3D.
- **Exports.** Flight data, the latest plan and reliability results as CSV.

## 6-DOF guidance in C++

`cpp/` adds full rigid-body landing guidance: attitude, body rates and gimbal limits, solved with successive convexification (Szmuk and Açıkmeşe, 2018) using ECOS.
- It matches an independent Python implementation to 0.0001% of propellant.
- It solves in 0.29 s.
- In closed loop with model errors and re-planning, it lands upright on target in 91 of 100 runs.

Method, checks and failure analysis are in [docs/SIXDOF.md](docs/SIXDOF.md).

## Validation

Each layer was checked against an independent reference. Full method and numbers are in [docs/VALIDATION.md](docs/VALIDATION.md).

- **Solver vs Clarabel and ECOS:** on 200 problems, fuel agrees within 0.0025%. It never calls an infeasible problem feasible.
- **Planner vs G-FOLD** (Açıkmeşe and Ploen 2007; Blackmore et al. 2010), on the published Mars lander case and 19 boosters: within 0.14% of the optimum propellant. The Mars case lands in 79.2 s using 362.5 kg, against 79.4 s and 362.5 kg for the reference. This check found and fixed a mass-model flaw that had been costing up to 3.6% extra propellant.
- **400 Monte Carlo flights** with wind, gusts and engine faults: 96.3% landed on the pad, and every flight without a fault did. Faults that left 87% thrust or more were all survived; deeper faults exceed the planner's 14% thrust reserve.

## Run it

Open `index.html` in a browser. Nothing to install.

```
python3 build.py        # rebuild index.html from src/
node test/landings.js   # fly the presets headless and check the results
python3 test/reference_check.py   # solver vs Clarabel and ECOS (needs: pip install cvxpy clarabel ecos)
```

`src/solver.js` is the optimiser, `src/sim.js` the flight model and guidance loop, `src/index.src.html` the interface (three.js).

## Simplifications

Point mass with ideal attitude control, one engine, flat ground, simple drag. Mass loss is handled by successive approximation inside each solve.
