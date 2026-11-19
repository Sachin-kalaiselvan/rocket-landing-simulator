# Rocket Landing Simulator

Fuel-optimal powered-descent guidance, solved live in the browser. A booster (or a hopper, Mars lander or lunar lander) comes in high and fast; the planner finds the thrust profile that lands it on the pad using the least propellant, then re-solves about once a second as wind, gusts, shoves and engine faults push it off course.

## What is in it

- **Convex guidance (lossless convexification, G-FOLD style).** The non-convex minimum-throttle limit is relaxed with a slack variable; the relaxed problem is convex and its optimum sits on the real limits. Constraints: exact dynamics with gravity, thrust bounds, engine tilt, glide slope.
- **A conic solver written from scratch.** No optimisation library. ADMM over 120 thrust variables with cached factorisation and adaptive step size. Flight time is found by golden-section search.
- **Closed-loop flight.** 120 Hz point-mass model with drag, Ornstein–Uhlenbeck gusts and thrust faults. Re-plans keep the touchdown time if they can and slip it otherwise. If no safe landing exists, it says so and brakes as hard as it can.
- **Your own vehicle.** Enter mass, thrust range, Isp and area; land on Earth, Mars or the Moon.
- **Reliability testing.** Flies a scenario 25 to 250 times with random start state, wind and engine faults across parallel Web Workers. Reports success rate, touchdown scatter, touchdown-speed histogram and failure causes, with replay of any failure in 3D.
- **Exports.** Flight data, the latest plan and reliability results as CSV.

## Run it

Open `index.html` in a browser. Nothing to install.

```
python3 build.py        # rebuild index.html from src/
node test/landings.js   # fly the presets headless and check the results
```

`src/solver.js` is the optimiser, `src/sim.js` the flight model and guidance loop, `src/index.src.html` the interface (three.js).

## Simplifications

Point mass with ideal attitude control, one engine, flat ground, simple drag. Mass loss is handled by successive approximation inside each solve.
