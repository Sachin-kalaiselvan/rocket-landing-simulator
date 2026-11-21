# 6-DOF landing guidance in C++

`cpp/` plans full rigid-body landings. The vehicle is described by its mass, position, velocity, attitude (a quaternion) and body rates. The engine is gimballed, so steering the thrust is also what turns the vehicle. The method is **successive convexification** (SCvx; Mao, Szmuk and Açıkmeşe, 2016) with free final time, as in Szmuk and Açıkmeşe, *Successive Convexification for 6-DoF Mars Rocket Powered Landing with Free-Final-Time* (AIAA 2018-0617). Each iteration:
1. linearises the nonlinear dynamics about the last trajectory;
2. discretises them exactly with a first-order hold, integrating the state-transition matrices with RK4;
3. solves one second-order cone program with ECOS through the Epigraph modelling library;
4. repeats until the trust region and the dynamics slack ("virtual control") vanish.

The constraints:
- thrust between minimum and maximum (the non-convex minimum is linearised each iteration);
- 20° gimbal;
- 75° tilt;
- 60°/s body rate;
- 20° glide slope;
- dry-mass floor;
- touchdown upright, with no pitch or yaw rate, at a set speed.

The horizontal landing point is a penalty, not a hard constraint. Without that, small disturbances made the problem infeasible (see below).

![6-DOF results](sixdof.png)

## Checks

| Check | Result |
|---|---|
| Independent implementation: Python, CVXPY + Clarabel (a different solver), adaptive RK45 (a different integrator) | propellant within 0.0001%, flight time within 1×10⁻⁵, every state within 4×10⁻⁵ |
| C++ controls flown through the nonlinear model with DOP853 at 10⁻¹⁰ tolerance | lands within 3.4×10⁻⁷ of the planned state |
| Lower bound: the same landing as a point mass, solved exactly by lossless convexification | 6-DOF uses 22% more propellant, which is the cost of having to rotate the vehicle at a limited rate |
| Solve time, nominal landing | C++ 0.29 s for 5 iterations; Python 17 s (59× slower) |

All values are non-dimensional (gravity 1, wet mass 2, dry mass 1).

## Closed loop: 100 landings with model errors

The truth model differs from the planner's:
- start position ±0.4 to ±1;
- start velocity ±0.3 to ±0.5;
- mass error 2%;
- inertia error 10%;
- thrust error 2%;
- engine misalignment 0.005 rad (a torque the planner does not know about);
- a steady side force of 3% of gravity.

Every 0.2 time units the planner re-solves from the measured state, warm-started from its last plan. An attitude tracking loop runs on top. It uses quaternion PD on the planned attitude, applied through the gimbal, plus thrust-magnitude feedback on the position and velocity error.

| Version | Landed upright on target* |
|---|---|
| Plan followed open loop | almost none: errors grow without feedback |
| + attitude tracking loop | 74 |
| + soft landing point, so re-plans stay feasible | 82 |
| + re-planning every 0.2 instead of 0.4, until 0.15 before touchdown | **91** |

\* Success means all of: within 0.05 of target, sink rate under 0.3, sideways speed under 0.15, tilt under 5°.

For landed runs:
- **Miss distance:** median 0.003.
- **Tilt at touchdown:** median 0.7°.

**Why runs fail.** The nine failures come in with large sideways errors late in the descent. Sideways position is corrected only by re-planning, and the 0.15 time units after the last re-plan are flown on the tracking loop alone. A faster solver (see below) would allow re-planning all the way down. So would a lateral position loop that tilts the vehicle, but adding one naively made results worse (63%), because it fought the attitude plan.

**Speed.** The median re-plan takes 0.5 s on this 2-core machine, against a 0.2 time-unit re-plan period. At a 100 m length scale one time unit is about 3.2 s, so the period is about 0.64 s and the solver keeps up on average but not in the worst case (5 s). Flight software uses customised interior-point solvers, for example Dueri et al., JGCD 2017, which solve this class of problem in milliseconds.

## Build and run

```
sudo apt install libeigen3-dev cmake g++
cmake -S cpp -B cpp/build && cmake --build cpp/build -j     # fetches Epigraph (MIT) with ECOS (GPLv3)
cpp/build/rocket6dof > landing.json                         # one landing, with its own nonlinear check
cpp/build/rocket6dof_mc 100 7 > mc.json                     # closed-loop Monte Carlo
python3 bench/scvx6_ref.py                                  # Python reference, nonlinear check, 3-DOF bound
```

**Licence note:** the C++ executables link ECOS, which is GPLv3, so binaries built from `cpp/` are GPLv3. The rest of the repository stays MIT.
