// Successive convexification (SCvx) for free-final-time 6-DOF powered descent,
// after Szmuk & Acikmese (AIAA 2018-0617) and Mao, Szmuk & Acikmese (2016).
// Each iteration: linearise the dynamics about the last trajectory, discretise with a
// first-order hold, solve one second-order cone program (ECOS via Epigraph), repeat until
// the trust region and the virtual control (dynamics slack) both go to zero.
#pragma once
#include "model.hpp"
#include <vector>
#include <string>

namespace r6 {
struct Trajectory {
    std::vector<Vx> X; std::vector<Vu> U; double sigma = 0;   // sigma = flight time
};
struct IterLog { int it; double defect, trust, nu, cost, ms; };
struct Result {
    Trajectory traj; bool converged = false; int iterations = 0;
    double solve_ms = 0;                                       // total wall time
    std::vector<IterLog> log;
};
struct Options {
    int K = 30, max_iter = 25;
    double w_nu = 1e5, w_trust = 1e-1, w_sigma = 1e-1, w_miss = 50.0;   // w_miss: propellant traded per unit of landing miss
    double tol_trust = 1e-3, tol_nu = 1e-6;
};

Trajectory initial_guess(const Params &p, int K);
// Discretise the linearised dynamics about `ref` (first-order hold, RK4 with 10 sub-steps per interval).
struct Discrete { std::vector<Mxx> A; std::vector<Mxu> Bm, Bp; std::vector<Vx> S, z; double defect = 0; };
Discrete discretise(const Params &p, const Trajectory &ref);
Vx propagate(const Params &p, const Vx &x, const Vu &u0, const Vu &u1, double dt, int steps = 10);
Result solve(const Params &p, const Options &o, const Trajectory *warm = nullptr);
}  // namespace r6
