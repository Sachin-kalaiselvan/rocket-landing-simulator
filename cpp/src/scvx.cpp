#include "scvx.hpp"
#include <epigraph.hpp>
#include <chrono>
#include <memory>
#include <cstdio>
#include <cstdlib>

namespace r6 {
using namespace cvx;

Trajectory initial_guess(const Params &p, int K) {
    Trajectory t; t.sigma = p.tf_guess;
    for (int k = 0; k < K; ++k) {
        const double a = double(k) / (K - 1);
        Vx x = (1 - a) * p.x0 + a * p.xf;
        x(0) = (1 - a) * p.m_wet + a * p.m_dry;
        x.segment<4>(7) = Eigen::Vector4d(1, 0, 0, 0);
        x.segment<3>(11).setZero();
        t.X.push_back(x);
        t.U.push_back(Vu(-p.g(0) * x(0), 0, 0));   // hover thrust along the body x axis
    }
    return t;
}

// Augmented state for the first-order-hold discretisation:
// [x, Phi (NX*NX), Phi^-1 B lambda- , Phi^-1 B lambda+, Phi^-1 f, Phi^-1 z]
namespace {
struct Aug { Vx x; Mxx Phi; Mxu Bm, Bp; Vx S, z; };
Aug deriv(const Params &p, const Aug &a, const Vu &u0, const Vu &u1, double s, double sigma) {
    const Vu u = (1 - s) * u0 + s * u1;
    Mxx A; Mxu B; jacobians(p, a.x, u, A, B);
    const Vx fx = f(p, a.x, u);
    Aug d;
    d.x = sigma * fx;
    d.Phi = sigma * A * a.Phi;
    const Eigen::PartialPivLU<Mxx> lu(a.Phi);
    d.Bm = lu.solve(sigma * B) * (1 - s);
    d.Bp = lu.solve(sigma * B) * s;
    d.S = lu.solve(fx);
    d.z = lu.solve(-sigma * A * a.x - sigma * B * u);
    return d;
}
Aug axpy(const Aug &a, const Aug &d, double h) {
    Aug r{a.x + h * d.x, a.Phi + h * d.Phi, a.Bm + h * d.Bm, a.Bp + h * d.Bp, a.S + h * d.S, a.z + h * d.z};
    return r;
}
}  // namespace

Discrete discretise(const Params &p, const Trajectory &ref) {
    const int K = ref.X.size(); const double dtau = 1.0 / (K - 1); const int n = 10;
    Discrete D;
    for (int k = 0; k < K - 1; ++k) {
        Aug a{ref.X[k], Mxx::Identity(), Mxu::Zero(), Mxu::Zero(), Vx::Zero(), Vx::Zero()};
        const double h = dtau / n;
        for (int i = 0; i < n; ++i) {
            const double s = double(i) / n, sh = s + 0.5 / n, s1 = double(i + 1) / n;
            const Aug k1 = deriv(p, a, ref.U[k], ref.U[k + 1], s, ref.sigma);
            const Aug k2 = deriv(p, axpy(a, k1, h / 2), ref.U[k], ref.U[k + 1], sh, ref.sigma);
            const Aug k3 = deriv(p, axpy(a, k2, h / 2), ref.U[k], ref.U[k + 1], sh, ref.sigma);
            const Aug k4 = deriv(p, axpy(a, k3, h), ref.U[k], ref.U[k + 1], s1, ref.sigma);
            Aug sum = axpy(axpy(axpy(k1, k2, 2.0), k3, 2.0), k4, 1.0);
            a = axpy(a, sum, h / 6);
        }
        D.A.push_back(a.Phi);
        D.Bm.push_back(a.Phi * a.Bm); D.Bp.push_back(a.Phi * a.Bp);
        D.S.push_back(a.Phi * a.S); D.z.push_back(a.Phi * a.z);
        D.defect = std::max(D.defect, (a.x - ref.X[k + 1]).norm());   // nonlinear propagation vs. the reference
    }
    return D;
}

Vx propagate(const Params &p, const Vx &x, const Vu &u0, const Vu &u1, double dt, int steps) {
    Vx s = x; const double h = dt / steps;
    for (int i = 0; i < steps; ++i) {
        const double a = double(i) / steps, b = (i + 0.5) / steps, c = double(i + 1) / steps;
        auto U = [&](double t) -> Vu { return (1 - t) * u0 + t * u1; };
        const Vx k1 = f(p, s, U(a)), k2 = f(p, s + h / 2 * k1, U(b)), k3 = f(p, s + h / 2 * k2, U(b)), k4 = f(p, s + h * k3, U(c));
        s += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    }
    return s;
}

Result solve(const Params &p, const Options &o, const Trajectory *warm) {
    const auto T0 = std::chrono::steady_clock::now();
    const int K = o.K;
    Result R; R.traj = warm ? *warm : initial_guess(p, K);
    Trajectory &ref = R.traj;

    // problem data that changes every iteration (epigraph dynamic parameters)
    std::vector<Mxx> A(K - 1); std::vector<Mxu> Bm(K - 1), Bp(K - 1); std::vector<Vx> S(K - 1), z(K - 1);
    std::vector<Vx> Xr(K); std::vector<Vu> Ur(K), Uhat(K); double sr = ref.sigma;
    for (int k = 0; k < K - 1; ++k) { A[k].setZero(); Bm[k].setZero(); Bp[k].setZero(); S[k].setZero(); z[k].setZero(); }

    OptimizationProblem P;
    MatrixX X = P.addVariable("X", NX, K), U = P.addVariable("U", NU, K), NU_ = P.addVariable("nu", NX, K - 1), NUabs = P.addVariable("nuabs", NX, K - 1);
    VectorX D = P.addVariable("trust", K); Scalar sigma = P.addVariable("sigma"), Ds = P.addVariable("dsigma");

    // boundary conditions: everything fixed at the start; at the end, all but mass
    P.addConstraint(equalTo(X.col(0), par(p.x0)));
    // final: position, velocity, upright (q_y = q_z = 0) and no pitch/yaw rate. Roll angle and roll rate are
    // left free: the gimballed engine cannot make a roll torque (a real vehicle would use roll thrusters).
    // Landing point: altitude and velocity are hard constraints; the horizontal miss is a penalty instead,
    // so a perturbed vehicle that can no longer reach the pad exactly still gets a safe, upright landing as
    // close to it as possible (the minimum-landing-error idea of Blackmore, Acikmese & Scharf, 2010).
    P.addConstraint(equalTo(X(1, K - 1), p.xf(1)));
    P.addConstraint(equalTo(X.col(K - 1).segment(4, 3), par(Vx(p.xf).segment(4, 3))));
    Scalar miss = P.addVariable("miss");
    P.addConstraint(lessThan((X.col(K - 1).segment(2, 2) - par(Vx(p.xf).segment(2, 2))).norm(), miss));
    P.addConstraint(equalTo(X.col(K - 1).segment(9, 2), par(Vx(p.xf).segment(9, 2))));
    P.addConstraint(equalTo(X.col(K - 1).segment(12, 2), par(Vx(p.xf).segment(12, 2))));
    for (int k = 0; k < K - 1; ++k) {
        P.addConstraint(equalTo(X.col(k + 1), dynpar(A[k]) * X.col(k) + dynpar(Bm[k]) * U.col(k) + dynpar(Bp[k]) * U.col(k + 1)
                                              + dynpar(S[k]) * sigma + dynpar(z[k]) + NU_.col(k)));
        P.addConstraint(box(-NUabs.col(k), NU_.col(k), NUabs.col(k)));
    }
    const double tan_glide = std::tan(p.glide), cos_gimbal = std::cos(p.gimbal_max), tilt = std::sqrt((1 - std::cos(p.tilt_max)) / 2);
    for (int k = 0; k < K; ++k) {
        if (k > 0) {   // node 0 is the measured state: path limits apply from the next node on
        P.addConstraint(greaterThan(X(0, k), p.m_dry));
        P.addConstraint(lessThan(X.block(2, k, 2, 1).norm(), (1.0 / tan_glide) * X(1, k)));  // glide slope
        P.addConstraint(lessThan(X.block(9, k, 2, 1).norm(), tilt));                     // tilt angle
        P.addConstraint(lessThan(X.block(11, k, 3, 1).norm(), p.w_max));                 // body rates
        }
        P.addConstraint(lessThan(U.col(k).norm(), (1.0 / cos_gimbal) * U(0, k)));        // gimbal angle
        P.addConstraint(lessThan(U.col(k).norm(), p.T_max));
        P.addConstraint(greaterThan(dynpar(Uhat[k]).dot(U.col(k)), p.T_min));             // T_min, linearised (Uhat = unit vector of last thrust)
        // trust region on state and control deviation from the reference
        VectorX dev(NX + NU); dev << X.col(k) - dynpar(Xr[k]), U.col(k) - dynpar(Ur[k]);
        P.addConstraint(lessThan(dev.norm(), D(k)));
    }
    P.addConstraint(box(-Ds, sigma - dynpar(sr), Ds));
    P.addConstraint(greaterThan(sigma, 0.1));
    P.addCostTerm(-X(0, K - 1) + o.w_miss * miss + o.w_nu * NUabs.sum() + o.w_trust * D.sum() + o.w_sigma * Ds);

    std::unique_ptr<ecos::ECOSSolver> solver;
    for (int it = 1; it <= o.max_iter; ++it) {
        const auto t0 = std::chrono::steady_clock::now();
        Discrete Dk = discretise(p, ref);
        for (int k = 0; k < K - 1; ++k) { A[k] = Dk.A[k]; Bm[k] = Dk.Bm[k]; Bp[k] = Dk.Bp[k]; S[k] = Dk.S[k]; z[k] = Dk.z[k]; }
        for (int k = 0; k < K; ++k) { Xr[k] = ref.X[k]; Ur[k] = ref.U[k]; Uhat[k] = ref.U[k].normalized(); }
        sr = ref.sigma;
        if (!solver) solver = std::make_unique<ecos::ECOSSolver>(P);
        solver->solve(false);
        if (std::getenv("R6_DEBUG")) std::fprintf(stderr, "    ecos exit %d %s\n", int(solver->getExitCode()), solver->getResultString().c_str());
        Trajectory next; next.sigma = eval(sigma);
        double trust = 0, nu = 0;
        for (int k = 0; k < K; ++k) {
            Vx x; for (int i = 0; i < NX; ++i) x(i) = eval(X(i, k));
            Vu u; for (int i = 0; i < NU; ++i) u(i) = eval(U(i, k));
            next.X.push_back(x); next.U.push_back(u);
            trust = std::max(trust, eval(D(k)));
            if (k < K - 1) for (int i = 0; i < NX; ++i) nu += std::abs(eval(NU_(i, k)));
        }
        const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        R.log.push_back({it, Dk.defect, trust, nu, -next.X.back()(0), ms});
        ref = next; R.iterations = it;
        if (trust < o.tol_trust && nu < o.tol_nu) { R.converged = true; break; }
    }
    R.solve_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - T0).count();
    return R;
}
}  // namespace r6
