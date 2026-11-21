// Closed-loop 6-DOF landings with re-planning (model-predictive guidance) and Monte Carlo dispersions.
//   rocket6dof_mc <runs> <seed> > results.json
// Truth model differs from the planner's model: mass and inertia errors, a steady side force (wind),
// thrust magnitude error and an engine misalignment torque. Every 0.4 time units the planner re-solves
// from the measured state, warm-started from its previous plan; it stops re-planning 0.6 before landing.
#include "scvx.hpp"
#include <cstdio>
#include <random>
#include <cstdlib>
#include <algorithm>
#include <string>

using namespace r6;

static Trajectory shift(const Trajectory &t, double elapsed, const Vx &now) {
    // resample the remaining part of the previous plan onto K nodes; first node = measured state
    const int K = t.X.size(); Trajectory out; out.sigma = std::max(0.3, t.sigma - elapsed);
    for (int k = 0; k < K; ++k) {
        const double tau = (elapsed + double(k) / (K - 1) * out.sigma) / t.sigma * (K - 1);
        const int i = std::min(K - 2, int(tau)); const double a = std::min(1.0, tau - i);
        out.X.push_back((1 - a) * t.X[i] + a * t.X[i + 1]); out.U.push_back((1 - a) * t.U[i] + a * t.U[i + 1]);
    }
    out.X[0] = now; return out;
}

int main(int argc, char **argv) {
    const int runs = argc > 1 ? std::atoi(argv[1]) : 50; const unsigned seed = argc > 2 ? std::atoi(argv[2]) : 1;
    const double D = (argc > 3 && std::string(argv[3]) == "nodisp") ? 0.0 : 1.0;   // scale of all dispersions
    std::mt19937 rng(seed); std::normal_distribution<double> N01(0, 1); std::uniform_real_distribution<double> U01(0, 1);
    std::printf("[");
    for (int run = 0; run < runs; ++run) {
        Params plan;                                     // what the guidance believes
        plan.x0(1) += D * 0.4 * (2 * U01(rng) - 1); plan.x0(2) += D * 1.0 * (2 * U01(rng) - 1); plan.x0(3) += D * 1.0 * (2 * U01(rng) - 1);
        plan.x0(4) += D * 0.3 * (2 * U01(rng) - 1); plan.x0(5) += D * 0.5 * (2 * U01(rng) - 1); plan.x0(6) += D * 0.5 * (2 * U01(rng) - 1);
        Params truth = plan;                             // what actually flies
        truth.x0(0) *= 1 + D * 0.02 * N01(rng);              // mass error 2% (1 sigma)
        truth.J *= 1 + D * 0.10 * N01(rng);                  // inertia error 10%
        const V3 wind(0, D * 0.03 * N01(rng), D * 0.03 * N01(rng));            // steady side acceleration (3% of g)
        const double thrust_err = 1 + D * 0.02 * N01(rng);                 // thrust magnitude error 2%
        const V3 misalign(0, D * 0.005 * N01(rng), D * 0.005 * N01(rng));     // engine misalignment (rad) -> torque
        Options o; Trajectory traj; double t = 0, solve_ms = 0, worst_ms = 0; int solves = 0, fails = 0; Vx x = truth.x0;
        bool first = true; double next_plan = 0; Trajectory cur; double cur_t0 = 0;
        const double dt = 0.005; double tilt_max_seen = 0, rate_max_seen = 0;
        while (true) {
            if (t >= next_plan && (first || cur.sigma - (t - cur_t0) > 0.15)) {
                Params p = plan; p.x0 = x; p.x0(0) = x(0) * (1 + 0.0);   // the planner measures the full state
                Trajectory warm = first ? initial_guess(p, o.K) : shift(cur, t - cur_t0, x);
                if (first) warm.sigma = p.tf_guess;
                Result R = solve(p, o, &warm);
                solves++; solve_ms += R.solve_ms; worst_ms = std::max(worst_ms, R.solve_ms);
                // accept a nearly converged plan (dynamics slack tiny) rather than keep flying a stale one
                const bool usable = R.converged || (!R.log.empty() && R.log.back().nu < 1e-4 && R.log.back().defect < 1e-3);
                if (usable || first) { cur = R.traj; cur_t0 = t; } else { fails++; if (std::getenv("R6_DEBUG")) { std::fprintf(stderr, "fail t=%.2f rem=%.2f\n", t, cur.sigma - (t - cur_t0)); for (auto &l : R.log) std::fprintf(stderr, "  %d defect %.2e trust %.2e nu %.2e cost %.4f\n", l.it, l.defect, l.trust, l.nu, l.cost); } }
                first = false; next_plan = t + 0.2;
            }
            // control from the current plan (first-order hold)
            const double tn = std::min(1.0, (t - cur_t0) / cur.sigma) * (o.K - 1);
            const int i = std::min(o.K - 2, int(tn)); const double a = tn - i;
            // tracking controller around the plan (runs every step; the planner closes the slow loop by re-solving):
            // attitude: PD on the quaternion error to the planned attitude, executed with the gimbal (pitch and yaw only;
            // the engine cannot roll the vehicle). Thrust magnitude: planned value plus a correction on the position
            // and velocity error along the body axis. Sideways position errors are left to the next re-plan.
            const Vx xr = (1 - a) * cur.X[i] + a * cur.X[i + 1]; const Vu ur = (1 - a) * cur.U[i] + a * cur.U[i + 1];
            const Eigen::Quaterniond qr(xr(7), xr(8), xr(9), xr(10)), qn(x(7), x(8), x(9), x(10));
            Eigen::Quaterniond qe = qr.normalized().conjugate() * qn.normalized(); if (qe.w() < 0) qe.coeffs() *= -1;
            const V3 e_att = 2 * qe.vec(), w = x.segment<3>(11), wr = xr.segment<3>(11);
            const V3 tau = plan.J * (-25.0 * e_att - 8.0 * (w - wr));                   // natural frequency 5, damping 0.8
            const M3 C = dcm(x.segment<4>(7)); const V3 bx = C.col(0);
            const V3 acc_corr = 0.8 * (xr.segment<3>(1) - x.segment<3>(1)) + 1.6 * (xr.segment<3>(4) - x.segment<3>(4));
            const double Tmag = std::clamp(ur(0) + x(0) * acc_corr.dot(bx), truth.T_min, truth.T_max);
            // torque = r_T x u with r_T = (-0.25, 0, 0): tau_y = 0.25 u_z, tau_z = -0.25 u_y; add to the planned gimbal
            Vu u(Tmag, ur(1) - tau(2) / 0.25, ur(2) + tau(1) / 0.25);
            const double lat = u.tail<2>().norm(), lat_max = Tmag * std::tan(truth.gimbal_max);
            if (lat > lat_max) u.tail<2>() *= lat_max / lat;
            u *= thrust_err;
            // truth dynamics: planner model + wind + misalignment torque (rotate thrust slightly in the body)
            Vu ua = u; ua(1) += misalign(1) * u(0); ua(2) += misalign(2) * u(0);
            Vx k1 = f(truth, x, ua); k1.segment<3>(4) += wind;
            x += dt * k1; x.segment<4>(7).normalize(); t += dt;
            const double tilt = std::acos(std::min(1.0, 1 - 2 * (x(9) * x(9) + x(10) * x(10))));
            tilt_max_seen = std::max(tilt_max_seen, tilt); rate_max_seen = std::max(rate_max_seen, x.segment<3>(11).norm());
            if (x(1) <= 0 || t > 3 * cur.sigma + cur_t0 || x(0) <= truth.m_dry) break;
        }
        const double miss = x.segment<2>(2).norm(), vs = -x(4), hs = x.segment<2>(5).norm();
        const double tilt = std::acos(std::min(1.0, 1 - 2 * (x(9) * x(9) + x(10) * x(10)))) * 180 / M_PI;
        std::printf("%s{\"run\": %d, \"miss\": %.5f, \"v_down\": %.5f, \"v_side\": %.5f, \"tilt_deg\": %.3f, \"rate\": %.4f, \"fuel\": %.5f, \"t\": %.4f, \"solves\": %d, \"failed_solves\": %d, \"mean_ms\": %.2f, \"worst_ms\": %.2f, \"max_tilt_deg\": %.2f, \"max_rate\": %.4f}\n",
                    run ? "," : "", run, miss, vs, hs, tilt, x.segment<3>(11).norm(), truth.x0(0) - x(0), t, solves, fails, solve_ms / solves, worst_ms, tilt_max_seen * 180 / M_PI, rate_max_seen);
        std::fflush(stdout);
    }
    std::printf("]\n");
}
