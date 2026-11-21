// rocket6dof: solve one landing and write the trajectory as JSON.
//   rocket6dof [--K 30] [--x0 m rx ry rz vx vy vz qw qx qy qz wx wy wz] > out.json
#include "scvx.hpp"
#include <cstdio>
#include <cstring>
#include <iostream>

int main(int argc, char **argv) {
    r6::Params p; r6::Options o;
    for (int i = 1; i < argc; ++i) {
        if (!std::strcmp(argv[i], "--K")) o.K = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--tf")) p.tf_guess = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--x0")) for (int j = 0; j < r6::NX; ++j) p.x0(j) = std::atof(argv[++i]);
    }
    const r6::Result R = r6::solve(p, o);
    // independent check: fly the solved controls through the nonlinear model with fine RK4 steps
    const double dt = R.traj.sigma / (o.K - 1);
    r6::Vx x = p.x0; double worst = 0;
    for (int k = 0; k < o.K - 1; ++k) { x = r6::propagate(p, x, R.traj.U[k], R.traj.U[k + 1], dt, 50); worst = std::max(worst, (x - R.traj.X[k + 1]).norm()); }
    std::printf("{\"converged\": %s, \"iterations\": %d, \"solve_ms\": %.2f, \"tf\": %.6f, \"m_final\": %.6f, \"open_loop_defect\": %.3e, \"final_state_error\": %.3e,\n",
                R.converged ? "true" : "false", R.iterations, R.solve_ms, R.traj.sigma, R.traj.X.back()(0), worst, (x - R.traj.X.back()).norm());
    std::printf("\"log\": [");
    for (size_t i = 0; i < R.log.size(); ++i) { const auto &l = R.log[i]; std::printf("%s[%d, %.3e, %.3e, %.3e, %.6f, %.2f]", i ? ", " : "", l.it, l.defect, l.trust, l.nu, l.cost, l.ms); }
    std::printf("],\n\"X\": [");
    for (int k = 0; k < o.K; ++k) { std::printf("%s[", k ? ", " : ""); for (int i = 0; i < r6::NX; ++i) std::printf("%s%.8f", i ? ", " : "", R.traj.X[k](i)); std::printf("]"); }
    std::printf("],\n\"U\": [");
    for (int k = 0; k < o.K; ++k) std::printf("%s[%.8f, %.8f, %.8f]", k ? ", " : "", R.traj.U[k](0), R.traj.U[k](1), R.traj.U[k](2));
    std::printf("]}\n");
    return R.converged ? 0 : 1;
}
