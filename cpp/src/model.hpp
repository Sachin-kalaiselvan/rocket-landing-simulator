// 6-DOF rocket: variable mass, thrust vector in the body frame (gimballed engine), quaternion attitude.
// Non-dimensional units in the style of Szmuk & Acikmese (2018): gravity 1, wet mass 2.
// State x = [m, r(3), v(3), q(4, scalar first), w(3)]  (x axis up), control u = thrust in body frame.
#pragma once
#include <Eigen/Dense>
#include <cmath>

namespace r6 {
constexpr int NX = 14, NU = 3;
using Vx = Eigen::Matrix<double, NX, 1>;
using Vu = Eigen::Matrix<double, NU, 1>;
using Mxx = Eigen::Matrix<double, NX, NX>;
using Mxu = Eigen::Matrix<double, NX, NU>;
using V3 = Eigen::Vector3d;
using M3 = Eigen::Matrix3d;

struct Params {
    double m_wet = 2.0, m_dry = 1.0;
    double alpha = 0.05;                 // fuel use per unit thrust per unit time
    double T_min = 0.3, T_max = 5.0;     // thrust magnitude limits
    double gimbal_max = 20.0 * M_PI / 180;
    double tilt_max = 75.0 * M_PI / 180; // body x axis from vertical
    double glide = 20.0 * M_PI / 180;    // minimum elevation of the vehicle seen from the pad
    double w_max = 60.0 * M_PI / 180;    // body rate limit
    V3 g{-1, 0, 0};
    V3 r_T{-0.25, 0, 0};                 // engine gimbal point, below the centre of mass
    M3 J = Eigen::Vector3d(0.02, 0.02, 0.02).asDiagonal();
    Vx x0, xf;                           // initial and final state (final mass is free)
    double tf_guess = 5.0;
    Params() {
        x0 << m_wet, 4, 4, 0, -1, -2, 0, 1, 0, 0, 0, 0, 0, 0;
        xf << m_dry, 0, 0, 0, -0.1, 0, 0, 1, 0, 0, 0, 0, 0, 0;
    }
};

inline M3 dcm(const Eigen::Vector4d &q) {   // body -> inertial
    const double w = q(0), x = q(1), y = q(2), z = q(3);
    M3 C;
    C << 1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y),
         2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x),
         2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y);
    return C;
}

inline Eigen::Matrix4d omega(const V3 &w) {
    Eigen::Matrix4d O;
    O << 0, -w(0), -w(1), -w(2),
         w(0), 0, w(2), -w(1),
         w(1), -w(2), 0, w(0),
         w(2), w(1), -w(0), 0;
    return O;
}

inline Vx f(const Params &p, const Vx &x, const Vu &u) {
    Vx dx;
    const double m = x(0);
    const Eigen::Vector4d q = x.segment<4>(7);
    const V3 w = x.segment<3>(11);
    dx(0) = -p.alpha * u.norm();
    dx.segment<3>(1) = x.segment<3>(4);
    dx.segment<3>(4) = dcm(q) * u / m + p.g;
    dx.segment<4>(7) = 0.5 * omega(w) * q;
    dx.segment<3>(11) = p.J.inverse() * (p.r_T.cross(u) - w.cross(p.J * w));
    return dx;
}

// Jacobians by central differences (step 1e-6; truncation error ~1e-12, well below solver tolerance)
inline void jacobians(const Params &p, const Vx &x, const Vu &u, Mxx &A, Mxu &B) {
    const double h = 1e-6;
    for (int i = 0; i < NX; ++i) { Vx a = x, b = x; a(i) += h; b(i) -= h; A.col(i) = (f(p, a, u) - f(p, b, u)) / (2 * h); }
    for (int i = 0; i < NU; ++i) { Vu a = u, b = u; a(i) += h; b(i) -= h; B.col(i) = (f(p, x, a) - f(p, x, b)) / (2 * h); }
}
}  // namespace r6
