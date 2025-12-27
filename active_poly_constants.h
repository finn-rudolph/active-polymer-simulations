#include <array>

namespace AP {
constexpr int d = 3;
constexpr int N = 2;  // #atoms per molecule

constexpr std::array<std::array<double, N - 1>, N - 1> Phi =
    {{-200.0}};

// constexpr std::array<std::array<double, N - 1>, N - 1> Phi =
//     {{{{-100.0}}}};

constexpr auto position_force_matrix() {
    std::array<std::array<double, N>, N> F = {};
    for (int i = 0; i < N - 1; ++i)
        for (int j = 0; j < N - 1; ++j) {
            // The force applied to each bead is only half the "force" we want
            // on the connection vector.
            F[i + 1][j + 1] = Phi[i][j] / 2;
            F[i + 1][0] -= Phi[i][j] / 2;
            F[0][j + 1] -= Phi[i][j] / 2;
            F[0][0] += Phi[i][j] / 2;
        }
    return F;
}

// The force matrix with respect to positions (Phi is wrt position differences).
constexpr auto F = position_force_matrix();

enum ParticleType {
    Linear,

    // Triangle with E = K * |a|^2 / 2 potential.
    PassiveTriangle,

    // Triangle with 3 springs + a rotation force orthogonal to the plane of
    // the triangle.
    ActiveTriangleOrthogonal,
};

constexpr ParticleType particle_type = ParticleType::Linear;
}  // namespace AP
