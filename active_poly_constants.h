#include <array>

namespace AP {
constexpr int d = 3;
constexpr int N = 2;  // #atoms per molecule

constexpr std::array<std::array<double, N - 1>, N - 1> Phi = {{-100.0}};

constexpr auto position_force_matrix() {
    std::array<std::array<double, N>, N> F = {};
    for (int i = 0; i < N - 1; ++i)
        for (int j = 0; j < N - 1; ++j) {
            F[i + 1][j + 1] = Phi[i][j];
            F[i + 1][0] -= Phi[i][j];
            F[0][j + 1] -= Phi[i][j];
            F[0][0] += Phi[i][j];
        }
    return F;
}

// The force matrix with respect to positions (Phi is wrt position differences).
constexpr auto F = position_force_matrix();
}  // namespace AP
