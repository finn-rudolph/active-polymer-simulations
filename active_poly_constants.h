#include <array>

namespace AP {
constexpr int d = 3;
constexpr int N = 2;  // #atoms per molecule

constexpr std::array<std::array<double, N - 1>, N - 1> Phi = {{-1.0}};

constexpr auto compute_origin_force() {
    std::array<double, N - 1> F_origin = {};
    for (int i = 0; i < N - 1; ++i)
        for (int j = 0; j < N - 1; ++j)
            F_origin[i] -= Phi[j][i];
    return F_origin;
}

// This is the force on the reference atom in the molecule, which is the origin 
// of the local coordinate system. They are computed from the other forces such 
// that the total force on the molecule is 0.
constexpr auto F_origin = compute_origin_force();
}  // namespace AP
