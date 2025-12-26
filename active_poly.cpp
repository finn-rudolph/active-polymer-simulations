#include <mpi.h>

#include <format>
#include <fstream>
#include <iostream>
#include <string>

#include "active_poly_constants.h"

using namespace std;

#define LAMMPS_LIB_MPI
#include "lammps/library.h"

#define cmd(s, ...) lammps_command(lmp, format(s, ##__VA_ARGS__).c_str())

constexpr std::array<char, 3> variables = {{'x', 'y', 'z'}};

constexpr double rho = 0.02;
constexpr double l = 40;
constexpr double T = 1.0;
constexpr double damp_coeff = 0.1;
constexpr double sigma = 1.22;  // for intermolecular Lennard-Jones potential

constexpr uint64_t equilibration_timesteps = 10000;
constexpr uint64_t run_timesteps = 100000;
constexpr double timestep = 0.001;

constexpr double shear_rate = 0.003;

// TODO: compare equilibrium second moments with theory
//       compare signs of stress differences with theory
//       check the random force derivation

void compute_cross_product_second_moments(void* lmp);

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    void* lmp = lammps_open(0, 0, MPI_COMM_WORLD, 0);

    // Create molecule file. The molecule structure is not really used for
    // anything, but ensures the total number of atoms is a multiple of N.
    // The forces are set via the atom tags.
    ofstream molecule_file("active_poly.txt");

    molecule_file << format(
        "\n"
        "{} atoms\n"
        "{} bonds\n",
        AP::N, AP::N - 1);

    molecule_file << "\nCoords\n\n";
    for (int i = 1; i <= AP::N; ++i)
        molecule_file << format("{}   {} 0.0 0.0\n", i, ((float)i - 1) / 10);

    molecule_file << "\nTypes\n\n";
    for (int i = 1; i <= AP::N; ++i)
        molecule_file << format("{}   1\n", i);

    // The only reason we have bonds is to ensure that all atoms in the molecule
    // appear as ghost atoms.
    molecule_file << "\nBonds\n\n";
    for (int i = 2; i <= AP::N; ++i)
        molecule_file << format("{}   {} {} {}\n", i - 1, 1, i - 1, i);

    // molecule_file << "\nSpecial Bond Counts\n\n";
    // for (int i = 1; i <= AP::N; ++i)
    //     molecule_file << format("{}   0 0 0\n", i);

    // molecule_file << "\nSpecial Bonds\n\n";
    // for (int i = 1; i <= AP::N; ++i)
    //     molecule_file << format("{}\n", i);

    molecule_file.close();

    // basic setup
    cmd("dimension {}", AP::d);
    cmd("units lj");

    cmd("atom_style bond");
    cmd("bond_style harmonic");
    cmd("comm_modify mode single cutoff {}", 4.0);
    cmd("newton on off");  // Try changing this.

    cmd("lattice sc {}", rho);
    cmd("region R prism 0 {} 0 {} 0 {} 0 0 0", l, l, l);

    cmd("create_box 1 R bond/types 1 extra/bond/per/atom 2");
    cmd("molecule active_poly active_poly.txt");
    cmd("create_atoms 0 region R mol active_poly 42");

    cmd("bond_coeff * 100 0");
    cmd("mass * 1");
    // inter-molecular interactions via Lennard-Jones potential
    // cmd("pair_style lj/cut {}", 2.5);
    // cmd("pair_coeff * * 1 {}", sigma);

    cmd("fix 1 all nve");
    cmd("fix 2 all langevin {} {} {} 42", T, T, damp_coeff);
    // cmd("fix 3 all active_poly_force");

    cmd("compute T all pressure NULL virial");
    cmd("variable Txx equal c_T[1]");
    cmd("variable Tyy equal c_T[2]");
    cmd("variable Tzz equal c_T[3]");
    cmd("variable Txy equal c_T[4]");
    cmd("variable Txz equal c_T[5]");
    cmd("variable Tyz equal c_T[6]");


    cmd("compute 1x1x all config_moment 1x 1x");
    cmd("compute 1y1y all config_moment 1y 1y");
    cmd("compute 1z1z all config_moment 1z 1z");
    cmd("compute 1x1y all config_moment 1x 1y");
    cmd("compute 1x1z all config_moment 1x 1z");
    cmd("compute 1y1z all config_moment 1y 1z");

    // cmd("compute 1x2x all config_moment 1x 2x");
    // cmd("compute 2x2x all config_moment 2x 2x");

    // cmd("compute S all active_stress");
    // for (int i = 0; i < 3; ++i)
    //     for (int j = 0; j < 3; ++j)
    //         cmd("variable S{}{} equal c_S[{}]", variables[i], variables[j], 3 * i + j);

    cmd("compute diam all particle_diameter");
    // compute_cross_product_second_moments(lmp);

    // for (auto& coord : {"x", "y", "z"}) {
    //     cmd("variable v{} equal vcm(all, {})", coord, coord);
    // }

    cmd("velocity all create {} 196883", T);
    cmd("thermo 1000");
    // cmd("thermo_style cus^tom step temp v_axx v_ayy v_azz v_axy v_ayz v_azx ke c_diam");
    // cmd("thermo_style custom step temp c_1x1x c_1y1y c_1x1y c_1x2x c_2x2x v_axx ke c_diam");
    // cmd("thermo_style custom step temp v_Txx v_Tyy v_Tzz v_Txy v_Txz v_Tyz c_diam");
    cmd("thermo_style custom step temp c_1x1x c_1y1y c_1z1z c_1x1y c_1x1z c_1y1z c_diam");

    cmd("timestep {}", timestep);
    cmd("run {}", equilibration_timesteps);

    if (shear_rate > 0.0) {
        cmd("fix 4 all deform 1 xy erate {} remap v", shear_rate);

        // This calculates temperature correctly under deformation.
        cmd("compute temp_deform all temp/deform");
        cmd("fix_modify 2 temp temp_deform");
        // cmd("velocity all ramp vx 0.0 1.0 y 0.0 128.0 temp temp_deform");
    }

    cmd("run {}", run_timesteps);

    lammps_close(lmp);
    MPI_Finalize();
}

void moment(void* lmp, std::string s) {
    std::string args = std::string(s.begin(), s.begin() + 2);
    for (int i = 2; i < s.size(); i += 2) {
        args.push_back(' ');
        args.append(s.begin() + i, s.begin() + i + 2);
    }
    std::replace(args.begin(), args.end(), '_', ' ');
    cmd("compute {} all config_moment {}", s, args);
}

void compute_cross_product_second_moments(void* lmp) {
    for (auto const& a : variables)
        for (auto const& b : variables) {
            if (a != b) {
                moment(lmp, format("1{}1{}2{}2{}", a, a, b, b));
                moment(lmp, format("1{}1{}2{}2{}", a, b, a, b));

                for (auto const& c : variables)
                    if (a != c && b != c) {
                        moment(lmp, format("1{}1{}2{}2{}", a, b, c, c));
                        moment(lmp, format("1{}1{}2{}2{}", a, a, b, c));
                        moment(lmp, format("1{}1{}2{}2{}", b, a, a, c));
                        moment(lmp, format("1{}1{}2{}2{}", a, b, c, a));
                    }
            }
        }

    // squares
    cmd("variable axx equal c_1y1y2z2z + c_1z1z2y2y - 2 * c_1y1z2y2z");
    cmd("variable ayy equal c_1z1z2x2x + c_1x1x2z2z - 2 * c_1z1x2z2x");
    cmd("variable azz equal c_1x1x2y2y + c_1y1y2x2x - 2 * c_1x1y2x2y");

    // mixed products
    cmd("variable axy equal c_1y1z2z2x - c_1y1x2z2z - c_1z1z2y2x + c_1z1x2y2z");
    cmd("variable ayz equal c_1z1x2x2y - c_1z1y2x2x - c_1x1x2z2y + c_1x1y2z2x");
    cmd("variable azx equal c_1x1y2y2z - c_1x1z2y2y - c_1y1y2x2z + c_1y1z2x2y");
}
