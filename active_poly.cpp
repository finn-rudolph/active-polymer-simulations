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

constexpr double rho = 0.01;
constexpr double l = 35;
constexpr double T = 1.0;
constexpr double damp_coeff = 0.1;
constexpr double sigma = 1.22;  // for intermolecular Lennard-Jones potential

constexpr uint64_t equilibration_timesteps = 10000;
constexpr uint64_t run_timesteps = 50000;
constexpr double timestep = 0.001;

constexpr double shear_rate = 0.0002;
// TODO: monitor molecule diameter

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
        molecule_file << format("{}   {}.0 0.0 0.0\n", i, i - 1);

    molecule_file << "\nTypes\n\n";
    for (int i = 1; i <= AP::N; ++i)
        molecule_file << format("{}   1\n", i);

    // The only reason we have bonds is to ensure that all atoms in the molecule
    // appear as ghost atoms.
    molecule_file << "\nBonds\n\n";
    for (int i = 2; i <= AP::N; ++i)
        molecule_file << format("{}   {} {} {}\n", i - 1, 1, i - 1, i);

    molecule_file << "\nSpecial Bond Counts\n\n";
    for (int i = 1; i <= AP::N; ++i)
        molecule_file << format("{}   0 0 0\n", i);

    molecule_file << "\nSpecial Bonds\n\n";
    for (int i = 1; i <= AP::N; ++i)
        molecule_file << format("{}\n", i);

    molecule_file.close();

    // basic setup
    cmd("dimension {}", AP::d);
    cmd("units lj");

    cmd("atom_style bond");
    cmd("bond_style zero");
    cmd("comm_modify mode single cutoff {}", AP::N + 2.0);
    cmd("newton on off");  // Try changing this.

    cmd("lattice sc {}", rho);
    cmd("region R prism 0 {} 0 {} 0 {} 0 0 0", l, l, l);

    cmd("create_box 1 R bond/types 1 extra/bond/per/atom 2");
    cmd("molecule active_poly active_poly.txt");
    cmd("create_atoms 0 region R mol active_poly 42");

    cmd("bond_coeff *");
    cmd("mass * 1");
    // inter-molecular interactions via Lennard-Jones potential
    // cmd("pair_style lj/cut {}", sigma);
    // cmd("pair_coeff * * 1 1");

    cmd("fix 1 all nve");
    cmd("fix 2 all langevin {} {} {} 42", T, T, damp_coeff);
    cmd("fix 3 all active_poly_force");

    cmd("compute stress all pressure NULL virial");
    cmd("variable Txx equal c_stress[1]");
    cmd("variable Tyy equal c_stress[2]");
    cmd("variable Tzz equal c_stress[3]");

    cmd("compute 1x1x all config_moment 1x 1x");
    cmd("compute 1y1y all config_moment 1y 1y");
    cmd("compute 1z1z all config_moment 1z 1z");
    cmd("compute 1x1y all config_moment 1x 1y");
    cmd("compute 1x1z all config_moment 1x 1z");
    cmd("compute 1y1z all config_moment 1y 1z");


    // cmd("compute 1x2x all config_moment 1x 2x");
    // cmd("compute 1x2y all config_moment 1x 2y");
    // cmd("compute 2x2x all config_moment 2x 2x");
    // cmd("compute 2x2y all config_moment 2x 2y");
    // cmd("compute 2y2y all config_moment 2y 2y");

    cmd("compute diam all particle_diameter");

    for (auto& coord : {"x", "y", "z"}) {
        cmd("variable v{} equal vcm(all, {})", coord, coord);
    }

    cmd("velocity all create {} 196883", T);
    cmd("thermo 1000");
    cmd("thermo_style custom step temp c_1x1x c_1y1y c_1z1z c_1x1y c_1x1z c_1y1z ke c_diam");  // c_1x2x c_1x2y c_2x2x c_2x2y c_2y2y
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
