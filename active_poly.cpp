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

constexpr double rho = 0.05;
constexpr double l = 50;
constexpr double T = 1.0;
constexpr double damp_coeff = 0.5;
constexpr double sigma = 1.22;  // for intermolecular Lennard-Jones potential

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
    cmd("comm_modify mode single cutoff {}", AP::N + 1.0);
    cmd("newton on off");  // Try changing this.

    cmd("lattice sc {}", rho);
    cmd("region R block 0 {} 0 {} 0 {}", l, l, l);

    cmd("create_box 1 R bond/types 1 extra/bond/per/atom 2");
    cmd("molecule active_poly active_poly.txt");
    cmd("create_atoms 0 region R mol active_poly 42");

    cmd("bond_coeff *");
    cmd("mass * 1");
    // inter-molecular interactions via Lennard-Jones potential
    cmd("pair_style lj/cut {}", sigma);
    cmd("pair_coeff * * 1 1");

    cmd("fix 1 all nve");
    cmd("fix 2 all langevin {} {} {} 42", T, T, damp_coeff);
    cmd("fix 3 all active_poly_force");

    // This will be necessary when the active force causes a nonzero movement of
    // the center of mass. Then we must tell the langevin fix not to cancel this.
    // cmd("compute temp_compute all temp/partial");
    // cmd("fix_modify langevin temp temp_compute [x y z]");

    cmd("compute stress all pressure NULL virial");
    cmd("variable Txx equal c_stress[1]");
    cmd("variable Tyy equal c_stress[2]");
    cmd("variable Tzz equal c_stress[3]");

    // for (auto& coord : {"x", "y", "z"}) {
    //     cmd("variable v{} equal vcm(all, {})", coord, coord);
    // }

    cmd("velocity all create {} 196883", T);
    cmd("thermo 1000");
    cmd("thermo_style custom step temp press etotal v_Txx v_Tyy v_Tzz");
    cmd("timestep 0.001");
    cmd("run 10000");

    lammps_close(lmp);
    MPI_Finalize();
}

// https://matsci.org/t/lammps-users-number-of-atoms-in-molecules/39599
