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

constexpr double rho = 0.1;
constexpr double l = 20;
constexpr double T = 1.0;
constexpr double damp_coeff = 0.5;
constexpr double sigma = 1.22;  // for intermolecular Lennard-Jones potential

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    void* lmp = lammps_open(0, 0, MPI_COMM_WORLD, 0);

    // create molecule file
    ofstream molecule_file("active_poly.txt");

    molecule_file << format(
        "\n"
        "{} atoms\n"
        "{} bonds\n",
        AP::N, AP::N - 1);

    molecule_file << "\nCoords\n\n";
    for (int i = 1; i <= AP::N; ++i)
        molecule_file << format("{}   {}.0 0.0 0.0\n", i, i - 1);

    // Atom 1 is the "reference atom" in the molecule, which is used as the
    // origin of the local coordinate system. It has type 2, all other atoms
    // have type 1.
    molecule_file << "\nTypes\n\n";
    molecule_file << "1   2\n";
    for (int i = 2; i <= AP::N; ++i)
        molecule_file << format("{}   1\n", i);

    // The molecule is a chain with the type 2 atom at one end. The chain
    // structure is useful to keep track of the indices, but there are no 
    // forces associated to the bonds.
    // [specific bond ID] [bond type] [atom 1 ID] [atom 2 ID]
    molecule_file << "\nBonds\n\n";
    for (int i = 2; i <= AP::N; ++i)
        molecule_file << format("{}   {} {} {}\n", i - 1, 1, i - 1, i);

    molecule_file.close();

    // basic setup
    cmd("dimension {}", AP::d);
    cmd("units lj");

    cmd("atom_style bond");
    cmd("bond_style zero");
    cmd("comm_modify vel yes cutoff 2");  // idk whether this is relevant.
    cmd("newton on off"); // Try changing this.
    
    cmd("lattice sc {}", rho);
    cmd("region R block 0 {} 0 {} 0 {}", l, l, l);
    
    cmd("create_box 2 R bond/types 1 extra/bond/per/atom {}", AP::N-1);
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

    // This will be necessary when the active force causes a nonzero movement of
    // the center of mass. Then we must tell the langevin fix not to cancel this.
    // cmd("compute temp_compute all temp/partial");
    // cmd("fix_modify langevin temp temp_compute [x y z]");

    for (auto& coord : {"x", "y", "z"}) {
        cmd("variable v{} equal vcm(all, {})", coord, coord);
    }

    cmd("velocity all create {} 196883", T);
    cmd("thermo 1000");
    cmd("thermo_style custom step temp press etotal v_vx v_vy v_vz");
    cmd("timestep 0.001");
    cmd("run 10000");

    lammps_close(lmp);
    MPI_Finalize();
}

// https://matsci.org/t/lammps-users-number-of-atoms-in-molecules/39599
