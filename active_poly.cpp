#include <mpi.h>

#include <format>
#include <fstream>
#include <iostream>
#include <string>

#include "active_poly_constants.h"

using namespace std;

#define LAMMPS_LIB_MPI
#include "lammps/library.h"

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    void* lmp = lammps_open(0, 0, MPI_COMM_WORLD, 0);

    auto cmd = [&lmp](string const& s) {
        lammps_command(lmp, s.c_str());
    };

    // create molecule file
    ofstream molecule_file("active_poly.txt");

    molecule_file << format("{} atoms\n\n", AP::N);
    molecule_file << "Coords\n\n";
    for (int i = 1; i <= AP::N; ++i)
        molecule_file << format("{}   {}.0 0.0 0.0\n", i, i - 1);

    molecule_file << "\nTypes\n\n";
    for (int i = 1; i <= AP::N; ++i)
        molecule_file << format("{}   1\n", i);

    molecule_file.close();

    // basic setup
    cmd(format("dimension {}", AP::d));
    cmd("units lj");

    cmd("molecule ");
    // cmd("atom_style template");
    cmd("create_atoms 0 region R mol active_poly 42");

    lammps_commands_string(lmp,
                           "# Parameters\n"
                           "variable rho      equal 0.2\n"
                           "variable T        equal 1.0\n"
                           "variable length   equal 42\n"
                           "variable sigma    equal 1.122462048309373 # for LJ potential\n"
                           "variable damp     equal 0.5\n"
                           "\n"
                           "# Molecule setup\n"
                           "atom_style bond\n"
                           "bond_style harmonic\n"
                           "comm_modify vel yes cutoff 2.5\n"
                           "# newton bonds thing?\n"
                           "\n"
                           "lattice sc ${rho}\n"
                           "region r block 0 ${length} 0 ${length} 0 ${length}\n"
                           "create_box 1 r bond/types 1 extra/bond/per/atom 1\n"
                           "\n"
                           "molecule dumbbell dumbbell.txt\n"
                           "create_atoms 0 region r mol dumbbell 42\n"
                           "\n"
                           "bond_coeff 1 100 1 # the first number is the bond_type (made up by us), which also appears in the molecule file.\n"
                           "mass 1 1\n"
                           "# pair_style      lj/cut ${sigma}\n"
                           "# pair_coeff      * * 1 1\n"
                           "\n"
                           "fix 1 all nve\n"
                           "fix 2 all langevin $T $T ${damp} 42\n"
                           "# compute temp_compute all temp/partial \n"
                           "# fix_modify 2 temp temp_compute [x y z] # <-- necessary if we have some center of mass movement by the forces, or maybe also in shear flow.\n"
                           "\n"
                           "variable vx equal vcm(all, x)\n"
                           "variable vy equal vcm(all, y)\n"
                           "variable vz equal vcm(all, z)\n"
                           "\n"
                           "velocity all create $T 196883\n"
                           "thermo 1000\n"
                           "thermo_style custom step temp press etotal v_vx v_vy v_vz\n"
                           "timestep 0.001\n"
                           "run 10000");

    lammps_close(lmp);
    MPI_Finalize();
}
