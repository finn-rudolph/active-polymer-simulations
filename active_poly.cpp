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

constexpr double rho = 0.2;
constexpr double l = 40;
constexpr double T = 1.0;
constexpr double damp_coeff = 0.5;
constexpr double sigma = 1.22;  // for intermolecular Lennard-Jones potential

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    void* lmp = lammps_open(0, 0, MPI_COMM_WORLD, 0);

    // auto cmd = [&lmp](string const& s) {
    //     lammps_command(lmp, s.c_str());
    // };

    // create molecule file
    ofstream molecule_file("active_poly.txt");

    molecule_file << format("\n{} atoms\n\n", AP::N);
    molecule_file << "Coords\n\n";
    for (int i = 1; i <= AP::N; ++i)
        molecule_file << format("{}   {}.0 0.0 0.0\n", i, i - 1);

    molecule_file << "\nTypes\n\n";
    for (int i = 1; i <= AP::N; ++i)
        molecule_file << format("{}   1\n", i);

    molecule_file.close();

    // basic setup
    cmd("dimension {}", AP::d);
    cmd("units lj");

    cmd("lattice sc {}", rho);
    cmd("region R block 0 {} 0 {} 0 {}", l, l, l);

    // cmd("comm_modify vel yes cutoff 2");  // idk whether this is relevant. There is also some newton bonds thing
    cmd("molecule active_poly active_poly.txt");
    // cmd("atom_style template");
    cmd("create_box 1 R");
    cmd("create_atoms 0 region R mol active_poly 42");

    cmd("mass 1 1");
    // inter-molecular interactions via Lennard-Jones potential
    // cmd("pair_style lj/cut {}", sigma);
    // cmd("pair_coeff * * 1 1");

    cmd("fix nve all nve");
    cmd("fix langevin all langevin {} {} {} 42", T, T, damp_coeff);

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
