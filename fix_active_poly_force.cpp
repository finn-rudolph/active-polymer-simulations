#include "fix_active_poly_force.h"

#include <cassert>
#include <vector>

#include "active_poly_constants.h"
#include "atom.h"
#include "error.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;

FixActivePolyForce::FixActivePolyForce(LAMMPS* lmp, int argc, char** argv) : Fix(lmp, argc, argv) {
    if (argc != 3) error->all(FLERR, "invalid number of args to `fix active_poly_force`");
}

int FixActivePolyForce::setmask() {
    return FixConst::POST_FORCE;
}

void FixActivePolyForce::post_force(int) {
    double** x = atom->x;
    double** f = atom->f;
    int nlocal = atom->nlocal;

    std::vector<bool> processed(atom->natoms);
    int unique_bonds = 0;
    printf("nbond %d\n", neighbor->nbondlist);
    for (int i = 0; i < neighbor->nbondlist; ++i) {
        int o = neighbor->bondlist[i][0];  // reference atom in the molecule
        int j = neighbor->bondlist[i][1];

        // if (o >= nlocal && j >= nlocal) exit(42);
        // if (o < nlocal && processed[o]) exit(42); processed[o] = 1;
        // if (j < nlocal && processed[j]) exit(42); processed[j] = 1;

        if (atom->type[o] != 2) continue;
        ++unique_bonds;

        if (j < nlocal) {
            for (int k = 0; k < AP::N - 1; ++k) {
                int m = atom->map(atom->bond_atom[o][k]);
                f[j][0] += AP::Phi[j][k] * (x[m][0] - x[o][0]);
                f[j][1] += AP::Phi[j][k] * (x[m][1] - x[o][1]);
                f[j][2] += AP::Phi[j][k] * (x[m][2] - x[o][2]);
            }
        }

        // We have to make a choice where to update the reference atom. The
        // choice is the first bond in the bond list of `o`.
        if (o < nlocal && atom->bond_atom[o][0] == atom->tag[j]) {
            for (int k = 0; k < AP::N - 1; ++k) {
                int m = atom->map(atom->bond_atom[o][k]);
                f[o][0] += AP::F_origin[k] * (x[m][0] - x[o][0]);
                f[o][1] += AP::F_origin[k] * (x[m][1] - x[o][1]);
                f[o][2] += AP::F_origin[k] * (x[m][2] - x[o][2]);
            }
        }
    }
    int nproc = 0;
    for (int i = 0; i < processed.size(); ++i) nproc += processed[i];
    // printf("ub %d  \n", unique_bonds);

    double avgforce = 0.0, avgdist = 0.0, avgdist2 = 0.0;
    for (int i = 0; i < atom->natoms; ++i) {
        // int j = neighbor->bondlist[i][0];
        // int k = neighbor->bondlist[i][1];

        double dx[3];
        for (int d = 0; d < 3; ++d)
            dx[d] = f[i][d];
        avgforce += sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
    }
    int z = 0;
    for (int i = 0; i < atom->natoms; ++i) {
        if (atom->num_bond[i] == 0) continue;
        ++z;
        int j = atom->bond_atom[i][0];
        // int k = neighbor->bondlist[i][1];

        double dx[3];
        for (int d = 0; d < 3; ++d)
            dx[d] = x[i][d] - x[j][d];
        avgdist2 += sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
    }
    printf("z %d\n", z);
    for (int i = 0; i < neighbor->nbondlist; ++i) {
        int j = neighbor->bondlist[i][0];
        int k = neighbor->bondlist[i][1];

        double dx[3];
        for (int d = 0; d < 3; ++d)
            dx[d] = x[j][d] - x[k][d];
        avgdist += sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
    }
    printf("avgforce: %lf  avgdist %lf  avgdist2 %lf", avgforce / atom->natoms, avgdist / neighbor->nbondlist, avgdist2 / atom->natoms);
}
