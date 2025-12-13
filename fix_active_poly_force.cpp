#include "fix_active_poly_force.h"

#include <cassert>

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

    for (int i = 0; i < neighbor->nbondlist; ++i) {
        int o = neighbor->bondlist[i][0];  // reference atom in the molecule
        int j = neighbor->bondlist[i][1];

        assert(atom->type[o] == 2);
        assert(atom->type[j] == 1);
        assert(atom->num_bond[o] == AP::N - 1);

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
}
