#include "fix_active_poly_force.h"

#include <cassert>
#include <utility>
#include <vector>

#include "active_poly_constants.h"
#include "active_poly_util.h"
#include "atom.h"
#include "domain.h"
#include "error.h"
#include "neighbor.h"
#include "update.h"

using namespace LAMMPS_NS;

FixActivePolyForce::FixActivePolyForce(LAMMPS* lmp, int argc, char** argv)
    : Fix(lmp, argc, argv) {
    if (argc != 3)
        error->all(FLERR, "invalid number of args to `fix active_poly_force`");
}

int FixActivePolyForce::setmask() {
    return FixConst::POST_FORCE;
}

void FixActivePolyForce::post_force(int) {
    double** x = atom->x;
    double** f = atom->f;
    int nlocal = atom->nlocal;

    double box_len[3];
    for (int d = 0; d < AP::d; ++d)
        box_len[d] = domain->boxhi[d] - domain->boxlo[d];

    for (int i = 0; i < atom->nlocal; ++i) {
        int tag = atom->tag[i];
        int index_in_molecule = (tag - 1) % AP::N;
        int molecule_begin = tag - index_in_molecule;

        // --- PASSIVE TRIANGLE ---

        if (AP::N != 3) exit(42);
        constexpr double spring_const = 10000.0;
        int j[3];
        for (int t = 0; t < 3; ++t) j[t] = atom->map(molecule_begin + t);

        double l1[3], l2[3];
        for (int d = 0; d < AP::d; ++d) {
            l1[d] = correct_coord_diff(atom->x[j[1]][d] - atom->x[j[0]][d], box_len[d]);
            l2[d] = correct_coord_diff(atom->x[j[2]][d] - atom->x[j[0]][d], box_len[d]);
        }
        double l1_sq = l1[0] * l1[0] + l1[1] * l1[1] + l1[2] * l1[2];
        double l2_sq = l2[0] * l2[0] + l2[1] * l2[1] + l2[2] * l2[2];
        double l1_dot_l2 = l1[0] * l2[0] + l1[1] * l2[1] + l1[2] * l2[2];

        if (index_in_molecule == 0) {
            for (int d = 0; d < AP::d; ++d)
                f[i][d] -= spring_const * (l1[d] * (l1_dot_l2 - l2_sq) + l2[d] * (l1_dot_l2 - l1_sq));
        } else if (index_in_molecule == 1) {
            for (int d = 0; d < AP::d; ++d)
                f[i][d] -= spring_const * (l1[d] * l2_sq - l2[d] * l1_dot_l2);
        } else {
            for (int d = 0; d < AP::d; ++d)
                f[i][d] -= spring_const * (l2[d] * l2_sq - l1[d] * l1_dot_l2);
        }

        continue;

        // --- LINEAR FORCES ---

        for (int t = 0; t < AP::N; ++t) {
            int j = atom->map(molecule_begin + t);

            // TODO: for numerical stability it may be better to use relative
            // positions + Phi
            for (int d = 0; d < AP::d; ++d) {
                f[i][d] += AP::F[index_in_molecule][t] *
                           correct_coord(x[j][d], x[i][d], box_len[d]);
            }
        }
    }
}
