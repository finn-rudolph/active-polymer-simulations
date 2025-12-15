#include "fix_active_poly_force.h"

#include <cassert>
#include <utility>
#include <vector>

#include "active_poly_constants.h"
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

inline double correct_coord(double coord, double reference, double len) {
    double difference = coord - reference;
    if (difference > len / 2)
        return coord - len;
    else if (difference < -len / 2)
        return coord + len;
    return coord;
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
