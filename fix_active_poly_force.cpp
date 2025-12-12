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

    int** molecule = neighbor->bondlist;
    assert(neighbor->nbondlist == AP::N);
}
