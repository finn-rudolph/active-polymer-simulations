#include "compute_active_stress_atom.h"

#include <array>

#include "active_poly_constants.h"
#include "active_poly_util.h"
#include "atom.h"
#include "domain.h"
#include "fix_active_poly_force.h"
#include "memory.h"

using namespace LAMMPS_NS;

ComputeActiveStressAtom::ComputeActiveStressAtom(LAMMPS* lmp, int narg, char** arg)
    : Compute(lmp, narg, arg) {
    peratom_flag = 1;
    size_peratom_cols = 9;
    pressatomflag = 1;
    timeflag = 1;
    nmax = 0;
    stress_tensor = 0;
}

ComputeActiveStressAtom::~ComputeActiveStressAtom() {
    memory->destroy(stress_tensor);
}

void ComputeActiveStressAtom::init() {}

void ComputeActiveStressAtom::compute_peratom() {
    if (atom->nmax > nmax) {
        memory->destroy(stress_tensor);
        nmax = atom->nmax;
        memory->create(stress_tensor, nmax, 9, "active_stress_atom");
        array_atom = stress_tensor;
    }

    double box_len[3];
    for (int d = 0; d < AP::d; ++d)
        box_len[d] = domain->boxhi[d] - domain->boxlo[d];

    for (int i = 0; i < atom->nlocal; ++i) {
        if ((atom->mask[i] & groupbit) && atom->tag[i] % AP::N == 1) {
            memset(stress_tensor[i], 0, 9 * sizeof(double));

            int molecule_begin = atom->tag[i];

            for (int t = 0; t < AP::N; ++t) {
                int j = atom->map(molecule_begin + t);
                auto intramolecular_force = active_poly_force(j, atom, box_len);

                for (int d = 0; d < AP::d; ++d)
                    for (int e = 0; e < AP::d; ++e)  // using -F \otimes r
                        stress_tensor[i][3 * d + e] -=
                            intramolecular_force[d] *
                            correct_coord(atom->x[j][e], atom->x[i][e], box_len[e]);
            }
        }
    }
}

double ComputeActiveStressAtom::memory_usage() {
    return nmax * 9 * sizeof(double);
}
