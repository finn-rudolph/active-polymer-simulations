#include "compute_active_stress.h"

#include "active_poly_constants.h"
#include "active_poly_util.h"
#include "atom.h"
#include "domain.h"
#include "fix_active_poly_force.h"

using namespace LAMMPS_NS;

ComputeActiveStress::ComputeActiveStress(LAMMPS* lmp, int argc, char** argv)
    : Compute(lmp, argc, argv) {
    vector_flag = 1;
    extvector = 0;
    size_vector = 9;
    vector = new double[size_vector];
}

ComputeActiveStress::~ComputeActiveStress() {
    delete[] vector;
}

void ComputeActiveStress::init() {}

void ComputeActiveStress::compute_vector() {
    double box_len[3];
    for (int d = 0; d < AP::d; ++d)
        box_len[d] = domain->boxhi[d] - domain->boxlo[d];

    double stress_tensor[9];  // T_ij is at index 3 * i + j.
    memset(stress_tensor, 0, sizeof stress_tensor);

    for (int i = 0; i < atom->nlocal; ++i) {
        // The first atom computes the contribution of the whole molecule.
        // For the Kirkwood formula to work, we need to make sure that the
        // positions of atoms in a molecule crossing the cyclic boundary have
        // correct coordinates relative to each other.

        if (atom->tag[i] % AP::N == 1) {
            int molecule_begin = atom->tag[i];

            for (int t = 0; t < AP::N; ++t) {
                int j = atom->map(molecule_begin + t);
                auto intramolecular_force = active_poly_force(j, atom, box_len);

                for (int d = 0; d < AP::d; ++d)
                    for (int e = 0; e < AP::d; ++e)  // using -F \otimes r
                        stress_tensor[3 * d + e] -=
                            intramolecular_force[d] *
                            correct_coord(atom->x[j][e], atom->x[i][e], box_len[e]);
            }
        }
    }

    MPI_Allreduce(stress_tensor, vector, size_vector, MPI_DOUBLE, MPI_SUM, world);

    for (int i = 0; i < size_vector; ++i)
        vector[i] /= atom->natoms / AP::N;  // maybe this is wrong
}