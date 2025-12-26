#include "compute_active_stress.h"

#include "active_poly_util.h"
#include "active_poly_constants.h"
#include "fix_active_poly_force.h"

using namespace LAMMPS_NS;

ComputeActiveStress::ComputeActiveStress(LAMMPS* lmp, int argc, char** argv)
    : Compute(lmp, argc, argv) {
    vector_flag = 1;
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

    double stress_tensor[9]; // T_ij is at index 3 * i + j.

    for (int i = 0; i < atom->nlocal; ++i) {
        auto intramolecular_force = active_poly_force(i, atom, box_len);
        for (int d = 0; d < AP::d; ++d)
            for (int e = 0; e < AP::d; ++e) // using -F \otimes r
                stress_tensor[3 * d + e] += intramolecular_force[d] * atom->x[]
    }
}