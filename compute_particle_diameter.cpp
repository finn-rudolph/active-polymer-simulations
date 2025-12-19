#include "compute_particle_diameter.h"

#include "active_poly_constants.h"
#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

ComputeParticleDiameter::ComputeParticleDiameter(LAMMPS* lmp, int argc, char** argv)
    : Compute(lmp, argc, argv) {
    if (argc >= 4)
        error->all(FLERR, "`compute config_particle_diameter` has no arguments");
    scalar_flag = 1;
    extscalar = 0;
}

void ComputeParticleDiameter::init() {}

inline double correct_coord_diff(double diff, double len) {
    if (diff > len / 2)
        return diff - len;
    else if (diff < -len / 2)
        return diff + len;
    return diff;
}

double ComputeParticleDiameter::compute_scalar() {
    double box_len[3];
    for (int d = 0; d < AP::d; ++d)
        box_len[d] = domain->boxhi[d] - domain->boxlo[d];

    double total = 0.0;
    for (int i = 0; i < atom->nlocal; ++i) {
        // The first atom computes the diameter of its molecule.
        if (atom->tag[i] % AP::N == 1) {
            int molecule_begin = atom->tag[i];
            double diameter = 0.0;
            for (int t = 0; t < AP::N; ++t) {
                int j = atom->map(molecule_begin + t);
                for (int s = t + 1; s < AP::N; ++s) {
                    int k = atom->map(molecule_begin + s);
                    double distance = 0.0;
                    for (int d = 0; d < AP::d; ++d) {
                        double diff = correct_coord_diff(atom->x[j][d] - atom->x[k][d], box_len[d]);
                        distance += diff * diff;
                    }
                    diameter = std::max(diameter, sqrt(distance));
                }
            }
            total += diameter;
        }
    }

    MPI_Allreduce(&total, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);

    scalar = scalar / (atom->natoms / AP::N);  // * force->nktv2p / (domain->xprd * domain->yprd * domain->zprd);
    return scalar;
}