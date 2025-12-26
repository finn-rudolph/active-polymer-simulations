#include "compute_config_moment.h"

#include <charconv>

#include "active_poly_constants.h"
#include "active_poly_util.h"
#include "atom.h"
#include "domain.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;

// syntax: compute config_moment [i_1][d_1] [i_2][d_2] ... (d_j = x, y, z)
// computes < l_{i_1,d_1} l_{i_2,d_2} ... >
// l_i is the connection vector from molecule 0 to i, i_j are 1-based
ComputeConfigMoment::ComputeConfigMoment(LAMMPS* lmp, int argc, char** argv)
    : Compute(lmp, argc, argv) {
    if (argc < 4)
        error->all(FLERR,
                   "`compute config_moment` requires at least one variable");

    for (int i = 3; i < argc; ++i) {
        int len = std::strlen(argv[i]), j;
        auto res = std::from_chars(argv[i], argv[i] + len - 1, j);
        int d = argv[i][len - 1] == 'x'
                    ? 0
                    : (argv[i][len - 1] == 'y' ? 1 : 2);
        if (d == 2 && argv[i][len - 1] != 'z')
            error->all(FLERR,
                       "the last character of each factor in "
                       "`compute config_moment` must be "
                       "`x`, `y` or `z`");
        vars.emplace_back(j, d);
    }

    scalar_flag = 1;
    extscalar = 0;
}

void ComputeConfigMoment::init() {}

double ComputeConfigMoment::compute_scalar() {
    double box_len[3];
    for (int d = 0; d < AP::d; ++d)
        box_len[d] = domain->boxhi[d] - domain->boxlo[d];

    double total = 0.0;
    for (int i = 0; i < atom->nlocal; ++i) {
        // The first atom computes the contribution of the whole molecule.
        if (atom->tag[i] % AP::N == 1) {
            int molecule_begin = atom->tag[i];
            double moment = 1.0;
            for (auto const& [index_in_molecule, d] : vars) {
                // index in molecule is 1-based
                int j = atom->map(molecule_begin + index_in_molecule);
                moment *= correct_coord_diff(atom->x[j][d] - atom->x[i][d], box_len[d]);
            }
            total += moment;
        }
    }

    // TODO: Is it really a good idea to do a sum here?? This could become huge...
    MPI_Allreduce(&total, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);

    scalar = scalar / (atom->natoms / AP::N);  // * force->nktv2p / (domain->xprd * domain->yprd * domain->zprd);
    return scalar;
}