#include "fix_active_poly_force.h"

#include <cassert>
#include <utility>
#include <vector>

#include "active_poly_constants.h"
#include "active_poly_util.h"
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

void cross_prod(double a[3], double b[3], double result[3]) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
}

void sub(double a[3], double b[3], double result[3]) {
    for (int d = 0; d < 3; ++d) result[d] = a[d] - b[d];
}

void correct_coords(double a[3], double box_len[3]) {
    for (int d = 0; d < 3; ++d)
        a[d] = correct_coord_diff(a[d], box_len[d]);
}

void linear_forces(
    Atom* atom, double box_len[3], int molecule_begin, int index_in_molecule, double result[3]) {
    int i = atom->map(molecule_begin + index_in_molecule);

    for (int t = 0; t < AP::N; ++t) {
        int j = atom->map(molecule_begin + t);

        // TODO: for numerical stability it may be better to use relative
        // positions + Phi
        for (int d = 0; d < AP::d; ++d) {
            result[d] += AP::F[index_in_molecule][t] *
                         correct_coord(atom->x[j][d], atom->x[i][d], box_len[d]);
        }
    }
}

void passive_triangle_forces(
    Atom* atom, double box_len[3], int molecule_begin, int index_in_molecule, double result[3]) {
    int i = atom->map(molecule_begin + index_in_molecule);

    constexpr double spring_const = 10000.0;
    int j[3];
    for (int t = 0; t < 3; ++t) j[t] = atom->map(molecule_begin + t);

    double l1[3], l2[3];
    sub(atom->x[j[1]], atom->x[j[0]], l1);
    correct_coords(l1, box_len);
    sub(atom->x[j[2]], atom->x[j[0]], l2);
    correct_coords(l2, box_len);

    double l1_sq = l1[0] * l1[0] + l1[1] * l1[1] + l1[2] * l1[2];
    double l2_sq = l2[0] * l2[0] + l2[1] * l2[1] + l2[2] * l2[2];
    double l1_dot_l2 = l1[0] * l2[0] + l1[1] * l2[1] + l1[2] * l2[2];

    if (index_in_molecule == 0) {
        for (int d = 0; d < AP::d; ++d)
            result[d] -= spring_const * (l1[d] * (l1_dot_l2 - l2_sq) + l2[d] * (l1_dot_l2 - l1_sq));
    } else if (index_in_molecule == 1) {
        for (int d = 0; d < AP::d; ++d)
            result[d] -= spring_const * (l1[d] * l2_sq - l2[d] * l1_dot_l2);
    } else {
        for (int d = 0; d < AP::d; ++d)
            result[d] -= spring_const * (l2[d] * l2_sq - l1[d] * l1_dot_l2);
    }
}

std::array<double, 3> active_poly_force(int i, Atom* atom, double box_len[3]) {
    std::array<double, 3> result = {{0.0, 0.0, 0.0}};

    int tag = atom->tag[i];
    int index_in_molecule = (tag - 1) % AP::N;
    int molecule_begin = tag - index_in_molecule;

    if (AP::particle_type == AP::ParticleType::PassiveTriangle) {
        passive_triangle_forces(atom, box_len, molecule_begin, index_in_molecule, result.data());
    } else if (AP::particle_type == AP::ParticleType::Linear) {
        linear_forces(atom, box_len, molecule_begin, index_in_molecule, result.data());
    } else if (AP::particle_type == AP::ParticleType::ActiveTriangleOrthogonal) {
        double rot_const = 200.0;
        // {{{{-300.0, 0.0}}, {{0.0, -300.0}}}};
        linear_forces(atom, box_len, molecule_begin, index_in_molecule, result.data());

        if (index_in_molecule >= 1) {
            int j[3];
            for (int t = 0; t < 3; ++t) j[t] = atom->map(molecule_begin + t);

            double l1[3], l2[3], a[3];
            sub(atom->x[j[1]], atom->x[j[0]], l1);
            correct_coords(l1, box_len);
            sub(atom->x[j[2]], atom->x[j[0]], l2);
            correct_coords(l2, box_len);
            cross_prod(l1, l2, a);

            if (index_in_molecule == 1) {
                for (int d = 0; d < 3; ++d)
                    result[d] += rot_const * a[d];
            } else {
                for (int d = 0; d < 3; ++d)
                    result[d] -= rot_const * a[d];
            }
        }
    }

    return result;
}

void FixActivePolyForce::post_force(int) {
    double** x = atom->x;
    double** f = atom->f;
    int nlocal = atom->nlocal;

    double box_len[3];
    for (int d = 0; d < AP::d; ++d)
        box_len[d] = domain->boxhi[d] - domain->boxlo[d];

    for (int i = 0; i < atom->nlocal; ++i) {
        auto intramolecular_force = active_poly_force(i, atom, box_len);
        for (int d = 0; d < 3; ++d)
            f[i][d] += intramolecular_force[d];
    }
}
