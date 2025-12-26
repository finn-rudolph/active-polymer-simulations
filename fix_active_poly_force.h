#ifdef FIX_CLASS

FixStyle(active_poly_force, FixActivePolyForce)

#else
#ifndef FIX_ACTIVE_FORCE_H
#define FIX_ACTIVE_FORCE_H

#include "fix.h"
#include "atom.h"

namespace LAMMPS_NS {

class FixActivePolyForce : public Fix {
   public:
    FixActivePolyForce(LAMMPS* lmp, int argc, char** argv);

    int setmask();

    void post_force(int);
};

};  // namespace LAMMPS_NS

std::array<double, 3> active_poly_force(int i, LAMMPS_NS::Atom* atom, double box_len[3]);

#endif
#endif