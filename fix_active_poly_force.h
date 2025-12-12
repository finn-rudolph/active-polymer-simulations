#ifdef FIX_CLASS

FixStyle(active_poly_force, FixActivePolyForce)

#else
#ifndef FIX_ACTIVE_FORCE_H
#define FIX_ACTIVE_FORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixActivePolyForce : public Fix {
   public:
    FixActivePolyForce(LAMMPS* lmp, int argc, char** argv);

    int setmask();

    void post_force(int);
};

};  // namespace LAMMPS_NS

#endif
#endif