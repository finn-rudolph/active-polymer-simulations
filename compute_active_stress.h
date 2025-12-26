#ifdef COMPUTE_CLASS

ComputeStyle(active_stress, ComputeActiveStress)

#else

#ifndef COMPUTE_ACTIVE_STRESS_H
#define COMPUTE_ACTIVE_STRESS_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeActiveStress : public Compute {
   public:
    ComputeActiveStress(class LAMMPS*, int, char**);

    virtual ~ComputeActiveStress();

    void init();

    virtual void compute_vector();
};

}  // namespace LAMMPS_NS

#endif
#endif
