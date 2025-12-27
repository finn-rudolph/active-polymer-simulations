#ifdef COMPUTE_CLASS

ComputeStyle(active_stress_atom, ComputeActiveStressAtom)

#else

#ifndef COMPUTE_ACTIVE_STRESS_ATOM_H
#define COMPUTE_ACTIVE_STRESS_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeActiveStressAtom : public Compute {
   public:
    ComputeActiveStressAtom(class LAMMPS*, int, char**);

    virtual ~ComputeActiveStressAtom();

    void init();

    void compute_peratom();

    double memory_usage();

   private:
    int nmax;
    double** stress_tensor;
};

}  // namespace LAMMPS_NS

#endif
#endif
