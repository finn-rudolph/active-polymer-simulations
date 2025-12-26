#ifdef COMPUTE_CLASS

ComputeStyle(particle_diameter, ComputeParticleDiameter)

#else

#ifndef COMPUTE_PARTICLE_DIAMETER_H
#define COMPUTE_PARTICLE_DIAMETER_H

#include <vector>

#include "compute.h"

namespace LAMMPS_NS {

class ComputeParticleDiameter : public Compute {
   public:
    ComputeParticleDiameter(class LAMMPS*, int, char**);

    void init();

    double compute_scalar();
};

}  // namespace LAMMPS_NS

#endif
#endif