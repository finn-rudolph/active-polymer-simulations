#ifdef COMPUTE_CLASS

ComputeStyle(config_moment, ComputeConfigMoment)

#else

#ifndef LMP_COMPUTE_CONFIG_MOMENT_H
#define LMP_COMPUTE_CONFIG_MOMENT_H

#include <vector>

#include "compute.h"

namespace LAMMPS_NS {

class ComputeConfigMoment : public Compute {
   public:
    ComputeConfigMoment(class LAMMPS*, int, char**);

    void init();

    double compute_scalar();

   private:
    std::vector<std::pair<int, int> > vars;
};

}  // namespace LAMMPS_NS

#endif
#endif