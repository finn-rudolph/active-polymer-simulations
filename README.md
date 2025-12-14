# Simulations of active polymer models

## Setup

1. Set the environment variable ´LAMMPS_DIR´ to the directory of the LAMMPS repository, which should have subdirectories src, build, etc.

## Remarks about LAMMPS

- LAMMPS tries to keep the system at constant temperature. If we have active forces in the particles, they
  pump energy into the system. The langevin fix (I guess) tries to counteract that by weakening the random
  forces, so the system approaches an equilibrium temperature and total energy. This is not correct, the
  energy should steadily increase. Maybe one can fix this by subtracting the amount of energy pumped into 
  the system from the total energy. This can be done after the force computation in the fix. Of course,
  this causes the simulation to become unstable long-term, and a continuous supply of active energy is also
  not very realistic. Maybe the particles should also be able to "absorb" energy somehow.

- There are local atom indices and atom tags. For access of any atom data, e.g. position or force, the index is required. The   
  `bond_atoms` array stores the tags.

- Even with one processor, there are ghost atoms and more bonds than there should be. This is due to the periodicity of the simulation box. If a bond crosses the boundary, it does not make sense to use the coordinates lying inside the box for both atoms, because one will be very small and the other very large. Thus (this is how I think it works, not 100% sure) LAMMPS creates two bonds and ghost atoms for every bond crossing the boundary, so the coordinates at the indices given in the bound are correct (relative to each other). As a consequence, one should not use the `bond_atoms` array plus `atom->map` to get the coordinates of the bond partner.

- What I don't understand: With one processor, not every atom occurs in a bond with its _owned_ index.
    ```
    for (int i = 0; i < neighbor->nbondlist; ++i) {
            int o = neighbor->bondlist[i][0];  // reference atom in the molecule
            int j = neighbor->bondlist[i][1];

            if (o < nlocal) processed[o] = 1;
            if (j < nlocal) processed[j] = 1;
    }
    int nproc = 0;
    for (int i = 0; i < processed.size(); ++i) nproc += processed[i];
    printf("nprocessed %d  \n", nproc);
    ```
    For 16000 atoms, this prints numbers fluctuating between 15988 and 15998.

