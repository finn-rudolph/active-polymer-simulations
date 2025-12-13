```
cmake -S ../cmake -B . -D BUILD_SHARED_LIBS=yes -D PKG_MOLECULE=yes -D CMAKE_INSTALL_PREFIX=/usr/
cmake --build . --parallel 16
sudo make install
```

## Setup

1. Set the environment variable ´LAMMPS_DIR´ to the directory of the LAMMPS repository
   (i.e. it should have subdirectories src, build, etc.)

## Force computation

Do each atom individually. If we can get the ordered list of atoms of a molecule,
e.g. via the neighbor list of the reference atom, then just go over this and compute
the differences again each time. Via constant folding we can optimize this a bit.

## Remarks about LAMMPS

- There are local atom indices and atom tags. For access of any atom data, e.g. position or force, the index is required. The `bond_atoms` array stores the tags.

- Even with one processor, there are ghost atoms and more bonds than there should be. This is due to the periodicity of the simulation box. If a bond crosses the boundary, it does not make sense to use the coordinates lying inside the box for both atoms, because one will be very small and the other very large. Thus (this is how I think it works, not 100% sure) LAMMPS creates two bonds and ghost atoms for every bond crossing the boundary, so the coordinates at the indices given in the bound are correct (relative to each other). As a consequence, one should not use the `bond_atoms` array plus `atom->map` to get the coordinates of the bond partner.