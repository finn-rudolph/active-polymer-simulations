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