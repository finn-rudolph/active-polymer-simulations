```
cmake -S ../cmake -B . -D BUILD_SHARED_LIBS=yes -D PKG_MOLECULE=yes -D CMAKE_INSTALL_PREFIX=/usr/
cmake --build . --parallel 16
sudo make install
```

## Setup

1. Set the environment variable ´LAMMPS_DIR´ to the directory of the LAMMPS repository
   (i.e. it should have subdirectories src, build, etc.)
