declare -a files_to_copy=("fix_active_poly_force" "compute_config_moment" "compute_particle_diameter" "active_poly_util")

for file in "${files_to_copy[@]}"
do 
    cp $file.h  $LAMMPS_DIR/src/$file.h
    cp $file.cpp $LAMMPS_DIR/src/$file.cpp
done

cp active_poly_constants.h $LAMMPS_DIR/src/active_poly_constants.h

cd $LAMMPS_DIR/build
cmake -S ../cmake -B . -D BUILD_SHARED_LIBS=yes -D PKG_MOLECULE=yes -D CMAKE_INSTALL_PREFIX=/usr/
cmake --build . --parallel 16
sudo make install