```
cmake -S ../cmake -B . -D BUILD_SHARED_LIBS=yes -D PKG_MOLECULE=yes -D CMAKE_INSTALL_PREFIX=/usr/
cmake --build . --parallel 16
sudo make install
```
