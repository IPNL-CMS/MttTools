#! /bin/bash

# YAML

wget --no-check-certificate https://yaml-cpp.googlecode.com/files/yaml-cpp-0.5.1.tar.gz
tar xf yaml-cpp-0.5.1.tar.gz

cd yaml-cpp-0.5.1
mkdir build
cd build

cmake -DYAML_CPP_BUILD_TOOLS=OFF -DYAML_CPP_BUILD_CONTRIB=OFF -DCMAKE_INSTALL_PREFIX:PATH=../../ ..

make -j4
make install

cd ../..
rm yaml-cpp-0.5.1.tar.gz
