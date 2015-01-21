#! /bin/bash

# YAML

curl -O https://yaml-cpp.googlecode.com/files/yaml-cpp-0.5.1.tar.gz
tar xf yaml-cpp-0.5.1.tar.gz

cd yaml-cpp-0.5.1
mkdir build
cd build

cmake -DBoost_NO_BOOST_CMAKE=TRUE -DYAML_CPP_BUILD_TOOLS=OFF -DYAML_CPP_BUILD_CONTRIB=OFF -DCMAKE_INSTALL_PREFIX:PATH=../../ ..

make -j4
make install

cd ../..
rm yaml-cpp-0.5.1.tar.gz

# TCLAP
curl -O http://optimate.dl.sourceforge.net/project/tclap/tclap-1.2.1.tar.gz
tar xf tclap-1.2.1.tar.gz

cd tclap-1.2.1
./configure --prefix=$PWD/../

make -j4
make install

cd ..
rm tclap-1.2.1.tar.gz

# LHAPDF
curl -O -L http://www.hepforge.org/archive/lhapdf/LHAPDF-6.1.5.tar.gz
tar xf LHAPDF-6.1.5.tar.gz

cd LHAPDF-6.1.5
./configure --prefix=$PWD/..

make -j20
make install

cd ..
rm LHAPDF-6.1.5.tar.gz

cd share/LHAPDF/
curl -O -L http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10.tar.gz
tar xf CT10.tar.gz
rm CT10.tar.gz

curl -O -L http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10as.tar.gz
tar xf CT10as.tar.gz
rm CT10as.tar.gz

curl -O -L http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nlo.tar.gz
tar xf CT10nlo.tar.gz
rm CT10nlo.tar.gz

curl -O -L http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nlo_as_0117.tar.gz
tar xf CT10nlo_as_0117.tar.gz
rm CT10nlo_as_0117.tar.gz

curl -O -L http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nlo_as_0119.tar.gz
tar xf CT10nlo_as_0119.tar.gz
rm CT10nlo_as_0119.tar.gz

curl -O -L http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo.tar.gz
tar xf CT10nnlo.tar.gz
rm CT10nnlo.tar.gz

curl -O -L http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0117.tar.gz
tar xf CT10nnlo_as_0117.tar.gz
rm CT10nnlo_as_0117.tar.gz

curl -O -L http://www.hepforge.org/archive/lhapdf/pdfsets/6.1/CT10nnlo_as_0119.tar.gz
tar xf CT10nnlo_as_0119.tar.gz
rm CT10nnlo_as_0119.tar.gz

cd ../../

# Json cpp
svn checkout svn://svn.code.sf.net/p/jsoncpp/code/trunk jsoncpp-code

cd jsoncpp-code/jsoncpp
mkdir -p ../build/release
cd ../build/release

cmake -DCMAKE_BUILD_TYPE=release -DCMAKE_INSTALL_PREFIX:PATH=$PWD/../../../ -DJSONCPP_WITH_TESTS=OFF -G "Unix Makefiles" ../../jsoncpp
make -j4
make install

cd ../../..
