#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

set -e

# Install a ROOT binary that we custom-built in a 64-bit Ubuntu VM
# for the correct python / ROOT version
[ -f root_v${ROOT}_python_2.7.tar.gz ] || time wget --no-check-certificate https://copy.com/s3BcYu1drmZa/ci/root_builds/root_v${ROOT}_python_2.7.tar.gz
time tar zxf root_v${ROOT}_python_2.7.tar.gz
mv root_v${ROOT}_python_2.7 root
source root/bin/thisroot.sh

# test ROOT install
# Check if ROOT
root -l -q

# setup newer compilers ( we need gcc >= 4.7 for c++11
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-${GCC_VERSION} 50;
sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-${GCC_VERSION} 50;
sudo update-alternatives --set gcc /usr/bin/gcc-${GCC_VERSION}; 
sudo update-alternatives --set g++ /usr/bin/g++-${GCC_VERSION};

sudo update-alternatives --install /usr/bin/clang clang /usr/bin/clang-${CLANG_VERSION} 50;
sudo update-alternatives --install /usr/bin/clang++ clang++ /usr/bin/clang++-${CLANG_VERSION} 50;
sudo update-alternatives --set clang /usr/bin/clang-${CLANG_VERSION};
sudo update-alternatives --set clang++ /usr/bin/clang++-${CLANG_VERSION};

export CC=/usr/bin/clang-${CLANG_VERSION};
export CXX=/usr/bin/clang++-${CLANG_VERSION};
# check version
echo "g++ version:"
g++ --version
echo "clang version:"
clang --version
which clang
echo "using CXX: ${CXX}"

ls -lah /usr/bin/ | grep clang

make
