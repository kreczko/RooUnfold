#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.
set -e
if [[ $TRAVIS_PYTHON_VERSION == 3* ]]; then
    export PYTHON_SUFFIX="3"
fi
# add repositories for gcc 4.8 and clang 3.5
sudo add-apt-repository --yes ppa:ubuntu-toolchain-r/test
sudo add-apt-repository --yes 'deb http://llvm.org/apt/precise/ llvm-toolchain-precise main'
wget -O - http://llvm.org/apt/llvm-snapshot.gpg.key | sudo apt-key add -
# Needed because sometimes travis' repositories get out of date
time sudo apt-get update -qq

# Install the dependencies we need
time sudo apt-get install -qq python${PYTHON_SUFFIX}-numpy python${PYTHON_SUFFIX}-sphinx python${PYTHON_SUFFIX}-nose clang-3.5 libclang-3.5-dev gcc-4.8 g++-4.8
# matplotlib and PyTables are not available for Python 3 as packages from the main repo yet.
if [[ $TRAVIS_PYTHON_VERSION == '2.7' ]]; then
  time sudo apt-get install -qq python${PYTHON_SUFFIX}-matplotlib python${PYTHON_SUFFIX}-tables
fi

# Install a ROOT binary that we custom-built in a 64-bit Ubuntu VM
# for the correct python / ROOT version
time wget --no-check-certificate https://copy.com/s3BcYu1drmZa/ci/root_builds/root_v${ROOT}_python_${TRAVIS_PYTHON_VERSION}.tar.gz
time tar zxf root_v${ROOT}_python_${TRAVIS_PYTHON_VERSION}.tar.gz
mv root_v${ROOT}_python_${TRAVIS_PYTHON_VERSION} root
source root/bin/thisroot.sh
# setup newer compilers for ROOT 6
if [[ $ROOT == 6* ]]; then
  sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 50
  sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.8 50
  sudo update-alternatives --set gcc /usr/bin/gcc-4.8
  sudo update-alternatives --set g++ /usr/bin/g++-4.8
fi

# compile RooUnfold
make
# add base path from setup_standalone to PYTHONPATH
export PYTHONPATH=$PYTHONPATH:$base
