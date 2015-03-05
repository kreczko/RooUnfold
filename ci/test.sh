#!/bin/bash
# This script is meant to be called by the "script" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.
set -e
# Check if root and PyROOT work ok
echo "Checking ROOT installation"
time root -l -q
echo "Checking PyROOT installation"
time python -c 'import ROOT; ROOT.TBrowser()'

# Check that RooUnfold can be imported
echo "Checking RooUnfold library in PyROOT"
time python -c 'from ROOT import gSystem; gSystem.Load("libRooUnfold.so");from ROOT import RooUnfoldResponse'
time python test/test_pyroot.py
echo "Executing unit tests"
make test
