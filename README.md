# RooUnfold


Fork from http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html

[![build status](https://travis-ci.org/kreczko/RooUnfold.svg?branch=master)](https://travis-ci.org/kreczko/RooUnfold)


RooUnfold is a framework for unfolding (AKA "deconvolution" or "unsmearing").
It currently implements several methods: iterative (Bayes),
singular value decomposition (SVD), unregularised matrix inversion,
bin-by-bin (simple correction factors), and an interface to ROOT's TUnfold.
It can be used from the ROOT prompt, or linked against the ROOT libraries.
RooUnfold requires ROOT 4 or later.

For the latest version and documentation, see here

  http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html

RooUnfold was written by Tim Adye, Kerstin Tackmann, and Fergus Wilson.

## Building
Since version 2 RooUnfold is using cmake to find dependencies and create the MakeFile:
```shell
cmake CMakeLists.txt
make jobs=2
```

to build the RooUnfold shared library.

### Tests
In order to run the boost tests for RooUnfold simply
```shell
make test
```

If you want to run specific tests you can do this with
```shell
# a single test suite
./RooUnfold_test --log_level=message --run_test=TestRooUnfoldSvd
# a single test case
./RooUnfold_test --log_level=message --run_test=TestRooUnfoldSvd/test_get_tau_from_constructor
```

## Running

RooUnfoldExample.cxx makes a simple test of RooUnfold.

  % root
  root [0] .L examples/RooUnfoldExample.cxx
  root [1] gSystem->Load("libRooUnfold.so");
  root [2] RooUnfoldExample()

See the web page (http://hepunx.rl.ac.uk/~adye/software/unfold/RooUnfold.html) for more examples and documentation.
