# NetworkInference

[![Build Status](https://travis-ci.org/Tchanders/NetworkInference.jl.svg?branch=master)](https://travis-ci.org/Tchanders/NetworkInference.jl)

[![Coverage Status](https://coveralls.io/repos/Tchanders/NetworkInference.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/Tchanders/NetworkInference.jl?branch=master)

[![codecov.io](http://codecov.io/github/Tchanders/NetworkInference.jl/coverage.svg?branch=master)](http://codecov.io/github/Tchanders/NetworkInference.jl?branch=master)

NB This package is still under development and will probably change significantly. Only one network inference algorithm is currently implemented (PIDC, explained in http://biorxiv.org/content/early/2017/04/26/082099), but we plan to include more.

## Installation

`Pkg.clone("git@github.com:Tchanders/NetworkInference.jl.git")`

## Basic usage

First make an array of `Gene`s from your data:

`genes = get_genes(path_to_data_file)`

Currently the package assumes the file is of the format:
line 1: headers (these are discarded for now)
other lines: GeneName value1 value2 value3 ...

Then get a network analysis:

`network_analysis = NetworkAnalysis(PIDCNetworkAnalysis(), genes)`

A `NetworkAnalysis` has an array of genes and an array of edges (sorted in descending order of confidence).
