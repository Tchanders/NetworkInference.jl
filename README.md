# NetworkInference

[![Build Status](https://travis-ci.org/Tchanders/NetworkInference.jl.svg?branch=master)](https://travis-ci.org/Tchanders/NetworkInference.jl)

[![Coverage Status](https://coveralls.io/repos/Tchanders/NetworkInference.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/Tchanders/NetworkInference.jl?branch=master)

[![codecov.io](http://codecov.io/github/Tchanders/NetworkInference.jl/coverage.svg?branch=master)](http://codecov.io/github/Tchanders/NetworkInference.jl?branch=master)

NB This package is still under development and could change significantly.

## Description

NetworkInference is a package for inferring (undirected) networks, given a set of measurements for each node.

Some things to note:
* The package was originally written for inferring biological networks using gene expression data, hence type names such as `Gene` and the use of "network" instead of "graph". However, these methods could be applied to other types of data.
* Four network inference algorithms are currently implemented (MI, CLR, PUC and PIDC, explained in http://biorxiv.org/content/early/2017/04/26/082099), but we plan to include more.
* Networks are assumed to be __undirected__, since all the algorithms included so far infer undirected networks. Hence:
	* in the `Edge` type, the order of the genes is arbitrary
	* when a network is written to file, edges are written in both directions, becuase downstream analyses sometimes require this
	* `infer_network` returns an edge list where the edges are only written in one (arbitrary) direction, to save space
* Inferred networks consist of all possible pairs of genes, and an edge score for each pair. See also [Scope](#scope).

## Installation

`Pkg.clone("git@github.com:Tchanders/NetworkInference.jl.git")`

## Basic usage

### One step

Given a data file and an inference algorithm, you can infer a network with a single function call:

`infer_network(<path to data file>, PIDCNetworkInference())`

This will return an edge list (of type `Array{Tuple{String,String,Float64},1}`). You can also write the inferred network to file, using the `out_file` keyword argument. See also [Options](#options).

### Multiple steps

First make an array of `Gene`s from your data:

`genes = get_genes(<path to data file>)`

Currently the package assumes the file is of the format:
* line 1: headers (these are discarded for now)
* other lines: GeneName value1 value2 value3 ...

Then infer a network:

`network_analysis = NetworkAnalysis(PIDCNetworkInference(), genes)`

A `NetworkAnalysis` has an array of genes and an array of edges between all possible gene pairs (sorted in descending order of edge score, a.k.a. confidence).

You can write the network to file:

`write_network_file(<path to output file>, network_analysis)`

You can get the network as an edge list:

`get_edge_list(network_analysis)`

## Options

The following keyword arguments can be passed in to `infer_network`:

**delim** (`Union{Char,Bool}`) Column delimiter
* `false` (default) Delimiter is whitespace

**discretizer** (`String`) Method for discretizing
* `"bayesian_blocks"` (default) Adaptive discretizer with varibale number of bins
* `"uniform_width"` Use this if Bayesian blocks fails, or if constant number of bins is required
* `"uniform_count"`

**estimator** (`String`) Estimator for estimating the probability distribution
* `"maximum_likelihood"` (default) Highly recommended for PUC and PIDC. (See inline comments for more information.)
* `"dirichlet"`
* `"shrinkage"`

**number_of_bins** (`Integer`)
* `10` (default)
(This will be ignored when using Bayesian blocks discretization)

**base** (`Number`) Base of the logarithm, i.e. the units for entropy
* `2` (default)

**out_file** (`String`) Path to the output network file
* `""` (default) No file will be written

NB **discretizer** and **estimator** defaults are explained in http://biorxiv.org/content/early/2017/04/26/082099

## Scope

This package is not designed for analysing networks/graphs or calculating network/graph metrics. In order to do such analyses, another package should be used (e.g. [LightGraphs](https://github.com/JuliaGraphs/LightGraphs.jl)). Of course, the edge list or the `NetworkAnalysis` will need to be parsed into the appropriate data structure first.

## Contributing

Bug reports, pull requests and other contributions are welcome!
