# NetworkInference

[![Build Status](https://travis-ci.org/Tchanders/NetworkInference.jl.svg?branch=master)](https://travis-ci.org/Tchanders/NetworkInference.jl)
[![codecov.io](http://codecov.io/github/Tchanders/NetworkInference.jl/coverage.svg?branch=master)](http://codecov.io/github/Tchanders/NetworkInference.jl?branch=master)

## Description

NetworkInference is a package for inferring (undirected) networks, given a set of measurements for each node. The main output is the `InferredNetwork` type, which represents a fully connected, weighted network, where an edge's weight indicates the relative confidence of that edge existing in the true network. See also [Scope](#scope).

Some things to note:
* The package was originally written for inferring biological networks using gene expression data, hence the use of "network" instead of "graph". However, these methods could be applied to other types of data.
* Four network inference algorithms are currently implemented (MI, CLR, PUC and PIDC, explained in [[1]](#references)), but we plan to include more.
* Networks are assumed to be __undirected__, since all the algorithms included so far infer undirected networks. Hence:
	* in the `Edge` type, the order of the nodes is arbitrary
	* when a network is written to file, edges are written in both directions, becuase downstream analyses sometimes require this
	* the `InferredNetwork` type contains a list of edges, with one edge for each pair of genes, in which the order of the genes is arbitrary

## Installation

`Pkg.add("NetworkInference")`

## Basic usage

First include the package at the start of your script or interactive session:

`using NetworkInference`

### One step

Given a data file and an inference algorithm, you can infer a network with a single function call:

`infer_network(<path to data file>, PIDCNetworkInference())`

This will return an `InferredNetwork` type. You can also write the inferred network to file, using the `out_file_path` keyword argument. See also [Options](#options).

### Multiple steps

First make an array of `Node`s from your data:

`nodes = get_nodes(<path to data file>)`

Currently the package assumes the file is of the format:
* line 1: headers (these are discarded for now)
* other lines: NodeLabel value1 value2 value3 ...

Then infer a network:

`inferred_network = InferredNetwork(PIDCNetworkInference(), nodes)`

An `InferredNetwork` has an array of nodes and an array of edges between all possible node pairs (sorted in descending order of edge weight, i.e. confidence of the edge existing in the true network).

You can write the network to file:

`write_network_file(<path to output file>, inferred_network)`

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

**out_file_path** (`String`) Path to the output network file
* `""` (default) No file will be written

Defaults for **discretizer** and **estimator** are explained in [[1]](#references)

## Scope

This package is not designed for analysing networks/graphs or calculating network/graph metrics. In order to do such analyses, another package should be used (e.g. [LightGraphs](https://github.com/JuliaGraphs/LightGraphs.jl)). Of course, the edge list or the `InferredNetwork` will need to be parsed into the appropriate data structure first; the method `get_adjacency_matrix` may help with this.

Note that the `InferredNetwork` type contains a list of every possible edge, and the confidence of each edge existing in the true network. For analysing the properties of an inferred network, you may first want to define a partially connected, unweighted network by classifying each edge as "in the network" or "not in the network", based on the confidences. The simplest ways to do this are either to decide that the top x percent of edges are "in the network", or to define a threshold confidence, above which edges are "in the network".

You can pass a threshold into `get_adjacency_matrix` to get the adjacency matrix of a thresholded network (as well as dictionaries to map the node labels to their numerical IDs within the matrix, and vice versa):

`get_adjacency_matrix(inferred_network, 0.1) # Keeps top 10% edges with the largest weights`

`get_adjacency_matrix(inferred_network, 0.1, absolute = true) # Keeps all edges with weights >= 0.1`

## Performance

It may be possible to speed up an analysis, particularly for large datasets, by using [multiple processes](https://docs.julialang.org/en/v1/manual/distributed-computing/).

If multiple processes are available, NetworkInference will distribute the most costly calculations across the processes. (These are the for loops in `get_mi_scores` and `get_puc_scores`.)

**Example**

```
$ ./julia -p 3

julia> using NetworkInference

julia> infer_network(<path to data file>, PIDCNetworkInference())
```

This opens the Julia REPL with 3 extra processes (so 4 in total). NetworkInference may then be used as normal; it will handle distributing the calculations.

Note that the performance gain from distributing calculations is offset by communicating between the processes, so for small datasets it is more efficient to use one process. For the same reason, using too many processes will degrade performance, so it is a good idea to do some timing tests with different numbers of processes.

## Contributing

Bug reports, pull requests and other contributions are welcome!

## References

[1] Chan, Stumpf and Babtie (2017) [Gene Regulatory Network Inference from Single-Cell Data Using Multivariate Information Measures](http://www.cell.com/cell-systems/fulltext/S2405-4712(17)30386-1) Cell Systems
