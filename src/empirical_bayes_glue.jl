# Glue functions for natively using EmpiricalBayes with NetworkInference structures.
# Will not load if EmpiricalBayes does not exist.

# Check for EmpiricalBayes package
EB_EXISTS = in("EmpiricalBayes",keys(Pkg.installed())) ? true : false

# Only load if package exists
if EB_EXISTS

import EmpiricalBayes.empirical_bayes # for overloading

"""
    to_index(label1::String, label2::String)

Convert any ordering of two strings into a unique tuple for that pair.
"""
function to_index(label1::AbstractString, label2::AbstractString)
    if label1 > label2
        return (label1, label2)
    else
        return (label2, label1)
    end
end

"""
    to_index(node1::Node, node2::Node)

Convert any ordering of two Node objects into a unique tuple for that pair.
"""
function to_index(node1::Node, node2::Node)
    return to_index(node1.label, node2.label)
end

"""
    to_index(nodes::Array{Node, 1})

Convert any ordering of two Node objects into a unique tuple for that pair.
"""
function to_index(nodes::Array{Node, 1})
    n1, n2 = nodes
    return to_index(n1.label, n2.label)
end

"""
    make_priors(filepath::String, weight = 1.0, skiplines = 1)

Convert a file containing prior data into a dictionary of (id, prior) pairs.

# Arguments
* `filepath` : path to prior file. File should contain data in three columns as
  follows: node1, node2, prior.
* `weight` : weight to apply to the prior.
* `skiplines = 1` : number of initial lines to skip in the file. Defaults to 1, to
  skip a header line.
"""
function make_priors(filepath::String, weight = 1.0, skiplines = 1)
    prior_mat = readdlm(filepath, skipstart = skiplines)

    num_priors = size(prior_mat, 1)

    prior_dict = Dict()
    for i in 1:num_priors
        n1, n2, p = prior_mat[i, :]
        n1 = string(n1)
        n2 = string(n2)
        index = to_index(n1, n2)
        if haskey(prior_dict, index)
            prior_dict[index] = max(p * weight, prior_dict[index])
        else
            prior_dict[index] = p * weight
        end
    end

    return prior_dict
end

"""
    make_priors(filepaths::Array{String}, weights = ones(length(filepaths)), skiplines = 1)

Convert a file containing prior data into a dictionary of (id, prior) pairs.

# Arguments
* `filepaths` : paths to prior files. Each file should contain data in three columns as
  follows: node1, node2, prior.
* `weights` : weights to apply to each prior, in the same order. Defaults to array of ones.
* `skiplines = 1` : number of initial lines to skip in the file. Defaults to 1, to
  skip a header line.
"""
function make_priors(filepaths::Array{String}, weights = ones(length(filepaths)), skiplines = 1)
    num_filepaths = length(filepaths)

    # Check there is a weight for each prior file; if not, fall back to default
    weights = length(weights) == num_filepaths ? weights : ones(num_filepaths)

    prior_dicts = make_priors.(filepaths, weights, skiplines)
    prior_dict = merge(+, prior_dicts...)

    return prior_dict
end

"""
    empirical_bayes(network::InferredNetwork, priors::Dict, num_bins, distr::Symbol, proportion_to_keep = 1.0, key_func = to_index)

Calculate the empirical Bayes posteriors of the input statistics using the priors.

# Arguments
* `network::InferredNetwork` : network to apply priors to.
* `priors::Dict` : dictionary of priors such that looking up a pair of nodes
  returns the prior value for the edge between them.
* `num_bins` : number of uniform width bins to discretize into.
* `distr` : form of the null distribution to be fitted.
* `proportion_to_keep = 1.0` : Proportion of lowest test statistics to
   keep when calculating null distribution.
* `key_func = to_index` : a function mapping an Edge object to a key useable in
   the `priors` dictionary.
* `tail` : Whether the test is two-tailed (:two) or one-tailed (:lower or :upper).
- `w0` : Default constant for the prior calculation
"""
function empirical_bayes(network::InferredNetwork, priors::Dict, num_bins, distr::Symbol; proportion_to_keep = 1.0, key_func = to_index,
    tail = :upper, w0 = 2.2)

    edge_list = network.edges
    test_statistics = [e.weight for e in edge_list]
    prior_list = [ get(priors, key_func(e.nodes), 0) for e in edge_list ]

    eb_edges = Array{Edge}(undef, length(edge_list))

    posteriors = empirical_bayes(test_statistics, prior_list, num_bins, distr, proportion_to_keep = proportion_to_keep,
        tail = tail, w0 = w0)

    for i in 1:length(edge_list)
        nodes = edge_list[i].nodes
        eb_edges[i] = Edge(nodes, posteriors[i])
    end

    # Remove infinite values
    eb_edges = filter(x->isfinite(x.weight), eb_edges)

    sort!(eb_edges, rev = true, by = x->x.weight)

    return InferredNetwork(network.nodes, eb_edges)
end

"""
  empirical_bayes(network::InferredNetwork num_bins, distr::Symbol, proportion_to_keep = 1.0, key_func = to_index)

Calculate the empirical Bayes posteriors of the input statistics with no priors.

# Arguments
* `network::InferredNetwork` : network to apply priors to.
* `num_bins` : number of uniform width bins to discretize into.
* `distr` : form of the null distribution to be fitted.
* `proportion_to_keep = 1.0` : proportion of lowest test statistics to
   keep when calculating null distribution.
* `key_func = to_index` : a function mapping an Edge object to a key useable in
   the `priors` dictionary.
* `tail` : whether the test is two-tailed (:two) or one-tailed (:lower or :upper).
- `w0` : Default constant for the prior calculation
"""
function empirical_bayes(network::InferredNetwork, num_bins, distr::Symbol; proportion_to_keep = 1.0, key_func = to_index, tail = :upper, w0 = 2.2)
    return empirical_bayes(network, Dict(), num_bins, distr, proportion_to_keep = proportion_to_keep, key_func = key_func, tail = tail, w0 = w0)
end

end
