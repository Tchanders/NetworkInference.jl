# Glue functions for natively using EmpiricalBayes with NetworkInference structures. 
# Will not load if EmpiricalBayes does not exist. 


# Check for EmpiricalBayes package 
try 
  ver = Pkg.installed("EmpiricalBayes") # Will throw errors if does not exist 
  using EmpiricalBayes
  eb_exists = true
catch 
  eb_exists = false 
end  

# Only load if package exists
if eb_exists


"""
  to_index(label1::String, label2::String)

Convert any ordering of two strings into a unique tuple for that pair.
"""
function to_index(label1::String, label2::String)
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
  make_priors(filepath::String, skiplines=1)

Convert a file containing prior data into a dictionary of (id, prior) pairs.

# Arguments
- `filepath` : path to prior file. File should contain data in three columns as
  follows: node1, node2, prior.
- `skiplines=1` : number of initial lines to skip in the file. Defaults to 1, to
  skip a header line.
"""
function make_priors(filepath::String, skiplines=1)
  prior_mat = readdlm(filepath, skipstart=skiplines)

  num_priors = size(prior_mat, 1)

  prior_dict = Dict()
  for i in 1:num_priors
    n1, n2, p = prior_mat[i, :]
    index = to_index(n1, n2)
    prior_dict[index] = p
  end

  return prior_dict
end


"""
  empirical_bayes(network::InferredNetwork, priors::Dict, key_func=to_index, num_bins::Int, proportion_to_keep::Float64=1.0)

Calculate the empirical Bayes posteriors of the input statistics using the priors.

# Arguments
- `network::InferredNetwork` : network to apply priors to.
- `priors::Dict` : dictionary of priors such that looking up a pair of nodes
  returns the prior value for the edge between them.
- `key_func=to_index` : a function mapping an Edge object to a key useable in
   the `priors` dictionary
- `num_bins::Int` : number of uniform width bins to discretize into.
- `proportion_to_keep::Float64=1.0` : Proportion of lowest test statistics to
   keep when calculating null distribution.
"""
function empirical_bayes(network::InferredNetwork, priors::Dict, key_func=to_index, num_bins::Int, proportion_to_keep::Float64=1.0)
  edge_list = get_edge_list(network)
  test_statistics = [e.weight for e in edge_list]
  prior_list = [ get(priors, key_func(e.nodes), 0) for e in edge_list]

  posteriors = empirical_bayes(test_statistics, priors, num_bins, proportion_to_keep)

  for i in 1:length(edge_list)
    nodes = edge_list[i].nodes
    edge_list[i] = Edge(nodes, posteriors[i])
  end

  return InferredNetwork(network.nodes, edge_list)
end


end 