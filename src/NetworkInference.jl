module NetworkInference

using InformationMeasures
using Distributions
using Distributed
using Pkg
using DelimitedFiles
using SharedArrays

export
    # Common types
    Node,
    Edge,
    InferredNetwork,
    # Network inference algorithms
    AbstractNetworkInference,
    MINetworkInference,
    CLRNetworkInference,
    PUCNetworkInference,
    PIDCNetworkInference,
    # Functions for inferring networks
    get_nodes,
    write_network_file,
    read_network_file,
    get_adjacency_matrix,
    infer_network

include("common.jl")
include("network_inference.jl")
include("infer_network.jl")
include("empirical_bayes_glue.jl")

# Optional exports
if EB_EXISTS
export
    # Empirical Bayes glue functions
    to_index,
    make_priors,
    empirical_bayes
end

end # module
