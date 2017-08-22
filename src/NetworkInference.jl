module NetworkInference

using InformationMeasures
using Distributions

export
    # Common types and functions
    Gene,
    Edge,
    NetworkAnalysis,
    # Network inference algorithms
    AbstractNetworkInference,
    MINetworkInference,
    CLRNetworkInference,
    PUCNetworkInference,
    PIDCNetworkInference,
    # Functions for inferring networks
    get_genes,
    write_network_file,
    get_edge_list,
    infer_network

include("common.jl")
include("network_inference.jl")
include("infer_network.jl")

end # module
