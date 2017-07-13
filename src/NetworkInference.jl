module NetworkInference

using InformationMeasures
using Distributions

export
    # Common types and functions
    Gene,
    Edge,
    Network,
    NetworkAnalysis,
    get_genes,
    write_network_file,
    # Network inference algorithms
    AbstractNetworkInference,
    MINetworkInference,
    CLRNetworkInference,
    PUCNetworkInference,
    PIDCNetworkInference

include("common.jl")
include("network_inference.jl")

end # module
