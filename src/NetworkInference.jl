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
    PIDNetworkInference,
    PIDCNetworkInference,
    MINetworkInference

include("common.jl")
include("network_inference.jl")

end # module
