using NetworkInference
using Base.Test

# These tests use a dataset generated from the 10-node Yeast1 network from http://gnw.sourceforge.net/
# GeneNetWeaver: In silico benchmark generation and performance profiling of network inference methods.
# Schaffter T, Marbach D, and Floreano D. Bioinformatics, 27(16):2263-70, 2011.

println("Getting nodes...")
nodes = get_nodes("data/yeast1_10_data.txt")

println("Inferring networks...")
mi_network = InferredNetwork(MINetworkInference(), nodes)
clr_network = InferredNetwork(CLRNetworkInference(), nodes)
puc_network = InferredNetwork(PUCNetworkInference(), nodes)
pidc_network = InferredNetwork(PIDCNetworkInference(), nodes)

# The benchmarks were inferred prior to the implementation of this package. MI, PUC and PIDC were
# inferred using scripts based on InformationMeasures.jl and CLR was inferred using MINET:
# https://bioconductor.org/packages/release/bioc/html/minet.html

mi_benchmark = readdlm("data/mi.txt")
clr_benchmark = readdlm("data/clr.txt")
puc_benchmark = readdlm("data/puc.txt")
pidc_benchmark = readdlm("data/pidc.txt")

# Compare a few selected edges throughout the inferred network
for i in (1, 5, 10, 20, 40)
    @test mi_network.edges[i].weight ≈ mi_benchmark[2*i, 3] atol = 0.0001
    println("MI network inference $i passed")
    @test clr_network.edges[i].weight ≈ clr_benchmark[2*i, 3] atol = 0.0001
    println("CLR network inference $i passed")
    @test puc_network.edges[i].weight ≈ puc_benchmark[2*i, 3] atol = 0.0001
    println("PUC network inference $i passed")
    @test pidc_network.edges[i].weight ≈ pidc_benchmark[2*i, 3] atol = 0.0001
    println("PIDC network inference $i passed")
end
