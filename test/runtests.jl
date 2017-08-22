using NetworkInference
using Base.Test

# These tests use a dataset generated from the 10-node Yeast1 network from http://gnw.sourceforge.net/
# GeneNetWeaver: In silico benchmark generation and performance profiling of network inference methods.
# Schaffter T, Marbach D, and Floreano D. Bioinformatics, 27(16):2263-70, 2011.

println("Getting genes...")
genes = get_genes("data/yeast1_10_data.txt")

println("Inferring networks...")
mi_network = NetworkAnalysis(MINetworkInference(), genes)
clr_network = NetworkAnalysis(CLRNetworkInference(), genes)
puc_network = NetworkAnalysis(PUCNetworkInference(), genes)
pidc_network = NetworkAnalysis(PIDCNetworkInference(), genes)

# The benchmarks were inferred prior to the implementation of this package. MI, PUC and PIDC were
# inferred using scripts based on InformationMeasures.jl and CLR was inferred using MINET:
# https://bioconductor.org/packages/release/bioc/html/minet.html

mi_benchmark = readdlm("data/mi.txt")
clr_benchmark = readdlm("data/clr.txt")
puc_benchmark = readdlm("data/puc.txt")
pidc_benchmark = readdlm("data/pidc.txt")

# Compare the edges with the top 5 highest confidences
# NB CLR differs more for lower-confidence edges, since the original implementation used a threshold
for i in (1, 2, 3, 4, 5)
    @test mi_network.edges[i].confidence ≈ mi_benchmark[2*i, 3] atol = eps()
    println("MI network inference $i passed")
    @test clr_network.edges[i].confidence ≈ clr_benchmark[2*i, 3] atol = eps()
    println("CLR network inference $i passed")
    @test puc_network.edges[i].confidence ≈ puc_benchmark[2*i, 3] atol = eps()
    println("PUC network inference $i passed")
    @test pidc_network.edges[i].confidence ≈ pidc_benchmark[2*i, 3] atol = eps()
    println("PIDC network inference $i passed")
end
