using NetworkInference
using Test
using DelimitedFiles

# These tests use a dataset generated from the 10-node Yeast1 network from http://gnw.sourceforge.net/
# GeneNetWeaver: In silico benchmark generation and performance profiling of network inference methods.
# Schaffter T, Marbach D, and Floreano D. Bioinformatics, 27(16):2263-70, 2011.

println("Getting nodes...")
data_path = joinpath(dirname(@__FILE__), "data")
data_file_path = joinpath(data_path, "yeast1_10_data.txt")
nodes = get_nodes(data_file_path)

println("Inferring networks...")
mi_network = InferredNetwork(MINetworkInference(), nodes)
clr_network = InferredNetwork(CLRNetworkInference(), nodes)
puc_network = InferredNetwork(PUCNetworkInference(), nodes)
pidc_network = InferredNetwork(PIDCNetworkInference(), nodes)

# The benchmarks were inferred prior to the implementation of this package. MI, PUC and PIDC were
# inferred using scripts based on InformationMeasures.jl and CLR was inferred using MINET:
# https://bioconductor.org/packages/release/bioc/html/minet.html

mi_benchmark = readdlm(joinpath(data_path, "mi.txt"))
clr_benchmark = readdlm(joinpath(data_path, "clr.txt"))
puc_benchmark = readdlm(joinpath(data_path, "puc.txt"))
pidc_benchmark = readdlm(joinpath(data_path, "pidc.txt"))

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

# Test infer_network
mi_network = infer_network(data_file_path, MINetworkInference())
@test mi_network.edges[40].weight ≈ mi_benchmark[2*40, 3] atol = 0.0001
println("infer_network passed")

# Test get_adjacency_matrix
present_node1, present_node2 = mi_network.edges[40].nodes
absent_node1, absent_node2 = mi_network.edges[41].nodes
adj_matrix_absolute, labels_to_ids, ids_to_labels = get_adjacency_matrix(mi_network, mi_network.edges[40].weight, absolute = true)
adj_matrix_percentage, labels_to_ids, ids_to_labels = get_adjacency_matrix(mi_network, 40 / length(mi_network.edges))
for adj_matrix in [adj_matrix_absolute, adj_matrix_percentage]
    @test adj_matrix[labels_to_ids[present_node1.label], labels_to_ids[present_node2.label]] === true
    @test adj_matrix[labels_to_ids[present_node2.label], labels_to_ids[present_node1.label]] === true
    @test adj_matrix[labels_to_ids[absent_node1.label], labels_to_ids[absent_node2.label]] === false
    @test adj_matrix[labels_to_ids[absent_node2.label], labels_to_ids[absent_node1.label]] === false
end
println("get_adjacency_matrix passed")

# Test read_network_file
test_network_filepath = joinpath(data_path, "mi.txt")
read_network = read_network_file(test_network_filepath)
@test length(read_network.edges) == 45
reference_nodes = ["G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "10"]
@test Set([x.label for x in read_network.nodes]) == Set(reference_nodes)
# Compare a few edges
@test read_network.edges[1].weight ≈ 0.2972651133683284 atol = 0.00001
@test Set([x.label for x in read_network.edges[1].nodes]) == Set(["G1", "G4"])
@test read_network.edges[15].weight ≈ 0.07202931598519413 atol = 0.00001
@test Set([x.label for x in read_network.edges[15].nodes]) == Set(["G3", "G9"])
@test read_network.edges[34].weight ≈ 0.02896110920075312 atol = 0.00001
@test Set([x.label for x in read_network.edges[34].nodes]) == Set(["G2", "G8"])
@test read_network.edges[35].weight ≈ 0.02848989653696199 atol = 0.00001
@test Set([x.label for x in read_network.edges[35].nodes]) == Set(["G9", "10"])
println("read_network_file passed")

# These tests will only run if the EmpiricalBayes package exists:
include("empirical_bayes_glue_tests.jl")
