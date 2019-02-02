using NetworkInference
using Test
using DelimitedFiles

# Only run tests if the EmpiricalBayes package exists
if NetworkInference.EB_EXISTS

using EmpiricalBayes

println("Testing empirical Bayes glue functions...")

data_folder_path = joinpath(dirname(@__FILE__), "data")

# Test to_index functions
@test to_index("bbb", "aaa") == to_index("aaa", "bbb")
@test to_index("aaa", "bbb") == ("bbb", "aaa")

n1 = Node("aaa", [], 0, [])
n2 = Node("bbb", [], 0, [])
@test to_index(n1, n2) == to_index(n2, n1)
@test to_index(n1, n2) == ("bbb", "aaa")

@test to_index([n1, n2]) == to_index([n2, n1])
@test to_index([n1, n2]) == ("bbb", "aaa")

# Test the make_priors function
prior_path = joinpath(data_folder_path, "test_priors.txt")
reference_priors = Dict( [ (("bbb", "aaa"), 1) ,
                     (("ccc", "aaa"), 0) ,
                     (("ccc", "bbb"), 1)
                   ])
@test make_priors(prior_path) == reference_priors
reference_priors = Dict( [ (("ccc", "aaa"), 0) ,
                     (("ccc", "bbb"), 1) ,
                     (("bbb", "aaa"), 0)
                   ])
@test make_priors(prior_path, 1.0, 2) == reference_priors
reference_priors = Dict( [ (("bbb", "aaa"), 2) ,
                     (("ccc", "aaa"), 0) ,
                     (("ccc", "bbb"), 2)
                   ])
@test make_priors(prior_path, 2.0) == reference_priors
@test make_priors([prior_path, prior_path]) == reference_priors
@test make_priors([prior_path, prior_path], [1.5, 0.5]) == reference_priors
reference_priors = Dict( [ (("ccc", "aaa"), 0) ,
                     (("ccc", "bbb"), 2) ,
                     (("bbb", "aaa"), 0)
                   ])
@test make_priors([prior_path, prior_path], [1.5, 0.5], 2) == reference_priors

# Test the empirical_bayes function
# Test with w0 = 2.2, because that is the default for empirical_bayes within NetworkInference
println("Inferring test empirical Bayes networks...")
mi_benchmark = readdlm(joinpath(data_folder_path, "mi.txt"))
mi_benchmark = mi_benchmark[1:2:end, :] # skip repeated edges

mi_priors_filepath = joinpath(data_folder_path, "mi_priors.txt")
mi_priors = readdlm(mi_priors_filepath)
mi_priors = mi_priors[1:2:end, :] # skip repeated edges
prior_dict = make_priors(mi_priors_filepath)

yeast_test_data = joinpath(data_folder_path, "yeast1_10_data.txt")
nodes = get_nodes(yeast_test_data)
mi_network = InferredNetwork(MINetworkInference(), nodes)

eb_network = empirical_bayes(mi_network, prior_dict, 5, :Gamma, tail = :two, w0 = 2.2)
eb_weights = [e.weight for e in eb_network.edges]

benchmark_stats = convert(Array{Float64}, mi_benchmark[:, 3])
benchmark_priors = convert(Array{Float64}, mi_priors[:, 3])
ref_weights = empirical_bayes(benchmark_stats, benchmark_priors, 5, :Gamma, tail = :two, w0 = 2.2)

@test eb_weights ≈ sort(ref_weights, rev = true) atol = 0.0001

# Test the empirical_bayes function with no priors
println("Inferring test empirical Bayes networks with no priors...")
eb_no_prior_network = empirical_bayes(mi_network, 5, :Gamma, tail = :two, w0 = 2.2)
eb_no_prior_weights = [e.weight for e in eb_no_prior_network.edges]

ref_no_prior_weights = empirical_bayes(benchmark_stats, 5, :Gamma, tail = :two, w0 = 2.2)

@test eb_no_prior_weights ≈ sort(ref_no_prior_weights, rev = true) atol = 0.0001

println("Empirical Bayes glue tests passed")

end
