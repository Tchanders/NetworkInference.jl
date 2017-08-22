immutable Gene
    name::String
    discretized_values::Array{Int64}
    number_of_bins::Int64
    probabilities::Array{Float64}
end

# Make a Gene from a data file line
# discretzed_values is expression_values mapped to bin ids, following discretization
function Gene(line::AbstractArray, discretizer, estimator, number_of_bins)

    name = string(line[1])
    expression_values = Array{Float64}(line[2:end])
    discretized_values = zeros(Int, length(expression_values))
    # If discretizer is Bayesian blocks, number_of_bins will be overwritten
    number_of_bins = get_bin_ids!(expression_values, discretizer, number_of_bins, discretized_values)
    probabilities = get_probabilities(estimator, get_frequencies_from_bin_ids(discretized_values, number_of_bins))

    return Gene(name, discretized_values, number_of_bins, probabilities)

end

immutable GenePair
    mi::Float64 # Mutual information
    si::Array{Float64} # Specific information
end

# Networks are assumed to be undirected for now: order of the genes is arbitrary
immutable Edge
    genes::Array{Gene}
    confidence::Float64
end
