immutable Gene
    name::String
    discretized_values::Array{Int64}
    number_of_bins::Int64 # Is this the best type?
    probabilities::Array{Float64}
end

# Make a Gene from a data file line
# TODO: allow choose discretizer?
function Gene(line::AbstractArray)

    name = string(line[1])
    expression_values = Array{Float64}(line[2:end])
    discretized_values = zeros(Int, length(expression_values))
    number_of_bins = get_bin_ids!(expression_values, "bayesian_blocks", 10, discretized_values)
    probabilities = get_probabilities("maximum_likelihood", get_frequencies_from_bin_ids(discretized_values, number_of_bins))

    return Gene(name, discretized_values, number_of_bins, probabilities)

end

immutable GenePair
    mi::Float64 # Mutual information
    si::Array{Float64} # Specific information
end

# Networks are assumed to be undirected for now: order of the genes is meaningless
immutable Edge
    genes::Array{Gene}
    confidence::Float64
end

immutable Network
    genes::Array{Gene}
    edges::Array{Edge}
end

immutable NetworkAnalysis
    genes::Array{Gene}
    edges::Array{Edge} # Edges in descending order of confidence
end

# TODO: Allow choose discretizer?
# TODO: Allow choose estimator?
function get_genes(data_file_path::String)

    lines = readdlm(open(data_file_path), skipstart = 1) # Assumes the first line is headers
    number_of_genes = size(lines, 1)
    genes = Array{Gene}(number_of_genes)

    for i in 1:number_of_genes
        genes[i] = Gene(lines[i:i, 1:end])
    end

    return genes

end

# Networks are assumed to be undirected for now: edges are written in both directions with the same confidence
function write_network_file(file_name::String, network_analysis::NetworkAnalysis)

    out_file = open(file_name, "w")

    for edge in network_analysis.edges
        genes = collect(edge.genes)
        write(out_file, string(
            genes[1].name, "\t", genes[2].name, "\t",
            edge.confidence, "\n",
            genes[2].name, "\t", genes[1].name, "\t",
            edge.confidence, "\n"
        ))
    end

    close(out_file)

end