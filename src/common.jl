immutable Gene
    name::String
    discretized_values::Array{Int64}
    number_of_bins::Int64 # Is this the best type?
end

# Make a Gene from a data file line
# TODO: allow choose discretizer
function Gene(line::Array{SubString{String}, 1})

    name = String(shift!(line))

    expression_values = map((x) -> parse(Float64, x), line)
    discretized_values = zeros(expression_values)
    number_of_bins = get_bin_ids!(expression_values, "bayesian_blocks", 10, discretized_values)

    return Gene(name, discretized_values, number_of_bins)

end

# TODO: Think about directedness; Set assumes the edge is undirected
immutable Edge
    genes::Set{Gene}
    confidence::Float64
end

immutable Network
    edges::Set{Edge}
    genes::Set{Gene}
end

immutable NetworkAnalysis
    genes::Array{Gene}
    edges_by_confidence::Array{Edge} # Edges in descending order of confidence
end

# TODO: Different discretization and estimation options
function get_genes(data_file_path::String)

    lines = readlines(open(data_file_path))
    shift!(lines) # In the future, every line should be GeneName expression1 expression2 ...
    genes = Array{Gene}(length(lines))    

    for (i, line) in enumerate(lines)
        genes[i] = Gene(split(line))
    end

    return genes

end

# TODO: This assumes the edges are undirected
function write_network_file(file_name::String, network_analysis::NetworkAnalysis)

    out_file = open(file_name, "w")

    for edge in network_analysis.edges_by_confidence
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