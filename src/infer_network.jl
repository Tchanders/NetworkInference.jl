# The maximum likelihood estimator is recommended for PUC and PIDC, because speedups
# are made here, based on the assumption that the marginal probability distribution for
# a gene, from the joint distribution with any two other genes is always the same. If
# the joint distributions are estimated using other estimators, this assumption is
# violated for PUC and PIDC in get_puc and get_joint_probabilities.
function get_genes(data_file_path::String; delim::Union{Char,Bool} = false, discretizer = "bayesian_blocks",
    estimator = "maximum_likelihood", number_of_bins = 10)

    if delim == false
        lines = readdlm(open(data_file_path); skipstart = 1) # Assumes the first line is headers
    else
        lines = readdlm(open(data_file_path), delim; skipstart = 1) # Assumes the first line is headers
    end
    number_of_genes = size(lines, 1)
    genes = Array{Gene}(number_of_genes)

    for i in 1:number_of_genes
        genes[i] = Gene(lines[i:i, 1:end], discretizer, estimator, number_of_bins)
    end

    return genes

end

# Networks are undirected; edges are written in both directions with the same confidence
function write_network_file(file_name::String, network_analysis::NetworkAnalysis)

    out_file = open(file_name, "w")

    for edge in network_analysis.edges
        genes = edge.genes
        write(out_file, string(
            genes[1].name, "\t", genes[2].name, "\t",
            edge.confidence, "\n",
            genes[2].name, "\t", genes[1].name, "\t",
            edge.confidence, "\n"
        ))
    end

    close(out_file)

end

# Networks are undirected; edges are only listed in one direction
function get_edge_list(network_analysis::NetworkAnalysis)
    edge_list = Tuple{String,String,Float64}[]

    for edge in network_analysis.edges
        genes = edge.genes
        push!(edge_list, (genes[1].name, genes[2].name, edge.confidence))
    end

    return edge_list
end

function infer_network(data_file_path::String, inference::AbstractNetworkInference; delim::Union{Char,Bool} = false,
    discretizer = "bayesian_blocks", estimator = "maximum_likelihood", number_of_bins = 10, base = 2,
    out_file = "")

    println("Getting genes...")
    genes = get_genes(
        data_file_path,
        delim = delim,
        discretizer = discretizer,
        estimator = estimator,
        number_of_bins = number_of_bins
    )

    println("Inferring network...")
    network_analysis = NetworkAnalysis(inference, genes, estimator = estimator, base = base)

    if length(out_file) > 1
        println("Writing network to file...")
        write_network_file(out_file, network_analysis)
    end

    println("Getting edge list...")
    return get_edge_list(network_analysis)

end
