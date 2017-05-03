abstract AbstractNetworkInference

immutable PIDNetworkInference <: AbstractNetworkInference end
immutable PIDCNetworkInference <: AbstractNetworkInference end
immutable MINetworkInference <: AbstractNetworkInference end

# TODO: allow choose estimator
# TODO: allow choose distribution


function NetworkAnalysis(::MINetworkInference, genes::Array{Gene})

 function get_mi(gene1, gene2)
   frequencies = get_frequencies_from_bin_ids(
       gene1.discretized_values,
       gene2.discretized_values,
       gene1.number_of_bins,
       gene2.number_of_bins,
   )
   probabilities = get_probabilities("maximum_likelihood", frequencies)
   mi = get_mutual_information(
     probabilities,
     base = 2,
     probabilities = true,
     estimator = "maximum_likelihood" #temporary, in future allow choice
   )
   return mi
 end

 function increment_mi_scores(gene_index_1, gene_index_2, unique)
     mi_score = unique
     mi_score = isnan(mi_score) ? 0 : mi_score
     mi_scores[gene_index_1, gene_index_2] += mi_score
     mi_scores[gene_index_2, gene_index_1] += mi_score

     return mi_scores
 end

 function get_mi_and_increment_mi_scores(i, j, k)
     mi = get_mi(genes[i], genes[j])

     increment_mi(i, j, mi) #output of MI is just a number
 end

function populate_edges_and_confidences(i, j, index)
        mi_score = mi_scores[i, j]
         # Ignore self-interaction zeros

       edges[index] = Edge(
           Set([genes[i], genes[j]]),
           mi_score
         )
       confidences[index] = mi_score
end

 number_of_genes = length(genes)
 mi_scores = zeros(number_of_genes, number_of_genes)
 edges = Array{Edge}(binomial(number_of_genes, 2))
 confidences = zeros(length(edges))

 # Get MI scores
 for i in 1 : number_of_genes
     println(i)
     for j in i+1 : number_of_genes #select a gene which has not been selected before
          get_mi_and_increment_mi_scores(i, j)
     end
 end

 # Get edges
 index = 0
 for i in 1 : number_of_genes
     for j in i+1 : number_of_genes
         index += 1
         populate_edges_and_confidences(i, j, index)
     end
 end

 # Sort edges by confidence
 indices = sortperm(confidences, rev = true)
 edges = edges[indices]

 return NetworkAnalysis(genes, edges)

end


function NetworkAnalysis(::PIDCNetworkInference, genes::Array{Gene})
    
    function get_pid(gene1, gene2, gene3)
        frequencies = get_frequencies_from_bin_ids(
            gene1.discretized_values,
            gene2.discretized_values,
            gene3.discretized_values,
            gene1.number_of_bins,
            gene2.number_of_bins,
            gene3.number_of_bins
        )
        probabilities = get_probabilities("maximum_likelihood", frequencies)
        pid = get_partial_information_decomposition(
            probabilities,
            probabilities = true,
            all_orientations = true,
            include_synergy = false
        )
        return pid
    end

    function increment_puc_scores(gene_index_1, gene_index_2, unique, redundancy)
        puc_score = unique / (unique + redundancy)
        puc_score = isnan(puc_score) ? 0 : puc_score
        puc_scores[gene_index_1, gene_index_2] += puc_score
        puc_scores[gene_index_2, gene_index_1] += puc_score
    end

    function get_pid_and_increment_puc_scores(i, j, k)
        pid = get_pid(genes[i], genes[j], genes[k])

        increment_puc_scores(i, k, pid["z"]["unique_1"], pid["z"]["redundancy"])
        increment_puc_scores(j, k, pid["z"]["unique_2"], pid["z"]["redundancy"])

        increment_puc_scores(i, j, pid["y"]["unique_1"], pid["y"]["redundancy"])
        increment_puc_scores(k, j, pid["y"]["unique_2"], pid["y"]["redundancy"])

        increment_puc_scores(j, i, pid["x"]["unique_1"], pid["x"]["redundancy"])
        increment_puc_scores(k, i, pid["x"]["unique_2"], pid["x"]["redundancy"])
    end

    function populate_edges_and_confidences(i, j, index)
        puc_score = puc_scores[i, j]
        # Ignore self-interaction zeros
        puc_scores_i = vcat(puc_scores[1:i-1, i], puc_scores[i+1:end, i])
        puc_scores_j = vcat(puc_scores[1:j-1, j], puc_scores[j+1:end, j])
        confidence =
            cdf(fit(Normal, puc_scores_i), puc_score) +
            cdf(fit(Normal, puc_scores_j), puc_score)
        edges[index] = Edge(
            Set([genes[i], genes[j]]),
            confidence
        )
        confidences[index] = confidence
    end

    number_of_genes = length(genes)
    puc_scores = zeros(number_of_genes, number_of_genes)
    edges = Array{Edge}(binomial(number_of_genes, 2))
    confidences = zeros(length(edges))

    # Get PUC scores
    for i in 1 : number_of_genes
        println(i)
        for j in i+1 : number_of_genes
            for k in j+1 : number_of_genes
                get_pid_and_increment_puc_scores(i, j, k)
            end
        end
    end

    # Get edges
    index = 0
    for i in 1 : number_of_genes
        for j in i+1 : number_of_genes
            index += 1
            populate_edges_and_confidences(i, j, index)
        end
    end

    # Sort edges by confidence
    indices = sortperm(confidences, rev = true)
    edges = edges[indices]

    return NetworkAnalysis(genes, edges)

end
