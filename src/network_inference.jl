abstract type AbstractNetworkInference end

immutable PIDNetworkInference <: AbstractNetworkInference end
immutable PIDCNetworkInference <: AbstractNetworkInference end
immutable MINetworkInference <: AbstractNetworkInference end

function NetworkAnalysis(::MINetworkInference, genes::Array{Gene}; print_status = false)

    function get_mi_scores(i, j)
        gene1 = genes[i]
        gene2 = genes[j]
        frequencies = get_frequencies_from_bin_ids(
            gene1.discretized_values,
            gene2.discretized_values,
            gene1.number_of_bins,
            gene2.number_of_bins,
        )
        probabilities = get_probabilities("maximum_likelihood", frequencies)
        probabilities_1 = sum(probabilities, 2)
        probabilities_2 = sum(probabilities, 1)
        mi = apply_mutual_information_formula(probabilities, probabilities_1, probabilities_2, b)
        mi_scores[i, j] += mi
        mi_scores[j, i] += mi
    end

    function populate_edges(i, j, index)
        mi_score = mi_scores[i, j]
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
    b = 2

    # Get MI scores
    for i in 1 : number_of_genes
        if print_status
            println(i)
        end
        for j in i+1 : number_of_genes
            get_mi_scores(i, j)
        end
    end

    # Get edges
    index = 0
    for i in 1 : number_of_genes
        for j in i+1 : number_of_genes
            index += 1
            populate_edges(i, j, index)
        end
    end

    # Sort edges by confidence
    indices = sortperm(confidences, rev = true)
    edges = edges[indices]

    return NetworkAnalysis(genes, edges)

end

# TODO: allow choose estimator
# TODO: allow choose distribution
function NetworkAnalysis(::PIDCNetworkInference, genes::Array{Gene}; print_status = false)

    function populate_gene_pairs(i, j)
        gene1 = genes[i]
        gene2 = genes[j]
        frequencies = get_frequencies_from_bin_ids(
            gene1.discretized_values,
            gene2.discretized_values,
            gene1.number_of_bins,
            gene2.number_of_bins
        )
        probabilities = get_probabilities("maximum_likelihood", frequencies)
        probabilities_1 = sum(probabilities, 2)
        probabilities_2 = sum(probabilities, 1)
        mi = apply_mutual_information_formula(probabilities, probabilities_1, probabilities_2, b)
        gene_pairs[i, j] = GenePair(mi, apply_specific_information_formula(probabilities, probabilities_1, probabilities_2, 1, b))
        gene_pairs[j, i] = GenePair(mi, apply_specific_information_formula(probabilities, probabilities_2, probabilities_1, 2, b))
    end

    function increment_puc_scores(gene_index_1, gene_index_2, mi, redundancy)
        unique = mi - redundancy
        puc_score = unique / mi
        puc_score = isfinite(puc_score) ? puc_score : zero(puc_score)
        puc_scores[gene_index_1, gene_index_2] += puc_score
        puc_scores[gene_index_2, gene_index_1] += puc_score
    end

    function get_pid_and_increment_puc_scores(i, j, k)
        redundancy_target_k = apply_redundancy_formula(genes[k].probabilities, gene_pairs[i, k].specific_information, gene_pairs[j, k].specific_information, b)
        redundancy_target_j = apply_redundancy_formula(genes[j].probabilities, gene_pairs[i, j].specific_information, gene_pairs[k, j].specific_information, b)
        redundancy_target_i = apply_redundancy_formula(genes[i].probabilities, gene_pairs[j, i].specific_information, gene_pairs[k, i].specific_information, b)

        increment_puc_scores(i, k, gene_pairs[i, k].mi, redundancy_target_k)
        increment_puc_scores(j, k, gene_pairs[j, k].mi, redundancy_target_k)

        increment_puc_scores(i, j, gene_pairs[i, j].mi, redundancy_target_j)
        increment_puc_scores(k, j, gene_pairs[k, j].mi, redundancy_target_j)

        increment_puc_scores(j, i, gene_pairs[j, i].mi, redundancy_target_i)
        increment_puc_scores(k, i, gene_pairs[k, i].mi, redundancy_target_i)
    end

    function populate_edges_and_confidences(i, j, index)
        puc_score = puc_scores[i, j]
        # Ignore self-interaction zeros
        puc_scores_i = vcat(puc_scores[1:i-1, i], puc_scores[i+1:end, i])
        puc_scores_j = vcat(puc_scores[1:j-1, j], puc_scores[j+1:end, j])

        confidence = puc_score
        confidence =
            cdf(fit(Gamma, puc_scores_i), puc_score) +
            cdf(fit(Gamma, puc_scores_j), puc_score)
        edges[index] = Edge(
            Set([genes[i], genes[j]]),
            confidence
        )
        confidences[index] = confidence
    end

    b = 2
    number_of_genes = length(genes)
    puc_scores = zeros(number_of_genes, number_of_genes)
    edges = Array{Edge}(binomial(number_of_genes, 2))
    confidences = zeros(length(edges))

    gene_pairs = Array{GenePair}(number_of_genes, number_of_genes)

    # Actions between pairs:
    # Make GenePairs and calculate MI and specific information between each pair 
    for i in 1 : number_of_genes
        for j in i+1 : number_of_genes
            populate_gene_pairs(i, j)
        end
    end

    # Actions between triples:
    # Get redundancy between each triple, and store PUC score for each gene
    for i in 1 : number_of_genes
        if print_status
            println(i)
        end
        for j in i+1 : number_of_genes
            for k in j+1 : number_of_genes
                get_pid_and_increment_puc_scores(i, j, k)
            end
        end
    end

    # Calculate confidences and get the edges
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
