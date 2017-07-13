abstract type AbstractNetworkInference end

immutable MINetworkInference <: AbstractNetworkInference end
immutable CLRNetworkInference <: AbstractNetworkInference end
immutable PUCNetworkInference <: AbstractNetworkInference end
immutable PIDCNetworkInference <: AbstractNetworkInference end

# Context trait
apply_context(::MINetworkInference) = false
apply_context(::CLRNetworkInference) = true
apply_context(::PUCNetworkInference) = false
apply_context(::PIDCNetworkInference) = true

# PUC trait
get_puc(::MINetworkInference) = false
get_puc(::CLRNetworkInference) = false
get_puc(::PUCNetworkInference) = true
get_puc(::PIDCNetworkInference) = true

# For sorting the edges
get_confidence(edge::Edge) = edge.confidence

function get_joint_probabilities(gene1, gene2)
    frequencies = get_frequencies_from_bin_ids(
        gene1.discretized_values,
        gene2.discretized_values,
        gene1.number_of_bins,
        gene2.number_of_bins
    )
    probabilities = get_probabilities("maximum_likelihood", frequencies)
    # probabilities is already property of a gene, but doing this gets correct array shapes
    probabilities1 = sum(probabilities, 2)
    probabilities2 = sum(probabilities, 1)
    return (probabilities, probabilities1, probabilities2)
end

# Mutual information
function get_mi(gene1, gene2, base)
    probabilities, probabilities1, probabilities2 = get_joint_probabilities(gene1, gene2)
    return apply_mutual_information_formula(probabilities, probabilities1, probabilities2, base)
end

# Mutual information and specific information
function get_mi_and_si(gene1, gene2, base)
    probabilities, probabilities1, probabilities2 = get_joint_probabilities(gene1, gene2)
    mi = apply_mutual_information_formula(probabilities, probabilities1, probabilities2, base)
    si1 = apply_specific_information_formula(probabilities, probabilities1, probabilities2, 1, base)
    si2 = apply_specific_information_formula(probabilities, probabilities2, probabilities1, 2, base)
    return (mi, si1, si2)
end

function get_mi_scores(genes, number_of_genes, base)

    mi_scores = zeros(number_of_genes, number_of_genes)
    for i in 1 : number_of_genes
        gene1 = genes[i]
        for j in i+1 : number_of_genes
            gene2 = genes[j]
            mi = get_mi(gene1, gene2, base)
            mi_scores[i, j] = mi
            mi_scores[j, i] = mi
        end
    end
    return mi_scores

end

function get_puc_scores(genes, number_of_genes, base)

    function increment_puc_scores(x, z, mi, redundancy)
        puc_score = (mi - redundancy) / mi
        puc_score = isfinite(puc_score) ? puc_score : zero(puc_score)
        puc_scores[x, z] += puc_score
        puc_scores[z, x] += puc_score
    end
    function get_puc(x, y, z)
        redundancy = apply_redundancy_formula(
            genes[z].probabilities,
            gene_pairs[x, z].si,
            gene_pairs[y, z].si,
            base
        )
        increment_puc_scores(x, z, gene_pairs[x, z].mi, redundancy)
        increment_puc_scores(y, z, gene_pairs[y, z].mi, redundancy)
    end

    gene_pairs = Array{GenePair}(number_of_genes, number_of_genes)
    puc_scores = zeros(number_of_genes, number_of_genes)
    for i in 1 : number_of_genes
        gene1 = genes[i]
        for j in i+1 : number_of_genes
            gene2 = genes[j]
            mi, si1, si2 = get_mi_and_si(gene1, gene2, base)
            gene_pairs[i, j] = GenePair(mi, si1)
            gene_pairs[j, i] = GenePair(mi, si2)
        end
    end
    for i in 1 : number_of_genes
        for j in i+1 : number_of_genes
            for k in j+1 : number_of_genes
                get_puc(i, j, k)
                get_puc(i, k, j)
                get_puc(j, k, i)
            end
        end
    end
    return puc_scores

end

# In their respective original implementations, CLR and PIDC applied network context in slightly
# different ways. Those differences are respected here; in informal tests, they have not been
# found to make much of a difference.
function get_confidence(::PIDCNetworkInference, i, j, scores, confidences)
    score = scores[i, j]
    scores_i = vcat(scores[1:i-1, i], scores[i+1:end, i])
    scores_j = vcat(scores[1:j-1, j], scores[j+1:end, j])
    confidences[i, j] = cdf(fit(Gamma, scores_i), score) + cdf(fit(Gamma, scores_j), score)
end
function get_confidence(::CLRNetworkInference, i, j, scores, confidences)
    score = scores[i, j]
    scores_i = vcat(scores[1:i-1, i], scores[i+1:end, i])
    scores_j = vcat(scores[1:j-1, j], scores[j+1:end, j])
    confidences[i, j] = sqrt(
        (score - mean(scores_i))^2 / var(scores_i) +
        (score - mean(scores_j))^2 / var(scores_j)
    )
end

function get_confidences(inference, scores, number_of_genes)

    confidences = zeros(number_of_genes, number_of_genes)
    for i in 1 : number_of_genes
        for j in i+1 : number_of_genes
            get_confidence(inference, i, j, scores, confidences)
        end
    end
    return confidences

end

# TODO: allow choose estimator?
# TODO: allow choose distribution?
function NetworkAnalysis(inference::AbstractNetworkInference, genes::Array{Gene}; print_status = false)

    # Constants and containers
    base = 2
    number_of_genes = length(genes)
    edges = Array{Edge}(binomial(number_of_genes, 2))

    # Get the raw scores
    if get_puc(inference)
        scores = get_puc_scores(genes, number_of_genes, base)
    else
        scores = get_mi_scores(genes, number_of_genes, base)
    end

    # Apply context if necessary
    if apply_context(inference)
        confidences = get_confidences(inference, scores, number_of_genes)
    else
        confidences = scores
    end

    # Get edges from scores
    index = 0
    for i in 1 : number_of_genes
        gene1 = genes[i]
        for j in i+1 : number_of_genes
            index += 1
            gene2 = genes[j]
            edges[index] = Edge(
                [gene1, gene2],
                confidences[i, j]
            )
        end
    end
    sort!(edges; by = get_confidence, rev = true)

    return NetworkAnalysis(genes, edges)

end
