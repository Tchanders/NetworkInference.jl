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

function get_joint_probabilities(gene1, gene2, estimator)
    frequencies = get_frequencies_from_bin_ids(
        gene1.discretized_values,
        gene2.discretized_values,
        gene1.number_of_bins,
        gene2.number_of_bins
    )
    probabilities = get_probabilities(estimator, frequencies)
    # probabilities is already property of a gene, but doing this gets correct array shapes.
    # Also, for MI and CLR, it means that we don't assume that the marginal probabilities for
    # a gene are always the same, no matter what the second gene is, meaning that we can use
    # estimators other than maximum likelihood. (We still can't do this for PUC and PIDC,
    # because we do make that assumption for 3-gene joint distributions, in get_puc.)
    probabilities1 = sum(probabilities, 2)
    probabilities2 = sum(probabilities, 1)
    return (probabilities, probabilities1, probabilities2)
end

function get_mi_scores(genes, number_of_genes, estimator, base)

    function get_mi(gene1, gene2, i, j, base, mi_scores)
        probabilities, probabilities1, probabilities2 = get_joint_probabilities(gene1, gene2, estimator)
        mi = apply_mutual_information_formula(probabilities, probabilities1, probabilities2, base)
        mi_scores[i, j] = mi
        mi_scores[j, i] = mi
    end

    mi_scores = SharedArray{Float64}(number_of_genes, number_of_genes)
    @sync @parallel for i in 1 : number_of_genes
        for j in i+1 : number_of_genes
            get_mi(genes[i], genes[j], i, j, base, mi_scores)
        end
    end
    return mi_scores

end

function get_puc_scores(genes, number_of_genes, estimator, base)

    function get_mi_and_si(gene1, gene2, base) # Mutual information and specific information
        probabilities, probabilities1, probabilities2 = get_joint_probabilities(gene1, gene2, estimator)
        mi = apply_mutual_information_formula(probabilities, probabilities1, probabilities2, base)
        si1 = apply_specific_information_formula(probabilities, probabilities1, probabilities2, 1, base)
        si2 = apply_specific_information_formula(probabilities, probabilities2, probabilities1, 2, base)
        return (mi, si1, si2)
    end

    function get_gene_pairs(gene1, gene2, i, j, base)
        mi, si1, si2 = get_mi_and_si(gene1, gene2, base)
        gene_pairs[i, j] = GenePair(mi, si1)
        gene_pairs[j, i] = GenePair(mi, si2)
    end

    function increment_puc_scores(x, z, mi, redundancy, puc_scores)
        puc_score = (mi - redundancy) / mi
        puc_score = isfinite(puc_score) ? puc_score : zero(puc_score)
        puc_scores[x, z] += puc_score
        puc_scores[z, x] += puc_score
    end

    function get_puc(target, source1_target, source2_target, x, y, z, puc_scores)
        redundancy = apply_redundancy_formula(
            target.probabilities,
            source1_target.si,
            source2_target.si,
            base
        )
        increment_puc_scores(x, z, source1_target.mi, redundancy, puc_scores)
        increment_puc_scores(y, z, source2_target.mi, redundancy, puc_scores)
    end

    gene_pairs = Array{GenePair}(number_of_genes, number_of_genes)
    puc_scores = SharedArray{Float64}(number_of_genes, number_of_genes)
    for i in 1 : number_of_genes
        for j in i+1 : number_of_genes
            get_gene_pairs(genes[i], genes[j], i, j, base)
        end
    end
    @sync @parallel for i in 1 : number_of_genes
        for j in i+1 : number_of_genes
            for k in j+1 : number_of_genes
                get_puc(genes[k], gene_pairs[i, k], gene_pairs[j, k], i, j, k, puc_scores)
                get_puc(genes[j], gene_pairs[i, j], gene_pairs[k, j], i, k, j, puc_scores)
                get_puc(genes[i], gene_pairs[j, i], gene_pairs[k, i], j, k, i, puc_scores)
            end
        end
    end
    return puc_scores

end

function get_confidences(inference, scores, number_of_genes, genes)

    # In their respective original implementations, CLR and PIDC applied network context in slightly
    # different ways. Those differences are respected here; in informal tests, they have not been
    # found to make much of a difference.
    function get_confidence(::PIDCNetworkInference, i, j, scores, confidences, genes)
        score = scores[i, j]
        scores_i = vcat(scores[1:i-1, i], scores[i+1:end, i])
        scores_j = vcat(scores[1:j-1, j], scores[j+1:end, j])
        try
            confidences[i, j] = cdf(fit(Gamma, scores_i), score) + cdf(fit(Gamma, scores_j), score)
        catch
            println(string("Gamma distribution failed for ", genes[i].name, " and ", genes[j].name, "; used normal instead."))
            get_clr_confidence(i, j, score, scores_i, scores_j, confidences)
        end
    end
    function get_confidence(::CLRNetworkInference, i, j, scores, confidences, genes)
        score = scores[i, j]
        scores_i = vcat(scores[1:i-1, i], scores[i+1:end, i])
        scores_j = vcat(scores[1:j-1, j], scores[j+1:end, j])
        get_clr_confidence(i, j, score, scores_i, scores_j, confidences)
    end
    function get_clr_confidence(i, j, score, scores_i, scores_j, confidences)
        confidences[i, j] = sqrt(
            (var(scores_i) == 0 ? 0 : (score - mean(scores_i))^2 / var(scores_i)) +
            (var(scores_j) == 0 ? 0 : (score - mean(scores_j))^2 / var(scores_j))
        )
    end

    confidences = SharedArray{Float64}(number_of_genes, number_of_genes)
    @sync @parallel for i in 1 : number_of_genes
        for j in i+1 : number_of_genes
            get_confidence(inference, i, j, scores, confidences, genes)
        end
    end
    return confidences

end

# The maximum likelihood estimator is recommended for PUC and PIDC, because speedups
# are made here, based on the assumption that the marginal probability distribution for
# a gene, from the joint distribution with any two other genes is always the same. If
# the joint distributions are estimated using other estimators, this assumption is
# violated for PUC and PIDC in get_puc and get_joint_probabilities.
function NetworkAnalysis(inference::AbstractNetworkInference, genes::Array{Gene}; estimator = "maximum_likelihood", base = 2)

    # Constants and containers
    number_of_genes = length(genes)
    edges = Array{Edge}(binomial(number_of_genes, base))

    # Get the raw scores
    if get_puc(inference)
        scores = get_puc_scores(genes, number_of_genes, estimator, base)
    else
        scores = get_mi_scores(genes, number_of_genes, estimator, base)
    end

    # Apply context if necessary
    if apply_context(inference)
        confidences = get_confidences(inference, scores, number_of_genes, genes)
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
