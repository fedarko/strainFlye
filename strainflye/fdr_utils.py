# Utilities for strainFlye fdr.


def run_estimate(
    vcf,
    diversity_indices,
    decoy_contig,
    decoy_context,
    high_p,
    high_r,
    output_fdr_info,
):
    # 1. Figure out range of p or r to use. Create a list, threshold_vals.
    # 2. Identify decoy genome
    # 3. For each value in threshold_vals, compute the decoy genome's mutation
    #    rate. Save this to a list, decoy_mut_rates -- this will have the same
    #    dimensions as threshold_vals.
    # 4. For each target genome...
    #    - For each value in threshold_vals...
    #      - Compute the mutation rate for this target genome at this
    #        threshold value.
    #      - Compute the FDR estimate for this pair of (target, threshold).
    #        Save to a list of target_fdr_ests, which has the same dimensions
    #        as threshold_vals.
    #    - Write out a new row to the FDR estimate file describing
    #      target_fdr_ests.
    pass


def run_fix(vcf, fdr_info, fdr, output_vcf):
    pass
