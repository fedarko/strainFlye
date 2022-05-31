import tempfile
import strainflye.fdr_utils as fu

# Header + first four mutations called on SheepGut using p = 0.15%
P_VCF = (
    "##fileformat=VCFv4.3\n"
    "##fileDate=20220526\n"
    '##source="strainFlye v0.0.1: p-mutation calling (--min-p = 0.15%)"\n'
    "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
    '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
    '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
    '##FILTER=<ID=strainflye_minp_15, Description="min p threshold (scaled up by 100)">\n'  # noqa: E501
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "edge_1\t43\t.\tA\tT\t.\t.\tMDP=380;AAD=2\n"
    "edge_1\t255\t.\tT\tC\t.\t.\tMDP=385;AAD=2\n"
    "edge_1\t356\t.\tA\tT\t.\t.\tMDP=403;AAD=2\n"
    "edge_1\t387\t.\tT\tC\t.\t.\tMDP=395;AAD=2\n"
)


def test_parse_vcf_good():
    fh = tempfile.NamedTemporaryFile(suffix=".vcf")
    with open(fh.name, "w") as f:
        f.write(P_VCF)

    vcf_obj, thresh_type, thresh_min = fu.parse_vcf(fh.name)

    # Let's get the easy checks out of the way first...
    assert thresh_type == "p"
    assert thresh_min == 15
    # This part of things -- parsing the VCF -- pysam handles, so we don't go
    # too in depth. Just verify it gets the positions, alleles, and info right.
    correct_pos = [43, 255, 356, 387]
    correct_mdp = [380, 385, 403, 395]
    correct_aad = [2, 2, 2, 2]
    correct_ref = ["A", "T", "A", "T"]
    correct_alt = ["T", "C", "T", "C"]
    for ri, rec in enumerate(vcf_obj.fetch()):
        assert rec.pos == correct_pos[ri]
        assert rec.info.get("MDP") == correct_mdp[ri]
        assert rec.ref == correct_ref[ri]

        # this part we'll have to change if we start calling multiple mutations
        # at a single position
        assert len(rec.alts) == 1
        assert rec.alts[0] == correct_alt[ri]

        # AAD is treated as a tuple containing one element by pysam, not as a
        # single value. This is because I've labelled the AAD INFO field with a
        # Number of "A", meaning that there is one of these per alternate
        # allele. Currently we don't call more than one mutation at a position,
        # so this doesn't do anything, but -- in case we modify this in the
        # future -- this leaves the door open for us to add this in.
        rec_aad = rec.info.get("AAD")
        assert len(rec_aad) == 1
        assert rec_aad[0] == correct_aad[ri]
