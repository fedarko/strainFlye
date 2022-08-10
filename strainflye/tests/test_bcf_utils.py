import os
import pytest
import strainflye.bcf_utils as bu
from strainflye.errors import ParameterError
from .utils_for_testing import write_vcf_tempfile


def test_parse_bcf_good():
    # Header + first four mutations called on SheepGut using p = 0.15%
    # ... just for reference, using StringIO didn't seem to work with pysam's
    # BCF reader, hence the use of tempfiles here
    # NOTE FROM MARCUS: pysam accepts either VCF or BCF files, and -- since we
    # just recently switched from working with VCF to BCF, mainly -- these
    # tests still output VCF files.
    fh = write_vcf_tempfile(
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
    bcf_obj, thresh_type, thresh_min = bu.parse_bcf(fh.name)

    # Let's get the easy checks out of the way first...
    assert thresh_type == "p"
    assert thresh_min == 15
    # This part of things -- parsing the BCF -- pysam handles, so we don't go
    # too in depth. Just verify it gets the positions, alleles, and info right.
    correct_pos = [43, 255, 356, 387]
    correct_mdp = [380, 385, 403, 395]
    correct_aad = [2, 2, 2, 2]
    correct_ref = ["A", "T", "A", "T"]
    correct_alt = ["T", "C", "T", "C"]
    for ri, rec in enumerate(bcf_obj.fetch()):
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


def test_parse_bcf_multi_sf_header():
    # Header + first four mutations called on SheepGut using p = 0.15%
    fh = write_vcf_tempfile(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220526\n"
        '##source="strainFlye v0.0.1: p-mutation calling (--min-p = 0.15%)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minp_15, Description="min p threshold (scaled up by 100)">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minp_20, Description="min p threshold (scaled up by 100)">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t43\t.\tA\tT\t.\t.\tMDP=380;AAD=2\n"
        "edge_1\t255\t.\tT\tC\t.\t.\tMDP=385;AAD=2\n"
        "edge_1\t356\t.\tA\tT\t.\t.\tMDP=403;AAD=2\n"
        "edge_1\t387\t.\tT\tC\t.\t.\tMDP=395;AAD=2\n"
    )
    with pytest.raises(ParameterError) as ei:
        bu.parse_bcf(fh.name)
    exp_patt = (
        f"BCF file {fh.name} has multiple strainFlye threshold filter headers."
    )
    assert str(ei.value) == exp_patt


def test_parse_bcf_no_sf_header():
    # Header + first four mutations called on SheepGut using p = 0.15%
    fh = write_vcf_tempfile(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220526\n"
        '##source="strainFlye v0.0.1: p-mutation calling (--min-p = 0.15%)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t43\t.\tA\tT\t.\t.\tMDP=380;AAD=2\n"
        "edge_1\t255\t.\tT\tC\t.\t.\tMDP=385;AAD=2\n"
        "edge_1\t356\t.\tA\tT\t.\t.\tMDP=403;AAD=2\n"
        "edge_1\t387\t.\tT\tC\t.\t.\tMDP=395;AAD=2\n"
    )
    with pytest.raises(ParameterError) as ei:
        bu.parse_bcf(fh.name)
    exp_patt = (
        f"BCF file {fh.name} doesn't seem to be from strainFlye: no threshold "
        "filter headers."
    )
    assert str(ei.value) == exp_patt


def test_parse_bcf_missing_info_fields():
    # Header + first four mutations called on SheepGut using p = 0.15%
    # kind of a gross way of doing this -- should ideally set this up as a list
    # from the start, theeen join it up, rather than going from string to list
    # to string again. but whatevs
    text = (
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
    split_text = text.splitlines()

    # remove just the MDP line
    fh = write_vcf_tempfile("\n".join(split_text[:4] + split_text[5:]))
    with pytest.raises(ParameterError) as ei:
        bu.parse_bcf(fh.name)
    exp_patt = f"BCF file {fh.name} needs to have MDP and AAD info fields."
    assert str(ei.value) == exp_patt

    # remove just the AAD line
    fh = write_vcf_tempfile("\n".join(split_text[:5] + split_text[6:]))
    with pytest.raises(ParameterError) as ei:
        bu.parse_bcf(fh.name)
    exp_patt = f"BCF file {fh.name} needs to have MDP and AAD info fields."
    assert str(ei.value) == exp_patt

    # remove both the MDP and AAD lines
    fh = write_vcf_tempfile("\n".join(split_text[:4] + split_text[6:]))
    with pytest.raises(ParameterError) as ei:
        bu.parse_bcf(fh.name)
    exp_patt = f"BCF file {fh.name} needs to have MDP and AAD info fields."
    assert str(ei.value) == exp_patt

    # finally, sanity check that putting back in all the lines works (doesn't
    # throw an error) -- we check this dataset in detail above, so no need to
    # do that again here
    fh = write_vcf_tempfile("\n".join(split_text))
    bu.parse_bcf(fh.name)


def test_get_mutated_positions_in_contig_good():
    test_bcf_path = os.path.join(
        "strainflye",
        "tests",
        "inputs",
        "small",
        "call-r-min3-di12345",
        "naive-calls.bcf",
    )
    bcf_obj, _, _ = bu.parse_bcf(test_bcf_path)

    mp_c1 = bu.get_mutated_positions_in_contig(bcf_obj, "c1")
    assert mp_c1 == set([3, 10, 12])

    mp_c2 = bu.get_mutated_positions_in_contig(bcf_obj, "c2")
    assert mp_c2 == set()

    mp_c3 = bu.get_mutated_positions_in_contig(bcf_obj, "c3")
    assert mp_c3 == set([6, 7])


def test_get_mutated_positions_in_contig_not_in_bcf():
    test_bcf_path = os.path.join(
        "strainflye",
        "tests",
        "inputs",
        "small",
        "call-r-min3-di12345",
        "naive-calls.bcf",
    )
    bcf_obj, _, _ = bu.parse_bcf(test_bcf_path)

    with pytest.raises(ValueError) as ei:
        bu.get_mutated_positions_in_contig(bcf_obj, "c4")

    assert str(ei.value) == (
        "Contig c4 is not described in the BCF object's header."
    )
