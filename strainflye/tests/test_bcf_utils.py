import os
import pytest
import pysam
import strainflye.bcf_utils as bu
from strainflye.errors import ParameterError
from .utils_for_testing import write_vcf_tempfile, write_indexed_bcf


BCF = os.path.join(
    "strainflye",
    "tests",
    "inputs",
    "small",
    "call-r-min3-di12345",
    "naive-calls.bcf",
)


def test_parse_sf_bcf_good():
    # Header + first four mutations called on SheepGut using p = 0.15%
    fp = write_indexed_bcf(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220526\n"
        '##source="strainFlye v0.0.1: p-mutation calling (--min-p = 0.15%)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        "##contig=<ID=edge_1,length=1000>\n"
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minp_15, Description="min p threshold (scaled up by 100)">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t43\t.\tA\tT\t.\t.\tMDP=380;AAD=2\n"
        "edge_1\t255\t.\tT\tC\t.\t.\tMDP=385;AAD=2\n"
        "edge_1\t356\t.\tA\tT\t.\t.\tMDP=403;AAD=2\n"
        "edge_1\t387\t.\tT\tC\t.\t.\tMDP=395;AAD=2\n"
    )
    bcf_obj, thresh_type, thresh_min = bu.parse_sf_bcf(fp)

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


def test_parse_sf_bcf_multi_sf_header():
    # Header + first four mutations called on SheepGut using p = 0.15%
    fp = write_indexed_bcf(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220526\n"
        '##source="strainFlye v0.0.1: p-mutation calling (--min-p = 0.15%)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        "##contig=<ID=edge_1,length=1000>\n"
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
        bu.parse_sf_bcf(fp)
    exp_patt = (
        f"BCF file {fp} has multiple strainFlye threshold filter headers."
    )
    assert str(ei.value) == exp_patt


def test_parse_sf_bcf_no_sf_header():
    # Header + first four mutations called on SheepGut using p = 0.15%
    fp = write_indexed_bcf(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220526\n"
        '##source="strainFlye v0.0.1: p-mutation calling (--min-p = 0.15%)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        "##contig=<ID=edge_1,length=1000>\n"
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t43\t.\tA\tT\t.\t.\tMDP=380;AAD=2\n"
        "edge_1\t255\t.\tT\tC\t.\t.\tMDP=385;AAD=2\n"
        "edge_1\t356\t.\tA\tT\t.\t.\tMDP=403;AAD=2\n"
        "edge_1\t387\t.\tT\tC\t.\t.\tMDP=395;AAD=2\n"
    )
    with pytest.raises(ParameterError) as ei:
        bu.parse_sf_bcf(fp)
    exp_patt = (
        f"BCF file {fp} doesn't seem to be from strainFlye: no threshold "
        "filter headers."
    )
    assert str(ei.value) == exp_patt


def test_parse_sf_bcf_missing_info_fields():
    # Header + first four mutations called on SheepGut using p = 0.15%
    # kind of a gross way of doing this -- should ideally set this up as a list
    # from the start, theeen join it up, rather than going from string to list
    # to string again. but whatevs
    text = (
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220526\n"
        '##source="strainFlye v0.0.1: p-mutation calling (--min-p = 0.15%)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        "##contig=<ID=edge_1,length=1000>\n"
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minp_15, Description="min p threshold (scaled up by 100)">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    split_text = text.splitlines()

    # There are multiple "versions" of the final line in this file we want to
    # check, so we create different strings for each version. (I guess this is
    # the easiest way to test this...) Note that the leading \n is beacuse
    # using "\n".join() will strip the final newline of the last line in
    # split_text. Gross, but this works at least.
    mut_line_all = "\nedge_1\t43\t.\tA\tT\t.\t.\tMDP=380;AAD=2\n"
    mut_line_no_mdp = "\nedge_1\t43\t.\tA\tT\t.\t.\tAAD=2\n"
    mut_line_no_aad = "\nedge_1\t43\t.\tA\tT\t.\t.\tMDP=380\n"
    mut_line_none = "\nedge_1\t43\t.\tA\tT\t.\t.\t.\n"

    # First, sanity check that all lines work (as in, no error is thrown).
    # We check this dataset in detail above, so no need to do that again here;
    # we just verify that we can parse it ok.
    fp = write_indexed_bcf("\n".join(split_text))
    bu.parse_sf_bcf(fp)

    # OK, now look for errors!
    ###########################################################################
    # If there is a mutation containing MDP/AAD info fields -- but if these
    # fields are not defined in the header -- then pysam itself will throw an
    # error. Let's just test this here real quick.
    #
    # These correspond to removing the MDP header, removing the AAD header,
    # then removing both
    for bounds in ((5, 6), (6, 7), (5, 7)):
        fp = write_indexed_bcf(
            "\n".join(split_text[: bounds[0]] + split_text[bounds[1] :])
            + mut_line_all
        )
        with pytest.raises(ValueError) as ei:
            bu.parse_sf_bcf(fp)
        assert "is it VCF/BCF format?" in str(ei.value)

    ###########################################################################
    # The more interesting case is when (for MDP, AAD, for both) there is no
    # info header field *and* no mutations have this info field.
    # In this case, the BCF is technically valid (so pysam won't complain) --
    # so we will need to complain.
    #
    # 1. Remove all evidence of MDP
    fp = write_indexed_bcf(
        "\n".join(split_text[:5] + split_text[6:]) + mut_line_no_mdp
    )
    with pytest.raises(ParameterError) as ei:
        bu.parse_sf_bcf(fp)
    # (sorry, the filename changes every time we write out a new temp. bcf
    # file, so we have to keep recreating these expected error messages. eesh.)
    assert str(ei.value) == (
        f"BCF file {fp} needs to have MDP and AAD info fields."
    )

    # 2. Remove all evidence of AAD
    fp = write_indexed_bcf(
        "\n".join(split_text[:6] + split_text[7:]) + mut_line_no_aad
    )
    with pytest.raises(ParameterError) as ei:
        bu.parse_sf_bcf(fp)
    assert str(ei.value) == (
        f"BCF file {fp} needs to have MDP and AAD info fields."
    )

    # 3. Remove evidence of both MDP and AAD
    fp = write_indexed_bcf(
        "\n".join(split_text[:5] + split_text[7:]) + mut_line_none
    )
    with pytest.raises(ParameterError) as ei:
        bu.parse_sf_bcf(fp)
    assert str(ei.value) == (
        f"BCF file {fp} needs to have MDP and AAD info fields."
    )


def test_parse_sf_bcf_no_contigs_in_header():
    # Header + first four mutations called on SheepGut using p = 0.15%
    vcf_text = (
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
    # BCF files *need* to have contig header tags, so trying to load a
    # BCF file without a contig header will make pysam complain:
    bcf_fp = write_indexed_bcf(vcf_text)
    with pytest.raises(ValueError) as ei:
        bu.parse_sf_bcf(bcf_fp)
    assert "is it VCF/BCF format?" in str(ei.value)

    # However, we can sneak past this by using a VCF instead of a VCF file. In
    # this case, our custom check triggers.
    vcf_fp = write_vcf_tempfile(vcf_text)
    with pytest.raises(ParameterError) as ei:
        bu.parse_sf_bcf(vcf_fp)
    # yeah yeah i know the error message says BCF instead of VCF but we assume
    # that these files are all BCF, like literally the parameter is --bcf or
    # whatever, not worth saying "BCF or VCF" esp when we assume that
    # everything is indexed BCF, this particular issue just comes up in testing
    # (... so, this check is kind of useless since it seems impossible to come
    # up in practice, but whatever).
    assert str(ei.value) == (
        f"BCF file {vcf_fp} doesn't describe any contigs in its header."
    )


def test_parse_sf_bcf_no_contig_length_in_header():
    # Header + first four mutations called on SheepGut using p = 0.15%
    fp = write_indexed_bcf(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220526\n"
        '##source="strainFlye v0.0.1: p-mutation calling (--min-p = 0.15%)"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        "##contig=<ID=edge_1>\n"
        '##INFO=<ID=MDP,Number=1,Type=Integer,Description="(Mis)match read depth">\n'  # noqa: E501
        '##INFO=<ID=AAD,Number=A,Type=Integer,Description="Alternate allele read depth">\n'  # noqa: E501
        '##FILTER=<ID=strainflye_minp_15, Description="min p threshold (scaled up by 100)">\n'  # noqa: E501
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t43\t.\tA\tT\t.\t.\tMDP=380;AAD=2\n"
        "edge_1\t255\t.\tT\tC\t.\t.\tMDP=385;AAD=2\n"
        "edge_1\t356\t.\tA\tT\t.\t.\tMDP=403;AAD=2\n"
        "edge_1\t387\t.\tT\tC\t.\t.\tMDP=395;AAD=2\n"
    )
    with pytest.raises(ParameterError) as ei:
        bu.parse_sf_bcf(fp)
    exp_patt = f"BCF file {fp} has no length given for contig edge_1."
    assert str(ei.value) == exp_patt


def test_get_mutated_positions_in_contig_good():
    bcf_obj = bu.parse_arbitrary_bcf(BCF)

    mp_c1 = bu.get_mutated_positions_in_contig(bcf_obj, "c1")
    assert mp_c1 == set([3, 10, 12])

    mp_c2 = bu.get_mutated_positions_in_contig(bcf_obj, "c2")
    assert mp_c2 == set()

    mp_c3 = bu.get_mutated_positions_in_contig(bcf_obj, "c3")
    assert mp_c3 == set([6, 7])


def test_verify_contig_in_bcf():
    bcf_obj = bu.parse_arbitrary_bcf(BCF)
    with pytest.raises(ParameterError) as ei:
        bu.verify_contig_in_bcf(bcf_obj, "c5")
    assert str(ei.value) == (
        "Contig c5 is not described in the BCF object's header."
    )


def test_get_mutated_positions_in_contig_not_in_bcf():
    bcf_obj = bu.parse_arbitrary_bcf(BCF)

    with pytest.raises(ParameterError) as ei:
        bu.get_mutated_positions_in_contig(bcf_obj, "c4")

    assert str(ei.value) == (
        "Contig c4 is not described in the BCF object's header."
    )


def test_get_mutated_position_details_in_contig_good():
    bcf_obj = bu.parse_arbitrary_bcf(BCF)
    assert bu.get_mutated_position_details_in_contig(bcf_obj, "c1") == {
        3: ("G", "T"),
        10: ("G", "A"),
        12: ("G", "A"),
    }
    assert bu.get_mutated_position_details_in_contig(bcf_obj, "c2") == {}
    assert bu.get_mutated_position_details_in_contig(bcf_obj, "c3") == {
        6: ("A", "T"),
        7: ("T", "C"),
    }


def test_get_mutated_position_details_in_contig_not_in_bcf():
    bcf_obj = bu.parse_arbitrary_bcf(BCF)
    with pytest.raises(ParameterError) as ei:
        bu.get_mutated_positions_in_contig(bcf_obj, "c4")
    assert str(ei.value) == (
        "Contig c4 is not described in the BCF object's header."
    )


def test_verify_no_multiallelic_mutations_in_bcf():
    fp = write_indexed_bcf(
        "##fileformat=VCFv4.3\n"
        "##fileDate=20220818\n"
        '##source="my source is i made it the heck up"\n'
        "##reference=/Poppy/mfedarko/sheepgut/main-workflow/output/all_edges.fasta\n"  # noqa: E501
        "##contig=<ID=edge_1>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "edge_1\t43\t.\tA\tT\t.\t.\t.\n"
        "edge_1\t43\t.\tA\tC\t.\t.\t.\n"
        "edge_1\t255\t.\tT\tC\t.\t.\t.\n"
        "edge_1\t356\t.\tA\tT\t.\t.\t.\n"
        "edge_1\t387\t.\tT\tC\t.\t.\t.\n"
    )
    bcf_obj = pysam.VariantFile(fp)
    with pytest.raises(ParameterError) as ei:
        bu.verify_no_multiallelic_mutations_in_bcf(bcf_obj, fp)
    assert str(ei.value) == (
        f"BCF file {fp} has multiple mutations at (1-indexed) "
        "position 43 on contig edge_1. strainFlye does not currently support "
        "BCF files containing multi-allelic mutations, sorry."
    )
