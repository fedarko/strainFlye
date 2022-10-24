import os
import pickle
import skbio
import pytest
import pandas as pd
import strainflye.matrix_utils as mu
from io import StringIO as sio
from collections import defaultdict
from strainflye import config
from strainflye.errors import ParameterError, WeirdError
from strainflye.tests.utils_for_testing import mock_log, mock_log_2


IN_DIR = os.path.join("strainflye", "tests", "inputs", "small")
FASTA = os.path.join(IN_DIR, "contigs.fasta")
BAM = os.path.join(IN_DIR, "alignment.bam")
GFF = os.path.join(IN_DIR, "genes.gff")
C1 = skbio.DNA("ACTGACACCCAAACCAAACCTAC")
C2 = skbio.DNA("AAAAAAGGGGGG")


def test_get_contig_cds_info_good(capsys):
    gff = "##gff-version 3\nc1	marcus	cds	5	19	.	+	0	ID=hi"
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    no_loops = True
    for contig, im in cim_tuples:
        no_loops = False
        cds_df, cc = mu.get_contig_cds_info(
            im, contig, C1, mock_log, mock_log_2
        )
        pd.testing.assert_frame_equal(
            cds_df,
            pd.DataFrame(
                {"LeftEnd": [5], "RightEnd": [19], "Strand": "+"},
                index=pd.Index(["hi"]),
            ),
        )
        assert cc.cds2left2counter == {
            "hi": {
                5: defaultdict(int),
                8: defaultdict(int),
                11: defaultdict(int),
                14: defaultdict(int),
                17: defaultdict(int),
            },
        }
        break
    # account for the silly case where skbio breaks and cim_tuples is empty, in
    # which case the for loop's body will never execute.
    if no_loops:
        raise WeirdError
    assert capsys.readouterr().out == (
        "MockLog2: Found 1 feature belonging to contig c1; inspecting...\n"
        f"MockLog2: Found 1 feature with a type in {config.CDS_TYPES} in "
        "contig c1. Going through alignments...\n"
    )


def test_get_contig_cds_info_multi_features(capsys):
    gff = """##gff-version 3
c1	marcus	cds	5	19	.	+	0	ID=hi
c1	marcus	gene	1	5	.	+	0	ID=2
c1	marcus	CDS	18	20	.	-	0	ID=3
c2	marcus	cds	1	6	50	+	0	ID=another_thing"""
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    seen_contigs = set()
    for contig, im in cim_tuples:
        seen_contigs.add(contig)
        if contig == "c1":
            seq = C1
        else:
            seq = C2
        cds_df, cc = mu.get_contig_cds_info(
            im, contig, seq, mock_log, mock_log_2
        )
        if contig == "c1":
            pd.testing.assert_frame_equal(
                cds_df,
                pd.DataFrame(
                    {
                        "LeftEnd": [5, 18],
                        "RightEnd": [19, 20],
                        "Strand": ["+", "-"],
                    },
                    index=pd.Index(["hi", "3"]),
                ),
            )
            assert cc.cds2left2counter == {
                "hi": {
                    5: defaultdict(int),
                    8: defaultdict(int),
                    11: defaultdict(int),
                    14: defaultdict(int),
                    17: defaultdict(int),
                },
                "3": {18: defaultdict(int)},
            }
        else:
            pd.testing.assert_frame_equal(
                cds_df,
                pd.DataFrame(
                    {"LeftEnd": [1], "RightEnd": [6], "Strand": ["+"]},
                    index=pd.Index(["another_thing"]),
                ),
            )
            assert cc.cds2left2counter == {
                "another_thing": {
                    1: defaultdict(int),
                    4: defaultdict(int),
                },
            }
    assert seen_contigs == set(["c1", "c2"])
    assert capsys.readouterr().out == (
        "MockLog2: Found 3 features belonging to contig c1; inspecting...\n"
        "MockLog2: Feature 2 on contig c1 has a type that is not in "
        f"{config.CDS_TYPES}; ignoring this feature.\n"
        f"MockLog2: Found 2 features with a type in {config.CDS_TYPES} in "
        "contig c1. Going through alignments...\n"
        "MockLog2: Found 1 feature belonging to contig c2; inspecting...\n"
        f"MockLog2: Found 1 feature with a type in {config.CDS_TYPES} in "
        "contig c2. Going through alignments...\n"
    )


def test_get_contig_cds_info_no_cds_features(capsys):
    gff = "##gff-version 3\nc1	marcus	gene	5	19	.	+	0	ID=hi"
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    no_loops = True
    for contig, im in cim_tuples:
        no_loops = False
        cds_df, cc = mu.get_contig_cds_info(
            im, contig, C1, mock_log, mock_log_2
        )
        assert cds_df is None
        assert cc is None
        assert capsys.readouterr().out == (
            "MockLog2: Found 1 feature belonging to contig c1; "
            "inspecting...\n"
            "MockLog2: Feature hi on contig c1 has a type that is not in "
            f"{config.CDS_TYPES}; ignoring this feature.\n"
            f"MockLog: Found 0 features with a type in {config.CDS_TYPES} in "
            "contig c1. Ignoring this contig.\n"
        )
        break
    if no_loops:
        raise WeirdError


def test_get_contig_cds_info_inconsistent_contig_lengths():
    # The code that makes this check is already tested elsewhere, but it makes
    # sense to verify that we are properly making this check from this point...
    gff = "##gff-version 3\nc1	marcus	cds	5	19	.	+	0	ID=hi"
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    no_loops = True
    for contig, im in cim_tuples:
        no_loops = False
        with pytest.raises(ParameterError) as ei:
            mu.get_contig_cds_info(
                im, contig, skbio.DNA("ACGTACGTACGTAC"), mock_log, mock_log_2
            )
        assert str(ei.value) == (
            "Feature hi on contig c1 has a (1-indexed) end coordinate of 19, "
            "which is greater than the contig's length of 14. We do not "
            "support 'circular' features yet."
        )
        break
    if no_loops:
        raise WeirdError


def check_c1_codoncounts(c1cc):
    # These are 1-indexed codon "left" positions. The "left" doesn't take
    # strand into account.
    # This should directly match the alignment pileup at these codons
    # (unless a gene is on the - strand, in which case each codon should
    # be RC'd here)
    assert c1cc.cds2left2counter == {
        "c1g1": {
            3: {"TGA": 6, "TTA": 3, "TCA": 2, "TAA": 1},
            6: {"CAC": 12},
            9: {"CCA": 5, "CCG": 7},
            12: {"AAC": 5, "AGC": 7},
        },
        # These are reverse complemented from AAA and CCT, respectively.
        "c1g2": {
            16: {"TTT": 12},
            19: {"AGG": 12},
        },
    }
    assert c1cc.cds2strand == {"c1g1": "+", "c1g2": "-"}
    assert c1cc.contig == "c1"
    # there can be extra metadata the skbio.DNA objects have that differ,
    # since C1 was created above by just converting from a str while
    # c1cc.contig_seq was loaded from a FASTA file. but that doesn't matter
    # for our purposes -- just verify that they record the same *sequence*.
    assert str(c1cc.contig_seq) == str(C1)


def test_run_count_good(capsys, tmp_path):
    cdir = tmp_path / "cdir"
    mu.run_count(FASTA, BAM, GFF, cdir, True, mock_log)
    assert sorted(os.listdir(cdir)) == ["c1_3mers.pickle", "c2_3mers.pickle"]

    with open(tmp_path / "cdir" / "c1_3mers.pickle", "rb") as f:
        c1cc = pickle.load(f)
        check_c1_codoncounts(c1cc)

    with open(tmp_path / "cdir" / "c2_3mers.pickle", "rb") as f:
        c2cc = pickle.load(f)
        assert c2cc.cds2left2counter == {"c2g1": {6: {"AGG": 11}}}
        assert c2cc.cds2strand == {"c2g1": "+"}
        assert c2cc.contig == "c2"
        assert str(c2cc.contig_seq) == str(C2)

    assert capsys.readouterr().out == (
        "PREFIX\nMockLog: Loading and checking FASTA and BAM files...\n"
        "MockLog: The FASTA file describes 3 contig(s).\n"
        "MockLog: All of these are included in the BAM file (which has 3 "
        "reference(s)), with the same lengths.\n"
        "PREFIX\nMockLog: Counting aligned 3-mers to codons in coding "
        "sequences in contigs...\n"
        "MockLog: Found 2 features belonging to contig c1; inspecting...\n"
        f"MockLog: Found 2 features with a type in {config.CDS_TYPES} in "
        "contig c1. Going through alignments...\n"
        "MockLog: Wrote out 3-mer count info for contig c1.\n"
        "MockLog: Found 1 feature belonging to contig c2; inspecting...\n"
        f"MockLog: Found 1 feature with a type in {config.CDS_TYPES} in "
        "contig c2. Going through alignments...\n"
        "MockLog: Wrote out 3-mer count info for contig c2.\n"
        "MockLog: Found 1 feature belonging to contig c3; inspecting...\n"
        "MockLog: Feature c3g1 on contig c3 has a type that is not in "
        f"{config.CDS_TYPES}; ignoring this feature.\n"
        f"MockLog: Found 0 features with a type in {config.CDS_TYPES} in "
        "contig c3. Ignoring this contig.\n"
        "MockLog: Done.\n"
    )


def test_run_count_contigs_in_gff_but_not_fasta(capsys, tmp_path):
    cdir = tmp_path / "cdir"
    # make a fasta without c2 or c3
    # (for some reason, skbio doesn't like reading the file when i create it
    # like tmp_path / "sus.fasta" -- so we use os.path.join() like chumps)
    tmp_fasta = os.path.join(tmp_path, "sus.fasta")
    with open(tmp_fasta, "w") as f:
        f.write(">c1\nACTGACACCCAAACCAAACCTAC\n")

    mu.run_count(tmp_fasta, BAM, GFF, cdir, True, mock_log)
    assert sorted(os.listdir(cdir)) == ["c1_3mers.pickle"]

    with open(tmp_path / "cdir" / "c1_3mers.pickle", "rb") as f:
        c1cc = pickle.load(f)
        check_c1_codoncounts(c1cc)

    assert capsys.readouterr().out == (
        "PREFIX\nMockLog: Loading and checking FASTA and BAM files...\n"
        "MockLog: The FASTA file describes 1 contig(s).\n"
        "MockLog: All of these are included in the BAM file (which has 3 "
        "reference(s)), with the same lengths.\n"
        "PREFIX\nMockLog: Counting aligned 3-mers to codons in coding "
        "sequences in contigs...\n"
        "MockLog: Found 2 features belonging to contig c1; inspecting...\n"
        f"MockLog: Found 2 features with a type in {config.CDS_TYPES} in "
        "contig c1. Going through alignments...\n"
        "MockLog: Wrote out 3-mer count info for contig c1.\n"
        "MockLog: Contig c2 is in the GFF3 file but not the FASTA file; "
        "ignoring it.\n"
        "MockLog: Contig c3 is in the GFF3 file but not the FASTA file; "
        "ignoring it.\n"
        "MockLog: Done.\n"
    )


def test_codon_counter_good():
    cc = mu.CodonCounter("c1", skbio.DNA("ACTGACACCCAAACCAAACCTAC"))
    assert str(cc) == "CodonCounter(c1, 23 bp, 0 CDSs)"
    cc.add_cds("cds1", 5, 7, "-")
    assert str(cc) == "CodonCounter(c1, 23 bp, 1 CDS)"
    assert cc.cds2left2counter == {"cds1": {5: defaultdict(int)}}
    assert cc.cds2strand == {"cds1": "-"}
    cc.add_cds("cds2", 6, 11, "+")
    assert str(cc) == "CodonCounter(c1, 23 bp, 2 CDSs)"
    assert cc.cds2left2counter == {
        "cds1": {5: defaultdict(int)},
        "cds2": {6: defaultdict(int), 9: defaultdict(int)},
    }
    assert cc.cds2strand == {"cds1": "-", "cds2": "+"}


def test_codon_counter_add_cds_twice():
    cc = mu.CodonCounter("c1", skbio.DNA("ACTGACACCCAAACCAAACCTAC"))
    assert str(cc) == "CodonCounter(c1, 23 bp, 0 CDSs)"
    cc.add_cds("cds1", 5, 7, "-")
    assert str(cc) == "CodonCounter(c1, 23 bp, 1 CDS)"
    with pytest.raises(WeirdError) as ei:
        cc.add_cds("cds1", 6, 8, "+")
    assert str(ei.value) == (
        "CDS cds1 has already been added to this CodonCounter."
    )


def test_codon_counter_add_count_good():
    cc = mu.CodonCounter("c1", skbio.DNA("ACTGACACCCAAACCAAACCTAC"))
    assert str(cc) == "CodonCounter(c1, 23 bp, 0 CDSs)"
    cc.add_cds("cds1", 5, 7, "-")
    assert str(cc) == "CodonCounter(c1, 23 bp, 1 CDS)"
    assert cc.cds2left2counter == {"cds1": {5: defaultdict(int)}}
    # reverse-complement of ACT is AGT
    cc.add_count("cds1", 5, "ACT")
    assert cc.cds2left2counter == {"cds1": {5: {"AGT": 1}}}
    # ... and the RC of AGT is ACT
    cc.add_count("cds1", 5, "AGT")
    assert cc.cds2left2counter == {"cds1": {5: {"AGT": 1, "ACT": 1}}}
    cc.add_count("cds1", 5, "AGT")
    assert cc.cds2left2counter == {"cds1": {5: {"AGT": 1, "ACT": 2}}}
    cc.add_count("cds1", 5, "AGT")
    assert cc.cds2left2counter == {"cds1": {5: {"AGT": 1, "ACT": 3}}}
    cc.add_count("cds1", 5, "TTA")
    assert cc.cds2left2counter == {"cds1": {5: {"AGT": 1, "ACT": 3, "TAA": 1}}}


def test_codon_counter_add_count_missing_cds_id_or_cpleft():
    cc = mu.CodonCounter("c1", skbio.DNA("ACTGACACCCAAACCAAACCTAC"))
    assert str(cc) == "CodonCounter(c1, 23 bp, 0 CDSs)"
    cc.add_cds("cds1", 5, 7, "-")

    # cds2 hasn't been added to this CodonCounter yet
    with pytest.raises(WeirdError) as ei:
        cc.add_count("cds2", 5, "AGT")
    assert (
        str(ei.value)
        == "No CDS named cds2 has been added to this CodonCounter."
    )

    # although position 6 is present in cds1, it isn't the leftmost position
    # of a CP in this CDS
    with pytest.raises(WeirdError) as ei:
        cc.add_count("cds1", 6, "AGT")
    assert (
        str(ei.value)
        == "6 is not a leftmost position in any codon in CDS cds1."
    )

    # check that both things going wrong at once doesn't "cancel out" lol
    with pytest.raises(WeirdError) as ei:
        cc.add_count("cds2", 6, "AGT")
    assert (
        str(ei.value)
        == "No CDS named cds2 has been added to this CodonCounter."
    )


def test_codon_counter_get_ref_codon_and_aa_good():
    cc = mu.CodonCounter("c1", skbio.DNA("ACTGACACCCAAACCAAACCTAC"))
    assert str(cc) == "CodonCounter(c1, 23 bp, 0 CDSs)"
    cc.add_cds("cds1", 5, 7, "-")
    # TGT codes for Cysteine (C)
    assert cc.get_ref_codon_and_aa("cds1", 5) == ("TGT", "C")

    cc.add_cds("cds1_but_+", 5, 7, "+")
    # ACA codes for Threonine (T)
    assert cc.get_ref_codon_and_aa("cds1_but_+", 5) == ("ACA", "T")


def test_codon_counter_get_ref_codon_and_aa_bad():
    # same error-handler as with add_count(), so we just test one case
    cc = mu.CodonCounter("c1", skbio.DNA("ACTGACACCCAAACCAAACCTAC"))
    assert str(cc) == "CodonCounter(c1, 23 bp, 0 CDSs)"
    cc.add_cds("cds1", 5, 7, "-")
    with pytest.raises(WeirdError) as ei:
        cc.get_ref_codon_and_aa("cds3", 5)
    assert (
        str(ei.value)
        == "No CDS named cds3 has been added to this CodonCounter."
    )
