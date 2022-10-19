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


def test_get_contig_cds_info_good():
    gff = "##gff-version 3\nc1	marcus	cds	5	19	.	+	0	ID=hi"
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    no_loops = True
    for contig, im in cim_tuples:
        no_loops = False
        cds_df, fid2codon2alignedcodons = mu.get_contig_cds_info(
            im, contig, {"c1": 23}, mock_log, mock_log_2
        )
        pd.testing.assert_frame_equal(
            cds_df,
            pd.DataFrame(
                {"LeftEnd": [5], "RightEnd": [19], "Strand": "+"},
                index=pd.Index(["hi"]),
            ),
        )
        assert fid2codon2alignedcodons == {
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


def test_get_contig_cds_info_multi_features():
    gff = """##gff-version 3
c1	marcus	cds	5	19	.	+	0	ID=hi
c1	marcus	gene	1	5	.	+	0	ID=2
c1	marcus	CDS	18	20	.	-	0	ID=3
c2	marcus	cds	1	6	50	+	0	ID=another_thing"""
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    seen_contigs = set()
    for contig, im in cim_tuples:
        seen_contigs.add(contig)
        cds_df, fid2codon2alignedcodons = mu.get_contig_cds_info(
            im, contig, {"c1": 23, "c2": 12}, mock_log, mock_log_2
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
            assert fid2codon2alignedcodons == {
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
            assert fid2codon2alignedcodons == {
                "another_thing": {
                    1: defaultdict(int),
                    4: defaultdict(int),
                },
            }
    assert seen_contigs == set(["c1", "c2"])


def test_get_contig_cds_info_contig_not_in_name2len(capsys):
    gff = "##gff-version 3\nc1	marcus	cds	5	19	.	+	0	ID=hi"
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    no_loops = True
    for contig, im in cim_tuples:
        no_loops = False
        cds_df, fid2codon2alignedcodons = mu.get_contig_cds_info(
            im, contig, {"c3": 16}, mock_log, mock_log_2
        )
        assert cds_df is None
        assert fid2codon2alignedcodons is None
        assert capsys.readouterr().out == (
            "MockLog2: Found 1 feature belonging to sequence c1 in the GFF3 "
            "file. Ignoring this sequence and its feature, since c1 isn't in "
            "the FASTA file.\n"
        )
        break
    if no_loops:
        raise WeirdError


def test_get_contig_cds_info_no_cds_features(capsys):
    gff = "##gff-version 3\nc1	marcus	gene	5	19	.	+	0	ID=hi"
    cim_tuples = skbio.io.read(sio(gff), format="gff3")
    no_loops = True
    for contig, im in cim_tuples:
        no_loops = False
        cds_df, fid2codon2alignedcodons = mu.get_contig_cds_info(
            im, contig, {"c1": 23}, mock_log, mock_log_2
        )
        assert cds_df is None
        assert fid2codon2alignedcodons is None
        assert capsys.readouterr().out == (
            "MockLog2: Found 1 feature belonging to contig c1; "
            "inspecting...\n"
            "MockLog2: Feature hi on contig c1 has a type that is not in "
            f"{config.CDS_TYPES}; ignoring it.\n"
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
                im, contig, {"c1": 14}, mock_log, mock_log_2
            )
        assert str(ei.value) == (
            "Feature hi on contig c1 has a (1-indexed) end coordinate of 19, "
            "which is greater than the contig's length of 14. We do not "
            "support 'circular' features yet."
        )
        break
    if no_loops:
        raise WeirdError


def test_run_count_good(capsys, tmp_path):
    cdir = tmp_path / "cdir"
    mu.run_count(FASTA, BAM, GFF, cdir, True, mock_log)
    assert sorted(os.listdir(cdir)) == ["c1_3mers.pickle", "c2_3mers.pickle"]

    with open(tmp_path / "cdir" / "c1_3mers.pickle", "rb") as f:
        c1data = pickle.load(f)
        # These are 1-indexed codon "left" positions. The "left" doesn't take
        # strand into account.
        # This should directly match the alignment pileup at these codons
        # (unless a gene is on the - strand, in which case each codon should
        # be RC'd here)
        assert c1data == {
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

    with open(tmp_path / "cdir" / "c2_3mers.pickle", "rb") as f:
        c1data = pickle.load(f)
        assert c1data == {"c2g1": {6: {"AGG": 11}}}
