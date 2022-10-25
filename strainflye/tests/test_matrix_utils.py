import os
import pickle
import skbio
import pytest
import numpy as np
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


def check_run_count_logs(logs):
    assert logs == (
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

    check_run_count_logs(capsys.readouterr().out)


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


def test_get_objs():
    assert mu.get_objs("codon") == config.CODONS
    assert mu.get_objs("aa") == config.AAS
    for bt in ("Codon", "CODON", "AA", "codons", "aas", "lol", "", 123, 0):
        with pytest.raises(WeirdError) as ei:
            mu.get_objs(bt)
        assert str(ei.value) == f"Unrecognized obj_type: {bt}"


def test_get_obj_type_hr():
    assert mu.get_obj_type_hr("codon") == "Codon"
    assert mu.get_obj_type_hr("aa") == "AminoAcid"
    for bt in ("Codon", "CODON", "AA", "codons", "aas", "lol", "", 123, 0):
        with pytest.raises(WeirdError) as ei:
            mu.get_obj_type_hr(bt)
        assert str(ei.value) == f"Unrecognized obj_type: {bt}"


def test_run_fill_no_counts(tmp_path):
    cdir = tmp_path / "empty_counts_dir"
    mdir = tmp_path / "matrix_dir"
    os.makedirs(cdir)
    with pytest.raises(FileNotFoundError) as ei:
        mu.run_fill(cdir, 50, 2, None, "tsv", mdir, True, mock_log)
    assert (
        str(ei.value) == f"Didn't find any 3-mer count information in {cdir}."
    )
    # creation of the output directory should be deferred until just before we
    # are ready to write out the first contig's matrix information
    assert not os.path.exists(mdir)


def check_matrix_fill_output(mdir):
    # check that output directories / files look as expected
    assert os.path.exists(mdir)
    # no c3 because no CDSs in it
    contigs = ["c1", "c2"]
    assert set(contigs).issubset(set(os.listdir(mdir)))
    for contig in contigs:
        assert sorted(os.listdir(mdir / contig)) == [
            f"{contig}_aa_matrix.tsv",
            f"{contig}_aa_refcounts.tsv",
            f"{contig}_codon_matrix.tsv",
            f"{contig}_codon_refcounts.tsv",
        ]

    # check c1's output info in depth
    # for reference, here's all reads aligned to positions 3 -- 14 in c1 (c1g1)
    # ---------------
    # TGA CAC CCA AAC <-- Reference c1 sequence
    # ---------------     (Note that c1g1 is a + strand CDS)
    # TGA CAC CCA AAC
    # TTA CAC CCA AAC
    # TCA CAC CCA AAC
    # TAA CAC CCA AAC
    # TCA CAC CCA AAC
    # TTA CAC CCG AGC
    # TTA CAC CCG AGC
    # TGA CAC CCG AGC
    # TGA CAC CCG AGC
    # TGA CAC CCG AGC
    # TGA CAC CCG AGC
    # TGA CAC CCG AGC
    #
    # ... and all aligned to positions 16 -- 21 in c1 (c1g2, - strand)
    # -------
    # AAA CCT
    # -------
    # AAA CCT
    # AAA CCT
    # AAA CCT
    # AAA CCT
    # AAA CCT
    # AAA CCT
    # AAA CCT
    # AAA CCT
    # AAA CCT
    # AAA CCT
    # AAA CCT
    # AAA CCT
    #
    # (... of course, these are fake genes -- obviously, real genes wouldn't
    # just start with a stop codon. we don't actually test the composition of
    # CDSs because we assume that whatever tool has produced them knows
    # better.)

    # 1. check ref counts -- the counts of codons and amino acids in c1's CDSs
    c1cr = pd.read_csv(
        mdir / "c1" / "c1_codon_refcounts.tsv", sep="\t", index_col=0
    )
    exp_codon_counts = [0] * len(config.CODONS)
    # c1g1 codons
    exp_codon_counts[config.CODONS.index("TGA")] = 1
    exp_codon_counts[config.CODONS.index("CAC")] = 1
    exp_codon_counts[config.CODONS.index("CCA")] = 1
    exp_codon_counts[config.CODONS.index("AAC")] = 1
    # c1g2 codons (reverse-complemented)
    exp_codon_counts[config.CODONS.index("AGG")] = 1
    exp_codon_counts[config.CODONS.index("TTT")] = 1
    pd.testing.assert_frame_equal(
        c1cr,
        pd.DataFrame(
            {"Count": exp_codon_counts},
            index=pd.Index(config.CODONS, name="Codon"),
        ),
    )
    c1ar = pd.read_csv(
        mdir / "c1" / "c1_aa_refcounts.tsv", sep="\t", index_col=0
    )
    exp_aa_counts = [0] * len(config.AAS)
    # c1g1 AAs (same order as above)
    exp_aa_counts[config.AAS.index("*")] = 1
    exp_aa_counts[config.AAS.index("H")] = 1
    exp_aa_counts[config.AAS.index("P")] = 1
    exp_aa_counts[config.AAS.index("N")] = 1
    # c1g2 AAs
    exp_aa_counts[config.AAS.index("R")] = 1
    exp_aa_counts[config.AAS.index("F")] = 1
    pd.testing.assert_frame_equal(
        c1ar,
        pd.DataFrame(
            {"Count": exp_aa_counts},
            index=pd.Index(config.AAS, name="AminoAcid"),
        ),
    )

    # 2. check actual mutation matrices
    c1cm = pd.read_csv(
        mdir / "c1" / "c1_codon_matrix.tsv", sep="\t", index_col=0
    )
    exp_matrix = pd.DataFrame(
        {c: [0] * len(config.CODONS) for c in config.CODONS},
        index=pd.Index(config.CODONS),
        columns=config.CODONS,
    )
    for codon in config.CODONS:
        exp_matrix[codon][codon] = np.nan
    # in order to modify a df, we need to access the column first -- so this is
    # backwards, we're making note of a mutation from TGA --> TTA
    exp_matrix["TTA"]["TGA"] = 1
    # Two of the codons in c1g1 are mutated, but both are "unreasonable" codon
    # mutations -- CCA and AAC are not the most common 3-mers in the alignment
    # to their codon, interestingly (i forgor about this lol). So there should
    # be just one entry in the codon mutation matrix.
    pd.testing.assert_frame_equal(c1cm, exp_matrix)

    c1am = pd.read_csv(mdir / "c1" / "c1_aa_matrix.tsv", sep="\t", index_col=0)
    exp_matrix = pd.DataFrame(
        {c: [0] * len(config.AAS) for c in config.AAS},
        index=pd.Index(config.AAS),
        columns=config.AAS,
    )
    # again, this is backwards, and really represents * --> L
    exp_matrix["L"]["*"] = 1
    pd.testing.assert_frame_equal(c1am, exp_matrix)

    # also check c2 stuff agoihsdoifaoidfj
    # it doesn't have any mutations so this is simpler
    c2cr = pd.read_csv(
        mdir / "c2" / "c2_codon_refcounts.tsv", sep="\t", index_col=0
    )
    exp_c2cr = pd.DataFrame(
        {"Count": [0] * len(config.CODONS)},
        index=pd.Index(config.CODONS, name="Codon"),
    )
    exp_c2cr["Count"]["AGG"] = 1
    pd.testing.assert_frame_equal(c2cr, exp_c2cr)

    c2ar = pd.read_csv(
        mdir / "c2" / "c2_aa_refcounts.tsv", sep="\t", index_col=0
    )
    exp_c2ar = pd.DataFrame(
        {"Count": [0] * len(config.AAS)},
        index=pd.Index(config.AAS, name="AminoAcid"),
    )
    exp_c2ar["Count"]["R"] = 1
    pd.testing.assert_frame_equal(c2ar, exp_c2ar)
    c2cm = pd.read_csv(
        mdir / "c2" / "c2_codon_matrix.tsv", sep="\t", index_col=0
    )
    exp_matrix = pd.DataFrame(
        {c: [0] * len(config.CODONS) for c in config.CODONS},
        index=pd.Index(config.CODONS),
        columns=config.CODONS,
    )
    for codon in config.CODONS:
        exp_matrix[codon][codon] = np.nan
    pd.testing.assert_frame_equal(c2cm, exp_matrix)

    c2am = pd.read_csv(mdir / "c2" / "c2_aa_matrix.tsv", sep="\t", index_col=0)
    exp_matrix = pd.DataFrame(
        {c: [0] * len(config.AAS) for c in config.AAS},
        index=pd.Index(config.AAS),
        columns=config.AAS,
    )
    pd.testing.assert_frame_equal(c2am, exp_matrix)


def test_matrix_integration(capsys, tmp_path):
    cdir = tmp_path / "cdir"
    mu.run_count(FASTA, BAM, GFF, cdir, True, mock_log)
    check_run_count_logs(capsys.readouterr().out)

    mdir = tmp_path / "mdir"
    # call matrix r-mutations for r = 1 (simplifies figuring out expected
    # output)
    mu.run_fill(cdir, None, 2, 1, "tsv", mdir, True, mock_log)
    check_matrix_fill_output(mdir)
    assert capsys.readouterr().out == (
        "PREFIX\nMockLog: Performing codon r-mutation calling (r = 1).\n"
        "PREFIX\nMockLog: Going through aligned 3-mer counts and creating "
        "matrices...\n"
        "MockLog: Creating matrices for contig c1...\n"
        "MockLog: Created an output directory for contig c1.\n"
        "MockLog: Creating matrices for contig c2...\n"
        "MockLog: Created an output directory for contig c2.\n"
        "MockLog: Done.\n"
    )


def test_matrix_integration_extra_file_in_matrix_dir_pmuts(capsys, tmp_path):
    # just testing other branches in the code...
    cdir = tmp_path / "cdir"
    mu.run_count(FASTA, BAM, GFF, cdir, True, mock_log)
    check_run_count_logs(capsys.readouterr().out)

    # add an extra file to the folder of counts -- make sure that we ignore it,
    # and warn the user about it
    with open(cdir / "sus_file.txt", "w") as f:
        f.write("I'm inconvenient!")

    mdir = tmp_path / "mdir"

    # call matrix p-mutations for p = 25% -- the lone codon mutation in c1 (TGA
    # --> TTA) occurs with frequency 3/12 = 25%, so this should pass (and the
    # mutation calls should thus be indistinguishable from
    # test_matrix_integration() above)
    mu.run_fill(cdir, 2500, 2, None, "tsv", mdir, True, mock_log)

    check_matrix_fill_output(mdir)
    # make sure sus_file didn't get overwritten or anything
    with open(cdir / "sus_file.txt", "r") as f:
        assert f.read() == "I'm inconvenient!"

    assert capsys.readouterr().out == (
        "PREFIX\nMockLog: Performing codon p-mutation calling (p = 25.00%).\n"
        "PREFIX\nMockLog: Going through aligned 3-mer counts and creating "
        "matrices...\n"
        "MockLog: Creating matrices for contig c1...\n"
        "MockLog: Created an output directory for contig c1.\n"
        "MockLog: Creating matrices for contig c2...\n"
        "MockLog: Created an output directory for contig c2.\n"
        "MockLog: Warning: found an unexpectedly named file (sus_file.txt) in "
        f"{cdir}. Ignoring it.\n"
        "MockLog: Done.\n"
    )


def test_matrix_integration_bad_output_format(tmp_path):
    cdir = tmp_path / "cdir"
    mu.run_count(FASTA, BAM, GFF, cdir, True, mock_log)
    mdir = tmp_path / "mdir"
    with pytest.raises(WeirdError) as ei:
        mu.run_fill(cdir, 2500, 2, None, "TSV", mdir, True, mock_log)
    assert str(ei.value) == 'Unrecognized output format: "TSV"'
