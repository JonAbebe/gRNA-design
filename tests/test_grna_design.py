"""Tests for grna_design.

``test_matches_golden_output`` is the regression guard: it pins the full
candidate list for the three example exons, so any change in behaviour shows
up as a diff.  ``test_reproduces_legacy_candidates`` additionally checks the
candidates the original script reported are all still found.
"""
import json
import sys
from dataclasses import asdict
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import grna_design as gd

DATA = Path(__file__).resolve().parent.parent / "data"
GOLDEN = Path(__file__).resolve().parent / "golden"
EXONS = ("exon44", "exon45", "exon53")


@pytest.fixture(scope="module")
def candidates():
    return {
        exon: gd.design(gd.read_sequence(DATA / f"{exon}.fasta")) for exon in EXONS
    }


# --------------------------------------------------------------------------
# Sequence preparation
# --------------------------------------------------------------------------


def test_clean_sequence_strips_numbers_and_whitespace():
    raw = "1\tacgt acgt\n101\tACGT\r\n"
    assert gd.clean_sequence(raw) == "acgtacgtACGT"


def test_clean_sequence_preserves_case():
    assert gd.clean_sequence("acgtACGT") == "acgtACGT"


def test_strand_sequences():
    seqs = gd.strand_sequences("acgtAAGG")
    assert seqs["top"] == "acgtAAGG"
    assert seqs["complement"] == "tgcaTTCC"
    assert seqs["reverse_complement"] == "CCTTacgt"


def test_exon_start_finds_first_uppercase():
    assert gd.exon_start("acgtACGT") == 4


def test_exon_start_without_exon_raises():
    with pytest.raises(ValueError):
        gd.exon_start("acgt")


def test_splice_site_positions():
    # intron, exon, intron -- the exon starts at index 6 on every strand.
    sites = gd.splice_site_positions(gd.strand_sequences("ttttagACGTgtaaaa"))
    assert sites["top"] == (4,)  # the 'a' of the acceptor 'ag'
    assert sites["complement"] == (5,)
    assert sites["reverse_complement"] == (5, 4)


# --------------------------------------------------------------------------
# PAM handling
# --------------------------------------------------------------------------


def test_expand_pam():
    assert gd.expand_pam("NGG") == ["AGG", "CGG", "GGG", "TGG"]
    assert len(gd.expand_pam("NGAN")) == 16


def test_pam_catalog_sizes_and_membership():
    catalog = gd.pam_catalog()
    assert catalog["spcas9"] == sorted(
        {*gd.expand_pam("NGG"), *gd.expand_pam("NAG")}
    )
    assert len(catalog["xcas9"]) == 16  # NGN
    # NGAN and NGNG each expand to 16 and overlap at the four NGAG sequences.
    assert len(catalog["VQR"]) == 28
    assert len(set(catalog["VQR"])) == len(catalog["VQR"])


def test_complement_pams_are_reversed_and_lowercase():
    reversed_pams = gd.complement_pams({"spcas9": ["AGG"]})
    assert reversed_pams["spcas9"] == ["gga"]


def test_pam_positions_finds_overlapping_matches():
    # "TGGG" holds NGG at index 0 (TGG) and index 1 (GGG).
    assert gd.pam_positions("TGGG", {"spcas9": ["TGG", "GGG"]}) == [0, 1]


def test_targeting_pam_positions_keeps_only_in_range():
    # A site 15 nt upstream is in range; one 5 nt upstream is not.
    assert gd.targeting_pam_positions([20], (5,), "top") == [20]
    assert gd.targeting_pam_positions([20], (15,), "top") == []


def test_targeting_pam_positions_complement_looks_downstream():
    assert gd.targeting_pam_positions([20], (35,), "complement") == [20]
    assert gd.targeting_pam_positions([20], (5,), "complement") == []


def test_nucleases_for_pam_matches_motifs_as_substrings():
    catalog = {"spcas9": ["AGG"], "xcas9": ["AGG", "TTT"]}
    assert gd.nucleases_for_pam("AGGC", catalog) == {
        "spcas9": ["AGG"],
        "xcas9": ["AGG"],
    }


def test_nucleases_for_pam_omits_nucleases_without_a_match():
    assert gd.nucleases_for_pam("AAAA", {"spcas9": ["AGG"]}) == {}


# --------------------------------------------------------------------------
# Guides and editors
# --------------------------------------------------------------------------


def test_build_guides_skips_truncated_guides():
    # A PAM at index 5 leaves too little sequence upstream for a 20 nt guide.
    sequence = "acgt" * 10
    assert gd.build_guides(sequence, "top", [5], {"x": ["cgt"]}) == []


def test_build_guides_returns_guide_and_pam():
    sequence = "a" * 20 + "TGGT" + "c" * 10
    guides = gd.build_guides(sequence, "top", [20], {"spcas9": ["TGG"]})
    assert len(guides) == 1
    assert guides[0].sequence == "a" * 20
    assert guides[0].pam == "TGGT"
    assert guides[0].nucleases == {"spcas9": ["TGG"]}


def test_editable_positions_finds_acceptor_bases():
    guide = "ccccccccccccccccccag"  # exon would start just past the guide
    assert gd.editable_positions(guide) == []
    # With an exonic base at index 18, the 'a' sits at 16 and 'c' at 17.
    assert gd.editable_positions("c" * 16 + "ac" + "GT") == [16, 17]


def test_editable_positions_without_exon_returns_empty():
    assert gd.editable_positions("acgt" * 5) == []


def test_evaluate_guide_picks_the_first_matching_editor():
    # The exon starts at index 7, putting the target 'a' at index 5 -- inside
    # ABE7.10's window, which is tried before ABE7.9's.
    guide = gd.Guide("top", 30, "TGGT", "ccccc" + "a" + "t" + "G" * 13)
    candidates = gd.evaluate_guide(guide)
    assert [c.editor for c in candidates] == ["ABE7.10"]
    assert candidates[0].target == "a6"
    assert candidates[0].offtargets == []


def test_evaluate_guide_reports_offtargets_in_window():
    # Target 'a' at index 3, with a second editable A at index 5.
    guide = gd.Guide("top", 30, "TGGT", "ccc" + "a" + "t" + "A" + "C" * 14)
    (candidate,) = gd.evaluate_guide(guide)
    assert candidate.target == "a4"
    assert candidate.offtargets == [5]


def test_editor_windows_are_documented_constants():
    assert gd.ABE7_10.window == (3, 7)
    assert gd.ABE7_9.window == (3, 9)
    assert gd.BE4.window == (3, 9)
    assert gd.BE4_COMPLEMENT.window == (11, 17)


# --------------------------------------------------------------------------
# Regression
# --------------------------------------------------------------------------


def test_matches_golden_output(candidates):
    expected = json.loads((GOLDEN / "candidates.json").read_text())
    actual = {exon: [asdict(c) for c in cs] for exon, cs in candidates.items()}
    assert actual == expected


def test_reproduces_legacy_candidates(candidates):
    """Every candidate the original script reported is still found.

    The original returned from inside its loop, so it reported only one
    candidate per strand, and which one varied with the hash seed.  The file
    pools what it reported over eight runs.
    """
    legacy = json.loads((GOLDEN / "legacy_candidates.json").read_text())
    for exon, entries in legacy.items():
        found = {
            (
                c.strand,
                c.guide,
                c.editor,
                c.target,
                tuple(f"{c.target_base}{i + 1}" for i in c.offtargets),
            )
            for c in candidates[exon]
        }
        for entry in entries:
            key = (
                entry["strand"],
                entry["guide"],
                entry["editor"],
                entry["target"],
                tuple(entry["offtargets"]),
            )
            assert key in found, f"{exon}: lost legacy candidate {key}"


def test_design_is_deterministic():
    sequence = gd.read_sequence(DATA / "exon45.fasta")
    first = [asdict(c) for c in gd.design(sequence)]
    second = [asdict(c) for c in gd.design(sequence)]
    assert first == second


def test_design_finds_more_than_the_legacy_script(candidates):
    legacy = json.loads((GOLDEN / "legacy_candidates.json").read_text())
    for exon, entries in legacy.items():
        assert len(candidates[exon]) >= len(entries)


# --------------------------------------------------------------------------
# Command line interface
# --------------------------------------------------------------------------


def test_cli_text_output(capsys):
    assert gd.main(["-f", str(DATA / "exon45.fasta")]) == 0
    out = capsys.readouterr().out
    assert "Base editor:" in out
    assert "ABE7.9" in out


def test_cli_tsv_output(capsys):
    assert gd.main(["-f", str(DATA / "exon45.fasta"), "--format", "tsv"]) == 0
    lines = capsys.readouterr().out.strip().splitlines()
    assert lines[0].split("\t")[0] == "strand"
    assert len(lines) == 7  # header plus six candidates


def test_cli_json_output(capsys):
    assert gd.main(["-f", str(DATA / "exon53.fasta"), "--format", "json"]) == 0
    parsed = json.loads(capsys.readouterr().out)
    assert len(parsed) == 4
    assert parsed[0]["editor"]


def test_cli_accepts_a_raw_sequence(capsys):
    sequence = gd.read_sequence(DATA / "exon44.fasta")
    assert gd.main(["-s", sequence, "--format", "tsv"]) == 0
    assert len(capsys.readouterr().out.strip().splitlines()) == 5


def test_cli_rejects_a_sequence_without_an_exon():
    with pytest.raises(SystemExit):
        gd.main(["-s", "acgtacgtacgt"])
