#!/usr/bin/env python3
"""Design base-editing guide RNAs that disrupt a splice-acceptor site.

The input is a genomic sequence in which exonic bases are uppercase and
intronic bases are lowercase (the usual Ensembl/UCSC convention).  The target
is the ``ag`` immediately 5' of the exon, so the editable base is located
relative to the first uppercase base of the sequence.

For every PAM that sits the right distance from that target the module builds
the adjacent 20 nt guide, records which nuclease recognises the PAM, and
reports the base editors whose editing window covers the target along with any
additional editable base inside the same window (a potential off-target edit).

The library functions are pure: they take and return data, and never print.
Formatting lives in :func:`format_candidates` and the command line interface.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass, field
from itertools import product
from pathlib import Path

from Bio.Seq import Seq

STRANDS = ("top", "complement", "reverse_complement")

GUIDE_LENGTH = 20

# The editable base must lie 12-17 nt upstream of the PAM, which places it at
# guide positions 3-8 -- the region the base editors below act on.
MIN_PAM_DISTANCE = 12
MAX_PAM_DISTANCE = 17

# PAM motifs per nuclease, written with N for any base.
NUCLEASE_PAMS = {
    "spcas9": ("NGG", "NAG"),
    "xcas9": ("NGN",),
    "VQR": ("NGAN", "NGNG"),
}

# Offsets from the first exonic base to the editable base of the splice
# acceptor, per strand.  The reverse complement is checked at both offsets.
SITE_OFFSETS = {
    "top": (-2,),
    "complement": (-1,),
    "reverse_complement": (-1, -2),
}

# Slice offsets relative to the PAM match position, as (start, stop) pairs.
# On the complement strand the sequence runs antiparallel to the index space,
# so the PAM reads one base earlier and the guide follows rather than precedes
# it.
PAM_SLICE = {
    "top": (0, 4),
    "complement": (-1, 2),
    "reverse_complement": (0, 4),
}
GUIDE_SLICE = {
    "top": (-GUIDE_LENGTH, 0),
    "complement": (2, 2 + GUIDE_LENGTH),
    "reverse_complement": (-GUIDE_LENGTH, 0),
}

# The editable base of the splice acceptor: 'a' from the AG dinucleotide, or
# the 'c' that pairs with its G.  Keyed by offset from the first exonic base.
TARGET_BASES = (("a", -2), ("c", -1))


@dataclass(frozen=True)
class BaseEditor:
    """A base editor and the guide positions it acts on.

    ``positions`` are the guide positions whose base this editor can change.
    ``window`` is the span reported as the editing window, and ``scan`` the
    span searched for additional editable bases.  ``scan`` reaches one base
    beyond ``window`` for the three editors below; that asymmetry is carried
    over from the original script and has not been verified against the
    source publications.
    """

    name: str
    target_base: str
    positions: range
    window: tuple[int, int]
    scan: range


ABE7_10 = BaseEditor("ABE7.10", "a", range(3, 7), (3, 7), range(3, 8))
ABE7_9 = BaseEditor("ABE7.9", "a", range(3, 10), (3, 9), range(3, 10))
BE4 = BaseEditor("BE4", "c", range(3, 9), (3, 9), range(3, 10))
# Guides on the complement strand run the other way, so the window sits
# further along the guide.
BE4_COMPLEMENT = BaseEditor("BE4", "c", range(11, 17), (11, 17), range(11, 17))

# Order matters: the first editor whose window covers the target base wins.
EDITORS = {
    "top": (ABE7_10, ABE7_9, BE4),
    "complement": (BE4_COMPLEMENT,),
    "reverse_complement": (ABE7_10, ABE7_9, BE4),
}


@dataclass
class Guide:
    """A 20 nt guide and the PAM it was built from."""

    strand: str
    pam_index: int
    pam: str
    sequence: str
    nucleases: dict[str, list[str]] = field(default_factory=dict)


@dataclass
class Candidate:
    """One guide, one editable base, one base editor."""

    strand: str
    pam_index: int
    pam: str
    nucleases: dict[str, list[str]]
    guide: str
    position: int
    target_base: str
    editor: str
    window: str
    offtargets: list[int]

    @property
    def target(self) -> str:
        """The target base with its 1-based position, e.g. ``a5``."""
        return f"{self.target_base}{self.position + 1}"


# --------------------------------------------------------------------------
# Sequence preparation
# --------------------------------------------------------------------------


def clean_sequence(text: str) -> str:
    """Strip the position numbers, whitespace and newlines out of a sequence.

    Case is preserved, since uppercase marks the exon.
    """
    return "".join(c for c in text if not c.isnumeric() and not c.isspace())


def strand_sequences(sequence: str) -> dict[str, str]:
    """Return the top, complement and reverse complement strands."""
    return {
        "top": sequence,
        "complement": str(Seq(sequence).complement()),
        "reverse_complement": str(Seq(sequence).reverse_complement()),
    }


def exon_start(sequence: str) -> int:
    """Index of the first exonic (uppercase) base.

    Raises:
        ValueError: if the sequence has no uppercase base.
    """
    for i, base in enumerate(sequence):
        if base.isupper():
            return i
    raise ValueError("sequence contains no uppercase (exonic) base")


def splice_site_positions(sequences: dict[str, str]) -> dict[str, tuple[int, ...]]:
    """Position of the editable splice-acceptor base on each strand."""
    return {
        strand: tuple(exon_start(seq) + off for off in SITE_OFFSETS[strand])
        for strand, seq in sequences.items()
    }


# --------------------------------------------------------------------------
# PAM discovery
# --------------------------------------------------------------------------


def expand_pam(pattern: str) -> list[str]:
    """Expand a PAM motif containing N into every matching sequence."""
    choices = ["ACGT" if base == "N" else base for base in pattern]
    return ["".join(combo) for combo in product(*choices)]


def pam_catalog() -> dict[str, list[str]]:
    """Every PAM sequence recognised by each nuclease."""
    catalog = {}
    for nuclease, patterns in NUCLEASE_PAMS.items():
        expanded = [pam for pattern in patterns for pam in expand_pam(pattern)]
        # NGAG matches two of VQR's motifs, so de-duplicate.
        catalog[nuclease] = sorted(dict.fromkeys(expanded))
    return catalog


def complement_pams(catalog: dict[str, list[str]]) -> dict[str, list[str]]:
    """The same PAMs as read on the complement strand.

    The complement strand is stored 3'->5' relative to the index space, so a
    PAM reads reversed.  Lowercasing restricts the match to intronic sequence,
    which is where a PAM upstream of the splice acceptor falls; this mirrors
    the uppercase (exonic) PAMs used on the other two strands.
    """
    return {
        nuclease: sorted(pam[::-1].lower() for pam in pams)
        for nuclease, pams in catalog.items()
    }


def pam_positions(sequence: str, catalog: dict[str, list[str]]) -> list[int]:
    """Every position in ``sequence`` where a catalogued PAM starts."""
    found = set()
    for pams in catalog.values():
        for pam in pams:
            start = sequence.find(pam)
            while start != -1:
                found.add(start)
                start = sequence.find(pam, start + 1)
    return sorted(found)


def targeting_pam_positions(
    positions: list[int], sites: tuple[int, ...], strand: str
) -> list[int]:
    """Keep the PAM positions whose guide would cover an editable base."""
    if strand == "complement":
        spans = [(pos + MIN_PAM_DISTANCE, pos + MAX_PAM_DISTANCE) for pos in positions]
    else:
        spans = [(pos - MAX_PAM_DISTANCE, pos - MIN_PAM_DISTANCE) for pos in positions]
    return [
        pos
        for pos, (low, high) in zip(positions, spans)
        if any(low <= site <= high for site in sites)
    ]


def nucleases_for_pam(pam: str, catalog: dict[str, list[str]]) -> dict[str, list[str]]:
    """Which nucleases recognise a PAM, and the motifs they matched.

    The slice taken from the sequence can be longer than a given nuclease's
    PAM, so motifs are matched as substrings.
    """
    matches = {}
    for nuclease, pams in catalog.items():
        hits = sorted({p for p in pams if p in pam})
        if hits:
            matches[nuclease] = hits
    return matches


# --------------------------------------------------------------------------
# Guide construction and scoring
# --------------------------------------------------------------------------


def build_guides(
    sequence: str, strand: str, positions: list[int], catalog: dict[str, list[str]]
) -> list[Guide]:
    """Build the guide adjacent to each PAM position."""
    pam_start, pam_stop = PAM_SLICE[strand]
    guide_start, guide_stop = GUIDE_SLICE[strand]
    guides = []
    for pos in positions:
        if pos + guide_start < 0:
            # Not enough sequence upstream for a full guide.
            continue
        pam = sequence[pos + pam_start : pos + pam_stop]
        guide = sequence[pos + guide_start : pos + guide_stop]
        if len(guide) != GUIDE_LENGTH:
            continue
        nucleases = nucleases_for_pam(pam, catalog)
        if not nucleases:
            continue
        guides.append(Guide(strand, pos, pam, guide, nucleases))
    return guides


def editable_positions(guide: str) -> list[int]:
    """Guide positions holding an editable splice-acceptor base."""
    try:
        start = exon_start(guide)
    except ValueError:
        return []
    positions = []
    for base, offset in TARGET_BASES:
        position = start + offset
        if position >= 0 and guide[position] == base:
            positions.append(position)
    return positions


def evaluate_guide(guide: Guide) -> list[Candidate]:
    """Every (editable base, base editor) pairing a guide supports."""
    candidates = []
    for position in editable_positions(guide.sequence):
        base = guide.sequence[position].lower()
        for editor in EDITORS[guide.strand]:
            if base != editor.target_base or position not in editor.positions:
                continue
            offtargets = [
                i
                for i in editor.scan
                if i != position and guide.sequence[i].lower() == editor.target_base
            ]
            window_start, window_stop = editor.window
            candidates.append(
                Candidate(
                    strand=guide.strand,
                    pam_index=guide.pam_index,
                    pam=guide.pam,
                    nucleases=guide.nucleases,
                    guide=guide.sequence,
                    position=position,
                    target_base=base,
                    editor=editor.name,
                    window=guide.sequence[window_start:window_stop],
                    offtargets=offtargets,
                )
            )
            # The first editor whose window covers the base wins.
            break
    return candidates


def design(sequence: str) -> list[Candidate]:
    """Design every guide/editor candidate for a cleaned sequence."""
    sequences = strand_sequences(sequence)
    sites = splice_site_positions(sequences)
    catalog = pam_catalog()
    reversed_catalog = complement_pams(catalog)

    candidates = []
    for strand in STRANDS:
        strand_catalog = reversed_catalog if strand == "complement" else catalog
        positions = pam_positions(sequences[strand], strand_catalog)
        positions = targeting_pam_positions(positions, sites[strand], strand)
        for guide in build_guides(sequences[strand], strand, positions, strand_catalog):
            candidates.extend(evaluate_guide(guide))
    return candidates


# --------------------------------------------------------------------------
# Input and output
# --------------------------------------------------------------------------


def read_sequence(path: Path) -> str:
    """Read a cleaned sequence from a FASTA or plain-text file."""
    lines = [
        line for line in path.read_text().splitlines() if not line.startswith(">")
    ]
    return clean_sequence("\n".join(lines))


def format_candidates(candidates: list[Candidate]) -> str:
    """Render candidates as a readable report."""
    if not candidates:
        return "No candidate guides found.\n"
    blocks = []
    for strand in STRANDS:
        for candidate in (c for c in candidates if c.strand == strand):
            offtargets = (
                ", ".join(
                    f"{candidate.target_base}{i + 1}" for i in candidate.offtargets
                )
                or "none"
            )
            nucleases = ", ".join(
                f"{name} ({'/'.join(pams)})"
                for name, pams in sorted(candidate.nucleases.items())
            )
            blocks.append(
                f"Strand:        {candidate.strand}\n"
                f"PAM:           {candidate.pam} at {candidate.pam_index}\n"
                f"Nuclease(s):   {nucleases}\n"
                f"Guide:         {candidate.guide}\n"
                f"Base editor:   {candidate.editor}\n"
                f"Target:        {candidate.target}\n"
                f"Window:        {candidate.window}\n"
                f"Off-target(s): {offtargets}\n"
            )
    return "\n".join(blocks)


def format_tsv(candidates: list[Candidate]) -> str:
    """Render candidates as a tab-separated table."""
    columns = [
        "strand",
        "pam_index",
        "pam",
        "nucleases",
        "guide",
        "editor",
        "target",
        "window",
        "offtargets",
    ]
    rows = ["\t".join(columns)]
    for c in candidates:
        rows.append(
            "\t".join(
                [
                    c.strand,
                    str(c.pam_index),
                    c.pam,
                    ";".join(sorted(c.nucleases)),
                    c.guide,
                    c.editor,
                    c.target,
                    c.window,
                    ";".join(f"{c.target_base}{i + 1}" for i in c.offtargets),
                ]
            )
        )
    return "\n".join(rows) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("-f", "--file", type=Path, help="FASTA or plain-text sequence")
    source.add_argument("-s", "--sequence", help="sequence given directly")
    parser.add_argument(
        "--format",
        choices=("text", "tsv", "json"),
        default="text",
        help="output format (default: text)",
    )
    args = parser.parse_args(argv)

    sequence = (
        read_sequence(args.file) if args.file else clean_sequence(args.sequence)
    )
    try:
        candidates = design(sequence)
    except ValueError as exc:
        parser.error(str(exc))

    if args.format == "json":
        print(json.dumps([asdict(c) for c in candidates], indent=2))
    elif args.format == "tsv":
        print(format_tsv(candidates), end="")
    else:
        print(format_candidates(candidates), end="")
    return 0


if __name__ == "__main__":
    sys.exit(main())
