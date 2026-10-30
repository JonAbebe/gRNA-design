# gRNA-design

Designs base-editing guide RNAs that disrupt a **splice-acceptor site**, so the
downstream exon is skipped. The example data are three exons from *DMD*
(exons 44, 45 and 53).

## How it works

The input is a genomic sequence in which **exonic bases are uppercase and
intronic bases are lowercase** (the usual Ensembl/UCSC convention). The target
is the `ag` immediately 5' of the exon, so the editable base is found relative
to the first uppercase base.

For each of the three strands (top, complement, reverse complement) the tool:

1. finds every PAM recognised by SpCas9 (NGG, NAG), xCas9 (NGN) or the VQR
   variant (NGAN, NGNG);
2. keeps the PAMs sitting 12–17 nt from the target, which places the target at
   guide positions 3–8;
3. builds the adjacent 20 nt guide;
4. reports each base editor (ABE7.10, ABE7.9, BE4) whose editing window covers
   the target, plus any other editable base in that window — a possible
   off-target edit.

## Installing

```bash
conda create -n grna-design -c conda-forge python=3.10 biopython pytest
conda activate grna-design
```

Or with pip: `pip install -r requirements.txt`.

## Usage

```bash
# A FASTA file (or plain text; '>' headers are ignored)
python grna_design.py -f data/exon45.fasta

# A sequence given directly
python grna_design.py -s "tttattctagTTGAAAGAATTCAGAATCAGTGG"

# Tab-separated or JSON output, for piping onward
python grna_design.py -f data/exon45.fasta --format tsv
python grna_design.py -f data/exon45.fasta --format json
```

Example output:

```
Strand:        top
PAM:           TGGC at 61
Nuclease(s):   spcas9 (TGG), xcas9 (GGC/TGG)
Guide:         atcttacagGAACTCCAGGA
Base editor:   ABE7.9
Target:        a8
Window:        ttacag
Off-target(s): a6
```

Importing the module gives the same results as data rather than text:

```python
from grna_design import design, read_sequence

for candidate in design(read_sequence("data/exon45.fasta")):
    print(candidate.strand, candidate.editor, candidate.target, candidate.guide)
```

## Tests

```bash
python -m pytest tests
```

`tests/golden/candidates.json` pins the full output for the three example
exons, so any change in behaviour shows up as a diff. Update it deliberately,
in its own commit, when a change is meant to alter results.

## Relation to the original script

The original script is kept at [legacy/v1_original.py](legacy/v1_original.py)
for comparison only; nothing imports it. It could not be run as it stood (its
pipeline was commented out), and it had three behavioural problems that the
rewrite fixes:

- **It reported only one candidate per strand.** Each analysis function
  returned from inside its loop, discarding the rest. The example exons yield
  14 candidates in total, where the original reported at most 3.
- **Its results were not reproducible.** It picked the candidate to report by
  looking up a value in a dictionary keyed by PAM, iterating a `set` of
  strings. Because Python randomises string hashing per process, the guide,
  PAM and nuclease it reported changed between runs. Running it under eight
  hash seeds produced 11 different candidates.
- **It under-reported nucleases.** Only the nucleases matching one PAM motif
  were listed, so a site recognised by both SpCas9 and xCas9 was credited to
  one of them.

`tests/golden/legacy_candidates.json` holds all 11 candidates the original
produced across those runs; `test_reproduces_legacy_candidates` checks each one
is still found, with the same off-target list.

Two behaviours were deliberately **kept** because they affect results and the
choice is a biological one, not a code one:

- PAM matching is case-sensitive, so the top and reverse-complement strands
  match exonic (uppercase) PAMs while the complement strand matches intronic
  (lowercase) ones. This is consistent with where a PAM falls relative to the
  splice acceptor, but it is implicit rather than stated.
- The span scanned for off-target bases reaches one base beyond the reported
  editing window for ABE7.10, ABE7.9 and BE4 on the top strand. See the
  `BaseEditor` docstring in [grna_design.py](grna_design.py).

Both are worth checking against the source publications before the tool is
used in anger.
