"""Capture the original v1 pipeline's output as a regression baseline.

Imports the legacy `v1` module (no .py extension) and replays the pipeline
that is commented out at the bottom of that file, for each embedded exon.
"""
import importlib.machinery
import importlib.util
import io
import json
import re
import sys
from contextlib import redirect_stdout
from pathlib import Path

REPO = str(Path(__file__).resolve().parent.parent)

spec = importlib.util.spec_from_loader(
    "v1", importlib.machinery.SourceFileLoader("v1", f"{REPO}/legacy/v1_original.py")
)
v1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v1)

raw = open(f"{REPO}/legacy/v1_original.py").read()
exons = dict(re.findall(r"'''(Exon\d+)\s*\n(.*?)'''", raw, re.S))

out = {}
for name, block in sorted(exons.items()):
    sequence = v1.input_edit(block)
    seqs = v1.all_seqs(sequence)
    sites = v1.all_sites(seqs)
    pam_seqs = v1.possible_pams()
    complement_pams = v1.complement_pams_func(pam_seqs)
    pam_indxs = v1.pam_indx(seqs, pam_seqs, complement_pams)
    pam_locations = v1.filter_location(sites, pam_indxs)

    grna_comp = v1.grna_design_compl(pam_locations, pam_seqs, complement_pams, seqs)
    grna_top = v1.grna_design(pam_locations, pam_seqs, complement_pams, seqs)[0]
    grna_rev = v1.grna_design(pam_locations, pam_seqs, complement_pams, seqs)[1]

    pam_top = v1.get_nuclease(grna_top, pam_seqs)
    pam_rev = v1.get_nuclease(grna_rev, pam_seqs)
    pam_comp = v1.get_nuclease(grna_comp, complement_pams)

    pam_grna_top = v1.fix_pam(grna_top, pam_top)
    pam_grna_comp = v1.fix_pam(grna_comp, pam_comp)
    pam_grna_rev = v1.fix_pam(grna_rev, pam_rev)

    site_top = v1.get_site(pam_grna_top)
    site_comp = v1.get_site(pam_grna_comp)
    site_rev = v1.get_site(pam_grna_rev)

    buf = io.StringIO()
    with redirect_stdout(buf):
        res_top = v1.off_target(site_top, pam_grna_top, pam_top, seqs["top"])
        res_comp = v1.off_target_comp(site_comp, pam_grna_comp, pam_comp, seqs["complement"])
        res_rev = v1.off_target(site_rev, pam_grna_rev, pam_rev, seqs["reverse_complement"])

    out[name] = {
        "sequence": sequence,
        "sites": {k: v for k, v in sites.items()},
        "pam_locations": {k: sorted(v) for k, v in pam_locations.items()},
        "nucleases": {"top": pam_top, "comp": pam_comp, "rev": pam_rev},
        "pam_grna": {"top": pam_grna_top, "comp": pam_grna_comp, "rev": pam_grna_rev},
        "edit_sites": {"top": site_top, "comp": site_comp, "rev": site_rev},
        "results": {"top": res_top, "comp": res_comp, "rev": res_rev},
        "stdout": buf.getvalue(),
    }

json.dump(out, sys.stdout, indent=2, sort_keys=True)
