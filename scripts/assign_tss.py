#!/usr/bin/env python3
"""Assign TSS positions according to STARR-seq template-switch data.

The code was developed with reference to the papers by Teemu Kivioja et al. (2021) and Sahu et al. (2022)

Usage:
    python3 assign_tss.py --conf example.conf
"""

from __future__ import annotations

import argparse
import gzip
import os
import pickle
import random
import re
import shlex
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


DEBUG = False


# ============================================================================
#                       Template-switch position model
#       (port of TemplateSwitchPosition.pm)
# ============================================================================

# Classification of template switch dinucleotide patterns.
#   C: certain TSS position
#   U: uncertain
#   E: error -- RNA and input DNA inconsistent under our assumptions
TS_CLASS: Dict[str, str] = {
    "HH-HH": "C", "HH-HG": "E", "HH-GH": "E", "HH-GG": "E",
    "HG-HH": "E", "HG-HG": "U", "HG-GH": "E", "HG-GG": "E",
    "GH-HH": "C", "GH-HG": "E", "GH-GH": "U", "GH-GG": "E",
    "GG-HH": "C", "GG-HG": "U", "GG-GH": "U", "GG-GG": "U",
}

# Collapse to four classes according to the first base only.
TS_FIRST_CLASS: Dict[str, str] = {
    "HH-HH": "H-H", "HH-HG": "H-H", "HH-GH": "H-G", "HH-GG": "H-G",
    "HG-HH": "H-H", "HG-HG": "H-H", "HG-GH": "H-G", "HG-GG": "H-G",
    "GH-HH": "G-H", "GH-HG": "G-H", "GH-GH": "G-G", "GH-GG": "G-G",
    "GG-HH": "G-H", "GG-HG": "G-H", "GG-GH": "G-G", "GG-GG": "G-G",
}

TS_FIRST_CLASS_TYPE: Dict[str, str] = {
    "H-H": "C",
    "H-G": "E",
    "G-H": "C",
    "G-G": "U",
}


def ts_pattern(b1: str, b2: str) -> str:
    """Return the dinucleotide pattern H/G for two bases."""
    return ("G" if b1 == "G" else "H") + ("G" if b2 == "G" else "H")


def ts_pattern_pair(rna_seq: str, input_seq: str) -> str:
    """Return the joined RNA-input dinucleotide pattern e.g. 'GH-HH'."""
    rp = ts_pattern(rna_seq[0], rna_seq[1])
    ip = ts_pattern(input_seq[0], input_seq[1])
    return f"{rp}-{ip}"


class TemplateSwitchPosition:
    """Counts template-switch patterns and turns them into TSS-shift calls."""

    def __init__(self) -> None:
        self.ts_counts: Dict[str, int] = {pat: 0 for pat in TS_CLASS}
        self.ts_first_counts: Dict[str, int] = {pat: 0 for pat in TS_FIRST_CLASS_TYPE}

        # Filled by estimate_*_probabilities().
        self.pcond: Dict[str, Dict[int, float]] = {}
        self.pcum: Dict[str, Dict[int, float]] = {}
        self.pcond_first: Dict[str, Dict[int, float]] = {}
        self.pcum_first: Dict[str, Dict[int, float]] = {}
        self.est_p0: float = 0.0
        self.est_p1: float = 0.0

    # ------------------------------------------------------------------ counts

    def add_pair(self, rna_seq: str, input_seq: str) -> str:
        """Record one RNA/input pair and return its pattern."""
        pair = ts_pattern_pair(rna_seq, input_seq)
        if pair not in self.ts_counts:
            raise ValueError(f"Unknown TS pair {pair}")
        self.ts_counts[pair] += 1
        return pair

    def print_counts(self, stream=sys.stderr) -> None:
        for pat in sorted(TS_CLASS):
            stream.write(f"{pat}\t{TS_CLASS[pat]}\t{self.ts_counts[pat]}\n")

    def _sum_ts_pair_counts(self) -> None:
        """Collapse 2-base counts into 1-base counts."""
        for pat in self.ts_first_counts:
            self.ts_first_counts[pat] = 0
        for pair_pat, count in self.ts_counts.items():
            self.ts_first_counts[TS_FIRST_CLASS[pair_pat]] += count

    # ------------------------------------------------------ probability estim.

    def estimate_2base_probabilities(self) -> Dict[str, Dict[int, float]]:
        """Estimate cumulative TSS-shift probs (0,1,2) using 2-base context."""
        c = self.ts_counts
        nc = sum(c[p] for p, k in TS_CLASS.items() if k == "C")
        nu = sum(c[p] for p, k in TS_CLASS.items() if k == "U")
        ne = sum(c[p] for p, k in TS_CLASS.items() if k == "E")
        if nc + nu + ne != sum(c.values()):
            raise RuntimeError(
                f"Class totals do not add up: {nc}+{nu}+{ne} != {sum(c.values())}"
            )

        p0 = c["HH-HH"] / nc
        p1 = c["GH-HH"] / nc
        p2 = c["GG-HH"] / nc
        print(f"Observed probabilities: p(0) = {p0}, p(1) = {p1}, p(2) = {p2}",
              file=sys.stderr)
        if abs((p0 + p1 + p2) - 1) > 1e-3:
            raise RuntimeError("Estimated TSS probabilities do not sum to one")

        p_uts: Dict[str, float] = {}
        p_cts: Dict[str, float] = {}
        for pat in sorted(TS_CLASS):
            if TS_CLASS[pat] == "U":
                p_uts[pat] = c[pat] / nu
                print(f"U: {pat} prob {p_uts[pat]}", file=sys.stderr)
            elif TS_CLASS[pat] == "C":
                p_cts[pat] = c[pat] / nc
                print(f"C: {pat} prob {p_cts[pat]}", file=sys.stderr)
        print(f"sum of uncertain classes: {sum(p_uts.values())}", file=sys.stderr)
        print(f"sum of certain classes:   {sum(p_cts.values())}", file=sys.stderr)

        a1 = p1 / p0
        a2 = p2 / p1
        x = 1 / (a1 + 1)
        y = 1 / (a2 + 1)
        z = (p1 - p_uts["GH-GH"] * (1 - x) - p_uts["GG-HG"] * y) / p_uts["GG-GG"]
        print(f"x = {x}, y = {y}, z = {z}", file=sys.stderr)

        p_cond: Dict[str, Dict[int, float]] = {
            "GH-GH": {0: x, 1: 1 - x},
            "GG-HG": {1: y, 2: 1 - y},
            "GG-GG": {
                0: (p0 - p_uts["HG-HG"] - p_uts["GH-GH"] * x) / p_uts["GG-GG"],
                1: z,
                2: (p2 - p_uts["GG-GH"] - p_uts["GG-HG"] * (1 - y)) / p_uts["GG-GG"],
            },
        }
        for pat in sorted(p_cond):
            for pos in sorted(p_cond[pat]):
                print(f"cond prob for {pat} at {pos} is {p_cond[pat][pos]}",
                      file=sys.stderr)
        self.pcond = p_cond

        # Cumulative: probability that the TSS is here OR earlier.
        p_cum: Dict[str, Dict[int, float]] = {
            "HH-HH": {0: 1.0},
            "HG-HG": {0: 1.0},
            "GH-HH": {0: 0.0, 1: 1.0},
            "GG-HH": {0: 0.0, 1: 0.0, 2: 1.0},
            "GG-GH": {0: 0.0, 1: 0.0, 2: 1.0},
            "GH-GH": {0: p_cond["GH-GH"][0], 1: 1.0},
            "GG-HG": {0: 0.0, 1: p_cond["GG-HG"][1], 2: 1.0},
            "GG-GG": {
                0: p_cond["GG-GG"][0],
                1: p_cond["GG-GG"][0] + p_cond["GG-GG"][1],
                2: 1.0,
            },
        }
        self.pcum = p_cum
        return p_cum

    def estimate_1base_probabilities(self) -> Dict[str, Dict[int, float]]:
        """Estimate cumulative TSS-shift probs (0,1) using 1-base context."""
        self._sum_ts_pair_counts()
        c = self.ts_first_counts

        nc = sum(c[p] for p, k in TS_FIRST_CLASS_TYPE.items() if k == "C")
        nu = sum(c[p] for p, k in TS_FIRST_CLASS_TYPE.items() if k == "U")
        ne = sum(c[p] for p, k in TS_FIRST_CLASS_TYPE.items() if k == "E")
        if nc + nu + ne != sum(c.values()):
            raise RuntimeError("Class totals do not add up")

        p0 = c["H-H"] / nc
        p1 = c["G-H"] / nc
        print(f"Observed probabilities: p(0) = {p0}, p(1) = {p1}", file=sys.stderr)
        self.est_p0 = p0
        self.est_p1 = p1
        if abs((p0 + p1) - 1) > 1e-3:
            raise RuntimeError("Estimated TSS probabilities do not sum to one")

        p_uts: Dict[str, float] = {}
        p_cts: Dict[str, float] = {}
        for pat, kind in TS_FIRST_CLASS_TYPE.items():
            if kind == "U":
                p_uts[pat] = c[pat] / nu
                print(f"U: {pat} prob {p_uts[pat]}", file=sys.stderr)
            elif kind == "C":
                p_cts[pat] = c[pat] / nc
                print(f"C: {pat} prob {p_cts[pat]}", file=sys.stderr)
        print(f"sum of uncertain classes: {sum(p_uts.values())}", file=sys.stderr)
        print(f"sum of certain classes:   {sum(p_cts.values())}", file=sys.stderr)

        a = p0 / p1
        x = 1 / (a + 1)

        p_cond = {"G-G": {0: 1 - x, 1: x}}
        for pat in sorted(p_cond):
            for pos in sorted(p_cond[pat]):
                print(f"cond prob for {pat} at {pos} is {p_cond[pat][pos]}",
                      file=sys.stderr)
        self.pcond_first = p_cond

        p_cum = {
            "H-H": {0: 1.0},
            "G-H": {0: 0.0, 1: 1.0},
            "G-G": {0: p_cond["G-G"][0], 1: 1.0},
        }
        self.pcum_first = p_cum
        return p_cum

    # ------------------------------------------------------------ assignment

    def assign_tss_shift_max_two(self, rna_seq: str, input_seq: str) -> int:
        pat = ts_pattern_pair(rna_seq, input_seq)
        kind = TS_CLASS[pat]
        if kind in ("C", "U"):
            r = random.random()
            cum = self.pcum[pat]
            if r < cum.get(0, 0):
                return 0
            if r < cum.get(1, 0):
                return 1
            if cum.get(2, 0) == 1:
                return 2
            raise RuntimeError(f"Could not assign TSS to pattern {pat}")
        if kind == "E":
            return -1
        raise RuntimeError(f"Could not assign TSS to pattern {pat}")

    def assign_tss_shift_max_one(self, rna_seq: str, input_seq: str) -> int:
        pair_pat = ts_pattern_pair(rna_seq, input_seq)
        kind = TS_CLASS[pair_pat]
        if kind in ("C", "U"):
            pat = TS_FIRST_CLASS[pair_pat]
            r = random.random()
            cum = self.pcum_first[pat]
            if r < cum.get(0, 0):
                return 0
            if r < cum.get(1, 0):
                return 1
            raise RuntimeError(f"Could not assign TSS to pattern {pat}")
        if kind == "E":
            return -1
        raise RuntimeError(f"Could not assign TSS to pattern {pair_pat}")

    def assign_tss_without_input_shift_max_one(self, rna_seq: str, pG: float) -> int:
        """Assign TSS when no matching input sequence is available."""
        if rna_seq[0] == "G":
            # P(0|G) = (P(0) - P(H)) / P(G)
            pH = 1 - pG
            p0G = (self.est_p0 - pH) / pG
            return 0 if random.random() < p0G else 1
        return 0


# ============================================================================
#                       Configuration loader
# ============================================================================

@dataclass
class Config:
    R1_rna_fastq_file: str = ""
    R1_input_fastq_file: str = ""
    R1_input_tab_sorted_file: str = ""
    output_dir: str = ""
    UMI_length: int = 10
    # promoter_length and utr_length are PAIRED lists. e.g.
    #     promoter_length = "60,100"
    #     utr_length      = "60,20"
    # produces two output sets, (prom=60,utr=60) and (prom=100,utr=20).
    promoter_length: List[int] = field(default_factory=lambda: [40])
    utr_length: List[int] = field(default_factory=lambda: [40])
    spliced_intron_seq: str = (
        "AAGCTTCTGCCTTCTCCCTCCTGTGAGTTTGGTTGGTGTACAGTAGCTTCCACC"
    )
    min_intron_overlap: int = 40
    run_revcoms_and_stats: int = 0


_INT_FIELDS = {
    "UMI_length",
    "min_intron_overlap",
    "run_revcoms_and_stats",
}

_INT_LIST_FIELDS = {
    "promoter_length",
    "utr_length",
}

# Lines look like: key = "value"  /  key = value  /  key="value"
_CONF_LINE = re.compile(r'^\s*([A-Za-z_]\w*)\s*=\s*"?([^"#]*?)"?\s*(?:#.*)?$')


def _parse_int_list(value: str) -> List[int]:
    """Parse '60', '60,100', '60 100', '60, 100' into a list of ints."""
    parts = [p for p in re.split(r"[,\s]+", value.strip()) if p]
    return [int(p) for p in parts]


def read_conf(path: str) -> Config:
    """Read an AppConfig-style file like example.conf into a Config object."""
    conf = Config()
    valid = {f for f in conf.__dataclass_fields__}
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            m = _CONF_LINE.match(line)
            if not m:
                continue
            key, value = m.group(1), m.group(2).strip()
            if key not in valid:
                # Unknown keys (anno_file, seq_out_file_template, ...) — ignore.
                continue
            if key in _INT_FIELDS:
                setattr(conf, key, int(value))
            elif key in _INT_LIST_FIELDS:
                setattr(conf, key, _parse_int_list(value))
            else:
                setattr(conf, key, value)
    return conf


# ============================================================================
#                       IO helpers
# ============================================================================

def _open_text(path: str):
    """Open plain text or gzip transparently, in text mode."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def run(cmd: str, *, shell_executable: str = "/bin/bash") -> None:
    """Run a shell command, echo to stderr, raise on non-zero exit."""
    print(f"\tcommand: {cmd}", file=sys.stderr)
    res = subprocess.run(cmd, shell=True, executable=shell_executable)
    if res.returncode != 0:
        raise RuntimeError(f"Command failed (exit {res.returncode}): {cmd}")


def _exists_non_empty(*paths: str) -> bool:
    """True only if every path exists and has nonzero size."""
    for p in paths:
        try:
            if not os.path.exists(p) or os.path.getsize(p) == 0:
                return False
        except OSError:
            return False
    return True


# ============================================================================
#                       Core analysis steps
# ============================================================================

@dataclass
class RNA:
    seq: str = ""
    umilabel: str = ""
    numg: int = 0
    rand_suffix: str = ""
    input_seq: Optional[str] = None
    input_id: Optional[str] = None
    match_pos: int = -1
    tss_shift: int = -2
    no_input_tss_shift: Optional[int] = None


def parse_template_switched_from_trim_info(
    info_file: str, min_rna_len: int, umi_len: int, suffix_file: str
) -> Dict[str, RNA]:
    """Parse cutadapt's --info-file output, extract template-switched RNAs."""
    rnas: Dict[str, RNA] = {}
    num_match = 0
    num_nogs = num3g = num4g = num5g = 0

    with open(info_file) as ifh, open(suffix_file, "w") as ofh:
        for raw in ifh:
            parts = raw.rstrip("\n").split("\t")
            # id, ne, pos, end_pos, left, match, right, adapter_name, qleft, qmatch, qright
            if len(parts) < 7:
                continue
            read_id = parts[0]
            try:
                ne = int(parts[1])
            except ValueError:
                continue
            left_flank_seq = parts[4]

            if ne < 0 or len(left_flank_seq) <= min_rna_len:
                continue

            umi_label = left_flank_seq[:umi_len]
            no_umi_seq = left_flank_seq[umi_len:]

            # strip anything after whitespace in the id
            read_id = re.sub(r"\s[\w:]+", "", read_id)

            numg = 0
            if no_umi_seq.startswith("GGGGG"):
                num5g += 1
                numg = 5
            elif no_umi_seq.startswith("GGGG"):
                num4g += 1
                numg = 4
            elif no_umi_seq.startswith("GGG"):
                num3g += 1
                numg = 3
            else:
                num_nogs += 1

            if numg == 0:
                continue

            rna_seq = no_umi_seq[3:]
            if len(rna_seq) <= min_rna_len:
                continue

            rand_suffix = left_flank_seq[-min_rna_len:]
            rnas[read_id] = RNA(
                seq=rna_seq,
                umilabel=umi_label,
                numg=numg,
                rand_suffix=rand_suffix,
            )
            ofh.write(f"{read_id}\t{rand_suffix}\n")
            num_match += 1

    print(
        f"\tparsing trimmed: starting with 5 Gs {num5g}, 4 Gs {num4g}; "
        f"3 Gs {num3g}; no 3 G {num_nogs}",
        file=sys.stderr,
    )
    print(
        f"\tparsing trimmed: parsed {num_match} sequences flanking the trimmed "
        f"sequence and at least {min_rna_len} long",
        file=sys.stderr,
    )
    return rnas


def generate_r1_match_index(
    rna_input_matches_file: str, rna_in_input_r1_file: str
) -> Dict[str, List[SeqRecord]]:
    """For every RNA id, list the input R1 SeqRecord(s) it matches."""
    input_to_rna_ids: Dict[str, List[str]] = defaultdict(list)
    with _open_text(rna_input_matches_file) as fh:
        for raw in fh:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            _r1_match_seq, input_id, rna_id = parts[0], parts[1], parts[2]
            input_to_rna_ids[input_id].append(rna_id)

    r1_match_index: Dict[str, List[SeqRecord]] = defaultdict(list)
    num_indexed = 0
    with _open_text(rna_in_input_r1_file) as fh:
        for rec in SeqIO.parse(fh, "fastq"):
            input_id = rec.id
            rna_ids = input_to_rna_ids.get(input_id)
            if rna_ids is None:
                raise RuntimeError(
                    f"Could not find matching RNA ids for input id {input_id}"
                )
            for rna_id in rna_ids:
                r1_match_index[rna_id].append(rec)
                num_indexed += 1

    print(f"Indexed {num_indexed} RNA input sequence pairs", file=sys.stderr)
    return r1_match_index


def match_to_index(
    pattern: str, rna_id: str, match_index: Dict[str, List[SeqRecord]]
) -> Tuple[int, Optional[SeqRecord]]:
    """Find the pattern in candidate input reads, accept only unique hits."""
    candidates = match_index.get(rna_id)
    if not candidates:
        return -1, None
    found = 0
    match_pos = -1
    match_rec: Optional[SeqRecord] = None
    for rec in candidates:
        pos = str(rec.seq).find(pattern)
        if pos != -1:
            found += 1
            match_pos = pos
            match_rec = rec
    if found == 1:
        return match_pos, match_rec
    return -1, None


def estimate_tss_probs(
    rnas: Dict[str, RNA],
    rna_input_matches_file: str,
    input_file: str,
    trim_min_len: int,
) -> Dict[str, RNA]:
    """Estimate TSS shift probabilities and tag each RNA with its shift."""
    tss_pos = TemplateSwitchPosition()
    rna2inputseq = generate_r1_match_index(rna_input_matches_file, input_file)

    num_matching = 0
    num_not_matching = 0
    num_not_matching_G = 0

    for rna_id, rna in rnas.items():
        # match a few bases downstream in case of extra Gs
        match_seq = rna.seq[2:2 + trim_min_len]
        match_pos, input_rec = match_to_index(match_seq, rna_id, rna2inputseq)

        if match_pos >= 0 and input_rec is not None:
            match_pos -= 2
            num_matching += 1

            input_seq = str(input_rec.seq)
            rna.input_seq = input_seq
            rna.input_id = input_rec.id
            rna.match_pos = match_pos

            input_match_seq = input_seq[match_pos:match_pos + trim_min_len]
            pair = tss_pos.add_pair(rna.seq, input_match_seq)
            if DEBUG:
                print(f"{pair}\nrna {rna.seq}\nin  {input_match_seq}", file=sys.stderr)
        else:
            num_not_matching += 1
            if rna.seq[:1] == "G":
                num_not_matching_G += 1

    num_tot = len(rnas)
    print(f"Number of matching {num_matching}, total {num_tot}", file=sys.stderr)
    pG_utr = num_not_matching_G / num_not_matching if num_not_matching else 0.0
    print(
        f"Number of non_matching {num_not_matching}, of which starts with G "
        f"{num_not_matching_G}, {pG_utr:.3f}",
        file=sys.stderr,
    )
    tss_pos.print_counts()
    tss_pos.estimate_1base_probabilities()

    # Now assign every RNA its TSS shift.
    at = ae = a0 = a1 = 0
    at_utr = a0_utr = a1_utr = g_utr = 0
    for rna in rnas.values():
        if rna.input_seq is not None:
            input_match_seq = rna.input_seq[rna.match_pos:rna.match_pos + trim_min_len]
            shift = tss_pos.assign_tss_shift_max_one(rna.seq, input_match_seq)
            if shift == 0:
                rna.tss_shift = 0
                a0 += 1
                at += 1
            elif shift == 1:
                rna.tss_shift = 1
                a1 += 1
                at += 1
            elif shift == -1:
                rna.tss_shift = -1
                ae += 1
            else:
                raise RuntimeError(f"Could not handle TSS shift {shift}")
        else:
            rna.tss_shift = -2  # no input match
            shift = tss_pos.assign_tss_without_input_shift_max_one(rna.seq, pG_utr)
            if rna.seq[:1] == "G":
                g_utr += 1
            if shift == 0:
                a0_utr += 1
                at_utr += 1
                rna.no_input_tss_shift = 0
            elif shift == 1:
                a1_utr += 1
                at_utr += 1
                rna.no_input_tss_shift = 1
            else:
                raise RuntimeError(f"Could not handle no input TSS shift {shift}")

    if at:
        print(
            f"Assigned probabilities:  p(0) = {a0/at:.3f}, p(1) = {a1/at:.3f}, "
            f"assigned {at}, rejected {ae}",
            file=sys.stderr,
        )
    if at_utr:
        print(
            f"Assigned probabilities for those lacking input:  "
            f"p(0) = {a0_utr/at_utr:.3f}, p(1) = {a1_utr/at_utr:.3f}, "
            f"assigned {at_utr}, p(G) = {g_utr/at_utr:.3f}",
            file=sys.stderr,
        )
    return rnas


def get_tss_seqs(
    rnas: Dict[str, RNA],
    prom_out_file: str,
    tss_out_file: str,
    prom_tss_out_file: str,
    umi_counts_out_file: str,
    anno_file: str,
    trim_min_len: int,
    prom_len: int,
    tss_len: int,
) -> None:
    """Write the promoter, UTR, both-flanks fasta files and the annotation tsv."""
    num_matching = num_not_matching = num_dup = num_written = 0
    umi_counts: Dict[str, Dict[int, Dict[str, int]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(int))
    )
    seen_input: Dict[str, bool] = {}

    with open(prom_out_file, "w") as prom_fh, open(tss_out_file, "w") as tss_fh, \
            open(prom_tss_out_file, "w") as both_fh, open(anno_file, "w") as anno_fh:

        anno_fh.write(
            "\t".join(
                ["input.id", "rna.id", "tss.shift", "tss.pos", "numg",
                 "tss.seq", "tss.flank.seq", "rna.seq", "inputstr", "umi"]
            ) + "\n"
        )

        for rna_id, rna in rnas.items():
            tss_shift = rna.tss_shift
            if tss_shift < 0:
                num_not_matching += 1
                continue

            num_matching += 1
            assert rna.input_seq is not None and rna.input_id is not None
            input_str = rna.input_seq
            input_id = rna.input_id
            tss_pos = rna.match_pos + tss_shift
            start_pos = tss_pos - prom_len  # 0-based
            if start_pos < 0:
                continue

            if input_id not in seen_input:
                seen_input[input_id] = True
                prom_str = input_str[start_pos:start_pos + prom_len]
                prom_id = f"{input_id}:pos_{tss_pos}"
                SeqIO.write(
                    [SeqRecord(Seq(prom_str), id=prom_id, description="")],
                    prom_fh, "fasta",
                )
                num_written += 1

                tss_str = input_str[tss_pos:tss_pos + tss_len]
                if tss_len > 0 and len(tss_str) >= tss_len:
                    tss_id = f"{input_id}:pos_{tss_pos}"
                    SeqIO.write(
                        [SeqRecord(Seq(tss_str), id=tss_id, description="")],
                        tss_fh, "fasta",
                    )
                    prom_tss_str = prom_str + tss_str
                    SeqIO.write(
                        [SeqRecord(Seq(prom_tss_str), id=tss_id, description="")],
                        both_fh, "fasta",
                    )
                    anno_fh.write(
                        "\t".join(str(x) for x in [
                            input_id, rna_id, tss_shift, tss_pos, rna.numg,
                            tss_str, prom_tss_str, rna.seq, input_str, rna.umilabel,
                        ]) + "\n"
                    )
                else:
                    anno_fh.write(
                        "\t".join(str(x) for x in [
                            input_id, rna_id, tss_shift, tss_pos, rna.numg,
                            ".", ".", rna.seq, input_str, rna.umilabel,
                        ]) + "\n"
                    )
            else:
                num_dup += 1

            umi_counts[input_id][tss_pos][rna.umilabel] += 1

    with open(umi_counts_out_file, "w") as ufh:
        for input_id in sorted(umi_counts):
            for pos in sorted(umi_counts[input_id]):
                for umi_label in sorted(umi_counts[input_id][pos]):
                    count = umi_counts[input_id][pos][umi_label]
                    ufh.write(f"{input_id}\t{pos}\t{umi_label}\t{count}\n")

    print(
        f"\tnumber of matched sequences {num_matching}, number of sequences "
        f"not matched (or inconsistent) {num_not_matching}",
        file=sys.stderr,
    )
    print(f"\twrote in total {num_written} sequences", file=sys.stderr)
    print(f"\tNumber of RNA:s matching to the same input: {num_dup}", file=sys.stderr)


# ============================================================================
#                       Main
# ============================================================================

def base_name_no_fastq_gz(path: str) -> Tuple[str, str, str]:
    """Mimic Perl's fileparse($p, '.fastq.gz')."""
    p = Path(path)
    dirs = str(p.parent) + os.sep if str(p.parent) not in ("", ".") else ""
    name = p.name
    if name.endswith(".fastq.gz"):
        return name[: -len(".fastq.gz")], dirs, ".fastq.gz"
    return p.stem, dirs, p.suffix


def _prom_utr_arg(value: str) -> Tuple[int, int]:
    """Parse a --prom-utr P,U argument into (P, U)."""
    parts = re.split(r"[,\s:]+", value.strip())
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(
            f"expected PROM,UTR (e.g. 60,60), got: {value!r}"
        )
    return int(parts[0]), int(parts[1])


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Assign TSS positions from STARR-seq template-switch data."
    )
    parser.add_argument("--conf", required=True, help="Configuration file")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument(
        "--prom-utr",
        action="append",
        type=_prom_utr_arg,
        metavar="PROM,UTR",
        help=("Promoter,UTR flank length pair. Repeatable, e.g. "
              "`--prom-utr 60,60 --prom-utr 100,20`. Overrides promoter_length "
              "/ utr_length from the config file when given."),
    )
    parser.add_argument(
        "--skip-existing",
        action="store_true",
        help=("Reuse intermediate files from a previous run when they exist. "
              "Skips cutadapt (step 1), the join/seqkit pipeline (step 2), "
              "and the TSS-shift assignment (step 3a) if their cached outputs "
              "are present and non-empty. Only the per-combo fasta/anno "
              "writeout is repeated."),
    )
    args = parser.parse_args(argv)

    conf = read_conf(args.conf)
    random.seed(args.seed)

    # Build the list of (promoter_length, utr_length) combos to emit.
    # Track whether they came from CLI: CLI combos always get a name suffix
    # so they don't collide with the default config-driven output filenames.
    if args.prom_utr:
        combos: List[Tuple[int, int]] = list(args.prom_utr)
        combos_from_cli = True
    else:
        if len(conf.promoter_length) != len(conf.utr_length):
            raise SystemExit(
                f"promoter_length ({conf.promoter_length}) and utr_length "
                f"({conf.utr_length}) must have the same number of entries"
            )
        combos = list(zip(conf.promoter_length, conf.utr_length))
        combos_from_cli = False
    if not combos:
        raise SystemExit("No promoter/utr length combinations to run")
    print(f"Running {len(combos)} promoter/utr combo(s): {combos}", file=sys.stderr)

    min_len = 20
    umi_len = conf.UMI_length
    out_dir = conf.output_dir
    if out_dir and not out_dir.endswith(os.sep):
        out_dir = out_dir + os.sep
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    base, dirs, suffix = base_name_no_fastq_gz(conf.R1_rna_fastq_file)
    print(f"Script start file: {base}, {dirs}, {suffix}", file=sys.stderr)

    r1_trimmed_fastq = f"{out_dir}{base}.intron_trimmed.fastq.gz"
    r1_trim_info = f"{out_dir}{base}.intron_info.txt"
    rna_rand_suffix_file = f"{out_dir}{base}.rand_suffix.tsv"
    rnas_pickle_file = f"{out_dir}{base}.rnas.pkl"

    def out_base(prom: int, utr: int) -> str:
        """Insert a `.prom{p}utr{u}` tag into the basename.

        Tagged when either:
          * the user passed --prom-utr on the CLI (so they don't clobber
            files from an earlier default run), OR
          * there is more than one combo this invocation.
        Otherwise (single combo from the config) we keep the original
        un-suffixed names for backward compatibility.
        """
        if combos_from_cli or len(combos) > 1:
            return f"{out_dir}{base}.prom{prom}utr{utr}"
        return f"{out_dir}{base}"

    rna_input_matches_file = f"{out_dir}{base}.rna_input_matches.tsv.gz"
    rna_input_matches_id_file = f"{out_dir}{base}.rna_input_matches_ids.txt"
    rna_in_input_r1_file = f"{out_dir}{base}.rna_in_input_R1.fastq"

    # ------------------------------------------------------------------
    # Cache fast-path: if we have a pickled `rnas` from a previous run,
    # skip ALL of steps 1, 2, 3a and just reload it.
    # ------------------------------------------------------------------
    if args.skip_existing and _exists_non_empty(rnas_pickle_file):
        print(
            f"Cache: loading TSS-assigned RNAs from {rnas_pickle_file} "
            f"(skipping steps 1-3a)",
            file=sys.stderr,
        )
        with open(rnas_pickle_file, "rb") as fh:
            rnas = pickle.load(fh)
    else:
        # Step 1 — trim spliced intron with cutadapt
        if args.skip_existing and _exists_non_empty(r1_trimmed_fastq, r1_trim_info):
            print(
                "Step 1: skipping cutadapt (trimmed fastq + info file already exist)",
                file=sys.stderr,
            )
        else:
            print("Step 1: Recognizing and trimming the spliced intron from R1...")
            cmd = (
                f"cutadapt --discard-untrimmed -e 0.05 "
                f"--overlap={conf.min_intron_overlap} "
                f"--minimum-length={min_len} -a {conf.spliced_intron_seq} "
                f"-o {shlex.quote(r1_trimmed_fastq)} "
                f"{shlex.quote(conf.R1_rna_fastq_file)} "
                f"--info-file {shlex.quote(r1_trim_info)}"
            )
            run(cmd)

        # Step 2 — parse trim info, then match RNA suffix to input
        rnas = parse_template_switched_from_trim_info(
            r1_trim_info, min_len, umi_len, rna_rand_suffix_file
        )

        if args.skip_existing and _exists_non_empty(
            rna_input_matches_file, rna_in_input_r1_file
        ):
            print(
                "Step 2: skipping join/seqkit pipeline "
                "(matches + input-R1 files already exist)",
                file=sys.stderr,
            )
        else:
            cmd2 = (
                f"join -1 2 -2 2 "
                f"<(zcat {shlex.quote(conf.R1_input_tab_sorted_file)}) "
                f"<(sort -k2 -b {shlex.quote(rna_rand_suffix_file)}) "
                f"| cut -f 1,2,3 -d ' ' | tr ' ' '\\t' "
                f"| gzip -c > {shlex.quote(rna_input_matches_file)}"
            )
            run(cmd2)

            cmd3 = (
                f"zcat {shlex.quote(rna_input_matches_file)} | cut -f 2 | sort -u "
                f"> {shlex.quote(rna_input_matches_id_file)}"
            )
            run(cmd3)

            cmd4 = (
                f"seqkit grep -f {shlex.quote(rna_input_matches_id_file)} "
                f"-o {shlex.quote(rna_in_input_r1_file)} "
                f"{shlex.quote(conf.R1_input_fastq_file)}"
            )
            run(cmd4)

        # Step 3a — assign TSS shifts once (random calls happen here)
        print("Step 3: Extracting the sequences around TSS ...")
        rnas = estimate_tss_probs(
            rnas, rna_input_matches_file, rna_in_input_r1_file, min_len
        )

        # Cache the fully-assigned rnas dict for future --skip-existing runs.
        try:
            with open(rnas_pickle_file, "wb") as fh:
                pickle.dump(rnas, fh, protocol=pickle.HIGHEST_PROTOCOL)
            print(f"Cache: wrote TSS-assigned RNAs to {rnas_pickle_file}",
                  file=sys.stderr)
        except OSError as exc:
            print(f"Warning: could not write rnas cache: {exc}", file=sys.stderr)

    # ...then emit one set of fasta + anno files per (prom, utr) combo.
    for prom_len, utr_len in combos:
        b = out_base(prom_len, utr_len)
        promoter_output_file = f"{b}.promoters.fasta"
        utr_output_file = f"{b}.utrs.fasta"
        both_flanks_output_file = f"{b}.tss_flanks.fasta"
        umi_count_file = f"{b}.umi_counts.tsv"
        anno_file = f"{b}.tss_flank_anno.tsv"

        print(
            f"  -> writing combo prom={prom_len}, utr={utr_len} "
            f"to {b}.*", file=sys.stderr,
        )
        get_tss_seqs(
            rnas,
            promoter_output_file,
            utr_output_file,
            both_flanks_output_file,
            umi_count_file,
            anno_file,
            min_len,
            prom_len,
            utr_len,
        )

        # Step 4 — optional reverse-complement + nucleotide stats (per combo)
        if conf.run_revcoms_and_stats:
            print(
                f"Step 4: Reverse complement and stats for prom={prom_len}, "
                f"utr={utr_len}", file=sys.stderr,
            )
            revcom_file = f"{b}.tss_flanks.revcom.fasta"
            run(f"seqkit seq -r -p {shlex.quote(both_flanks_output_file)} "
                f"-o {shlex.quote(revcom_file)}")

            out_nuc_dist_file = f"{b}.tss_flanks_nuc_dist.txt"
            run(f"cat {shlex.quote(both_flanks_output_file)} "
                f"| fasta-get-markov -norc -m 5 > {shlex.quote(out_nuc_dist_file)}")

            prom_revcom_file = f"{b}.promoters.revcom.fasta"
            run(f"seqkit seq -r -p {shlex.quote(promoter_output_file)} "
                f"-o {shlex.quote(prom_revcom_file)}")

    print("All done!", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
