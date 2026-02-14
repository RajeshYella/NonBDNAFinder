"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Triplex DNA Detector - Seed-and-extend mirror repeats (H-DNA)               │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.2           │
│ References: Frank-Kamenetskii 1995; Soyfer & Potaman 1995                   │
└──────────────────────────────────────────────────────────────────────────────┘
"""

import re
from typing import List, Dict, Any, Tuple
from collections import defaultdict
from ..base.base_detector import BaseMotifDetector
from Utilities.core.motif_normalizer import normalize_class_subclass

# Literature constraints
MIN_ARM = 10
MAX_ARM = 100
MAX_LOOP = 8
PURITY_THRESHOLD = 0.90
SEED_SIZE = 6
SCORE_THRESHOLD = 0.25


class TriplexDetector(BaseMotifDetector):

    MIN_ARM = MIN_ARM
    MAX_ARM = MAX_ARM
    MAX_LOOP = MAX_LOOP
    PURITY_THRESHOLD = PURITY_THRESHOLD
    SEED_SIZE = SEED_SIZE
    SCORE_THRESHOLD = SCORE_THRESHOLD

    # ------------------------------------------------------------------ #
    # Required Interface
    # ------------------------------------------------------------------ #

    def get_motif_class_name(self) -> str:
        return "Triplex"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        return {
            "mirror_triplex": [
                (r"", "TRX_MIRROR", "Mirror repeat triplex",
                 "Triplex_Motif", self.MIN_ARM,
                 "structural_triplex_score", self.SCORE_THRESHOLD,
                 "Seed-based mirror repeat detection",
                 "Frank-Kamenetskii 1995")
            ]
        }

    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        return 1.0 if sequence else 0.0

    # ------------------------------------------------------------------ #
    # Mirror Repeat Core (C-style logic reproduction)
    # ------------------------------------------------------------------ #

    def _find_mirror_repeats(self, sequence: str) -> List[Dict[str, Any]]:

        seq = sequence.upper()
        n = len(seq)
        k = self.SEED_SIZE
        hits = []

        seed_index = defaultdict(list)

        # Build seed index
        for i in range(n - k + 1):
            seed_index[seq[i:i+k]].append(i)

        for i in range(n - k + 1):

            seed = seq[i:i+k]
            mirror = seed[::-1]

            if mirror not in seed_index:
                continue

            for j in seed_index[mirror]:

                if j <= i:
                    continue

                loop_len = j - (i + k)

                # Strict spacer constraint
                if loop_len < 0 or loop_len > self.MAX_LOOP:
                    continue

                left_start = i
                right_start = j
                arm_len = k

                # Extend symmetrically (exact mirror)
                while (
                    left_start > 0 and
                    right_start + arm_len < n and
                    seq[left_start - 1] == seq[right_start + arm_len] and
                    arm_len < self.MAX_ARM
                ):
                    left_start -= 1
                    arm_len += 1

                if arm_len < self.MIN_ARM:
                    continue

                left_seq = seq[left_start:left_start + arm_len]

                # Purine / Pyrimidine 90% rule
                purine_fraction = sum(b in "AG" for b in left_seq) / arm_len
                pyr_fraction = sum(b in "CT" for b in left_seq) / arm_len

                if max(purine_fraction, pyr_fraction) < self.PURITY_THRESHOLD:
                    continue

                hits.append({
                    "start": left_start,
                    "end": right_start + arm_len,
                    "arm_length": arm_len,
                    "loop_length": loop_len,
                    "score": 1.0,
                    "sequence": seq[left_start:right_start + arm_len]
                })

        # Sort longest arms first (like C maximizing stem)
        hits.sort(key=lambda x: (-x["arm_length"], x["loop_length"], x["start"]))
        return hits

    # ------------------------------------------------------------------ #
    # Annotation with overlap control
    # ------------------------------------------------------------------ #

    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:

        seq = sequence.upper()
        results = []
        used = [False] * len(seq)

        mirrors = self._find_mirror_repeats(seq)

        for m in mirrors:

            s = m["start"]
            e = m["end"]

            if any(used[s:e]):
                continue

            for i in range(s, e):
                used[i] = True

            results.append({
                "class_name": "Triplex",
                "pattern_id": "TRX_MIRROR",
                "start": s,
                "end": e,
                "length": e - s,
                "score": m["score"],
                "matched_seq": m["sequence"]
            })

        results.sort(key=lambda r: r["start"])
        return results

    # ------------------------------------------------------------------ #
    # Final Detection Interface
    # ------------------------------------------------------------------ #

    def detect_motifs(self,
                      sequence: str,
                      sequence_name: str = "sequence") -> List[Dict[str, Any]]:

        sequence = sequence.upper().strip()
        motifs = []
        annotations = self.annotate_sequence(sequence)

        for annotation in annotations:

            canonical_class, canonical_subclass = normalize_class_subclass(
                self.get_motif_class_name(),
                "Triplex_Motif",
                strict=False,
                auto_correct=True
            )

            motifs.append({
                "ID": f"{sequence_name}_TRX_{annotation['start']+1}",
                "Sequence_Name": sequence_name,
                "Class": canonical_class,
                "Subclass": canonical_subclass,
                "Start": annotation["start"] + 1,
                "End": annotation["end"],
                "Length": annotation["length"],
                "Sequence": annotation["matched_seq"],
                "Score": 1.0,
                "Strand": "+",
                "Method": "Triplex_seed_mirror_detection",
                "Pattern_ID": annotation["pattern_id"]
            })

        return motifs
