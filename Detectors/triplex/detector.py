"""Triplex DNA detector: seed-and-extend mirror repeats (H-DNA)
Frank-Kamenetskii 1995; Soyfer & Potaman 1995
Sticky DNA handled separately.
"""

import re
from typing import List, Dict, Any, Tuple
from collections import defaultdict
from ..base.base_detector import BaseMotifDetector
from Utilities.core.motif_normalizer import normalize_class_subclass


def revcomp(seq: str) -> str:
    comp = str.maketrans('ATGC', 'TACG')
    return seq.translate(comp)[::-1]


class TriplexDetector(BaseMotifDetector):

    # Literature parameters
    MIN_ARM = 10
    MAX_ARM = 100
    MAX_LOOP = 8
    PURITY_THRESHOLD = 0.90
    SEED_SIZE = 6
    SCORE_THRESHOLD = 0.25

    def get_motif_class_name(self) -> str:
        return "Triplex"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return patterns for Triplex/H-DNA detection."""
        return {
            'mirror_triplex': [
                (r'', 'TRX_MIRROR', 'Mirror repeat triplex', 'H-DNA', self.MIN_ARM, 'structural_triplex_score', 0.25, 'Seed-based mirror repeat detection', 'Frank-Kamenetskii 1995')
            ],
            'sticky_dna': [
                (r'(?:GAA){4,}', 'TRX_STICKY_GAA', 'GAA repeat (Sticky DNA)', 'Sticky_DNA', 12, 'sticky_dna_score', 0.20, 'Friedreich ataxia-associated GAA repeat', 'Soyfer & Potaman 1995'),
                (r'(?:TTC){4,}', 'TRX_STICKY_TTC', 'TTC repeat (Sticky DNA)', 'Sticky_DNA', 12, 'sticky_dna_score', 0.20, 'TTC repeat (complement of GAA)', 'Soyfer & Potaman 1995'),
            ]
        }

    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        """Calculate triplex score based on pattern type."""
        if not sequence:
            return 0.0
        
        seq = sequence.upper()
        
        # Check if this is a sticky DNA pattern
        if pattern_info and len(pattern_info) > 1:
            pattern_id = pattern_info[1]
            if 'STICKY' in pattern_id:
                return self._sticky_dna_score(seq)
        
        # For mirror-based triplex, use structural score with default loop length
        # When called from base class detect_motifs, we don't have arm/loop info
        # so we estimate based on sequence length
        estimated_loop = min(self.MAX_LOOP, max(1, len(seq) // 10))
        return self._structural_triplex_score(seq, estimated_loop)

    # ---------------------------------------------------
    # Structural Stability Model (Non-thermodynamic)
    # ---------------------------------------------------

    def _structural_triplex_score(self, arm_seq: str, loop_len: int) -> float:

        if not arm_seq:
            return 0.0

        length_term = len(arm_seq) / self.MAX_ARM

        purine_fraction = sum(1 for b in arm_seq if b in "AG") / len(arm_seq)
        pyr_fraction = sum(1 for b in arm_seq if b in "CT") / len(arm_seq)

        purity_term = max(purine_fraction, pyr_fraction)

        optimal_loop = 5
        loop_term = max(0.0, 1 - abs(loop_len - optimal_loop) / optimal_loop)

        score = length_term * purity_term * loop_term

        return round(min(1.0, max(0.0, score)), 6)

    # ---------------------------------------------------
    # Seed-Based Mirror Detection
    # ---------------------------------------------------

    def _find_mirror_repeats(self, sequence: str) -> List[Dict[str, Any]]:

        seq = sequence.upper()
        n = len(seq)
        k = self.SEED_SIZE

        hits = []

        seed_index = defaultdict(list)

        # Index k-mers
        for i in range(n - k + 1):
            seed = seq[i:i+k]
            seed_index[seed].append(i)

        # Scan seeds
        for i in range(n - k + 1):
            seed = seq[i:i+k]
            mirror = seed[::-1]  # mirror, not reverse complement

            if mirror not in seed_index:
                continue

            for j in seed_index[mirror]:

                if j <= i:
                    continue

                loop_len = j - (i + k)

                if loop_len < 0 or loop_len > self.MAX_LOOP:
                    continue

                left_start = i
                right_start = j
                arm_len = k

                # Extend symmetrically
                while True:

                    if left_start == 0 or right_start + arm_len >= n:
                        break

                    next_left = seq[left_start - 1]
                    next_right = seq[right_start + arm_len]

                    if next_left != next_right:
                        break

                    left_start -= 1
                    arm_len += 1

                    if arm_len >= self.MAX_ARM:
                        break

                if arm_len < self.MIN_ARM:
                    continue

                left_seq = seq[left_start:left_start + arm_len]

                # Purity requirement
                purine_fraction = sum(1 for b in left_seq if b in "AG") / len(left_seq)
                pyr_fraction = sum(1 for b in left_seq if b in "CT") / len(left_seq)

                if max(purine_fraction, pyr_fraction) < self.PURITY_THRESHOLD:
                    continue

                score = self._structural_triplex_score(left_seq, loop_len)

                hits.append({
                    "start": left_start,
                    "end": right_start + arm_len,
                    "arm_length": arm_len,
                    "loop_length": loop_len,
                    "score": score,
                    "sequence": seq[left_start:right_start + arm_len]
                })

        hits.sort(key=lambda x: (-x["score"], x["start"]))
        return hits

    # ---------------------------------------------------
    # Sticky DNA
    # ---------------------------------------------------

    def _sticky_dna_score(self, sequence: str) -> float:
        repeat_units = len(sequence) // 3
        return round(min(1.0, repeat_units / 20.0), 6)

    # ---------------------------------------------------
    # Annotation
    # ---------------------------------------------------

    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:

        seq = sequence.upper()
        results = []
        used = [False] * len(seq)

        # Mirror-based triplex
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

        # Sticky DNA (GAA/TTC)
        for pat, pid in [(r'(?:GAA){4,}', 'TRX_STICKY_GAA'),
                         (r'(?:TTC){4,}', 'TRX_STICKY_TTC')]:

            for m in re.finditer(pat, seq):
                s, e = m.span()

                if any(used[s:e]):
                    continue

                for i in range(s, e):
                    used[i] = True

                results.append({
                    "class_name": "Sticky_DNA",
                    "pattern_id": pid,
                    "start": s,
                    "end": e,
                    "length": e - s,
                    "score": self._sticky_dna_score(seq[s:e]),
                    "matched_seq": seq[s:e]
                })

        results.sort(key=lambda r: r["start"])
        return results
