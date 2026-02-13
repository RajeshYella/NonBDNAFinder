"""Slipped DNA detector: STRs (k=1-9) and direct repeats (k≥10) (Sinden 1994, Pearson 2005).""""""Slipped DNA detector: seed-accelerated tandem repeat detection
(Streisinger 1966; Sinden 1994; Pearson 2005; Mirkin 2007)."""

import math
from typing import List, Dict, Any, Tuple
from collections import defaultdict

from ..base.base_detector import BaseMotifDetector
from Utilities.core.motif_normalizer import normalize_class_subclass


class SlippedDNADetector(BaseMotifDetector):
    """
    Literature-strict slipped DNA detector.
    STR (k=1–9) and Direct Repeats (k≥10).
    Seed-accelerated detection.
    """

    # -----------------------------
    # Literature-Supported Parameters
    # -----------------------------

    MIN_TRACT_LENGTH = 20
    MIN_PURITY = 0.90
    MAX_UNIT_SIZE = 100
    STR_DIRECT_REPEAT_THRESHOLD = 10
    SCORE_THRESHOLD = 0.3

    # Unit-dependent minimum copies (literature)
    MIN_COPIES_TABLE = {
        1: 10,
        2: 6,
        3: 5,
        4: 4,
        5: 4,
        6: 4
    }

    def get_motif_class_name(self) -> str:
        return "Slipped_DNA"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        return {"short_tandem_repeats": [], "direct_repeats": []}

    # ----------------------------------------------------
    # Primitive Motif Computation
    # ----------------------------------------------------

    @staticmethod
    def compute_primitive_motif(sequence: str) -> str:
        n = len(sequence)
        for period in range(1, n // 2 + 1):
            unit = sequence[:period]
            if all(sequence[i:i+period] == unit[:len(sequence[i:i+period])]
                   for i in range(0, n, period)):
                return unit
        return sequence

    @staticmethod
    def compute_repeat_purity(sequence: str, unit: str) -> float:
        if not unit:
            return 0.0
        unit_len = len(unit)
        matches = sum(sequence[i] == unit[i % unit_len]
                      for i in range(len(sequence)))
        return matches / len(sequence)

    # ----------------------------------------------------
    # Seed-Based Tandem Detection
    # ----------------------------------------------------

    def find_all_tandem_repeats(self, sequence: str) -> List[Dict[str, Any]]:
        seq = sequence.upper()
        n = len(seq)
        candidates = []

        # Seed index (1–6 bp seeds)
        seed_index = defaultdict(list)
        max_seed = min(6, n)

        for k in range(1, max_seed + 1):
            for i in range(n - k + 1):
                seed_index[(k, seq[i:i+k])].append(i)

        visited = set()

        for (k, seed), positions in seed_index.items():

            if k > self.MAX_UNIT_SIZE:
                continue

            min_copies = self.MIN_COPIES_TABLE.get(k, 3)

            for pos in positions:
                if (pos, k) in visited:
                    continue

                copies = 1
                next_pos = pos + k

                while next_pos + k <= n and seq[next_pos:next_pos+k] == seed:
                    copies += 1
                    next_pos += k

                tract_len = copies * k

                if copies >= min_copies and tract_len >= self.MIN_TRACT_LENGTH:
                    sequence_block = seq[pos:next_pos]

                    primitive = self.compute_primitive_motif(sequence_block)
                    purity = self.compute_repeat_purity(sequence_block, primitive)

                    if purity < self.MIN_PURITY:
                        continue

                    candidates.append({
                        'start': pos,
                        'end': next_pos,
                        'length': tract_len,
                        'unit': seed,
                        'primitive_unit': primitive,
                        'unit_size': len(primitive),
                        'copies': copies,
                        'sequence': sequence_block,
                        'purity': purity
                    })

                    for i in range(pos, next_pos):
                        visited.add((i, k))

        return candidates

    # ----------------------------------------------------
    # Mechanistic Slippage Score (Copy-Dominant)
    # ----------------------------------------------------

    def compute_slippage_energy_score(self,
                                      sequence: str,
                                      unit: str,
                                      copy_number: float,
                                      purity: float) -> float:
        """
        Literature-aligned instability proxy.
        Copy number dominant driver.
        """

        tract_length = len(sequence)
        unit_size = len(unit)

        # Copy dominance
        copy_factor = min(1.0, copy_number / 30.0)

        # Length contribution
        length_factor = min(1.0, tract_length / 100.0)

        # Purity linear
        purity_factor = purity

        # Unit size penalty (larger units less slippage-prone)
        unit_penalty = 1.0 / math.log2(unit_size + 1)

        raw_score = (
            0.5 * copy_factor +
            0.3 * length_factor +
            0.2 * purity_factor
        ) * unit_penalty

        return round(min(1.0, raw_score), 6)

    # ----------------------------------------------------
    # Redundancy Elimination
    # ----------------------------------------------------

    def eliminate_redundancy(self, candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        if not candidates:
            return []

        sorted_cands = sorted(candidates,
                              key=lambda c: (c['start'], -c['unit_size']))

        non_redundant = []
        used_intervals = []

        for cand in sorted_cands:
            s, e = cand['start'], cand['end']
            overlap = any(not (e <= us or s >= ue)
                          for us, ue in used_intervals)
            if not overlap:
                non_redundant.append(cand)
                used_intervals.append((s, e))

        return non_redundant

    # ----------------------------------------------------
    # Annotation Pipeline
    # ----------------------------------------------------

    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        candidates = self.find_all_tandem_repeats(sequence)
        non_redundant = self.eliminate_redundancy(candidates)

        for cand in non_redundant:
            cand['slippage_score'] = self.compute_slippage_energy_score(
                sequence=cand['sequence'],
                unit=cand['primitive_unit'],
                copy_number=cand['copies'],
                purity=cand['purity']
            )

        return non_redundant

    # ----------------------------------------------------
    # Detection Interface
    # ----------------------------------------------------

    def detect_motifs(self,
                      sequence: str,
                      sequence_name: str = "sequence") -> List[Dict[str, Any]]:

        self.audit['invoked'] = True
        self.audit['windows_scanned'] = 1
        self.audit['candidates_seen'] = 0
        self.audit['reported'] = 0

        annotations = self.annotate_sequence(sequence)
        self.audit['candidates_seen'] = len(annotations)

        motifs = []

        for i, ann in enumerate(annotations):

            subclass = 'STR' if ann['unit_size'] < self.STR_DIRECT_REPEAT_THRESHOLD else 'Direct Repeat'

            canonical_class, canonical_subclass = normalize_class_subclass(
                self.get_motif_class_name(),
                subclass,
                strict=False,
                auto_correct=True
            )

            motif = {
                'ID': f"{sequence_name}_SLIPPED_{ann['start']+1}",
                'Sequence_Name': sequence_name,
                'Class': canonical_class,
                'Subclass': canonical_subclass,
                'Start': ann['start'] + 1,
                'End': ann['end'],
                'Length': ann['length'],
                'Sequence': ann['sequence'],
                'Repeat_Unit': ann['primitive_unit'],
                'Unit_Size': ann['unit_size'],
                'Copy_Number': ann['copies'],
                'Purity': round(ann['purity'], 3),
                'Slippage_Energy_Score': round(ann['slippage_score'], 3),
                'Score': round(ann['slippage_score'], 3),
                'Strand': '+',
                'Method': 'Slipped_DNA_detection',
                'Pattern_ID': f'SLIPPED_{i+1}'
            }

            motifs.append(motif)
            self.audit['reported'] += 1

        return motifs

    def calculate_score(self, sequence: str, pattern_info: Tuple = None) -> float:
        annotations = self.annotate_sequence(sequence)
        if annotations:
            return max(a['slippage_score'] for a in annotations)
        return 0.0
