"""
Slipped DNA detector: seed-accelerated tandem repeat detection
(Streisinger 1966; Sinden 1994; Pearson 2005; Mirkin 2007)

Subclasses:
- Slipped STRs (k = 1–9)
- Slipped DRs (k ≥ 10)
"""

import math
from typing import List, Dict, Any, Tuple
from collections import defaultdict

from ..base.base_detector import BaseMotifDetector
from Utilities.core.motif_normalizer import normalize_class_subclass


class SlippedDNADetector(BaseMotifDetector):

    # -----------------------------
    # Strict Literature Parameters
    # -----------------------------

    MIN_TRACT_LENGTH = 20
    MIN_PURITY = 0.90
    MAX_UNIT_SIZE = 100
    STR_DIRECT_REPEAT_THRESHOLD = 10
    SCORE_THRESHOLD = 0.30

    # Literature-informed minimum copies
    MIN_COPIES_TABLE = {
        1: 10,   # Homopolymers unstable below ~10
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
    # Primitive Motif (fast version)
    # ----------------------------------------------------

    @staticmethod
    def compute_primitive_motif(sequence: str) -> str:
        n = len(sequence)
        for period in range(1, min(50, n // 2 + 1)):
            if n % period != 0:
                continue
            unit = sequence[:period]
            if unit * (n // period) == sequence:
                return unit
        return sequence

    @staticmethod
    def compute_repeat_purity(sequence: str, unit: str) -> float:
        if not unit:
            return 0.0
        ulen = len(unit)
        matches = 0
        for i, base in enumerate(sequence):
            if base == unit[i % ulen]:
                matches += 1
        return matches / len(sequence)

    # ----------------------------------------------------
    # 🚀 Seed-Accelerated Detection (O(n))
    # ----------------------------------------------------

    def find_all_tandem_repeats(self, sequence: str) -> List[Dict[str, Any]]:

        seq = sequence.upper()
        n = len(seq)
        candidates = []

        # iterate possible unit sizes
        for k in range(1, min(self.MAX_UNIT_SIZE, n // 2) + 1):

            min_copies = self.MIN_COPIES_TABLE.get(k, 3)

            i = 0
            while i <= n - k * min_copies:

                unit = seq[i:i+k]

                # fast seed check:
                # require next copy to match before full extension
                if seq[i+k:i+2*k] != unit:
                    i += 1
                    continue

                # extend copies
                j = i
                copies = 0
                while j + k <= n and seq[j:j+k] == unit:
                    copies += 1
                    j += k

                tract_length = copies * k

                if copies >= min_copies and tract_length >= self.MIN_TRACT_LENGTH:

                    sequence_block = seq[i:j]
                    primitive = self.compute_primitive_motif(sequence_block)
                    purity = self.compute_repeat_purity(sequence_block, primitive)

                    if purity >= self.MIN_PURITY:

                        candidates.append({
                            'start': i,
                            'end': j,
                            'length': tract_length,
                            'primitive_unit': primitive,
                            'unit_size': len(primitive),
                            'copies': copies,
                            'sequence': sequence_block,
                            'purity': purity
                        })

                    i = j  # skip whole tract (performance critical)
                else:
                    i += 1

        return candidates

    # ----------------------------------------------------
    # Copy-Dominant Slippage Score
    # ----------------------------------------------------

    def compute_slippage_energy_score(self,
                                      sequence: str,
                                      unit: str,
                                      copy_number: float,
                                      purity: float) -> float:

        tract_length = len(sequence)
        unit_size = len(unit)

        copy_factor = min(1.0, copy_number / 30.0)
        length_factor = min(1.0, tract_length / 100.0)
        purity_factor = purity
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

    def eliminate_redundancy(self, candidates):

        if not candidates:
            return []

        candidates.sort(key=lambda x: (x['start'], -x['unit_size']))

        final = []
        occupied = []

        for cand in candidates:
            s, e = cand['start'], cand['end']
            overlap = any(not (e <= os or s >= oe) for os, oe in occupied)
            if not overlap:
                final.append(cand)
                occupied.append((s, e))

        return final

    # ----------------------------------------------------
    # Annotation Pipeline
    # ----------------------------------------------------

    def annotate_sequence(self, sequence: str):

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

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence"):

        self.audit['invoked'] = True
        self.audit['windows_scanned'] = 1
        self.audit['candidates_seen'] = 0
        self.audit['reported'] = 0

        annotations = self.annotate_sequence(sequence)
        self.audit['candidates_seen'] = len(annotations)

        motifs = []

        for i, ann in enumerate(annotations):

            subclass = (
                'Slipped STRs'
                if ann['unit_size'] < self.STR_DIRECT_REPEAT_THRESHOLD
                else 'Slipped DRs'
            )

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
