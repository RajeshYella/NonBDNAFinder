"""
Ultra-fast G-Quadruplex detector:
Seeded scanning + G-only G4Hunter scoring
(Huppert 2005, Bedrat 2016 inspired)
"""

import re
from typing import Dict, List, Tuple, Any
from ..base.base_detector import BaseMotifDetector
from Utilities.core.motif_normalizer import normalize_class_subclass

# Tunable constants
WINDOW_SIZE_DEFAULT = 25
MIN_REGION_LEN = 8

CLASS_PRIORITY = [
    "telomeric_g4",
    "stacked_canonical_g4s",
    "stacked_g4s_linker",
    "canonical_g4",
    "extended_loop_g4",
    "higher_order_g4",
    "g_triplex",
    "weak_pqs"
]

class GQuadruplexDetector(BaseMotifDetector):
    """Ultra-fast seeded G4 detector with priority logic retained."""

    # -------------------------
    # Core Interface
    # -------------------------

    def get_motif_class_name(self) -> str:
        return "G-Quadruplex"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        return {
            'telomeric_g4': [(r'(?:TTAGGG){4,}', 'G4_TEL', 'Telomeric G4', 'Telomeric G4')],
            'stacked_canonical_g4s': [(r'(?:(?:G{3,}[ACGT]{1,7}){3}G{3,}){2,}', 'G4_STK_CAN', 'Stacked canonical G4s', 'Stacked canonical G4s')],
            'stacked_g4s_linker': [(r'(?:(?:G{3,}[ACGT]{1,7}){3}G{3,})(?:[ACGT]{0,20}(?:(?:G{3,}[ACGT]{1,7}){3}G{3,})){1,}', 'G4_STK_LNK', 'Stacked G4s with linker', 'Stacked G4s with linker')],
            'canonical_g4': [(r'G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}', 'G4_CAN', 'Canonical intramolecular G4', 'Canonical intramolecular G4')],
            'extended_loop_g4': [(r'G{3,}[ACGT]{1,12}G{3,}[ACGT]{1,12}G{3,}[ACGT]{1,12}G{3,}', 'G4_EXT', 'Extended-loop canonical', 'Extended-loop canonical')],
            'higher_order_g4': [(r'(?:G{3,}[ACGT]{1,7}){7,}', 'G4_HIGH', 'Higher-order G4 array/G4-wire', 'Higher-order G4 array/G4-wire')],
            'g_triplex': [(r'G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}', 'G4_TRX', 'Intramolecular G-triplex', 'Intramolecular G-triplex')],
            'weak_pqs': [(r'G{2,}[ACGT]{1,7}G{2,}[ACGT]{1,7}G{2,}[ACGT]{1,7}G{2,}', 'G4_WEAK', 'Two-tetrad weak PQS', 'Two-tetrad weak PQS')],
        }

    # -------------------------
    # Public API
    # -------------------------

    def calculate_score(self, sequence: str, pattern_info: Tuple = None) -> float:
        annotations = self.annotate_sequence(sequence)
        return float(sum(a['score'] for a in annotations))

    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        seq = sequence.upper()
        candidates = self._seed_and_scan(seq)
        scored = [self._score_candidate(c, seq) for c in candidates]
        accepted = self._resolve_overlaps(scored)
        return accepted

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        sequence = sequence.upper().strip()
        motifs = []
        annotations = self.annotate_sequence(sequence)

        subclass_map = {
            'telomeric_g4': 'Telomeric G4',
            'stacked_canonical_g4s': 'Stacked canonical G4s',
            'stacked_g4s_linker': 'Stacked G4s with linker',
            'canonical_g4': 'Canonical intramolecular G4',
            'extended_loop_g4': 'Extended-loop canonical',
            'higher_order_g4': 'Higher-order G4 array/G4-wire',
            'g_triplex': 'Intramolecular G-triplex',
            'weak_pqs': 'Two-tetrad weak PQS'
        }

        for ann in annotations:
            subclass = subclass_map.get(ann['class_name'], ann['class_name'])
            canonical_class, canonical_subclass = normalize_class_subclass(
                self.get_motif_class_name(),
                subclass,
                strict=False,
                auto_correct=True
            )

            motif = {
                'ID': f"{sequence_name}_{ann['pattern_id']}_{ann['start']+1}",
                'Sequence_Name': sequence_name,
                'Class': canonical_class,
                'Subclass': canonical_subclass,
                'Start': ann['start'] + 1,
                'End': ann['end'],
                'Length': ann['end'] - ann['start'],
                'Sequence': sequence[ann['start']:ann['end']],
                'Score': round(ann['score'], 4),
                'Strand': '+',
                'Method': 'Seeded_G4Hunter',
                'Pattern_ID': ann['pattern_id']
            }

            motifs.append(motif)

        return motifs

    # -------------------------
    # Ultra-Fast Seeding
    # -------------------------

    def _seed_and_scan(self, seq: str) -> List[Dict[str, Any]]:
        """Seed on G3+ tracts, then local regex refinement."""
        candidates = []
        seed_positions = [m.start() for m in re.finditer(r'G{3,}', seq)]

        if not seed_positions:
            return []

        patterns = self.get_patterns()
        scanned_windows = set()

        for seed in seed_positions:
            window_start = max(0, seed - 50)
            window_end = min(len(seq), seed + 200)
            window_key = (window_start, window_end)

            if window_key in scanned_windows:
                continue
            scanned_windows.add(window_key)

            window_seq = seq[window_start:window_end]

            for class_name, pattern_list in patterns.items():
                for pat in pattern_list:
                    regex = pat[0]
                    pattern_id = pat[1] if len(pat) > 1 else f"{class_name}_pat"

                    for m in re.finditer(regex, window_seq):
                        s = window_start + m.start()
                        e = window_start + m.end()
                        if (e - s) >= MIN_REGION_LEN:
                            candidates.append({
                                'class_name': class_name,
                                'pattern_id': pattern_id,
                                'start': s,
                                'end': e
                            })

        return candidates

    # -------------------------
    # G-only G4Hunter Scoring
    # -------------------------

    def _score_candidate(self, candidate: Dict[str, Any], seq: str,
                         window_size: int = WINDOW_SIZE_DEFAULT) -> Dict[str, Any]:

        s, e = candidate['start'], candidate['end']
        region = seq[s:e]
        L = len(region)

        # G-only scoring (C ignored completely)
        vals = [1 if ch == 'G' else 0 for ch in region]
        ws = min(window_size, L)
        max_sum = 0

        if ws > 0:
            cur = sum(vals[:ws])
            max_sum = cur
            for i in range(1, L - ws + 1):
                cur += vals[i + ws - 1] - vals[i - 1]
                if cur > max_sum:
                    max_sum = cur

        normalized_score = max_sum / float(ws) if ws > 0 else 0.0
        region_score = normalized_score * (L / float(ws)) if ws > 0 else 0.0

        out = candidate.copy()
        out['score'] = float(region_score)
        return out

    # -------------------------
    # Priority-based Overlap Resolution
    # -------------------------

    def _resolve_overlaps(self, scored_candidates: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        if not scored_candidates:
            return []

        def class_prio_idx(class_name):
            try:
                return CLASS_PRIORITY.index(class_name)
            except ValueError:
                return len(CLASS_PRIORITY)

        scored_sorted = sorted(
            scored_candidates,
            key=lambda x: (
                -x['score'],
                class_prio_idx(x['class_name']),
                -(x['end'] - x['start'])
            )
        )

        accepted = []
        occupied = []

        for cand in scored_sorted:
            s, e = cand['start'], cand['end']
            conflict = any(not (e <= os or s >= oe) for (os, oe) in occupied)
            if not conflict:
                accepted.append(cand)
                occupied.append((s, e))

        accepted.sort(key=lambda x: x['start'])
        return accepted
