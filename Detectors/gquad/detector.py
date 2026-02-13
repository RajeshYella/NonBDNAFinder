"""
Ultra-fast G-Quadruplex detector:
Seeded G4Hunter scoring with priority-aware overlap resolution
(Huppert 2005; Bedrat 2016 optimized for genome-scale scanning)
"""

import re
from typing import Dict, List, Tuple, Any
from ..base.base_detector import BaseMotifDetector
from Utilities.core.motif_normalizer import normalize_class_subclass

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

    SEED_PATTERN = re.compile(r'GGG')  # Core seed

    def get_motif_class_name(self) -> str:
        return "G-Quadruplex"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        return {
            'telomeric_g4': [(r'(?:TTAGGG){4,}', 'G4_TEL')],
            'stacked_canonical_g4s': [(r'(?:(?:G{3,}[ACGT]{1,7}){3}G{3,}){2,}', 'G4_STK_CAN')],
            'stacked_g4s_linker': [(r'(?:(?:G{3,}[ACGT]{1,7}){3}G{3,})(?:[ACGT]{0,20}(?:(?:G{3,}[ACGT]{1,7}){3}G{3,})){1,}', 'G4_STK_LNK')],
            'canonical_g4': [(r'G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}', 'G4_CAN')],
            'extended_loop_g4': [(r'G{3,}[ACGT]{1,12}G{3,}[ACGT]{1,12}G{3,}[ACGT]{1,12}G{3,}', 'G4_EXT')],
            'higher_order_g4': [(r'(?:G{3,}[ACGT]{1,7}){7,}', 'G4_HIGH')],
            'g_triplex': [(r'G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}', 'G4_TRX')],
            'weak_pqs': [(r'G{2,}[ACGT]{1,7}G{2,}[ACGT]{1,7}G{2,}[ACGT]{1,7}G{2,}', 'G4_WEAK')],
        }

    # ============================================================
    # 🚀 Seeded Candidate Detection (Genome Optimized)
    # ============================================================

    def _find_all_candidates(self, seq: str) -> List[Dict[str, Any]]:

        patterns = self.get_patterns()
        candidates = []
        seq_len = len(seq)

        for seed_match in self.SEED_PATTERN.finditer(seq):

            seed_pos = seed_match.start()

            # Local window around seed (bounded extension)
            local_start = max(0, seed_pos - 150)
            local_end = min(seq_len, seed_pos + 150)
            local_seq = seq[local_start:local_end]

            for class_name, pattern_list in patterns.items():
                for regex, pattern_id in pattern_list:
                    for m in re.finditer(regex, local_seq):
                        s = local_start + m.start()
                        e = local_start + m.end()

                        if (e - s) >= MIN_REGION_LEN:
                            candidates.append({
                                'class_name': class_name,
                                'pattern_id': pattern_id,
                                'start': s,
                                'end': e
                            })

        return candidates

    # ============================================================
    # G-Only G4Hunter Scoring
    # ============================================================

    def _score_candidate(self, candidate: Dict[str, Any], seq: str) -> Dict[str, Any]:

        s, e = candidate['start'], candidate['end']
        region = seq[s:e]
        L = len(region)

        # G-only scoring (no C penalty)
        vals = [1 if ch == 'G' else 0 for ch in region]

        ws = min(WINDOW_SIZE_DEFAULT, L)
        max_sum = 0

        if ws > 0:
            cur = sum(vals[:ws])
            max_sum = cur

            for i in range(1, L - ws + 1):
                cur += vals[i + ws - 1] - vals[i - 1]
                if cur > max_sum:
                    max_sum = cur

        normalized_window = (max_sum / ws) if ws > 0 else 0.0

        # Tract bonus
        g_tracts = re.findall(r'G{2,}', region)
        n_g = len(g_tracts)
        total_g_len = sum(len(t) for t in g_tracts)

        tract_bonus = 0.0
        if n_g >= 3:
            tract_bonus = min(0.5, 0.08 * (n_g - 2))

        normalized_score = min(1.0, normalized_window + tract_bonus)
        region_score = normalized_score * (L / float(ws)) if ws > 0 else 0.0

        out = candidate.copy()
        out['score'] = float(region_score)
        out['details'] = {
            'n_g_tracts': n_g,
            'total_g_len': total_g_len
        }

        return out

    # ============================================================
    # Priority-Aware Overlap Resolution (UNCHANGED LOGIC)
    # ============================================================

    def _resolve_overlaps(self, scored_candidates: List[Dict[str, Any]]):

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
            conflict = any(not (e <= os or s >= oe) for os, oe in occupied)

            if not conflict:
                accepted.append(cand)
                occupied.append((s, e))

        accepted.sort(key=lambda x: x['start'])
        return accepted

    # ============================================================
    # Public API
    # ============================================================

    def annotate_sequence(self, sequence: str):

        seq = sequence.upper()
        candidates = self._find_all_candidates(seq)
        scored = [self._score_candidate(c, seq) for c in candidates]
        accepted = self._resolve_overlaps(scored)

        return accepted

    def calculate_score(self, sequence: str, pattern_info: Tuple = None) -> float:

        annotations = self.annotate_sequence(sequence)
        if not annotations:
            return 0.0

        return float(sum(a['score'] for a in annotations))

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence"):

        sequence = sequence.upper().strip()
        motifs = []
        annotations = self.annotate_sequence(sequence)

        for ann in annotations:

            canonical_class, canonical_subclass = normalize_class_subclass(
                self.get_motif_class_name(),
                ann['class_name'],
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
                'Score': round(ann['score'], 3),
                'Strand': '+',
                'Method': 'G4Hunter_seeded',
                'Pattern_ID': ann['pattern_id']
            }

            motifs.append(motif)

        return motifs
