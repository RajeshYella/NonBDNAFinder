"""Triplex DNA detector: mirror repeats (12-100 nt arms, ≤8 nt spacer) 
with structural triplet scoring and separate Sticky DNA class.
Frank-Kamenetskii 1995; Soyfer & Potaman 1995; Lexa 2011; Sakamoto 1999.
"""

import re
from typing import List, Dict, Any, Tuple
from ..base.base_detector import BaseMotifDetector
from Utilities.core.motif_normalizer import normalize_class_subclass

try:
    from scanner import find_mirror_repeats_optimized as _find_mirror_repeats_optimized
except ImportError:
    _find_mirror_repeats_optimized = None

try:
    from motif_patterns import TRIPLEX_PATTERNS
except ImportError:
    TRIPLEX_PATTERNS = {}


def revcomp(seq: str) -> str:
    comp = str.maketrans('ATGCRYSWKMBDHVN', 'TACGYRSWMKVHDBN')
    return seq.translate(comp)[::-1]


class TriplexDetector(BaseMotifDetector):
    """
    Triplex DNA detector:
    - Canonical intramolecular triplex (H-DNA) via mirror repeats
    - Sticky DNA (GAA/TTC expansions) as separate subclass
    """

    # Literature-supported parameters
    MIN_ARM = 10
    MAX_ARM = 100
    MAX_LOOP = 8
    PURINE_PYRIMIDINE_THRESHOLD = 0.90

    def get_motif_class_name(self) -> str:
        return "Triplex"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        return self._load_patterns(TRIPLEX_PATTERNS, lambda: {
            'triplex_forming_sequences': [
                (r'(?:GAA){4,}', 'TRX_STICKY_GAA', 'GAA repeat', 'Sticky_DNA', 12,
                 'sticky_dna_score', 0.95, 'Disease-associated repeats', 'Sakamoto 1999'),
                (r'(?:TTC){4,}', 'TRX_STICKY_TTC', 'TTC repeat', 'Sticky_DNA', 12,
                 'sticky_dna_score', 0.95, 'Disease-associated repeats', 'Sakamoto 1999'),
            ]
        })

    # ------------------------------------------------------------------
    # Structural Triplet Scoring (Literature-Based)
    # ------------------------------------------------------------------

    def _triplet_score(self, base: str) -> int:
        """
        Literature-supported canonical triplets.
        Strong = 2
        Moderate = 1
        Others = 0
        """
        strong = {'A', 'T'}   # A·A:T , T·A:T (stable)
        moderate = {'G', 'C'} # G·G:C , C+·G:C (context dependent)

        if base in strong:
            return 2
        elif base in moderate:
            return 1
        return 0

    def _structural_triplex_score(self, arm_seq: str, loop_len: int) -> float:
        """
        Structural stability score (0–1 normalized).
        Based on canonical triplet contribution and loop penalty.
        """
        if not arm_seq:
            return 0.0

        total_score = sum(self._triplet_score(b) for b in arm_seq)
        max_score = 2 * len(arm_seq)

        if max_score == 0:
            return 0.0

        triplet_component = total_score / max_score

        # Optimal loop ≈ 5 nt (Haasnoot 1986)
        optimal_loop = 5
        loop_penalty = max(0.0, 1 - abs(loop_len - optimal_loop) / optimal_loop)

        return round(triplet_component * loop_penalty, 6)

    # ------------------------------------------------------------------
    # Sticky DNA Scoring (Expansion-Based)
    # ------------------------------------------------------------------

    def _sticky_dna_score(self, sequence: str) -> float:
        """
        Expansion-driven scoring for GAA/TTC.
        Length-dependent instability model.
        """
        if len(sequence) < 12:
            return 0.0

        repeat_units = len(sequence) // 3
        uninterrupted_bonus = min(repeat_units / 20.0, 1.0)

        return round(min(0.5 + uninterrupted_bonus * 0.5, 1.0), 6)

    # ------------------------------------------------------------------
    # Core Annotation
    # ------------------------------------------------------------------

    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        seq = sequence.upper()
        results = []
        used = [False] * len(seq)

        patterns = self.get_patterns()['triplex_forming_sequences']

        # -------- Mirror Repeat Triplex --------
        if _find_mirror_repeats_optimized is not None:
            mirror_results = _find_mirror_repeats_optimized(
                seq,
                min_arm=self.MIN_ARM,
                max_arm=self.MAX_ARM,
                max_loop=self.MAX_LOOP,
                purine_pyrimidine_threshold=self.PURINE_PYRIMIDINE_THRESHOLD
            )

            for mr in mirror_results:
                if not mr.get('Is_Triplex', False):
                    continue

                s = mr['Start'] - 1
                e = mr['End']
                if any(used[s:e]):
                    continue

                for i in range(s, e):
                    used[i] = True

                left_arm = mr.get('Left_Arm', '')
                loop_len = mr.get('Loop', 0)

                score = self._structural_triplex_score(left_arm, loop_len)

                results.append({
                    'class_name': 'Triplex',
                    'pattern_id': 'TRX_MIRROR',
                    'start': s,
                    'end': e,
                    'length': e - s,
                    'score': score,
                    'matched_seq': seq[s:e],
                    'details': {
                        'type': 'Mirror repeat (H-DNA)',
                        'reference': 'Frank-Kamenetskii 1995',
                        'arm_length': mr.get('Arm_Length', 0),
                        'loop_length': loop_len,
                        'left_arm': left_arm,
                        'right_arm': mr.get('Right_Arm', ''),
                        'loop_seq': mr.get('Loop_Seq', '')
                    }
                })

        # -------- Sticky DNA --------
        for patinfo in patterns:
            pat, pid, name, cname, minlen, scoretype, cutoff, desc, ref = patinfo
            for m in re.finditer(pat, seq):
                s, e = m.span()
                if any(used[s:e]):
                    continue

                for i in range(s, e):
                    used[i] = True

                score = self._sticky_dna_score(seq[s:e])

                results.append({
                    'class_name': cname,
                    'pattern_id': pid,
                    'start': s,
                    'end': e,
                    'length': e - s,
                    'score': score,
                    'matched_seq': seq[s:e],
                    'details': {
                        'type': name,
                        'reference': ref,
                        'description': desc
                    }
                })

        results.sort(key=lambda r: r['start'])
        return results

    # ------------------------------------------------------------------
    # Integration Layer
    # ------------------------------------------------------------------

    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        scoring_method = pattern_info[5] if len(pattern_info) > 5 else 'structural'
        if scoring_method == 'sticky_dna_score':
            return self._sticky_dna_score(sequence)
        return 0.0

    def passes_quality_threshold(self, sequence: str, score: float, pattern_info: Tuple) -> bool:
        return score >= 0.25

    # ------------------------------------------------------------------
    # Detection (Both Strands)
    # ------------------------------------------------------------------

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        self.audit['invoked'] = True
        self.audit['windows_scanned'] = 0
        self.audit['candidates_seen'] = 0
        self.audit['candidates_filtered'] = 0
        self.audit['reported'] = 0
        self.audit['both_strands_scanned'] = True

        sequence = sequence.upper().strip()
        motifs = []

        # Forward strand
        results_fwd = self.annotate_sequence(sequence)
        self.audit['candidates_seen'] += len(results_fwd)

        for result in results_fwd:
            canonical_class, canonical_subclass = normalize_class_subclass(
                self.get_motif_class_name(),
                'Sticky DNA' if result['class_name'] == 'Sticky_DNA' else 'Triplex',
                strict=False,
                auto_correct=True
            )

            motif_dict = {
                'ID': f"{sequence_name}_{result['pattern_id']}_{result['start']+1}",
                'Sequence_Name': sequence_name,
                'Class': canonical_class,
                'Subclass': canonical_subclass,
                'Start': result['start'] + 1,
                'End': result['end'],
                'Length': result['length'],
                'Sequence': result['matched_seq'],
                'Score': round(result['score'], 3),
                'Strand': '+',
                'Method': 'Triplex_detection',
                'Pattern_ID': result['pattern_id']
            }

            motifs.append(motif_dict)
            self.audit['reported'] += 1

        # Reverse strand
        seq_rc = revcomp(sequence)
        results_rc = self.annotate_sequence(seq_rc)
        seq_len = len(sequence)

        for result in results_rc:
            fwd_start = seq_len - result['end']
            fwd_end = seq_len - result['start']

            canonical_class, canonical_subclass = normalize_class_subclass(
                self.get_motif_class_name(),
                'Sticky DNA' if result['class_name'] == 'Sticky_DNA' else 'Triplex',
                strict=False,
                auto_correct=True
            )

            motif_dict = {
                'ID': f"{sequence_name}_{result['pattern_id']}_RC_{fwd_start+1}",
                'Sequence_Name': sequence_name,
                'Class': canonical_class,
                'Subclass': canonical_subclass,
                'Start': fwd_start + 1,
                'End': fwd_end,
                'Length': result['length'],
                'Sequence': sequence[fwd_start:fwd_end],
                'Score': round(result['score'], 3),
                'Strand': '-',
                'Method': 'Triplex_detection',
                'Pattern_ID': result['pattern_id'] + '_RC'
            }

            motifs.append(motif_dict)
            self.audit['reported'] += 1

        return motifs
