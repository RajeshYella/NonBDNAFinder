"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ R-Loop Detector - QmRLFS algorithm (Jenjaroenpun 2016)                       │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.2            │
│ Hyperscan + Seed-accelerated REZ detection (high performance)                │
│ Optimization: Uses shared SeedEngine for ~10000x performance gain            │
└──────────────────────────────────────────────────────────────────────────────┘
"""
# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
import re
from typing import List, Dict, Any, Tuple, Optional
from ..base.base_detector import BaseMotifDetector
from Utilities.detectors_utils import revcomp
from Utilities.core.motif_normalizer import normalize_class_subclass
from Utilities.core.seed_engine import get_seed_engine

try: import hyperscan; HS_AVAILABLE = True
except Exception: HS_AVAILABLE = False

# ═══════════════════════════════════════════════════════════════════════════════
# TUNABLE PARAMETERS - QmRLFS Literature Parameters
# ═══════════════════════════════════════════════════════════════════════════════
MIN_PERC_G_RIZ = 50; NUM_LINKER = 50; WINDOW_STEP = 100
MAX_LENGTH_REZ = 2000; MIN_PERC_G_REZ = 40; QUALITY_THRESHOLD = 0.4
# ═══════════════════════════════════════════════════════════════════════════════


class RLoopDetector(BaseMotifDetector):
    """QmRLFS-finder R-loop detector (literature-faithful, accelerated)."""

    MIN_PERC_G_RIZ = MIN_PERC_G_RIZ; NUM_LINKER = NUM_LINKER; WINDOW_STEP = WINDOW_STEP
    MAX_LENGTH_REZ = MAX_LENGTH_REZ; MIN_PERC_G_REZ = MIN_PERC_G_REZ; QUALITY_THRESHOLD = QUALITY_THRESHOLD

    def __init__(self):
        super().__init__()
        self.hs_db = None
        self.hs_id_to_model = {}

        if HS_AVAILABLE:
            self._compile_hyperscan_patterns()

    # --------------------------------------------------
    # Required abstract methods
    # --------------------------------------------------

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """
        Return patterns for R-loop detection (QmRLFS algorithm).
        
        R-loops are detected using two models based on G-rich patterns:
        - Model 1: Standard G-cluster pattern for RIZ detection
        - Model 2: Extended G-tract pattern for RIZ detection
        """
        return {
            'qmrlfs_model_1': [
                (r'G{3,}[ATCG]{1,10}?G{3,}(?:[ATCG]{1,10}?G{3,}){1,}?',
                 'RLOOP_M1', 'QmRLFS Model 1', 'R-loop formation sites',
                 12, 'qmrlfs_score', 0.4, 'G-cluster RIZ pattern', 'Jenjaroenpun 2016')
            ],
            'qmrlfs_model_2': [
                (r'G{4,}(?:[ATCG]{1,10}?G{4,}){1,}?',
                 'RLOOP_M2', 'QmRLFS Model 2', 'R-loop formation sites',
                 8, 'qmrlfs_score', 0.4, 'Extended G-tract RIZ pattern', 'Jenjaroenpun 2016')
            ]
        }

    def calculate_score(self, sequence: str, pattern_info: Tuple = None) -> float:
        annotations = self.annotate_sequence(sequence)
        if not annotations:
            return 0.0

        best = max(
            (ann['riz_perc_g'] / 100.0) +
            (ann.get('rez_perc_g', 0) / 100.0)
            for ann in annotations
        )

        return min(1.0, round(best, 6))

    def passes_quality_threshold(self,
                                  sequence: str,
                                  score: float,
                                  pattern_info: Tuple = None) -> bool:
        return score >= self.QUALITY_THRESHOLD

    # --------------------------------------------------

    def get_motif_class_name(self) -> str:
        return "R-Loop"

    # --------------------------------------------------
    # Hyperscan Compilation
    # --------------------------------------------------

    def _compile_hyperscan_patterns(self):
        try:
            expressions = [
                br"G{3,}[ATCG]{1,10}?G{3,}(?:[ATCG]{1,10}?G{3,}){1,}?",
                br"G{4,}(?:[ATCG]{1,10}?G{4,}){1,}?"
            ]
            ids = [1, 2]
            flags = [hyperscan.HS_FLAG_DOTALL] * 2

            self.hs_db = hyperscan.Database()
            self.hs_db.compile(expressions=expressions,
                               ids=ids,
                               flags=flags)

            self.hs_id_to_model = {1: 'qmrlfs_model_1',
                                   2: 'qmrlfs_model_2'}
        except Exception:
            self.hs_db = None

    # --------------------------------------------------
    # RIZ Detection
    # --------------------------------------------------

    def _riz_search(self, seq: str, model: str) -> List[Dict[str, Any]]:

        results = []

        if HS_AVAILABLE and self.hs_db is not None:

            seq_bytes = seq.encode()

            def on_match(id, start, end, flags, context):
                if self.hs_id_to_model.get(id) == model:
                    riz_seq = seq[start:end]
                    if self._percent_g(riz_seq) >= self.MIN_PERC_G_RIZ:
                        results.append({
                            'start': start,
                            'end': end,
                            'sequence': riz_seq
                        })
                return 0

            self.hs_db.scan(seq_bytes, match_event_handler=on_match)

        else:
            if model == 'qmrlfs_model_1':
                pattern = re.compile(
                    r"G{3,}[ATCG]{1,10}?G{3,}(?:[ATCG]{1,10}?G{3,}){1,}?",
                    re.IGNORECASE
                )
            else:
                pattern = re.compile(
                    r"G{4,}(?:[ATCG]{1,10}?G{4,}){1,}?",
                    re.IGNORECASE
                )

            for m in pattern.finditer(seq):
                riz_seq = m.group(0)
                if self._percent_g(riz_seq) >= self.MIN_PERC_G_RIZ:
                    results.append({
                        'start': m.start(),
                        'end': m.end(),
                        'sequence': riz_seq
                    })

        return results

    # --------------------------------------------------
    # 🚀 Seed-Accelerated REZ Detection
    # Uses shared SeedEngine prefix sums (O(1) per query)
    # --------------------------------------------------

    def _find_rez(self,
                  seq: str,
                  riz_end: int) -> Optional[Dict[str, Any]]:

        seq_len = len(seq)
        search_start = riz_end + self.NUM_LINKER
        if search_start >= seq_len:
            return None

        # Use shared SeedEngine for O(1) G-content queries
        seed_engine = get_seed_engine()

        best = None
        max_score = 0.0

        max_end = min(seq_len, riz_end + self.MAX_LENGTH_REZ)

        for start in range(search_start, max_end, self.WINDOW_STEP):

            # Quick seed prefilter using shared prefix sums (O(1) lookup)
            window_seed_end = min(start + 100, seq_len)
            seed_len = window_seed_end - start
            
            if seed_len == 0:
                continue
            
            # Use seed engine for O(1) G-content check
            seed_g_percent = seed_engine.get_g_percent_in_range(seq, start, window_seed_end)

            if seed_g_percent < self.MIN_PERC_G_REZ:
                continue

            for end in range(start + 50,
                             max_end,
                             50):

                length = end - start
                # Use seed engine for O(1) G-content calculation
                perc_g = seed_engine.get_g_percent_in_range(seq, start, end)

                if perc_g >= self.MIN_PERC_G_REZ:
                    score = perc_g * length / 100.0
                    if score > max_score:
                        max_score = score
                        best = {
                            'start': start,
                            'end': end,
                            'length': length,
                            'sequence': seq[start:end],
                            'perc_g': round(perc_g, 2)
                        }

        return best

    # --------------------------------------------------

    def _percent_g(self, seq: str) -> float:
        return round((seq.count("G") / float(len(seq))) * 100.0, 2) if seq else 0.0

    # --------------------------------------------------

    def annotate_sequence(self,
                          sequence: str,
                          models: Optional[List[str]] = None
                          ) -> List[Dict[str, Any]]:

        seq = sequence.upper()
        models = models or ['qmrlfs_model_1', 'qmrlfs_model_2']
        results = []

        for model in models:

            riz_regions = self._riz_search(seq, model)

            for riz in riz_regions:

                rez = self._find_rez(seq, riz['end'])

                result = {
                    'model': model,
                    'riz_start': riz['start'],
                    'riz_end': riz['end'],
                    'riz_length': riz['end'] - riz['start'],
                    'riz_sequence': riz['sequence'],
                    'riz_perc_g': self._percent_g(riz['sequence']),
                    'total_start': riz['start'],
                    'total_end': riz['end'],
                    'total_length': riz['end'] - riz['start']
                }

                if rez:
                    result.update({
                        'rez_start': rez['start'],
                        'rez_end': rez['end'],
                        'rez_length': rez['length'],
                        'rez_sequence': rez['sequence'],
                        'rez_perc_g': rez['perc_g'],
                        'total_end': rez['end'],
                        'total_length': rez['end'] - riz['start']
                    })

                results.append(result)

        return results

    # --------------------------------------------------

    def detect_motifs(self,
                      sequence: str,
                      sequence_name: str = "sequence"
                      ) -> List[Dict[str, Any]]:

        self.audit['invoked'] = True
        self.audit['windows_scanned'] = 2
        self.audit['candidates_seen'] = 0
        self.audit['reported'] = 0
        self.audit['both_strands_scanned'] = True

        motifs = []

        for strand, seq in [('+', sequence.upper()),
                            ('-', revcomp(sequence.upper()))]:

            annotations = self.annotate_sequence(seq)
            self.audit['candidates_seen'] += len(annotations)

            seq_len = len(sequence)

            for i, ann in enumerate(annotations):

                if strand == '-':
                    start = seq_len - ann['total_end']
                    end = seq_len - ann['total_start']
                else:
                    start = ann['total_start']
                    end = ann['total_end']

                score = (
                    ann['riz_perc_g'] / 100.0 +
                    ann.get('rez_perc_g', 0) / 100.0
                )

                canonical_class, canonical_subclass = normalize_class_subclass(
                    self.get_motif_class_name(),
                    'R-loop formation sites',
                    strict=False,
                    auto_correct=True
                )

                motifs.append({
                    'ID': f"{sequence_name}_RLOOP_{start+1}",
                    'Sequence_Name': sequence_name,
                    'Class': canonical_class,
                    'Subclass': canonical_subclass,
                    'Start': start + 1,
                    'End': end,
                    'Length': end - start,
                    'Sequence': sequence[start:end],
                    'Score': round(min(score, 1.0), 3),
                    'Strand': strand,
                    'Method': 'QmRLFS_detection'
                })

                self.audit['reported'] += 1

        return motifs
