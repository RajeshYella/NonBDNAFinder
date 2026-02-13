"""A-philic DNA detector: 10-mer scoring table with hyperscan acceleration."""

import logging
from typing import Any, Dict, List, Tuple

try: from Detectors.base.base_detector import BaseMotifDetector
except ImportError:
    import sys; from pathlib import Path
    parent_dir = str(Path(__file__).parent.parent.parent)
    if parent_dir not in sys.path: sys.path.insert(0, parent_dir)
    from Detectors import BaseMotifDetector

from Utilities.core.motif_normalizer import normalize_class_subclass
from .tenmer_table import TENMER_LOG2

try: from Detectors.zdna import hyperscan_backend; _HYPERSCAN_AVAILABLE = hyperscan_backend.is_hyperscan_available()
except ImportError: _HYPERSCAN_AVAILABLE = False; hyperscan_backend = None

try: from motif_patterns import APHILIC_DNA_PATTERNS
except ImportError: APHILIC_DNA_PATTERNS = {}

logger = logging.getLogger(__name__)


class APhilicDetector(BaseMotifDetector):
    """A-philic DNA detector using 10-mer scoring table."""

    # Tunable parameters
    MIN_SUM_LOG2 = 0.5  # Minimum sum_log2 for A-philic regions

    def get_motif_class_name(self) -> str: return "A-philic_DNA"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        return {"a_philic_10mers": [(r"", "APH_10MER", "A-philic 10-mer table", "A-philic DNA", 10, "a_philic_10mer_score", 0.9, "A-philic 10mer motif", "user_table")]}

    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        seq = sequence.upper(); merged_regions = self._find_and_merge_10mer_matches(seq)
        if not merged_regions: return 0.0
        contrib = self._build_per_base_contrib(seq)
        return float(sum(sum(contrib[s:e]) for s, e in merged_regions))

    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        """Return merged A-philic region annotations."""
        seq = sequence.upper(); matches = self._find_10mer_matches(seq)
        if not matches: return []
        merged = self._merge_matches(matches); contrib = self._build_per_base_contrib(seq); annotations = []
        for region in merged:
            s, e, region_matches = region; sum_log2 = sum(contrib[s:e]); n_10 = len(region_matches)
            mean_per10 = (sum(m[2] for m in region_matches) / n_10) if n_10 > 0 else 0.0
            annotations.append({"start": s, "end": e, "length": e - s, "sum_log2": round(sum_log2, 6),
                                "mean_log2_per10mer": round(mean_per10, 6), "n_10mers": n_10,
                                "contributing_10mers": [{"tenmer": m[1], "start": m[0], "log2": m[2]} for m in region_matches]})
        return annotations

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Detect A-philic regions using 10-mer scoring."""
        self.audit['invoked'] = True
        self.audit['windows_scanned'] = 1
        self.audit['candidates_seen'] = 0
        self.audit['reported'] = 0
        
        sequence = sequence.upper().strip(); motifs = []; annotations = self.annotate_sequence(sequence)
        self.audit['candidates_seen'] = len(annotations)
        
        for i, region in enumerate(annotations):
            if region.get('sum_log2', 0) > self.MIN_SUM_LOG2 and region.get('n_10mers', 0) >= 1:
                start_pos, end_pos = region['start'], region['end']
                canonical_class, canonical_subclass = normalize_class_subclass(self.get_motif_class_name(), 'A-philic DNA', strict=False, auto_correct=True)
                motifs.append({
                    'ID': f"{sequence_name}_APHIL_{start_pos+1}", 'Sequence_Name': sequence_name, 'Class': canonical_class,
                    'Subclass': canonical_subclass, 'Start': start_pos + 1, 'End': end_pos, 'Length': region['length'],
                    'Sequence': sequence[start_pos:end_pos], 'Score': round(region['sum_log2'], 3), 'Strand': '+',
                    'Method': 'A-philic_detection', 'Pattern_ID': f'APHIL_{i+1}'
                })
                self.audit['reported'] += 1
        return motifs

    def _find_10mer_matches(self, seq: str) -> List[Tuple[int, str, float]]:
        """Find all exact 10-mer matches. Uses Hyperscan if available, else pure-Python."""
        if _HYPERSCAN_AVAILABLE:
            try: return self._hs_find_matches(seq)
            except Exception as e:
                logger.warning(f"Hyperscan failed for A-philic, falling back: {e}"); return self._py_find_matches(seq)
        return self._py_find_matches(seq)

    def _hs_find_matches(self, seq: str) -> List[Tuple[int, str, float]]:
        try: return hyperscan_backend.hs_find_matches(seq, TENMER_LOG2)
        except Exception as e: logger.error(f"Hyperscan matching failed: {e}"); raise

    def _py_find_matches(self, seq: str) -> List[Tuple[int, str, float]]:
        try: return hyperscan_backend.py_find_matches(seq, TENMER_LOG2)
        except Exception as e: logger.warning(f"Optimized matching failed: {e}"); return self._py_find_matches_loop(seq)

    def _py_find_matches_loop(self, seq: str) -> List[Tuple[int, str, float]]:
        n = len(seq); matches: List[Tuple[int, str, float]] = []
        for i in range(0, n - 10 + 1):
            ten = seq[i:i + 10]; log2 = TENMER_LOG2.get(ten)
            if log2 is not None: matches.append((i, ten, float(log2)))
        return matches

    def _merge_matches(self, matches: List[Tuple[int, str, float]], merge_gap: int = 0) -> List[Tuple[int, int, List[Tuple[int, str, float]]]]:
        """Merge overlapping/adjacent 10-mer matches."""
        if not matches: return []
        merged = []; cur_start, cur_end = matches[0][0], matches[0][0] + 10; cur_matches = [matches[0]]
        for m in matches[1:]:
            s = m[0]; m_end = s + 10
            if s <= cur_end + merge_gap: cur_end = max(cur_end, m_end); cur_matches.append(m)
            else: merged.append((cur_start, cur_end, cur_matches)); cur_start, cur_end = s, m_end; cur_matches = [m]
        merged.append((cur_start, cur_end, cur_matches))
        return merged

    def _find_and_merge_10mer_matches(self, seq: str, merge_gap: int = 0) -> List[Tuple[int, int]]:
        matches = self._find_10mer_matches(seq); merged = self._merge_matches(matches, merge_gap=merge_gap)
        return [(s, e) for (s, e, _) in merged]

    def _build_per_base_contrib(self, seq: str) -> List[float]:
        """Build per-base contribution array. Each 10-mer's log2 distributed as L/10 to each base."""
        n = len(seq); contrib = [0.0] * n; matches = self._find_10mer_matches(seq)
        for (start, ten, log2) in matches:
            per_base = float(log2) / 10.0
            for k in range(start, min(start + 10, n)): contrib[k] += per_base
        return contrib
