"""
R-Loop detector: QmRLFS algorithm (Jenjaroenpun 2016)
Seed-accelerated + Hyperscan optimized
Literature-faithful implementation
"""

import re
from typing import List, Dict, Any, Tuple, Optional
from ..base.base_detector import BaseMotifDetector
from Utilities.detectors_utils import revcomp
from Utilities.core.motif_normalizer import normalize_class_subclass

try:
    import hyperscan
    HS_AVAILABLE = True
except Exception:
    HS_AVAILABLE = False

try:
    from motif_patterns import RLOOP_PATTERNS
except ImportError:
    RLOOP_PATTERNS = {}


class RLoopDetector(BaseMotifDetector):

    # Literature parameters (UNCHANGED)
    MIN_PERC_G_RIZ = 50
    NUM_LINKER = 50
    WINDOW_STEP = 100
    MAX_LENGTH_REZ = 2000
    MIN_PERC_G_REZ = 40
    QUALITY_THRESHOLD = 0.4

    def __init__(self):
        super().__init__()
        self.hs_db = None
        self.hs_id_to_model = {}
        if HS_AVAILABLE:
            self._compile_hyperscan_patterns()

    def get_motif_class_name(self) -> str:
        return "R-Loop"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        return self._load_patterns(RLOOP_PATTERNS, lambda: {
            'qmrlfs_model_1': [(r'G{3,}[ATCG]{1,10}?G{3,}(?:[ATCG]{1,10}?G{3,}){1,}?',)],
            'qmrlfs_model_2': [(r'G{4,}(?:[ATCG]{1,10}?G{4,}){1,}?',)]
        })

    # ------------------------------------------------------------
    # Hyperscan Compilation
    # ------------------------------------------------------------

    def _compile_hyperscan_patterns(self):
        try:
            expressions = []
            ids = []
            flags = []
            next_id = 1

            for model, plist in self.get_patterns().items():
                for p in plist:
                    expressions.append(p[0].encode())
                    ids.append(next_id)
                    flags.append(hyperscan.HS_FLAG_DOTALL)
                    self.hs_id_to_model[next_id] = model
                    next_id += 1

            if expressions:
                self.hs_db = hyperscan.Database()
                self.hs_db.compile(expressions=expressions, ids=ids, flags=flags)
        except Exception:
            self.hs_db = None

    # ------------------------------------------------------------
    # FAST GC via prefix sum (O(1))
    # ------------------------------------------------------------

    def _build_g_prefix(self, seq: str):
        prefix = [0]
        for base in seq:
            prefix.append(prefix[-1] + (1 if base == "G" else 0))
        return prefix

    def _percent_g_fast(self, prefix, start, end):
        length = end - start
        if length == 0:
            return 0.0
        g_count = prefix[end] - prefix[start]
        return (g_count / length) * 100.0

    # ------------------------------------------------------------
    # Seed-Accelerated RIZ Search
    # ------------------------------------------------------------

    def _riz_search(self, seq: str, model: str) -> List[Dict[str, Any]]:

        if HS_AVAILABLE and self.hs_db:
            results = []
            seq_bytes = seq.encode()

            def on_match(id, from_, to, flags, context):
                if self.hs_id_to_model.get(id) == model:
                    results.append((from_, to))
                return 0

            self.hs_db.scan(seq_bytes, match_event_handler=on_match)
        else:
            pattern = re.compile(self.get_patterns()[model][0][0])
            results = [(m.start(), m.end()) for m in pattern.finditer(seq)]

        prefix = self._build_g_prefix(seq)
        riz_list = []

        for start, end in results:
            perc_g = self._percent_g_fast(prefix, start, end)
            if perc_g >= self.MIN_PERC_G_RIZ:
                riz_list.append({
                    "start": start,
                    "end": end,
                    "length": end - start,
                    "perc_g": perc_g,
                    "seq": seq[start:end]
                })

        return riz_list

    # ------------------------------------------------------------
    # Seed-Accelerated REZ Search (Optimized)
    # ------------------------------------------------------------

    def _find_rez(self, seq: str, riz_end: int, prefix) -> Optional[Dict[str, Any]]:

        search_start = riz_end + self.NUM_LINKER
        if search_start >= len(seq):
            return None

        best_rez = None
        max_len = 0

        # seed prefilter: must contain at least one GGG
        for i in range(search_start,
                       min(len(seq), riz_end + self.MAX_LENGTH_REZ),
                       self.WINDOW_STEP):

            if "GGG" not in seq[i:i+100]:
                continue

            for j in range(i + 100,
                           min(len(seq), i + self.MAX_LENGTH_REZ),
                           100):

                perc_g = self._percent_g_fast(prefix, i, j)
                if perc_g >= self.MIN_PERC_G_REZ:
                    length = j - i
                    if length > max_len:
                        max_len = length
                        best_rez = {
                            "start": i,
                            "end": j,
                            "length": length,
                            "perc_g": perc_g,
                            "seq": seq[i:j]
                        }

        return best_rez

    # ------------------------------------------------------------
    # Core Annotation
    # ------------------------------------------------------------

    def annotate_sequence(self, sequence: str,
                          models: Optional[List[str]] = None) -> List[Dict[str, Any]]:

        seq = sequence.upper()
        prefix = self._build_g_prefix(seq)
        models = models or ['qmrlfs_model_1', 'qmrlfs_model_2']
        results = []

        for model in models:
            riz_regions = self._riz_search(seq, model)

            for riz in riz_regions:
                rez = self._find_rez(seq, riz["end"], prefix)

                total_end = rez["end"] if rez else riz["end"]

                results.append({
                    "model": model,
                    "riz_start": riz["start"],
                    "riz_end": riz["end"],
                    "riz_length": riz["length"],
                    "riz_perc_g": riz["perc_g"],
                    "riz_sequence": riz["seq"],
                    "rez_start": rez["start"] if rez else None,
                    "rez_end": rez["end"] if rez else None,
                    "rez_length": rez["length"] if rez else 0,
                    "rez_perc_g": rez["perc_g"] if rez else 0.0,
                    "rez_sequence": rez["seq"] if rez else "",
                    "total_start": riz["start"],
                    "total_end": total_end,
                    "total_length": total_end - riz["start"]
                })

        return results

    # ------------------------------------------------------------
    # Detection Interface (UNCHANGED STRUCTURE)
    # ------------------------------------------------------------

    def detect_motifs(self, sequence: str,
                      sequence_name: str = "sequence"):

        self.audit['invoked'] = True
        self.audit['both_strands_scanned'] = True

        motifs = []

        for strand, seq in [("+", sequence.upper()),
                            ("-", revcomp(sequence.upper()))]:

            annotations = self.annotate_sequence(seq)

            seq_len = len(sequence)

            for i, ann in enumerate(annotations):

                if strand == "-":
                    start = seq_len - ann["total_end"]
                    end = seq_len - ann["total_start"]
                else:
                    start = ann["total_start"]
                    end = ann["total_end"]

                canonical_class, canonical_subclass = normalize_class_subclass(
                    self.get_motif_class_name(),
                    "R-loop formation sites",
                    strict=False,
                    auto_correct=True
                )

                motifs.append({
                    "ID": f"{sequence_name}_RLOOP_{start+1}",
                    "Sequence_Name": sequence_name,
                    "Class": canonical_class,
                    "Subclass": canonical_subclass,
                    "Start": start + 1,
                    "End": end,
                    "Length": end - start,
                    "Sequence": sequence[start:end],
                    "Score": round(
                        (ann["riz_perc_g"]/100) +
                        (ann["rez_perc_g"]/100),
                        3
                    ),
                    "Strand": strand,
                    "Method": "QmRLFS_detection",
                    "Pattern_ID": f"QmRLFS_{i+1}"
                })

        return motifs
