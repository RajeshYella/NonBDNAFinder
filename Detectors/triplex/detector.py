"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Triplex DNA Detector - Seed-and-extend mirror repeats (H-DNA) & Sticky DNA  │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.3           │
│ References: Frank-Kamenetskii 1995; Soyfer & Potaman 1995; Sakamoto 1999    │
│ Subclasses: Triplex (mirror repeats), Sticky DNA (GAA/TTC disease repeats)  │
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

# Sticky DNA specific constants (GAA/TTC repeats - disease-associated)
MIN_STICKY_COPIES = 4  # Minimum 4 copies for Sticky DNA
MIN_STICKY_LENGTH = 12  # Minimum tract length (4 copies * 3 bp)


class TriplexDetector(BaseMotifDetector):

    MIN_ARM = MIN_ARM
    MAX_ARM = MAX_ARM
    MAX_LOOP = MAX_LOOP
    PURITY_THRESHOLD = PURITY_THRESHOLD
    SEED_SIZE = SEED_SIZE
    SCORE_THRESHOLD = SCORE_THRESHOLD
    MIN_STICKY_COPIES = MIN_STICKY_COPIES
    MIN_STICKY_LENGTH = MIN_STICKY_LENGTH

    # Sticky DNA patterns (GAA/TTC repeats associated with Friedreich ataxia)
    STICKY_PATTERNS = [
        (re.compile(r'(?:GAA){4,}', re.IGNORECASE), 'TRX_STICKY_GAA', 'GAA repeat'),
        (re.compile(r'(?:TTC){4,}', re.IGNORECASE), 'TRX_STICKY_TTC', 'TTC repeat'),
    ]

    # ------------------------------------------------------------------ #
    # Required Interface
    # ------------------------------------------------------------------ #

    def get_motif_class_name(self) -> str:
        return "Triplex"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        return {
            "mirror_triplex": [
                (r"", "TRX_MIRROR", "Mirror repeat triplex",
                 "Triplex", self.MIN_ARM,
                 "structural_triplex_score", self.SCORE_THRESHOLD,
                 "Seed-based mirror repeat detection",
                 "Frank-Kamenetskii 1995")
            ],
            "sticky_dna": [
                (r'(?:GAA){4,}', 'TRX_STICKY_GAA', 'GAA repeat',
                 'Sticky DNA', 12, 'sticky_dna_score', 0.95,
                 'Disease-associated repeats', 'Sakamoto 1999'),
                (r'(?:TTC){4,}', 'TRX_STICKY_TTC', 'TTC repeat',
                 'Sticky DNA', 12, 'sticky_dna_score', 0.95,
                 'Disease-associated repeats', 'Sakamoto 1999'),
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
    # Sticky DNA Detection (GAA/TTC repeats)
    # ------------------------------------------------------------------ #

    def _find_sticky_dna(self, sequence: str) -> List[Dict[str, Any]]:
        """
        Find Sticky DNA regions (GAA/TTC repeats).
        
        GAA and TTC repeats are associated with Friedreich ataxia and form 
        unusual sticky DNA structures. Reference: Sakamoto 1999.
        """
        seq = sequence.upper()
        hits = []

        for pattern, pattern_id, description in self.STICKY_PATTERNS:
            for match in pattern.finditer(seq):
                start, end = match.start(), match.end()
                tract_length = end - start
                
                # Calculate copy number
                if 'GAA' in pattern_id:
                    unit = 'GAA'
                else:
                    unit = 'TTC'
                copy_number = tract_length // 3
                
                # Calculate score based on copy number (more copies = higher pathogenic potential)
                # Scoring formula: score = 1.0 + (copies - 4) * 0.2, capped at 3.0
                # - Base score 1.0 for minimum 4 copies (detection threshold)
                # - Increment 0.2 per additional copy to reflect increasing instability
                # - Cap at 3.0 (maximum pathogenic potential)
                # Score interpretation by copy number:
                #   - Weak (1.0-1.5):   4-6 copies
                #   - Moderate (1.6-2.5): 7-11 copies
                #   - Strong (2.6-3.0): 12+ copies
                # Note: Pathogenic threshold for Friedreich ataxia is typically >66 copies
                score = min(3.0, 1.0 + (copy_number - 4) * 0.2)
                
                hits.append({
                    "start": start,
                    "end": end,
                    "length": tract_length,
                    "score": score,
                    "sequence": seq[start:end],
                    "pattern_id": pattern_id,
                    "subclass": "Sticky DNA",
                    "repeat_unit": unit,
                    "copy_number": copy_number,
                    "description": description
                })

        # Sort by length (longer tracts more significant)
        hits.sort(key=lambda x: (-x["length"], x["start"]))
        return hits

    # ------------------------------------------------------------------ #
    # Annotation with overlap control
    # ------------------------------------------------------------------ #

    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        """
        Annotate sequence with both mirror triplex and sticky DNA motifs.
        
        Note: Sticky DNA and mirror triplex motifs are allowed to overlap since
        they represent different structural features. However, overlaps within
        the same subclass are removed.
        """
        seq = sequence.upper()
        results = []
        
        # Detect mirror repeats (Triplex subclass)
        mirror_used = [False] * len(seq)
        mirrors = self._find_mirror_repeats(seq)

        for m in mirrors:
            s = m["start"]
            e = m["end"]

            if any(mirror_used[s:e]):
                continue

            for i in range(s, e):
                mirror_used[i] = True

            results.append({
                "class_name": "Triplex",
                "subclass": "Triplex",
                "pattern_id": "TRX_MIRROR",
                "start": s,
                "end": e,
                "length": e - s,
                "score": m["score"],
                "matched_seq": m["sequence"]
            })

        # Detect sticky DNA (Sticky DNA subclass) - separate overlap tracking
        sticky_used = [False] * len(seq)
        sticky_hits = self._find_sticky_dna(seq)
        
        for hit in sticky_hits:
            s = hit["start"]
            e = hit["end"]
            
            if any(sticky_used[s:e]):
                continue
            
            for i in range(s, e):
                sticky_used[i] = True
            
            results.append({
                "class_name": "Triplex",
                "subclass": "Sticky DNA",
                "pattern_id": hit["pattern_id"],
                "start": s,
                "end": e,
                "length": hit["length"],
                "score": hit["score"],
                "matched_seq": hit["sequence"],
                "repeat_unit": hit.get("repeat_unit"),
                "copy_number": hit.get("copy_number")
            })

        results.sort(key=lambda r: r["start"])
        return results

    # ------------------------------------------------------------------ #
    # Final Detection Interface
    # ------------------------------------------------------------------ #

    def detect_motifs(self,
                      sequence: str,
                      sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """
        Detect triplex DNA motifs including mirror repeats and sticky DNA.
        
        Returns motifs with proper canonical subclasses:
        - "Triplex" for H-DNA forming mirror repeats
        - "Sticky DNA" for GAA/TTC disease-associated repeats
        
        Overlaps between different subclasses are allowed (e.g., a region can be
        both a Triplex and overlap with Slipped DNA since they represent different
        structural properties).
        """
        sequence = sequence.upper().strip()
        motifs = []
        annotations = self.annotate_sequence(sequence)

        for annotation in annotations:
            # Use the subclass from annotation (Triplex or Sticky DNA)
            subclass = annotation.get("subclass", "Triplex")
            
            canonical_class, canonical_subclass = normalize_class_subclass(
                self.get_motif_class_name(),
                subclass,
                strict=False,
                auto_correct=True
            )
            
            # Determine method and ID prefix based on subclass
            if canonical_subclass == "Sticky DNA":
                method = "Sticky_DNA_detection"
                id_prefix = "STICKY"
            else:
                method = "Triplex_seed_mirror_detection"
                id_prefix = "TRX"

            motif = {
                "ID": f"{sequence_name}_{id_prefix}_{annotation['start']+1}",
                "Sequence_Name": sequence_name,
                "Class": canonical_class,
                "Subclass": canonical_subclass,
                "Start": annotation["start"] + 1,
                "End": annotation["end"],
                "Length": annotation["length"],
                "Sequence": annotation["matched_seq"],
                "Score": annotation.get("score", 1.0),
                "Strand": "+",
                "Method": method,
                "Pattern_ID": annotation["pattern_id"]
            }
            
            # Add extra fields for Sticky DNA
            if canonical_subclass == "Sticky DNA":
                if annotation.get("repeat_unit"):
                    motif["Repeat_Unit"] = annotation["repeat_unit"]
                if annotation.get("copy_number"):
                    motif["Copy_Number"] = annotation["copy_number"]

            motifs.append(motif)

        return motifs
