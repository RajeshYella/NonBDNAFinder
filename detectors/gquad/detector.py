"""G-Quadruplex detector: G4Hunter scoring with overlap resolution (Huppert 2005, Bedrat 2016)."""

import re
from typing import Dict, List, Tuple, Any
from ..base.base_detector import BaseMotifDetector
from core.motif_normalizer import normalize_class_subclass

# Tunable constants
WINDOW_SIZE_DEFAULT = 25; MIN_REGION_LEN = 8
CLASS_PRIORITY = ["telomeric_g4", "stacked_canonical_g4s", "stacked_g4s_linker", "canonical_g4",
                  "extended_loop_g4", "higher_order_g4", "g_triplex", "weak_pqs"]

class GQuadruplexDetector(BaseMotifDetector):
    """G-Quadruplex detector: G4Hunter scoring with overlap resolution."""

    def get_motif_class_name(self) -> str: return "G-Quadruplex"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return G4 motif patterns with G4Hunter scoring."""
        return {
            'telomeric_g4': [(r'(?:TTAGGG){4,}', 'G4_TEL', 'Telomeric G4', 'Telomeric G4', 24, 'g4hunter_score', 0.95, 'Human telomeric G4', 'Parkinson 2002')],
            'stacked_canonical_g4s': [(r'(?:(?:G{3,}[ACGT]{1,7}){3}G{3,}){2,}', 'G4_STK_CAN', 'Stacked canonical G4s', 'Stacked canonical G4s', 30, 'g4hunter_score', 0.90, 'Stacking polymorphism', 'Phan 2007')],
            'stacked_g4s_linker': [(r'(?:(?:G{3,}[ACGT]{1,7}){3}G{3,})(?:[ACGT]{0,20}(?:(?:G{3,}[ACGT]{1,7}){3}G{3,})){1,}', 'G4_STK_LNK', 'Stacked G4s with linker', 'Stacked G4s with linker', 30, 'g4hunter_score', 0.85, 'Clustered G4s', 'Hänsel-Hertsch 2016')],
            'canonical_g4': [(r'G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}', 'G4_CAN', 'Canonical intramolecular G4', 'Canonical intramolecular G4', 15, 'g4hunter_score', 0.80, 'Canonical G4', 'Huppert 2005')],
            'extended_loop_g4': [(r'G{3,}[ACGT]{1,12}G{3,}[ACGT]{1,12}G{3,}[ACGT]{1,12}G{3,}', 'G4_EXT', 'Extended-loop canonical', 'Extended-loop canonical', 15, 'g4hunter_score', 0.70, 'Long-loop G4s', 'Chambers 2015')],
            'higher_order_g4': [(r'(?:G{3,}[ACGT]{1,7}){7,}', 'G4_HIGH', 'Higher-order G4 array/G4-wire', 'Higher-order G4 array/G4-wire', 49, 'g4hunter_score', 0.65, 'Higher-order G4s', 'Wong 2003')],
            'g_triplex': [(r'G{3,}[ACGT]{1,7}G{3,}[ACGT]{1,7}G{3,}', 'G4_TRX', 'Intramolecular G-triplex', 'Intramolecular G-triplex', 12, 'g_triplex_score', 0.45, 'G-triplex structures', 'Lim 2013')],
            'weak_pqs': [(r'G{2,}[ACGT]{1,7}G{2,}[ACGT]{1,7}G{2,}[ACGT]{1,7}G{2,}', 'G4_WEAK', 'Two-tetrad weak PQS', 'Two-tetrad weak PQS', 11, 'g4hunter_score', 0.25, 'QGRS Mapper weak PQS', 'Kikin 2006')],
        }

    def calculate_score(self, sequence: str, pattern_info: Tuple = None) -> float:
        """Compute total score for all accepted G4 regions after overlap resolution."""
        seq = sequence.upper(); candidates = self._find_all_candidates(seq)
        scored = [self._score_candidate(c, seq) for c in candidates]; accepted = self._resolve_overlaps(scored)
        return float(sum(a['score'] for a in accepted))

    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        """Annotate all accepted G4 regions after overlap resolution."""
        seq = sequence.upper(); candidates = self._find_all_candidates(seq)
        scored = [self._score_candidate(c, seq) for c in candidates]; accepted = self._resolve_overlaps(scored)
        return [{'class_name': a['class_name'], 'pattern_id': a['pattern_id'], 'start': a['start'], 'end': a['end'],
                 'length': a['end'] - a['start'], 'score': round(a['score'], 6), 'matched_seq': seq[a['start']:a['end']],
                 'details': a['details']} for a in accepted]

    def _find_all_candidates(self, seq: str) -> List[Dict[str, Any]]:
        """Find all regions matching any G4 motif."""
        patt_groups = self.get_patterns(); candidates = []
        for class_name, patterns in patt_groups.items():
            for pat in patterns:
                regex = pat[0]; pattern_id = pat[1] if len(pat) > 1 else f"{class_name}_pat"
                for m in re.finditer(regex, seq):
                    s, e = m.start(), m.end()
                    if (e - s) >= MIN_REGION_LEN:
                        candidates.append({'class_name': class_name, 'pattern_id': pattern_id, 'start': s, 'end': e, 'match_text': seq[s:e]})
        return candidates

    def _score_candidate(self, candidate: Dict[str, Any], seq: str, window_size: int = WINDOW_SIZE_DEFAULT) -> Dict[str, Any]:
        """Calculate per-region G4Hunter score plus tract/GC penalties."""
        s, e = candidate['start'], candidate['end']; region = seq[s:e]; L = len(region)
        vals = [1 if ch == 'G' else -1 if ch == 'C' else 0 for ch in region]; ws = min(window_size, L); max_abs = 0
        if ws > 0:
            cur = sum(vals[0:ws]); max_abs = abs(cur)
            for i in range(1, L - ws + 1): cur += vals[i + ws - 1] - vals[i - 1]; max_abs = max(max_abs, abs(cur))
        normalized_window = (max_abs / ws) if ws > 0 else 0.0
        g_tracts = re.findall(r'G{2,}', region); n_g = len(g_tracts); total_g_len = sum(len(t) for t in g_tracts)
        tract_bonus = min(0.5, 0.08 * (n_g - 2) * ((total_g_len / n_g) / 4.0)) if n_g >= 3 else 0.0
        total_c, total_g = region.count('C'), region.count('G'); gc_balance = (total_g - total_c) / (L if L > 0 else 1)
        gc_penalty = 0.2 if gc_balance < -0.3 else 0.1 if gc_balance < -0.1 else 0.0
        normalized_score = max(0.0, min(1.0, normalized_window + tract_bonus - gc_penalty))
        region_score = normalized_score * (L / float(ws)) if ws > 0 else 0.0
        stems, loops = self._extract_g4_components(region)
        gc_total = self._calc_gc(region); gc_stems = self._calc_gc(''.join(stems)) if stems else 0
        details = {'n_g_tracts': n_g, 'total_g_len': total_g_len, 'gc_balance': round(gc_balance, 4),
                   'max_window_abs': float(max_abs), 'normalized_window': round(normalized_window, 6),
                   'tract_bonus': round(tract_bonus, 6), 'gc_penalty': round(gc_penalty, 6),
                   'normalized_score': round(normalized_score, 6), 'region_score': round(region_score, 6),
                   'stems': stems, 'loops': loops, 'num_stems': len(stems), 'num_loops': len(loops),
                   'stem_lengths': [len(s) for s in stems], 'loop_lengths': [len(l) for l in loops],
                   'GC_Total': round(gc_total, 2), 'GC_Stems': round(gc_stems, 2)}
        out = candidate.copy(); out['score'] = float(region_score); out['details'] = details
        return out

    def _extract_g4_components(self, sequence: str) -> Tuple[List[str], List[str]]:
        """Extract stems (G-tracts) and loops from G-quadruplex sequence."""
        matches = list(re.finditer(r'G{2,}', sequence)); stems = []; loops = []
        if len(matches) >= 2:
            for i, match in enumerate(matches):
                stems.append(match.group())
                if i < len(matches) - 1:
                    loop_start, loop_end = match.end(), matches[i + 1].start()
                    if loop_end > loop_start: loops.append(sequence[loop_start:loop_end])
        return stems, loops

    def _resolve_overlaps(self, scored_candidates: List[Dict[str, Any]], merge_gap: int = 0) -> List[Dict[str, Any]]:
        """Select non-overlapping candidates by score and class priority."""
        if not scored_candidates: return []
        def class_prio_idx(class_name):
            try: return CLASS_PRIORITY.index(class_name)
            except ValueError: return len(CLASS_PRIORITY)
        scored_sorted = sorted(scored_candidates, key=lambda x: (-x['score'], class_prio_idx(x['class_name']), -(x['end'] - x['start'])))
        accepted = []; occupied = []
        for cand in scored_sorted:
            s, e = cand['start'], cand['end']; conflict = any(not (e <= as_ - merge_gap or s >= ae + merge_gap) for (as_, ae) in occupied)
            if not conflict: accepted.append(cand); occupied.append((s, e))
        accepted.sort(key=lambda x: x['start'])
        return accepted

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Detect G4 motifs with overlap resolution and component info."""
        sequence = sequence.upper().strip(); motifs = []; annotations = self.annotate_sequence(sequence)
        subclass_map = {'telomeric_g4': 'Telomeric G4', 'stacked_canonical_g4s': 'Stacked canonical G4s',
                        'stacked_g4s_linker': 'Stacked G4s with linker', 'canonical_g4': 'Canonical intramolecular G4',
                        'extended_loop_g4': 'Extended-loop canonical', 'higher_order_g4': 'Higher-order G4 array/G4-wire',
                        'g_triplex': 'Intramolecular G-triplex', 'weak_pqs': 'Two-tetrad weak PQS'}
        for annotation in annotations:
            class_name = annotation['class_name']; subclass = subclass_map.get(class_name, class_name)
            canonical_class, canonical_subclass = normalize_class_subclass(self.get_motif_class_name(), subclass, strict=False, auto_correct=True)
            start_pos, end_pos = annotation['start'], annotation['end']; details = annotation.get('details', {})
            motif = {'ID': f"{sequence_name}_{annotation['pattern_id']}_{start_pos+1}", 'Sequence_Name': sequence_name,
                     'Class': canonical_class, 'Subclass': canonical_subclass, 'Start': start_pos + 1, 'End': end_pos,
                     'Length': annotation['length'], 'Sequence': annotation['matched_seq'], 'Score': round(annotation['score'], 3),
                     'Strand': '+', 'Method': 'G4Hunter_detection', 'Pattern_ID': annotation['pattern_id']}
            if details:
                motif.update({'Stems': details.get('stems', []), 'Loops': details.get('loops', []),
                              'Num_Stems': details.get('num_stems', 0), 'Num_Loops': details.get('num_loops', 0),
                              'Stem_Lengths': details.get('stem_lengths', []), 'Loop_Lengths': details.get('loop_lengths', []),
                              'GC_Total': details.get('GC_Total', 0), 'GC_Stems': details.get('GC_Stems', 0)})
                stem_lengths, loop_lengths = details.get('stem_lengths', []), details.get('loop_lengths', [])
                if stem_lengths: motif['Stem_Length'] = sum(stem_lengths) / len(stem_lengths)
                if loop_lengths: motif['Loop_Length'] = sum(loop_lengths) / len(loop_lengths)
            motifs.append(motif)
        return motifs

