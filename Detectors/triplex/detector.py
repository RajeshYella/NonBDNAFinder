"""Triplex DNA detector: mirror repeats (10-100 nt arms) and sticky DNA (Frank-Kamenetskii 1995, Sakamoto 1999)."""

import re
from typing import List, Dict, Any, Tuple
from ..base.base_detector import BaseMotifDetector
from Utilities.core.motif_normalizer import normalize_class_subclass

try: from scanner import find_mirror_repeats_optimized as _find_mirror_repeats_optimized
except ImportError: _find_mirror_repeats_optimized = None

try: from motif_patterns import TRIPLEX_PATTERNS
except ImportError: TRIPLEX_PATTERNS = {}

def revcomp(seq: str) -> str:
    comp = str.maketrans('ATGCRYSWKMBDHVN', 'TACGYRSWMKVHDBN'); return seq.translate(comp)[::-1]


class TriplexDetector(BaseMotifDetector):
    """Triplex DNA detector: mirror repeats (10-100 nt arms, ≤8 nt spacer, 90% purine/pyrimidine) and sticky DNA."""

    # Tunable parameters
    MIN_ARM = 10; MAX_ARM = 100; MAX_LOOP = 8; PURINE_PYRIMIDINE_THRESHOLD = 0.9

    def get_motif_class_name(self) -> str: return "Triplex"

    def get_patterns(self) -> Dict[str, List[Tuple]]:
        """Return sticky DNA patterns; mirror repeats use optimized k-mer scanner."""
        return self._load_patterns(TRIPLEX_PATTERNS, lambda: {
            'triplex_forming_sequences': [
                (r'(?:GAA){4,}', 'TRX_5_4', 'GAA repeat', 'Sticky_DNA', 12, 'sticky_dna_score', 0.95, 'Disease-associated repeats', 'Sakamoto 1999'),
                (r'(?:TTC){4,}', 'TRX_5_5', 'TTC repeat', 'Sticky_DNA', 12, 'sticky_dna_score', 0.95, 'Disease-associated repeats', 'Sakamoto 1999'),
            ]
        })
    
    def annotate_sequence(self, sequence: str) -> List[Dict[str, Any]]:
        seq = sequence.upper(); results = []; used = [False] * len(seq)
        patterns = self.get_patterns()['triplex_forming_sequences']

        if _find_mirror_repeats_optimized is not None:
            mirror_results = _find_mirror_repeats_optimized(seq, min_arm=self.MIN_ARM, max_arm=self.MAX_ARM,
                                                            max_loop=self.MAX_LOOP, purine_pyrimidine_threshold=self.PURINE_PYRIMIDINE_THRESHOLD)
            for mr_rec in mirror_results:
                if mr_rec.get('Is_Triplex', False):
                    s, e = mr_rec['Start'] - 1, mr_rec['End']
                    if any(used[s:e]): continue
                    for i in range(s, e): used[i] = True
                    pur_frac, pyr_frac = mr_rec['Purine_Fraction'], mr_rec['Pyrimidine_Fraction']
                    subtype = 'Homopurine mirror repeat' if pur_frac >= self.PURINE_PYRIMIDINE_THRESHOLD else 'Homopyrimidine mirror repeat'
                    pid = 'TRX_MR_PU' if pur_frac >= self.PURINE_PYRIMIDINE_THRESHOLD else 'TRX_MR_PY'
                    results.append({'class_name': 'Triplex', 'pattern_id': pid, 'start': s, 'end': e, 'length': e - s,
                                    'score': self._triplex_potential(seq[s:e]), 'matched_seq': seq[s:e],
                                    'details': {'type': subtype, 'reference': 'Frank-Kamenetskii 1995', 'description': 'H-DNA formation',
                                                'arm_length': mr_rec['Arm_Length'], 'loop_length': mr_rec['Loop'],
                                                'left_arm': mr_rec.get('Left_Arm', ''), 'right_arm': mr_rec.get('Right_Arm', ''),
                                                'loop_seq': mr_rec.get('Loop_Seq', ''), 'purine_fraction': pur_frac, 'pyrimidine_fraction': pyr_frac}})
        else:
            # Fallback regex for homopurine
            pat_pu = rf'([GA]{{{self.MIN_ARM},{self.MAX_ARM}}})([ATGC]{{0,{self.MAX_LOOP}}})([GA]{{{self.MIN_ARM},{self.MAX_ARM}}})'
            for m in re.finditer(pat_pu, seq):
                s, e = m.span()
                if any(used[s:e]): continue
                arm1, arm2, loop = m.group(1), m.group(3), m.group(2)
                if len(arm1) < self.MIN_ARM or len(arm1) > self.MAX_ARM or len(arm2) < self.MIN_ARM or len(arm2) > self.MAX_ARM: continue
                if len(loop) > self.MAX_LOOP: continue
                pur_ct = sum(1 for b in arm1+arm2 if b in 'AG') / max(1, len(arm1+arm2))
                if pur_ct < self.PURINE_PYRIMIDINE_THRESHOLD: continue
                for i in range(s, e): used[i] = True
                results.append({'class_name': 'Triplex', 'pattern_id': 'TRX_MR_PU', 'start': s, 'end': e, 'length': e-s,
                                'score': self._triplex_potential(seq[s:e]), 'matched_seq': seq[s:e],
                                'details': {'type': 'Homopurine mirror repeat', 'reference': 'Frank-Kamenetskii 1995', 'description': 'H-DNA formation (homopurine)'}})
            # Fallback regex for homopyrimidine
            pat_py = rf'([CT]{{{self.MIN_ARM},{self.MAX_ARM}}})([ATGC]{{0,{self.MAX_LOOP}}})([CT]{{{self.MIN_ARM},{self.MAX_ARM}}})'
            for m in re.finditer(pat_py, seq):
                s, e = m.span()
                if any(used[s:e]): continue
                arm1, arm2, loop = m.group(1), m.group(3), m.group(2)
                if len(arm1) < self.MIN_ARM or len(arm1) > self.MAX_ARM or len(arm2) < self.MIN_ARM or len(arm2) > self.MAX_ARM: continue
                if len(loop) > self.MAX_LOOP: continue
                pyr_ct = sum(1 for b in arm1+arm2 if b in 'CT') / max(1, len(arm1+arm2))
                if pyr_ct < self.PURINE_PYRIMIDINE_THRESHOLD: continue
                for i in range(s, e): used[i] = True
                results.append({'class_name': 'Triplex', 'pattern_id': 'TRX_MR_PY', 'start': s, 'end': e, 'length': e-s,
                                'score': self._triplex_potential(seq[s:e]), 'matched_seq': seq[s:e],
                                'details': {'type': 'Homopyrimidine mirror repeat', 'reference': 'Frank-Kamenetskii 1995', 'description': 'H-DNA formation (homopyrimidine)'}})

        # Sticky DNA patterns (GAA/TTC)
        for patinfo in patterns:
            pat, pid, name, cname, minlen, scoretype, cutoff, desc, ref = patinfo
            for m in re.finditer(pat, seq):
                s, e = m.span()
                if any(used[s:e]): continue
                for i in range(s, e): used[i] = True
                results.append({'class_name': cname, 'pattern_id': pid, 'start': s, 'end': e, 'length': e-s,
                                'score': self.calculate_score(seq[s:e], patinfo), 'matched_seq': seq[s:e],
                                'details': {'type': name, 'reference': ref, 'description': desc}})
        results.sort(key=lambda r: r['start'])
        return results

    def calculate_score(self, sequence: str, pattern_info: Tuple) -> float:
        scoring_method = pattern_info[5] if len(pattern_info) > 5 else 'triplex_potential'
        return self._triplex_potential(sequence) if scoring_method == 'triplex_potential' else self._sticky_dna_score(sequence) if scoring_method == 'sticky_dna_score' else 0.0

    def _triplex_potential(self, sequence: str) -> float:
        """Score: tract length and purine/pyrimidine content (≥90%)."""
        if len(sequence) < 20: return 0.0
        pur = sum(1 for b in sequence if b in "AG") / len(sequence); pyr = sum(1 for b in sequence if b in "CT") / len(sequence)
        score = (pur if pur > 0.9 else 0) + (pyr if pyr > 0.9 else 0)
        return min(score * len(sequence) / 150, 1.0)

    def _sticky_dna_score(self, sequence: str) -> float:
        """Score for sticky DNA: repeat density and length."""
        if len(sequence) < 12: return 0.0
        gaa_count, ttc_count = sequence.count("GAA"), sequence.count("TTC"); rep_total = gaa_count + ttc_count
        density = (rep_total * 3) / len(sequence)
        extras = sum(len(m.group(0)) for m in re.finditer(r'(?:GAA){2,}', sequence)) + sum(len(m.group(0)) for m in re.finditer(r'(?:TTC){2,}', sequence))
        cons_bonus = extras / len(sequence) if len(sequence) else 0
        return min(0.7 * density + 0.3 * cons_bonus, 1.0)

    def passes_quality_threshold(self, sequence: str, score: float, pattern_info: Tuple) -> bool: return score >= 0.2

    def detect_motifs(self, sequence: str, sequence_name: str = "sequence") -> List[Dict[str, Any]]:
        """Detect triplex motifs on BOTH strands (strand-agnostic)."""
        self.audit['invoked'] = True; self.audit['windows_scanned'] = 0; self.audit['candidates_seen'] = 0
        self.audit['candidates_filtered'] = 0; self.audit['reported'] = 0; self.audit['both_strands_scanned'] = True
        sequence = sequence.upper().strip(); motifs = []

        self.audit['windows_scanned'] += 1; results_fwd = self.annotate_sequence(sequence); self.audit['candidates_seen'] += len(results_fwd)
        for i, result in enumerate(results_fwd):
            subtype = result['details']['type']
            canonical_subclass = 'Sticky DNA' if 'Sticky' in subtype or 'GAA' in subtype or 'TTC' in subtype else 'Triplex'
            canonical_class, canonical_subclass = normalize_class_subclass(self.get_motif_class_name(), canonical_subclass, strict=False, auto_correct=True)
            motif_dict = {'ID': f"{sequence_name}_{result['pattern_id']}_{result['start']+1}", 'Sequence_Name': sequence_name,
                          'Class': canonical_class, 'Subclass': canonical_subclass, 'Start': result['start'] + 1, 'End': result['end'],
                          'Length': result['length'], 'Sequence': result['matched_seq'], 'Score': round(result['score'], 3),
                          'Strand': '+', 'Method': 'Triplex_detection', 'Pattern_ID': result['pattern_id']}
            if 'mirror repeat' in result['details']['type'].lower():
                details = result['details']
                if 'left_arm' in details and details['left_arm']:
                    motif_dict.update({'Left_Arm': details['left_arm'], 'Right_Arm': details['right_arm'], 'Loop_Seq': details['loop_seq'],
                                       'Arm_Length': details.get('arm_length', len(details['left_arm'])), 'Loop_Length': details.get('loop_length', len(details['loop_seq'])),
                                       'GC_Left_Arm': round(self._calc_gc(details['left_arm']), 2), 'GC_Right_Arm': round(self._calc_gc(details['right_arm']), 2),
                                       'GC_Loop': round(self._calc_gc(details['loop_seq']), 2)})
                elif details.get('arm_length', 0) > 0 and details.get('loop_length', 0) > 0:
                    arm_len, loop_len = details['arm_length'], details['loop_length']; matched_seq = result['matched_seq']
                    left_arm, loop_seq, right_arm = matched_seq[:arm_len], matched_seq[arm_len:arm_len+loop_len], matched_seq[arm_len+loop_len:arm_len+loop_len+arm_len]
                    motif_dict.update({'Left_Arm': left_arm, 'Right_Arm': right_arm, 'Loop_Seq': loop_seq, 'Arm_Length': arm_len, 'Loop_Length': loop_len,
                                       'GC_Left_Arm': round(self._calc_gc(left_arm), 2), 'GC_Right_Arm': round(self._calc_gc(right_arm), 2), 'GC_Loop': round(self._calc_gc(loop_seq), 2)})
            motifs.append(motif_dict); self.audit['reported'] += 1

        self.audit['windows_scanned'] += 1; seq_rc = revcomp(sequence); results_rc = self.annotate_sequence(seq_rc); self.audit['candidates_seen'] += len(results_rc)
        seq_len = len(sequence)
        for i, result in enumerate(results_rc):
            fwd_start, fwd_end = seq_len - result['end'], seq_len - result['start']
            subtype = result['details']['type']
            canonical_subclass = 'Sticky DNA' if 'Sticky' in subtype or 'GAA' in subtype or 'TTC' in subtype else 'Triplex'
            canonical_class, canonical_subclass = normalize_class_subclass(self.get_motif_class_name(), canonical_subclass, strict=False, auto_correct=True)
            motif_dict = {'ID': f"{sequence_name}_{result['pattern_id']}_RC_{fwd_start+1}", 'Sequence_Name': sequence_name,
                          'Class': canonical_class, 'Subclass': canonical_subclass, 'Start': fwd_start + 1, 'End': fwd_end,
                          'Length': result['length'], 'Sequence': sequence[fwd_start:fwd_end], 'Score': round(result['score'], 3),
                          'Strand': '-', 'Method': 'Triplex_detection', 'Pattern_ID': result['pattern_id'] + '_RC'}
            if 'mirror repeat' in result['details']['type'].lower():
                details = result['details']
                if 'left_arm' in details and details['left_arm']:
                    motif_dict.update({'Left_Arm': details['left_arm'], 'Right_Arm': details['right_arm'], 'Loop_Seq': details['loop_seq'],
                                       'Arm_Length': details.get('arm_length', len(details['left_arm'])), 'Loop_Length': details.get('loop_length', len(details['loop_seq'])),
                                       'GC_Left_Arm': round(self._calc_gc(details['left_arm']), 2), 'GC_Right_Arm': round(self._calc_gc(details['right_arm']), 2),
                                       'GC_Loop': round(self._calc_gc(details['loop_seq']), 2)})
                elif details.get('arm_length', 0) > 0 and details.get('loop_length', 0) > 0:
                    arm_len, loop_len = details['arm_length'], details['loop_length']; matched_seq = result['matched_seq']
                    left_arm, loop_seq, right_arm = matched_seq[:arm_len], matched_seq[arm_len:arm_len+loop_len], matched_seq[arm_len+loop_len:arm_len+loop_len+arm_len]
                    motif_dict.update({'Left_Arm': left_arm, 'Right_Arm': right_arm, 'Loop_Seq': loop_seq, 'Arm_Length': arm_len, 'Loop_Length': loop_len,
                                       'GC_Left_Arm': round(self._calc_gc(left_arm), 2), 'GC_Right_Arm': round(self._calc_gc(right_arm), 2), 'GC_Loop': round(self._calc_gc(loop_seq), 2)})
            motifs.append(motif_dict); self.audit['reported'] += 1
        return motifs
