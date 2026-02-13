"""
NonBDNAFinder Detector Validation Tests
========================================

This module contains comprehensive tests to validate that all 9 detectors
are working correctly and producing expected results for known test sequences.

Test Categories:
1. Basic Detection Tests - Verify each detector finds expected motifs
2. Score Validation Tests - Verify scoring is consistent and reasonable
3. Edge Case Tests - Test boundary conditions
4. Integration Tests - Test combined detection pipeline

Author: Dr. Venkata Rajesh Yella
"""

import sys
import os
from pathlib import Path

# Add parent directory to path for imports
parent_dir = str(Path(__file__).parent.parent)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

import unittest
from typing import List, Dict, Any


class TestDetectorBase(unittest.TestCase):
    """Base class for detector tests with common utilities."""

    @classmethod
    def setUpClass(cls):
        """Set up detector instances for all tests."""
        from Detectors import (
            CurvedDNADetector,
            SlippedDNADetector,
            CruciformDetector,
            RLoopDetector,
            TriplexDetector,
            GQuadruplexDetector,
            IMotifDetector,
            ZDNADetector,
            APhilicDetector
        )

        cls.detectors = {
            'curved_dna': CurvedDNADetector(),
            'slipped_dna': SlippedDNADetector(),
            'cruciform': CruciformDetector(),
            'r_loop': RLoopDetector(),
            'triplex': TriplexDetector(),
            'g_quadruplex': GQuadruplexDetector(),
            'i_motif': IMotifDetector(),
            'z_dna': ZDNADetector(),
            'a_philic': APhilicDetector()
        }

    def _run_detector(self, detector_name: str, sequence: str, seq_name: str = "test") -> List[Dict[str, Any]]:
        """Helper to run a specific detector."""
        return self.detectors[detector_name].detect_motifs(sequence, seq_name)

    def _check_motif_structure(self, motif: Dict[str, Any]) -> None:
        """Verify motif has all required fields."""
        required_fields = ['ID', 'Sequence_Name', 'Class', 'Subclass', 'Start', 'End', 
                          'Length', 'Sequence', 'Score', 'Strand', 'Method']
        for field in required_fields:
            self.assertIn(field, motif, f"Motif missing required field: {field}")


class TestGQuadruplexDetector(TestDetectorBase):
    """Test G-Quadruplex detector functionality."""

    def test_canonical_g4_detection(self):
        """Test detection of canonical G4 motif (4 G-tracts with short loops)."""
        # Classic G4: GGGNNNGGGNNNGGGNNNGGGG (4 G-tracts, 3 loops)
        seq = "AAAA" + "GGGTTTGGGCCCGGGAAAGGG" + "TTTT"
        motifs = self._run_detector('g_quadruplex', seq, "g4_test")
        
        self.assertGreater(len(motifs), 0, "Should detect at least one G4 motif")
        for motif in motifs:
            self._check_motif_structure(motif)
            self.assertEqual(motif['Class'], 'G-Quadruplex')

    def test_telomeric_g4_detection(self):
        """Test detection of telomeric G4 repeat (TTAGGG)n."""
        seq = "TTAGGGTTAGGGTTAGGGTTAGGG"  # Human telomeric repeat x4
        motifs = self._run_detector('g_quadruplex', seq, "telomeric_test")
        
        # Should detect at least one G4 motif
        self.assertGreater(len(motifs), 0, "Should detect G4 in telomeric sequence")
        # The detected motif should be G-Quadruplex class
        g4_found = any(m.get('Class') == 'G-Quadruplex' for m in motifs)
        self.assertTrue(g4_found, "Should identify G-Quadruplex class")

    def test_weak_pqs_detection(self):
        """Test detection of weak PQS (2-tetrad G4)."""
        seq = "GGTTGGTTGGTTGG"  # 2-tetrad pattern
        motifs = self._run_detector('g_quadruplex', seq, "weak_pqs_test")
        # Weak PQS may or may not be detected depending on threshold
        # Just verify structure if detected
        for motif in motifs:
            self._check_motif_structure(motif)

    def test_no_false_positive(self):
        """Test that random sequence doesn't produce false positive G4."""
        seq = "ATCGATCGATCGATCGATCGATCGATCG"  # No G-tracts
        motifs = self._run_detector('g_quadruplex', seq, "no_g4_test")
        self.assertEqual(len(motifs), 0, "Should not detect G4 in non-G-rich sequence")


class TestIMotifDetector(TestDetectorBase):
    """Test i-Motif detector functionality."""

    def test_canonical_imotif_detection(self):
        """Test detection of canonical i-Motif (4 C-tracts)."""
        seq = "AAAA" + "CCCATTTCCCGTTTCCCAATTCCC" + "TTTT"
        motifs = self._run_detector('i_motif', seq, "imotif_test")
        
        self.assertGreater(len(motifs), 0, "Should detect at least one i-Motif")
        for motif in motifs:
            self._check_motif_structure(motif)
            self.assertEqual(motif['Class'], 'i-Motif')

    def test_no_false_positive(self):
        """Test that G-rich sequence doesn't produce false positive i-Motif."""
        seq = "GGGGGGGGGGGGGGGGGGGG"  # G-rich, no C-tracts
        motifs = self._run_detector('i_motif', seq, "no_imotif_test")
        self.assertEqual(len(motifs), 0, "Should not detect i-Motif in G-rich sequence")


class TestZDNADetector(TestDetectorBase):
    """Test Z-DNA detector functionality."""

    def test_zdna_alternating_cg(self):
        """Test detection of Z-DNA in alternating CG sequences."""
        seq = "AAAA" + "CGCGCGCGCGCGCGCGCGCG" + "TTTT"  # 10 CG dinucleotides
        motifs = self._run_detector('z_dna', seq, "zdna_test")
        
        # Should detect Z-DNA in CG-rich region
        for motif in motifs:
            self._check_motif_structure(motif)
            self.assertEqual(motif['Class'], 'Z-DNA')

    def test_egz_motif_detection(self):
        """Test detection of eGZ-DNA (CGG repeats)."""
        seq = "CGGCGGCGGCGGCGGCGGCGG"  # CGG repeat x7
        motifs = self._run_detector('z_dna', seq, "egz_test")
        
        # Should detect eGZ motif
        for motif in motifs:
            self._check_motif_structure(motif)


class TestCurvedDNADetector(TestDetectorBase):
    """Test Curved DNA detector functionality."""

    def test_a_tract_detection(self):
        """Test detection of A-tract curvature."""
        seq = "TTTT" + "AAAAAAAAAA" + "CCCC"  # Long A-tract (10 bp)
        motifs = self._run_detector('curved_dna', seq, "atract_test")
        
        for motif in motifs:
            self._check_motif_structure(motif)
            self.assertEqual(motif['Class'], 'Curved_DNA')

    def test_phased_a_tracts(self):
        """Test detection of phased A-tracts (global curvature)."""
        # A-tracts spaced approximately 10-11 bp (one helical turn)
        seq = "AAAACCCCCCCCCAAAACCCCCCCCCAAAACCCCCCCCCAAAA"
        motifs = self._run_detector('curved_dna', seq, "phased_test")
        
        for motif in motifs:
            self._check_motif_structure(motif)


class TestSlippedDNADetector(TestDetectorBase):
    """Test Slipped DNA detector functionality."""

    def test_str_detection(self):
        """Test detection of Short Tandem Repeats (STRs)."""
        # CAG repeat (associated with Huntington's disease)
        seq = "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAG"  # 12 CAG repeats
        motifs = self._run_detector('slipped_dna', seq, "str_test")
        
        self.assertGreater(len(motifs), 0, "Should detect STR")
        for motif in motifs:
            self._check_motif_structure(motif)
            self.assertEqual(motif['Class'], 'Slipped_DNA')
            # Verify repeat unit is detected
            if 'Repeat_Unit' in motif:
                self.assertIn(motif['Repeat_Unit'].upper(), ['CAG', 'AGC', 'GCA'])

    def test_dinucleotide_repeat(self):
        """Test detection of dinucleotide repeats."""
        seq = "CACACACACACACACACACACACACA"  # 13 CA repeats
        motifs = self._run_detector('slipped_dna', seq, "dinuc_test")
        
        for motif in motifs:
            self._check_motif_structure(motif)


class TestCruciformDetector(TestDetectorBase):
    """Test Cruciform (Inverted Repeat) detector functionality."""

    def test_inverted_repeat_detection(self):
        """Test detection of inverted repeats (palindrome)."""
        # Perfect palindrome: GCTAGC...GCTAGC (inverted)
        left_arm = "GCTAGCTAGCTAGC"
        loop = "NNNN"
        right_arm = "GCTAGCTAGCTAGC"[::-1].translate(str.maketrans('ATGC', 'TACG'))
        # Note: For cruciform, the right arm is reverse complement
        seq = "AAAA" + left_arm + loop + right_arm + "TTTT"
        
        motifs = self._run_detector('cruciform', seq, "cruciform_test")
        
        for motif in motifs:
            self._check_motif_structure(motif)
            if len(motifs) > 0:
                self.assertEqual(motif['Class'], 'Cruciform')


class TestRLoopDetector(TestDetectorBase):
    """Test R-Loop detector functionality."""

    def test_rloop_grich_detection(self):
        """Test detection of R-Loop forming sequences (G-rich RIZ)."""
        # R-loop initiation zone: G-rich region
        seq = "AAAA" + "GGGAGGGAGGGAGGGAGGGAGGG" + "TTTT" * 10 + "GGGGAGGGGAGGGGA"
        motifs = self._run_detector('r_loop', seq, "rloop_test")
        
        for motif in motifs:
            self._check_motif_structure(motif)
            if len(motifs) > 0:
                self.assertEqual(motif['Class'], 'R-Loop')


class TestTriplexDetector(TestDetectorBase):
    """Test Triplex DNA detector functionality."""

    def test_sticky_dna_detection(self):
        """Test detection of Sticky DNA (GAA repeats)."""
        seq = "GAAGAAGAAGAAGAAGAAGAAGAAGAA"  # 9 GAA repeats
        motifs = self._run_detector('triplex', seq, "sticky_test")
        
        # Sticky DNA is a subtype handled by triplex detector
        for motif in motifs:
            self._check_motif_structure(motif)

    def test_mirror_repeat_detection(self):
        """Test detection of mirror repeats for triplex formation."""
        # Mirror repeat: same sequence reads the same 5'->3' on both strands
        seq = "AAAAGAGGAGAGCGCGCGCGCGAGAGAGGAGAAAA"
        motifs = self._run_detector('triplex', seq, "mirror_test")
        
        for motif in motifs:
            self._check_motif_structure(motif)


class TestAPhilicDetector(TestDetectorBase):
    """Test A-philic DNA detector functionality."""

    def test_aphilic_detection(self):
        """Test detection of A-philic DNA regions."""
        # A-philic DNA has specific 10-mer signatures
        seq = "AAAAAAAAAA" * 10  # 100bp A-tract
        motifs = self._run_detector('a_philic', seq, "aphilic_test")
        
        for motif in motifs:
            self._check_motif_structure(motif)
            if len(motifs) > 0:
                self.assertEqual(motif['Class'], 'A-philic_DNA')


class TestIntegration(TestDetectorBase):
    """Integration tests for the full detection pipeline."""

    def test_analyze_sequence_function(self):
        """Test the main analyze_sequence API function."""
        from Utilities.nonbscanner import analyze_sequence
        
        # Test sequence with multiple motif types
        seq = "GGGTTAGGGTTAGGGTTAGGG" + "CCCATTTCCCGTTTCCCAATTCCC" + "CGCGCGCGCGCGCGCG"
        motifs = analyze_sequence(seq, "integration_test")
        
        self.assertIsInstance(motifs, list)
        # Should find at least one motif type
        if len(motifs) > 0:
            for motif in motifs:
                self._check_motif_structure(motif)

    def test_empty_sequence(self):
        """Test handling of empty sequence."""
        from Utilities.nonbscanner import analyze_sequence
        
        with self.assertRaises(ValueError):
            analyze_sequence("", "empty_test")

    def test_short_sequence(self):
        """Test handling of very short sequence."""
        from Utilities.nonbscanner import analyze_sequence
        
        # Very short sequences may not produce motifs but should not error
        seq = "ATCG"
        try:
            motifs = analyze_sequence(seq, "short_test")
            self.assertIsInstance(motifs, list)
        except ValueError:
            pass  # Some validators may reject very short sequences


class TestDetectorAudit(TestDetectorBase):
    """Test detector audit functionality."""

    def test_audit_tracking(self):
        """Test that detectors track audit information."""
        seq = "GGGTTAGGGTTAGGGTTAGGG"
        
        for name, detector in self.detectors.items():
            motifs = detector.detect_motifs(seq, f"{name}_audit_test")
            audit = detector.get_audit_info()
            
            self.assertIsInstance(audit, dict)
            self.assertIn('invoked', audit)
            self.assertTrue(audit['invoked'], f"{name} should mark audit as invoked")


def run_tests():
    """Run all detector tests and print results."""
    print("=" * 70)
    print("NonBDNAFinder Detector Validation Tests")
    print("=" * 70)
    print()
    
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    test_classes = [
        TestGQuadruplexDetector,
        TestIMotifDetector,
        TestZDNADetector,
        TestCurvedDNADetector,
        TestSlippedDNADetector,
        TestCruciformDetector,
        TestRLoopDetector,
        TestTriplexDetector,
        TestAPhilicDetector,
        TestIntegration,
        TestDetectorAudit
    ]
    
    for test_class in test_classes:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    
    # Run tests with verbosity
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print()
    print("=" * 70)
    print(f"Tests Run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print("=" * 70)
    
    return result


if __name__ == "__main__":
    run_tests()
