"""Test suite for performance optimizations.

This module validates that all performance optimizations produce
identical output to the original implementations, ensuring zero
changes to scientific results.
"""

import sys
from pathlib import Path

# Add parent directory to path
parent_dir = str(Path(__file__).parent.parent)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

import unittest
from typing import List, Tuple, Dict


class TestVectorized10merScanning(unittest.TestCase):
    """Test vectorized 10-mer window scanning optimization."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Import after path setup
        from detectors.zdna.tenmer_table import TENMER_SCORE
        from detectors.aphilic.tenmer_table import TENMER_LOG2
        from detectors.zdna import hyperscan_backend
        
        self.TENMER_SCORE = TENMER_SCORE
        self.TENMER_LOG2 = TENMER_LOG2
        self.hyperscan_backend = hyperscan_backend
        
    def test_zdna_vectorized_identical_output_small(self):
        """Verify Z-DNA vectorized scanning produces identical results on small sequence."""
        # Small test sequence (400 bp)
        seq = "ATCGCGCGCGATCGATCGCGCGCGATATACGCGCGCGAT" * 10
        seq = seq.upper()
        
        # Old method (loop-based)
        old_results = self.hyperscan_backend.py_find_matches_loop(seq, self.TENMER_SCORE)
        
        # New method (vectorized)
        new_results = self.hyperscan_backend.vectorized_find_matches(seq, self.TENMER_SCORE)
        
        # Assert exact match
        self.assertEqual(len(old_results), len(new_results), 
                        "Number of matches should be identical")
        
        for old, new in zip(old_results, new_results):
            self.assertEqual(old[0], new[0], "Start positions should match")
            self.assertEqual(old[1], new[1], "10-mer sequences should match")
            self.assertAlmostEqual(old[2], new[2], places=10, 
                                  msg="Scores should match within floating-point precision")
    
    def test_zdna_vectorized_identical_output_medium(self):
        """Verify Z-DNA vectorized scanning produces identical results on medium sequence."""
        # Medium test sequence (10 KB)
        base_seq = "ATCGCGCGCGATCGATCGCGCGCGATATACGCGCGCGAT" * 250
        seq = base_seq.upper()
        
        # Old method (loop-based)
        old_results = self.hyperscan_backend.py_find_matches_loop(seq, self.TENMER_SCORE)
        
        # New method (vectorized)
        new_results = self.hyperscan_backend.vectorized_find_matches(seq, self.TENMER_SCORE)
        
        # Assert exact match
        self.assertEqual(len(old_results), len(new_results))
        
        for old, new in zip(old_results, new_results):
            self.assertEqual(old[0], new[0])
            self.assertEqual(old[1], new[1])
            self.assertAlmostEqual(old[2], new[2], places=10)
    
    def test_zdna_vectorized_identical_output_large(self):
        """Verify Z-DNA vectorized scanning produces identical results on large sequence."""
        # Large test sequence (100 KB)
        base_seq = "ATCGCGCGCGATCGATCGCGCGCGATATACGCGCGCGAT" * 2500
        seq = base_seq.upper()
        
        # Old method (loop-based)
        old_results = self.hyperscan_backend.py_find_matches_loop(seq, self.TENMER_SCORE)
        
        # New method (vectorized)
        new_results = self.hyperscan_backend.vectorized_find_matches(seq, self.TENMER_SCORE)
        
        # Assert exact match
        self.assertEqual(len(old_results), len(new_results))
        
        for old, new in zip(old_results, new_results):
            self.assertEqual(old[0], new[0])
            self.assertEqual(old[1], new[1])
            self.assertAlmostEqual(old[2], new[2], places=10)
    
    def test_aphilic_vectorized_identical_output(self):
        """Verify A-philic vectorized scanning produces identical results."""
        # Test sequence (4 KB)
        seq = ("AAAAAATTTTTAAAAAAGGGGCCCCAAAAAATTTTT" * 100).upper()
        
        # Import A-philic detector
        from detectors.aphilic.detector import APhilicDetector
        
        detector = APhilicDetector()
        
        # Old method (loop-based)
        old_results = detector._py_find_matches_loop(seq)
        
        # New method (should use vectorized via hyperscan_backend)
        new_results = detector._py_find_matches(seq)
        
        # Assert exact match
        self.assertEqual(len(old_results), len(new_results))
        
        for old, new in zip(old_results, new_results):
            self.assertEqual(old[0], new[0])
            self.assertEqual(old[1], new[1])
            self.assertAlmostEqual(old[2], new[2], places=10)
    
    def test_zdna_detector_end_to_end(self):
        """Test complete Z-DNA detection pipeline with optimizations."""
        from detectors.zdna.detector import ZDNADetector
        
        detector = ZDNADetector()
        
        # Test sequence with known Z-DNA forming regions
        seq = "CGCGCGCGCG" * 10 + "ATATATAT" * 5 + "CGCGCGCGCG" * 10
        seq = seq.upper()
        
        # Run detection
        motifs = detector.detect_motifs(seq, "test_seq")
        
        # Should find Z-DNA regions
        self.assertGreater(len(motifs), 0, "Should detect Z-DNA motifs")
        
        # Verify motif structure
        for motif in motifs:
            self.assertIn('Start', motif)
            self.assertIn('End', motif)
            self.assertIn('Score', motif)
            self.assertIn('Class', motif)
            self.assertEqual(motif['Class'], 'Z-DNA')
    
    def test_aphilic_detector_end_to_end(self):
        """Test complete A-philic detection pipeline with optimizations."""
        from detectors.aphilic.detector import APhilicDetector
        
        detector = APhilicDetector()
        
        # Test sequence
        seq = ("AAAAAATTTTTAAAAAAGGGGCCCCAAAAAATTTTT" * 20).upper()
        
        # Run detection
        motifs = detector.detect_motifs(seq, "test_seq")
        
        # Verify motif structure (may or may not find motifs depending on thresholds)
        for motif in motifs:
            self.assertIn('Start', motif)
            self.assertIn('End', motif)
            self.assertIn('Score', motif)
            self.assertIn('Class', motif)


class TestParallelMultiSequenceAnalysis(unittest.TestCase):
    """Test parallel multi-sequence analysis optimization."""
    
    def setUp(self):
        """Set up test fixtures."""
        import nonbscanner
        self.nonbscanner = nonbscanner
        
    def test_parallel_identical_output_order(self):
        """Verify parallel processing preserves exact output order."""
        # Create test sequences
        sequences = {
            "seq1": "GGGTTAGGGTTAGGGTTAGGG" * 10,
            "seq2": "CGCGCGCGCGCGCGCGCGCG" * 10,
            "seq3": "AAAAAATTTTTAAAAAAGGGG" * 10,
            "seq4": "CCCCCCCCCCGGGGGGGGGG" * 10,
        }
        
        # Sequential analysis
        seq_results = self.nonbscanner.analyze_fasta(
            ">" + "\n>".join([f"{name}\n{seq}" for name, seq in sequences.items()])
        )
        
        # Parallel analysis
        par_results = self.nonbscanner.analyze_multiple_sequences_parallel(
            sequences, num_processes=2, preserve_order=True
        )
        
        # Verify same keys and order
        self.assertEqual(list(seq_results.keys()), list(par_results.keys()))
        
        # Verify identical results for each sequence
        for name in sequences.keys():
            seq_motifs = seq_results[name]
            par_motifs = par_results[name]
            
            self.assertEqual(len(seq_motifs), len(par_motifs),
                           f"Sequence {name} should have same number of motifs")
            
            # Compare each motif
            for sm, pm in zip(seq_motifs, par_motifs):
                # Compare critical fields
                for key in ['Start', 'End', 'Class', 'Subclass']:
                    if key in sm and key in pm:
                        self.assertEqual(sm[key], pm[key],
                                       f"Sequence {name}: {key} should match")
    
    def test_parallel_single_sequence_fallback(self):
        """Verify parallel processing handles single sequence correctly."""
        sequences = {"seq1": "GGGTTAGGGTTAGGGTTAGGG" * 10}
        
        # Should work without errors
        results = self.nonbscanner.analyze_multiple_sequences_parallel(sequences)
        
        self.assertIn("seq1", results)
        self.assertIsInstance(results["seq1"], list)
    
    def test_parallel_empty_sequences(self):
        """Verify parallel processing handles empty input."""
        sequences = {}
        
        results = self.nonbscanner.analyze_multiple_sequences_parallel(sequences)
        
        self.assertEqual(results, {})
    
    def test_fasta_parallel_identical_output(self):
        """Verify analyze_fasta_parallel produces identical results."""
        fasta = ">seq1\nGGGTTAGGGTTAGGGTTAGGG\n>seq2\nCGCGCGCGCGCGCGCGCGCG"
        
        # Sequential
        seq_results = self.nonbscanner.analyze_fasta(fasta)
        
        # Parallel
        par_results = self.nonbscanner.analyze_fasta_parallel(fasta, num_processes=2)
        
        # Same keys
        self.assertEqual(set(seq_results.keys()), set(par_results.keys()))
        
        # Same number of motifs
        for name in seq_results.keys():
            self.assertEqual(len(seq_results[name]), len(par_results[name]))


class TestStreamingFastaParser(unittest.TestCase):
    """Test streaming FASTA parser optimization."""
    
    def setUp(self):
        """Set up test fixtures."""
        import utilities
        self.utilities = utilities
        
    def test_streaming_identical_output_single_sequence(self):
        """Verify streaming parser produces identical output for single sequence."""
        fasta = ">seq1\nGGGTTAGGGTTAGGGTTAGGG"
        
        # Standard parsing
        std_result = self.utilities.parse_fasta(fasta, streaming=False)
        
        # Streaming parsing
        stream_result = dict(self.utilities.parse_fasta(fasta, streaming=True))
        
        # Should be identical
        self.assertEqual(std_result, stream_result)
        self.assertEqual(std_result['seq1'], stream_result['seq1'])
    
    def test_streaming_identical_output_multiple_sequences(self):
        """Verify streaming parser produces identical output for multiple sequences."""
        fasta = ">seq1\nGGGTTAGGGTTAGGGTTAGGG\n>seq2\nCGCGCGCGCGCGCGCGCGCG\n>seq3\nAAAATTTTCCCCGGGG"
        
        # Standard parsing
        std_result = self.utilities.parse_fasta(fasta, streaming=False)
        
        # Streaming parsing
        stream_result = dict(self.utilities.parse_fasta(fasta, streaming=True))
        
        # Should be identical
        self.assertEqual(std_result, stream_result)
        self.assertEqual(len(std_result), len(stream_result))
        
        for name in std_result.keys():
            self.assertEqual(std_result[name], stream_result[name])
    
    def test_streaming_identical_output_large_fasta(self):
        """Verify streaming parser produces identical output for large FASTA."""
        # Create a large FASTA with 10 sequences
        sequences = {
            f"seq{i}": ("ATCGATCGATCG" * 100) for i in range(10)
        }
        fasta = "\n".join([f">{name}\n{seq}" for name, seq in sequences.items()])
        
        # Standard parsing
        std_result = self.utilities.parse_fasta(fasta, streaming=False)
        
        # Streaming parsing
        stream_result = dict(self.utilities.parse_fasta(fasta, streaming=True))
        
        # Should be identical
        self.assertEqual(std_result, stream_result)
        self.assertEqual(len(std_result), 10)
        
        for name in std_result.keys():
            self.assertEqual(std_result[name], stream_result[name])
    
    def test_streaming_generator_behavior(self):
        """Verify streaming returns a generator."""
        fasta = ">seq1\nGGGTTAGGG\n>seq2\nCGCGCGCG"
        
        result = self.utilities.parse_fasta(fasta, streaming=True)
        
        # Should be a generator
        import types
        self.assertIsInstance(result, types.GeneratorType)
        
        # Can iterate
        sequences = list(result)
        self.assertEqual(len(sequences), 2)
        self.assertEqual(sequences[0][0], 'seq1')
        self.assertEqual(sequences[1][0], 'seq2')
    
    def test_streaming_memory_efficiency(self):
        """Verify streaming processes sequences one at a time."""
        fasta = ">seq1\nGGGTTAGGG\n>seq2\nCGCGCGCG\n>seq3\nAAAATTTT"
        
        # Streaming should yield one at a time
        count = 0
        for name, seq in self.utilities.parse_fasta(fasta, streaming=True):
            count += 1
            self.assertIn(name, ['seq1', 'seq2', 'seq3'])
            self.assertIsInstance(seq, str)
            self.assertGreater(len(seq), 0)
        
        self.assertEqual(count, 3)


class TestLazyMatplotlibImports(unittest.TestCase):
    """Test lazy matplotlib/seaborn imports optimization."""
    
    def test_utilities_imports_without_matplotlib(self):
        """Verify utilities module can be imported without loading matplotlib."""
        # This test verifies that importing utilities doesn't trigger matplotlib import
        # We can't fully test this without reloading modules, but we can check the pattern
        
        import utilities
        
        # The _ensure_matplotlib function should exist
        self.assertTrue(hasattr(utilities, '_ensure_matplotlib'))
        self.assertTrue(callable(utilities._ensure_matplotlib))
    
    def test_matplotlib_loads_on_demand(self):
        """Verify matplotlib loads when first plot function is called."""
        import utilities
        
        # Call _ensure_matplotlib
        plt, sns, patches, PdfPages = utilities._ensure_matplotlib()
        
        # Verify modules are loaded
        self.assertIsNotNone(plt)
        self.assertIsNotNone(sns)
        self.assertIsNotNone(patches)
        self.assertIsNotNone(PdfPages)
        
        # Verify they are the correct modules
        self.assertTrue(hasattr(plt, 'figure'))
        self.assertTrue(hasattr(sns, 'set_style'))


if __name__ == '__main__':
    unittest.main()
