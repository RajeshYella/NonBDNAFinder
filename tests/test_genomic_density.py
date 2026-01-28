"""
Test suite for genomic density and GC content calculations.

This module validates that:
1. Hybrid and Cluster motifs are excluded from genomic density calculations
2. GC percentage calculations are consistent across the codebase
3. Coverage calculations handle overlaps correctly
"""

import sys
from pathlib import Path

# Add parent directory to path
parent_dir = str(Path(__file__).parent.parent)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

import unittest
from typing import List, Dict, Any
from utilities import calculate_genomic_density, gc_content


class TestGenomicDensityExclusions(unittest.TestCase):
    """Test that Hybrid and Cluster motifs are excluded from genomic density calculations."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.sequence_length = 1000
        
    def test_exclude_hybrid_motifs(self):
        """Verify Hybrid motifs are excluded from density calculations."""
        motifs = [
            {'Class': 'G-Quadruplex', 'Start': 1, 'End': 50, 'Subclass': 'Test'},  # 50 bases
            {'Class': 'Hybrid', 'Start': 25, 'End': 75, 'Subclass': 'Dynamic overlaps'},
            {'Class': 'Z-DNA', 'Start': 100, 'End': 150, 'Subclass': 'Test'},  # 51 bases
        ]
        
        # Calculate density
        density = calculate_genomic_density(motifs, self.sequence_length, by_class=False)
        
        # Expected coverage: 50 bp (G-Quadruplex) + 51 bp (Z-DNA) = 101 bp = 10.1%
        # Hybrid should NOT be counted
        self.assertAlmostEqual(density['Overall'], 10.1, places=1)
        
    def test_exclude_cluster_motifs(self):
        """Verify Cluster motifs are excluded from density calculations."""
        motifs = [
            {'Class': 'G-Quadruplex', 'Start': 1, 'End': 50, 'Subclass': 'Test'},  # 50 bases
            {'Class': 'Non-B_DNA_Clusters', 'Start': 25, 'End': 100, 'Subclass': 'Dynamic clusters'},
            {'Class': 'Z-DNA', 'Start': 200, 'End': 250, 'Subclass': 'Test'},  # 51 bases
        ]
        
        # Calculate density
        density = calculate_genomic_density(motifs, self.sequence_length, by_class=False)
        
        # Expected coverage: 50 bp (G-Quadruplex) + 51 bp (Z-DNA) = 101 bp = 10.1%
        # Cluster should NOT be counted
        self.assertAlmostEqual(density['Overall'], 10.1, places=1)
        
    def test_exclude_both_hybrid_and_cluster(self):
        """Verify both Hybrid and Cluster motifs are excluded."""
        motifs = [
            {'Class': 'G-Quadruplex', 'Start': 1, 'End': 100, 'Subclass': 'Test'},  # 100 bases
            {'Class': 'Hybrid', 'Start': 50, 'End': 150, 'Subclass': 'Dynamic overlaps'},
            {'Class': 'Non-B_DNA_Clusters', 'Start': 100, 'End': 200, 'Subclass': 'Dynamic clusters'},
            {'Class': 'Z-DNA', 'Start': 300, 'End': 350, 'Subclass': 'Test'},  # 51 bases
        ]
        
        # Calculate density
        density = calculate_genomic_density(motifs, self.sequence_length, by_class=False)
        
        # Expected coverage: 100 bp (G-Quadruplex) + 51 bp (Z-DNA) = 151 bp = 15.1%
        # Both Hybrid and Cluster should NOT be counted
        self.assertAlmostEqual(density['Overall'], 15.1, places=1)
        
    def test_only_hybrid_and_cluster_returns_zero(self):
        """Verify that if only Hybrid and Cluster motifs exist, density is 0."""
        motifs = [
            {'Class': 'Hybrid', 'Start': 1, 'End': 100, 'Subclass': 'Dynamic overlaps'},
            {'Class': 'Non-B_DNA_Clusters', 'Start': 200, 'End': 300, 'Subclass': 'Dynamic clusters'},
        ]
        
        # Calculate density
        density = calculate_genomic_density(motifs, self.sequence_length, by_class=False)
        
        # Expected: 0% because all motifs are excluded
        self.assertEqual(density['Overall'], 0.0)
        
    def test_exclude_hybrid_cluster_by_class(self):
        """Verify exclusion works with by_class=True."""
        motifs = [
            {'Class': 'G-Quadruplex', 'Start': 1, 'End': 100, 'Subclass': 'Test'},  # 100 bases
            {'Class': 'Hybrid', 'Start': 50, 'End': 150, 'Subclass': 'Dynamic overlaps'},
            {'Class': 'Z-DNA', 'Start': 200, 'End': 300, 'Subclass': 'Test'},  # 101 bases
        ]
        
        # Calculate density by class
        density = calculate_genomic_density(motifs, self.sequence_length, by_class=True)
        
        # Should have entries for G-Quadruplex and Z-DNA, but NOT Hybrid
        self.assertIn('G-Quadruplex', density)
        self.assertIn('Z-DNA', density)
        self.assertNotIn('Hybrid', density)
        self.assertNotIn('Non-B_DNA_Clusters', density)
        
        # Overall density should be 20.1% (100 + 101 = 201 bp)
        self.assertAlmostEqual(density['Overall'], 20.1, places=1)
        
    def test_exclude_hybrid_cluster_by_subclass(self):
        """Verify exclusion works with by_subclass=True."""
        motifs = [
            {'Class': 'G-Quadruplex', 'Start': 1, 'End': 100, 'Subclass': 'Canonical'},  # 100 bases
            {'Class': 'Hybrid', 'Start': 50, 'End': 150, 'Subclass': 'Dynamic overlaps'},
            {'Class': 'Non-B_DNA_Clusters', 'Start': 100, 'End': 200, 'Subclass': 'Dynamic clusters'},
            {'Class': 'Z-DNA', 'Start': 300, 'End': 400, 'Subclass': 'CG-rich'},  # 101 bases
        ]
        
        # Calculate density by subclass
        density = calculate_genomic_density(motifs, self.sequence_length, by_subclass=True)
        
        # Should have entries for G-Quadruplex and Z-DNA subclasses, but NOT Hybrid or Cluster
        self.assertIn('G-Quadruplex:Canonical', density)
        self.assertIn('Z-DNA:CG-rich', density)
        self.assertNotIn('Hybrid:Dynamic overlaps', density)
        self.assertNotIn('Non-B_DNA_Clusters:Dynamic clusters', density)
        
        # Overall density should be 20.1% (100 + 101 = 201 bp)
        self.assertAlmostEqual(density['Overall'], 20.1, places=1)


class TestGenomicDensityOverlapHandling(unittest.TestCase):
    """Test that overlapping motifs are handled correctly in density calculations."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.sequence_length = 1000
        
    def test_overlapping_motifs_same_class(self):
        """Verify overlapping motifs of same class don't double-count."""
        motifs = [
            {'Class': 'G-Quadruplex', 'Start': 1, 'End': 100, 'Subclass': 'Test'},
            {'Class': 'G-Quadruplex', 'Start': 50, 'End': 150, 'Subclass': 'Test'},
        ]
        
        # Calculate density
        density = calculate_genomic_density(motifs, self.sequence_length, by_class=False)
        
        # Expected coverage: unique positions 1-150 = 150 bp = 15%
        self.assertEqual(density['Overall'], 15.0)
        
    def test_overlapping_motifs_different_classes(self):
        """Verify overlapping motifs of different classes count each position once."""
        motifs = [
            {'Class': 'G-Quadruplex', 'Start': 1, 'End': 100, 'Subclass': 'Test'},
            {'Class': 'Z-DNA', 'Start': 50, 'End': 150, 'Subclass': 'Test'},
        ]
        
        # Calculate density
        density = calculate_genomic_density(motifs, self.sequence_length, by_class=False)
        
        # Expected coverage: unique positions 1-150 = 150 bp = 15%
        self.assertEqual(density['Overall'], 15.0)
        
    def test_multiple_overlapping_regions(self):
        """Verify complex overlapping regions are handled correctly."""
        motifs = [
            {'Class': 'G-Quadruplex', 'Start': 1, 'End': 50, 'Subclass': 'Test'},  # 50 bases
            {'Class': 'Z-DNA', 'Start': 25, 'End': 75, 'Subclass': 'Test'},  # Overlaps
            {'Class': 'Cruciform', 'Start': 60, 'End': 100, 'Subclass': 'Test'},  # Overlaps
            {'Class': 'Slipped', 'Start': 200, 'End': 250, 'Subclass': 'Test'},  # 51 bases
        ]
        
        # Calculate density
        density = calculate_genomic_density(motifs, self.sequence_length, by_class=False)
        
        # Expected coverage: 1-100 (100 bp) + 200-250 (51 bp) = 151 bp = 15.1%
        self.assertAlmostEqual(density['Overall'], 15.1, places=1)


class TestGCContentConsistency(unittest.TestCase):
    """Test that GC content calculations are consistent."""
    
    def test_gc_content_basic(self):
        """Test basic GC content calculation."""
        seq = "ATCGATCG"
        gc_pct = gc_content(seq)
        
        # 50% GC (4 G/C out of 8 total)
        self.assertEqual(gc_pct, 50.0)
        
    def test_gc_content_all_gc(self):
        """Test GC content with all GC bases."""
        seq = "GGGGCCCC"
        gc_pct = gc_content(seq)
        
        # 100% GC
        self.assertEqual(gc_pct, 100.0)
        
    def test_gc_content_all_at(self):
        """Test GC content with all AT bases."""
        seq = "AAAATTTT"
        gc_pct = gc_content(seq)
        
        # 0% GC
        self.assertEqual(gc_pct, 0.0)
        
    def test_gc_content_empty_sequence(self):
        """Test GC content with empty sequence."""
        seq = ""
        gc_pct = gc_content(seq)
        
        # Should return 0.0 for empty sequence
        self.assertEqual(gc_pct, 0.0)
        
    def test_gc_content_case_insensitive(self):
        """Test GC content is case-insensitive."""
        seq1 = "atcgatcg"
        seq2 = "ATCGATCG"
        seq3 = "AtCgAtCg"
        
        gc1 = gc_content(seq1)
        gc2 = gc_content(seq2)
        gc3 = gc_content(seq3)
        
        # All should return same result
        self.assertEqual(gc1, 50.0)
        self.assertEqual(gc2, 50.0)
        self.assertEqual(gc3, 50.0)
        
    def test_gc_content_with_ambiguous_bases(self):
        """Test GC content calculation with standard bases only."""
        # GC content should only count G and C
        seq = "ATCGATCGNNN"
        gc_pct = gc_content(seq)
        
        # 4 G/C out of 11 total = 36.36%
        self.assertAlmostEqual(gc_pct, 36.36, places=2)


class TestGenomicDensityEdgeCases(unittest.TestCase):
    """Test edge cases in genomic density calculations."""
    
    def test_empty_motifs_list(self):
        """Test density calculation with empty motifs list."""
        motifs = []
        density = calculate_genomic_density(motifs, 1000, by_class=False)
        
        self.assertEqual(density['Overall'], 0.0)
        
    def test_zero_sequence_length(self):
        """Test density calculation with zero sequence length."""
        motifs = [{'Class': 'G-Quadruplex', 'Start': 1, 'End': 100, 'Subclass': 'Test'}]
        density = calculate_genomic_density(motifs, 0, by_class=False)
        
        self.assertEqual(density['Overall'], 0.0)
        
    def test_motif_covers_entire_sequence(self):
        """Test density when motif covers entire sequence."""
        motifs = [{'Class': 'G-Quadruplex', 'Start': 1, 'End': 1000, 'Subclass': 'Test'}]
        density = calculate_genomic_density(motifs, 1000, by_class=False)
        
        # Should be 100%
        self.assertEqual(density['Overall'], 100.0)
        
    def test_density_capped_at_100(self):
        """Test that density is capped at 100% even with calculation errors."""
        # This tests the min(..., 100.0) cap in the function
        motifs = [
            {'Class': 'G-Quadruplex', 'Start': 1, 'End': 1000, 'Subclass': 'Test'},
            {'Class': 'Z-DNA', 'Start': 1, 'End': 1000, 'Subclass': 'Test'},
        ]
        density = calculate_genomic_density(motifs, 1000, by_class=False)
        
        # Should be capped at 100%
        self.assertEqual(density['Overall'], 100.0)


if __name__ == '__main__':
    unittest.main()
