"""
Chunk Size Comparison Analysis for E. coli Genome
=================================================

This script analyzes the E. coli genome with different chunk sizes
and compares the detection results to evaluate:
1. Total motifs detected at each chunk size
2. Processing time for each chunk size
3. Memory usage patterns
4. Motif distribution by class
5. Detection consistency across chunk sizes

Author: Dr. Venkata Rajesh Yella
"""

import sys
import os
import time
import json
from pathlib import Path
from datetime import datetime
from collections import defaultdict
from typing import List, Dict, Any, Tuple

# Add parent directory to path for imports
parent_dir = str(Path(__file__).parent.parent)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

import pandas as pd

# Import scanner with chunking capability
from Utilities.nonbscanner import NonBScanner, analyze_sequence


def generate_test_sequence(length: int, gc_content: float = 0.5) -> str:
    """Generate a random DNA sequence with specified GC content for testing."""
    import random
    
    gc_bases = ['G', 'C']
    at_bases = ['A', 'T']
    
    sequence = []
    for _ in range(length):
        if random.random() < gc_content:
            sequence.append(random.choice(gc_bases))
        else:
            sequence.append(random.choice(at_bases))
    
    return ''.join(sequence)


def generate_ecoli_like_sequence(length: int = 100000) -> str:
    """Generate a sequence with E. coli-like characteristics."""
    import random
    
    # E. coli has ~50.8% GC content
    gc_content = 0.508
    
    # Add some known motif regions
    sequence_parts = []
    position = 0
    
    # Add some G4-forming regions periodically
    while position < length:
        # Random segment
        seg_length = random.randint(1000, 5000)
        if position + seg_length > length:
            seg_length = length - position
        
        sequence_parts.append(generate_test_sequence(seg_length, gc_content))
        position += seg_length
        
        # Occasionally add G4 motif
        if position < length and random.random() < 0.1:
            g4_motif = "GGGTTTGGGCCCGGGAAAGGG"
            sequence_parts.append(g4_motif)
            position += len(g4_motif)
        
        # Occasionally add i-motif
        if position < length and random.random() < 0.1:
            imotif = "CCCATTTCCCGTTTCCCAATTCCC"
            sequence_parts.append(imotif)
            position += len(imotif)
        
        # Occasionally add A-tract
        if position < length and random.random() < 0.15:
            atract = "AAAAAAAAAA"
            sequence_parts.append(atract)
            position += len(atract)
    
    return ''.join(sequence_parts)[:length]


def analyze_with_chunking(sequence: str, sequence_name: str, chunk_size: int, 
                          overlap: int = 500) -> Tuple[List[Dict], float]:
    """
    Analyze a sequence with explicit chunking.
    
    Args:
        sequence: DNA sequence to analyze
        sequence_name: Name for the sequence
        chunk_size: Size of each chunk in bp
        overlap: Overlap between chunks in bp
        
    Returns:
        Tuple of (motifs_list, processing_time)
    """
    scanner = NonBScanner(enable_all_detectors=True)
    
    start_time = time.time()
    all_motifs = []
    seq_len = len(sequence)
    
    # Process in chunks
    chunk_start = 0
    chunk_num = 0
    
    while chunk_start < seq_len:
        chunk_end = min(chunk_start + chunk_size, seq_len)
        chunk_seq = sequence[chunk_start:chunk_end]
        
        # Analyze chunk
        chunk_motifs = scanner.analyze_sequence(chunk_seq, f"{sequence_name}_chunk{chunk_num}")
        
        # Adjust positions to be relative to full sequence
        for motif in chunk_motifs:
            motif['Start'] = motif.get('Start', 0) + chunk_start
            motif['End'] = motif.get('End', 0) + chunk_start
            motif['Chunk'] = chunk_num
        
        all_motifs.extend(chunk_motifs)
        
        # Move to next chunk with overlap
        chunk_start = chunk_end - overlap if chunk_end < seq_len else seq_len
        chunk_num += 1
    
    # Remove duplicate motifs from overlapping regions
    all_motifs = remove_chunk_duplicates(all_motifs)
    
    elapsed = time.time() - start_time
    return all_motifs, elapsed


def remove_chunk_duplicates(motifs: List[Dict]) -> List[Dict]:
    """Remove duplicate motifs that appear in overlapping chunk regions."""
    if not motifs:
        return motifs
    
    # Sort by position and score
    motifs.sort(key=lambda x: (x.get('Start', 0), -x.get('Score', 0)))
    
    unique = []
    seen_positions = set()
    
    for motif in motifs:
        key = (motif.get('Class'), motif.get('Start'), motif.get('End'))
        if key not in seen_positions:
            unique.append(motif)
            seen_positions.add(key)
    
    return unique


def run_chunk_comparison(test_sequence: str, chunk_sizes: List[int], 
                         sequence_name: str = "test_genome") -> Dict[str, Any]:
    """
    Run comparison analysis across different chunk sizes.
    
    Args:
        test_sequence: The DNA sequence to analyze
        chunk_sizes: List of chunk sizes to test (in bp)
        sequence_name: Name for the sequence
        
    Returns:
        Dictionary with comparison results
    """
    results = {
        'sequence_length': len(test_sequence),
        'sequence_name': sequence_name,
        'timestamp': datetime.now().isoformat(),
        'chunk_comparisons': [],
        'class_distributions': {}
    }
    
    print(f"\n{'='*70}")
    print(f"CHUNK SIZE COMPARISON ANALYSIS")
    print(f"Sequence: {sequence_name}")
    print(f"Length: {len(test_sequence):,} bp")
    print(f"{'='*70}")
    
    baseline_motifs = None
    
    for chunk_size in chunk_sizes:
        print(f"\n📊 Testing chunk size: {chunk_size:,} bp")
        print("-" * 50)
        
        # Calculate optimal overlap (10% of chunk size, min 500, max 2500)
        overlap = min(max(int(chunk_size * 0.1), 500), 2500)
        
        motifs, elapsed = analyze_with_chunking(
            test_sequence, sequence_name, chunk_size, overlap
        )
        
        # Count by class
        class_counts = defaultdict(int)
        for motif in motifs:
            class_counts[motif.get('Class', 'Unknown')] += 1
        
        # Calculate statistics
        num_chunks = (len(test_sequence) + chunk_size - 1) // chunk_size
        throughput = len(test_sequence) / elapsed if elapsed > 0 else 0
        
        comparison = {
            'chunk_size': chunk_size,
            'overlap': overlap,
            'num_chunks': num_chunks,
            'total_motifs': len(motifs),
            'processing_time': round(elapsed, 3),
            'throughput_bp_per_sec': round(throughput, 1),
            'class_counts': dict(class_counts)
        }
        
        # Calculate consistency with baseline
        if baseline_motifs is None:
            baseline_motifs = motifs
            comparison['vs_baseline'] = 'baseline'
        else:
            # Compare with baseline
            baseline_positions = {(m.get('Class'), m.get('Start'), m.get('End')) 
                                  for m in baseline_motifs}
            current_positions = {(m.get('Class'), m.get('Start'), m.get('End')) 
                                 for m in motifs}
            
            common = baseline_positions & current_positions
            only_baseline = baseline_positions - current_positions
            only_current = current_positions - baseline_positions
            
            comparison['vs_baseline'] = {
                'common_motifs': len(common),
                'only_in_baseline': len(only_baseline),
                'only_in_current': len(only_current),
                'jaccard_similarity': round(len(common) / len(baseline_positions | current_positions), 4) if baseline_positions | current_positions else 1.0
            }
        
        results['chunk_comparisons'].append(comparison)
        
        # Print summary
        print(f"  Total motifs found: {len(motifs)}")
        print(f"  Processing time: {elapsed:.3f}s")
        print(f"  Throughput: {throughput:,.0f} bp/s")
        print(f"  Number of chunks: {num_chunks}")
        print(f"\n  Motifs by class:")
        for cls, count in sorted(class_counts.items()):
            print(f"    {cls}: {count}")
        
        if comparison['vs_baseline'] != 'baseline':
            print(f"\n  Comparison with baseline (chunk_size={chunk_sizes[0]}bp):")
            print(f"    Common motifs: {comparison['vs_baseline']['common_motifs']}")
            print(f"    Only in baseline: {comparison['vs_baseline']['only_in_baseline']}")
            print(f"    Only in current: {comparison['vs_baseline']['only_in_current']}")
            print(f"    Jaccard similarity: {comparison['vs_baseline']['jaccard_similarity']:.4f}")
    
    return results


def generate_report(results: Dict[str, Any]) -> str:
    """Generate a markdown report from the comparison results."""
    report = []
    report.append("# NonBDNAFinder Chunk Size Comparison Report")
    report.append("")
    report.append(f"**Generated:** {results['timestamp']}")
    report.append(f"**Sequence:** {results['sequence_name']}")
    report.append(f"**Sequence Length:** {results['sequence_length']:,} bp")
    report.append("")
    
    report.append("## Summary Table")
    report.append("")
    report.append("| Chunk Size | Overlap | Chunks | Motifs | Time (s) | Throughput (bp/s) |")
    report.append("|------------|---------|--------|--------|----------|-------------------|")
    
    for comp in results['chunk_comparisons']:
        report.append(
            f"| {comp['chunk_size']:,} bp | {comp['overlap']} bp | "
            f"{comp['num_chunks']} | {comp['total_motifs']} | "
            f"{comp['processing_time']} | {comp['throughput_bp_per_sec']:,.0f} |"
        )
    
    report.append("")
    report.append("## Detection Consistency")
    report.append("")
    
    baseline_chunk = results['chunk_comparisons'][0]['chunk_size']
    report.append(f"Baseline chunk size: **{baseline_chunk:,} bp**")
    report.append("")
    
    report.append("| Chunk Size | Common | Only Baseline | Only Current | Jaccard |")
    report.append("|------------|--------|---------------|--------------|---------|")
    
    for comp in results['chunk_comparisons'][1:]:
        vs = comp['vs_baseline']
        report.append(
            f"| {comp['chunk_size']:,} bp | {vs['common_motifs']} | "
            f"{vs['only_in_baseline']} | {vs['only_in_current']} | "
            f"{vs['jaccard_similarity']:.4f} |"
        )
    
    report.append("")
    report.append("## Motif Distribution by Class")
    report.append("")
    
    # Get all classes
    all_classes = set()
    for comp in results['chunk_comparisons']:
        all_classes.update(comp['class_counts'].keys())
    
    classes = sorted(all_classes)
    
    header = "| Chunk Size | " + " | ".join(classes) + " |"
    separator = "|------------|" + "|".join(["-" * 15 for _ in classes]) + "|"
    
    report.append(header)
    report.append(separator)
    
    for comp in results['chunk_comparisons']:
        row = f"| {comp['chunk_size']:,} bp | "
        row += " | ".join(str(comp['class_counts'].get(c, 0)) for c in classes)
        row += " |"
        report.append(row)
    
    report.append("")
    report.append("## Analysis")
    report.append("")
    report.append("### Key Observations:")
    report.append("")
    
    # Analyze results
    chunk_data = results['chunk_comparisons']
    
    # Find optimal chunk size (best throughput)
    best_throughput = max(chunk_data, key=lambda x: x['throughput_bp_per_sec'])
    report.append(f"1. **Best throughput:** {best_throughput['chunk_size']:,} bp "
                  f"({best_throughput['throughput_bp_per_sec']:,.0f} bp/s)")
    
    # Find most consistent chunk size
    if len(chunk_data) > 1:
        consistencies = [(c['chunk_size'], c['vs_baseline'].get('jaccard_similarity', 1.0)) 
                         for c in chunk_data[1:]]
        best_consistency = max(consistencies, key=lambda x: x[1])
        report.append(f"2. **Most consistent:** {best_consistency[0]:,} bp "
                      f"(Jaccard: {best_consistency[1]:.4f})")
    
    # Motif count variation
    motif_counts = [c['total_motifs'] for c in chunk_data]
    report.append(f"3. **Motif count range:** {min(motif_counts)} - {max(motif_counts)}")
    
    report.append("")
    report.append("### Recommendations:")
    report.append("")
    report.append("Based on the analysis:")
    
    # Calculate optimal based on balance of speed and consistency
    if len(chunk_data) > 1:
        scores = []
        for c in chunk_data[1:]:
            # Normalize throughput and jaccard
            throughput_norm = c['throughput_bp_per_sec'] / best_throughput['throughput_bp_per_sec']
            jaccard = c['vs_baseline'].get('jaccard_similarity', 1.0)
            combined = 0.4 * throughput_norm + 0.6 * jaccard  # Weight consistency higher
            scores.append((c['chunk_size'], combined, c))
        
        optimal = max(scores, key=lambda x: x[1])
        report.append(f"- **Recommended chunk size:** {optimal[0]:,} bp")
        report.append(f"  - Provides good balance of speed and detection consistency")
    
    report.append(f"- **For maximum speed:** Use {best_throughput['chunk_size']:,} bp chunks")
    report.append(f"- **For genome-scale analysis:** Consider 25,000-50,000 bp chunks with 2,500 bp overlap")
    
    return "\n".join(report)


def main():
    """Main function to run the chunk size comparison analysis."""
    print("\n" + "="*70)
    print("NonBDNAFinder - Chunk Size Comparison Analysis")
    print("="*70)
    
    # Test with a moderate-sized synthetic sequence
    # (In production, you would use a real E. coli genome file)
    print("\nGenerating test sequence with E. coli-like characteristics...")
    
    # Generate a 100kb test sequence (smaller for faster testing)
    test_sequence = generate_ecoli_like_sequence(100000)  # 100kb
    
    print(f"Test sequence generated: {len(test_sequence):,} bp")
    
    # Define chunk sizes to test
    chunk_sizes = [5000, 10000, 25000, 50000, 100000]
    
    # Run comparison
    results = run_chunk_comparison(test_sequence, chunk_sizes, "ecoli_like_test")
    
    # Generate report
    report = generate_report(results)
    
    # Save report
    report_path = Path(__file__).parent / "chunk_comparison_report.md"
    with open(report_path, 'w') as f:
        f.write(report)
    
    print(f"\n{'='*70}")
    print(f"Report saved to: {report_path}")
    print("="*70)
    
    # Also save JSON results
    json_path = Path(__file__).parent / "chunk_comparison_results.json"
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"JSON results saved to: {json_path}")
    
    # Print report summary
    print("\n" + "="*70)
    print("REPORT SUMMARY")
    print("="*70)
    print(report)
    
    return results


if __name__ == "__main__":
    main()
