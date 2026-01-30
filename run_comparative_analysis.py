#!/usr/bin/env python3
"""
Comparative Genomics Analysis Script for Non-B DNA Motifs

This script analyzes all genomes in the Genomes folder using NonBDNAFinder
and generates comprehensive comparative statistics and results.

Author: NonBDNAFinder Team
Date: 2024
"""

import os
import sys
import time
import json
import pandas as pd
from collections import defaultdict
from datetime import datetime

# Add the parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import nonbscanner as nbs

# Genome metadata with full scientific names and characteristics
GENOME_METADATA = {
    'Buchnera aphidicola.fna': {
        'full_name': 'Buchnera aphidicola',
        'strain': 'str. APS',
        'category': 'Obligate endosymbiont',
        'phylum': 'Proteobacteria',
        'note': 'Smallest free-living bacterial genome'
    },
    'Candidatus Carsonella ruddii.fna': {
        'full_name': 'Candidatus Carsonella ruddii',
        'strain': 'PV',
        'category': 'Obligate endosymbiont',
        'phylum': 'Proteobacteria',
        'note': 'One of the smallest known bacterial genomes'
    },
    'ecoli.fna': {
        'full_name': 'Escherichia coli',
        'strain': 'K-12 substr. MG1655',
        'category': 'Free-living',
        'phylum': 'Proteobacteria',
        'note': 'Model organism'
    },
    'mtb.fasta': {
        'full_name': 'Mycobacterium tuberculosis',
        'strain': 'H37Rv',
        'category': 'Pathogen',
        'phylum': 'Actinobacteria',
        'note': 'High-GC content, causative agent of TB'
    },
    'saureus.fna': {
        'full_name': 'Staphylococcus aureus',
        'strain': 'NCTC 8325',
        'category': 'Pathogen',
        'phylum': 'Firmicutes',
        'note': 'Opportunistic human pathogen'
    },
    'Streptococcus pneumoniae.fna': {
        'full_name': 'Streptococcus pneumoniae',
        'strain': 'TIGR4',
        'category': 'Pathogen',
        'phylum': 'Firmicutes',
        'note': 'Respiratory pathogen'
    },
    'Cellulomonas shaoxiangyii.fna': {
        'full_name': 'Cellulomonas shaoxiangyii',
        'strain': 'Z28',
        'category': 'Free-living',
        'phylum': 'Actinobacteria',
        'note': 'High-GC, cellulose-degrading'
    },
    'Miltoncostaea marina.fna': {
        'full_name': 'Miltoncostaea marina',
        'strain': 'DSM 45384',
        'category': 'Free-living',
        'phylum': 'Actinobacteria',
        'note': 'High-GC marine bacterium'
    },
    'Scer.fna': {
        'full_name': 'Saccharomyces cerevisiae',
        'strain': 'S288C',
        'category': 'Model eukaryote',
        'phylum': 'Ascomycota',
        'note': 'Baker\'s yeast, model eukaryotic organism'
    }
}


def calculate_gc_content(sequence):
    """Calculate GC content percentage of a sequence."""
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total = len(sequence)
    if total == 0:
        return 0.0
    return (gc_count / total) * 100


def read_fasta_basic(filepath):
    """Read FASTA file and return sequences dictionary."""
    sequences = {}
    current_name = None
    current_seq = []
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    sequences[current_name] = ''.join(current_seq)
                current_name = line[1:].split()[0]  # Take first word as name
                current_seq = []
            elif current_name:
                current_seq.append(line)
    
    if current_name:
        sequences[current_name] = ''.join(current_seq)
    
    return sequences


def analyze_genome(filepath, genome_filename):
    """
    Analyze a single genome file and return results with statistics.
    """
    metadata = GENOME_METADATA.get(genome_filename, {'full_name': genome_filename.replace('.fna', '').replace('.fasta', '')})
    
    print(f"\n{'='*70}")
    print(f"Analyzing: {metadata.get('full_name', genome_filename)}")
    print(f"File: {genome_filename}")
    print(f"{'='*70}")
    
    # Read sequences
    start_time = time.time()
    sequences = read_fasta_basic(filepath)
    read_time = time.time() - start_time
    
    # Calculate total length and GC content
    total_length = sum(len(seq) for seq in sequences.values())
    full_sequence = ''.join(sequences.values())
    gc_content = calculate_gc_content(full_sequence)
    
    print(f"Sequences: {len(sequences)}")
    print(f"Total length: {total_length:,} bp ({total_length/1e6:.2f} Mb)")
    print(f"GC content: {gc_content:.1f}%")
    print(f"Read time: {read_time:.2f}s")
    
    # Analyze with NonBDNAFinder
    print(f"\nRunning Non-B DNA motif detection...")
    detection_start = time.time()
    
    all_motifs = []
    for seq_name, seq in sequences.items():
        short_name = seq_name[:50] + '...' if len(seq_name) > 50 else seq_name
        print(f"  Processing: {short_name} ({len(seq):,} bp)", end='', flush=True)
        
        seq_start = time.time()
        motifs = nbs.analyze_sequence(seq, seq_name)
        seq_time = time.time() - seq_start
        
        all_motifs.extend(motifs)
        print(f" -> {len(motifs)} motifs ({seq_time:.1f}s)")
    
    detection_time = time.time() - detection_start
    
    # Calculate motif class statistics
    class_counts = defaultdict(int)
    subclass_counts = defaultdict(int)
    for motif in all_motifs:
        class_counts[motif['Class']] += 1
        subclass_counts[f"{motif['Class']}:{motif['Subclass']}"] += 1
    
    # Build results dictionary
    results = {
        'genome_file': genome_filename,
        'organism': metadata.get('full_name', ''),
        'strain': metadata.get('strain', ''),
        'category': metadata.get('category', ''),
        'phylum': metadata.get('phylum', ''),
        'note': metadata.get('note', ''),
        'num_sequences': len(sequences),
        'total_length_bp': total_length,
        'total_length_mb': total_length / 1e6,
        'gc_content': gc_content,
        'total_motifs': len(all_motifs),
        'motifs_per_kb': len(all_motifs) / (total_length / 1000) if total_length > 0 else 0,
        'motifs_per_mb': len(all_motifs) / (total_length / 1e6) if total_length > 0 else 0,
        'class_counts': dict(class_counts),
        'subclass_counts': dict(subclass_counts),
        'detection_time_s': detection_time,
        'bp_per_second': total_length / detection_time if detection_time > 0 else 0,
        'motifs': all_motifs
    }
    
    # Summary
    print(f"\n--- Results Summary ---")
    print(f"Total motifs: {len(all_motifs):,}")
    print(f"Motif density: {results['motifs_per_kb']:.2f} motifs/kb")
    print(f"Detection time: {detection_time:.1f}s ({results['bp_per_second']:.0f} bp/s)")
    print(f"\nMotif classes found:")
    for cls, count in sorted(class_counts.items(), key=lambda x: -x[1]):
        print(f"  {cls}: {count:,}")
    
    return results


def run_comparative_analysis(output_dir='Genomes/analysis_results'):
    """
    Run comparative analysis on all genomes in the Genomes folder.
    """
    genomes_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Genomes')
    
    # Create output directory
    output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), output_dir)
    os.makedirs(output_path, exist_ok=True)
    
    # Find all genome files
    genome_files = []
    for f in os.listdir(genomes_dir):
        if f.endswith(('.fna', '.fasta')) and not f.startswith('.'):
            genome_files.append(f)
    
    print(f"\n{'#'*70}")
    print(f"# NON-B DNA COMPARATIVE GENOMICS ANALYSIS")
    print(f"# NonBDNAFinder v{nbs.__version__}")
    print(f"# Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"# Genomes to analyze: {len(genome_files)}")
    print(f"{'#'*70}")
    
    # Analyze each genome
    all_results = []
    for genome_file in sorted(genome_files):
        filepath = os.path.join(genomes_dir, genome_file)
        try:
            results = analyze_genome(filepath, genome_file)
            all_results.append(results)
        except Exception as e:
            print(f"ERROR analyzing {genome_file}: {e}")
            import traceback
            traceback.print_exc()
    
    # Generate comparative statistics CSV
    print(f"\n{'='*70}")
    print("Generating comparative statistics...")
    print(f"{'='*70}")
    
    stats_data = []
    for r in all_results:
        row = {
            'Organism': r['organism'],
            'Strain': r['strain'],
            'Category': r['category'],
            'Phylum': r['phylum'],
            'Genome_Size_bp': r['total_length_bp'],
            'Genome_Size_Mb': round(r['total_length_mb'], 2),
            'GC_Content_pct': round(r['gc_content'], 1),
            'Num_Sequences': r['num_sequences'],
            'Total_Motifs': r['total_motifs'],
            'Motifs_per_Kb': round(r['motifs_per_kb'], 3),
            'Motifs_per_Mb': round(r['motifs_per_mb'], 0),
            'Detection_Time_s': round(r['detection_time_s'], 1),
            'Throughput_bp_s': round(r['bp_per_second'], 0)
        }
        # Add class-specific counts
        for cls in ['Curved_DNA', 'Slipped_DNA', 'Cruciform', 'R-Loop', 'Triplex', 
                    'G-Quadruplex', 'i-Motif', 'Z-DNA', 'A-philic_DNA', 'Hybrid', 'Non-B_DNA_Clusters']:
            row[cls] = r['class_counts'].get(cls, 0)
        stats_data.append(row)
    
    # Create DataFrame and save
    stats_df = pd.DataFrame(stats_data)
    stats_csv = os.path.join(output_path, 'genome_statistics.csv')
    stats_df.to_csv(stats_csv, index=False)
    print(f"Saved: {stats_csv}")
    
    # Save detailed results as JSON (excluding raw motifs to keep file manageable)
    summary_results = []
    for r in all_results:
        summary = {k: v for k, v in r.items() if k != 'motifs'}
        summary_results.append(summary)
    
    summary_json = os.path.join(output_path, 'analysis_summary.json')
    with open(summary_json, 'w') as f:
        json.dump(summary_results, f, indent=2)
    print(f"Saved: {summary_json}")
    
    # Save individual genome motif files
    for r in all_results:
        organism_name = r['organism'].replace(' ', '_').replace('.', '')
        motifs_json = os.path.join(output_path, f'{organism_name}_motifs.json')
        with open(motifs_json, 'w') as f:
            json.dump({
                'organism': r['organism'],
                'total_motifs': r['total_motifs'],
                'class_counts': r['class_counts'],
                'motifs': r['motifs']
            }, f)
        print(f"Saved: {motifs_json}")
    
    # Print final comparative summary
    print(f"\n{'='*70}")
    print("COMPARATIVE ANALYSIS SUMMARY")
    print(f"{'='*70}")
    print(f"\n{'Organism':<35} {'Size (Mb)':<10} {'GC%':<8} {'Motifs':<10} {'Motifs/kb':<12}")
    print("-" * 75)
    for r in sorted(all_results, key=lambda x: -x['total_motifs']):
        print(f"{r['organism']:<35} {r['total_length_mb']:<10.2f} {r['gc_content']:<8.1f} {r['total_motifs']:<10,} {r['motifs_per_kb']:<12.3f}")
    
    print(f"\n{'='*70}")
    print("Analysis complete!")
    print(f"{'='*70}")
    
    return all_results


if __name__ == '__main__':
    results = run_comparative_analysis()
