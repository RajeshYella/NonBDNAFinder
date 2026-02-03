#!/usr/bin/env python3
"""
Subclass-Level Comparative Genomics Analysis Script for Non-B DNA Motifs

This script extends the class-level analysis to provide detailed subclass-level
statistics with special emphasis on clusters and hybrid motifs.

Author: NonBDNAFinder Team
Date: 2026
"""

import os
import sys
import json
import pandas as pd
from collections import defaultdict
from datetime import datetime

# Add the parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config.motif_taxonomy import MOTIF_CLASSIFICATION, CLASS_TO_SUBCLASSES


def load_analysis_results(results_dir='Genomes/analysis_results'):
    """Load existing analysis results from JSON files."""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    results_path = os.path.join(base_dir, results_dir)
    
    # Load summary
    summary_file = os.path.join(results_path, 'analysis_summary.json')
    with open(summary_file, 'r') as f:
        summary = json.load(f)
    
    return summary


def extract_subclass_statistics(summary):
    """Extract detailed subclass statistics from analysis summary."""
    
    subclass_data = []
    
    for genome in summary:
        organism = genome['organism']
        gc_content = genome['gc_content']
        total_length_mb = genome['total_length_mb']
        phylum = genome['phylum']
        category = genome['category']
        
        subclass_counts = genome.get('subclass_counts', {})
        
        for subclass_key, count in subclass_counts.items():
            # Parse "Class:Subclass" format
            if ':' in subclass_key:
                class_name, subclass_name = subclass_key.split(':', 1)
            else:
                class_name = subclass_key
                subclass_name = subclass_key
            
            # Calculate density
            density_per_mb = count / total_length_mb if total_length_mb > 0 else 0
            
            subclass_data.append({
                'Organism': organism,
                'Phylum': phylum,
                'Category': category,
                'GC_Content_pct': round(gc_content, 1),
                'Genome_Size_Mb': round(total_length_mb, 2),
                'Class': class_name,
                'Subclass': subclass_name,
                'Count': count,
                'Density_per_Mb': round(density_per_mb, 2)
            })
    
    return pd.DataFrame(subclass_data)


def generate_cluster_analysis(summary):
    """Generate detailed analysis of Non-B DNA clusters."""
    
    cluster_data = []
    
    for genome in summary:
        organism = genome['organism']
        gc_content = genome['gc_content']
        total_length_mb = genome['total_length_mb']
        phylum = genome['phylum']
        
        subclass_counts = genome.get('subclass_counts', {})
        
        # Extract cluster subclasses
        for subclass_key, count in subclass_counts.items():
            if 'Cluster' in subclass_key or 'cluster' in subclass_key:
                # Parse cluster type (e.g., "Mixed_Cluster_3_classes")
                if ':' in subclass_key:
                    _, cluster_type = subclass_key.split(':', 1)
                else:
                    cluster_type = subclass_key
                
                density_per_mb = count / total_length_mb if total_length_mb > 0 else 0
                
                cluster_data.append({
                    'Organism': organism,
                    'Phylum': phylum,
                    'GC_Content_pct': round(gc_content, 1),
                    'Cluster_Type': cluster_type,
                    'Count': count,
                    'Density_per_Mb': round(density_per_mb, 2)
                })
    
    return pd.DataFrame(cluster_data)


def generate_hybrid_analysis(summary):
    """Generate detailed analysis of Hybrid motifs (overlaps)."""
    
    hybrid_data = []
    
    for genome in summary:
        organism = genome['organism']
        gc_content = genome['gc_content']
        total_length_mb = genome['total_length_mb']
        phylum = genome['phylum']
        
        subclass_counts = genome.get('subclass_counts', {})
        
        # Extract hybrid subclasses (overlaps)
        for subclass_key, count in subclass_counts.items():
            if 'Overlap' in subclass_key or 'overlap' in subclass_key:
                # Parse hybrid type (e.g., "G-Quadruplex_R-Loop_Overlap")
                if ':' in subclass_key:
                    _, hybrid_type = subclass_key.split(':', 1)
                else:
                    hybrid_type = subclass_key
                
                # Extract the two motif classes involved
                hybrid_type_clean = hybrid_type.replace('_Overlap', '')
                parts = hybrid_type_clean.split('_')
                
                # Reconstruct class names (handle hyphenated names)
                motif_1 = parts[0] if len(parts) > 0 else 'Unknown'
                motif_2 = parts[1] if len(parts) > 1 else 'Unknown'
                
                density_per_mb = count / total_length_mb if total_length_mb > 0 else 0
                
                hybrid_data.append({
                    'Organism': organism,
                    'Phylum': phylum,
                    'GC_Content_pct': round(gc_content, 1),
                    'Hybrid_Type': hybrid_type,
                    'Motif_1': motif_1,
                    'Motif_2': motif_2,
                    'Count': count,
                    'Density_per_Mb': round(density_per_mb, 2)
                })
    
    return pd.DataFrame(hybrid_data)


def generate_gquadruplex_subclass_analysis(summary):
    """Generate detailed G-Quadruplex subclass analysis."""
    
    g4_data = []
    
    # Define G4 subclasses
    g4_subclasses = CLASS_TO_SUBCLASSES.get('G-Quadruplex', [])
    
    for genome in summary:
        organism = genome['organism']
        gc_content = genome['gc_content']
        total_length_mb = genome['total_length_mb']
        phylum = genome['phylum']
        
        subclass_counts = genome.get('subclass_counts', {})
        
        row = {
            'Organism': organism,
            'Phylum': phylum,
            'GC_Content_pct': round(gc_content, 1),
            'Genome_Size_Mb': round(total_length_mb, 2),
            'Total_G4': genome.get('class_counts', {}).get('G-Quadruplex', 0)
        }
        
        # Add each G4 subclass count
        for subclass in g4_subclasses:
            key = f'G-Quadruplex:{subclass}'
            count = subclass_counts.get(key, 0)
            # Create clean column name
            col_name = subclass.replace(' ', '_').replace('/', '_').replace('-', '_')
            row[col_name] = count
            row[f'{col_name}_per_Mb'] = round(count / total_length_mb, 2) if total_length_mb > 0 else 0
        
        g4_data.append(row)
    
    return pd.DataFrame(g4_data)


def generate_imotif_subclass_analysis(summary):
    """Generate detailed i-Motif subclass analysis."""
    
    imotif_data = []
    
    # Define i-Motif subclasses
    imotif_subclasses = CLASS_TO_SUBCLASSES.get('i-Motif', [])
    
    for genome in summary:
        organism = genome['organism']
        gc_content = genome['gc_content']
        total_length_mb = genome['total_length_mb']
        phylum = genome['phylum']
        
        subclass_counts = genome.get('subclass_counts', {})
        
        row = {
            'Organism': organism,
            'Phylum': phylum,
            'GC_Content_pct': round(gc_content, 1),
            'Genome_Size_Mb': round(total_length_mb, 2),
            'Total_iMotif': genome.get('class_counts', {}).get('i-Motif', 0)
        }
        
        # Add each i-Motif subclass count
        for subclass in imotif_subclasses:
            key = f'i-Motif:{subclass}'
            count = subclass_counts.get(key, 0)
            col_name = subclass.replace(' ', '_').replace('-', '_')
            row[col_name] = count
            row[f'{col_name}_per_Mb'] = round(count / total_length_mb, 2) if total_length_mb > 0 else 0
        
        imotif_data.append(row)
    
    return pd.DataFrame(imotif_data)


def generate_zdna_subclass_analysis(summary):
    """Generate detailed Z-DNA subclass analysis."""
    
    zdna_data = []
    
    # Define Z-DNA subclasses
    zdna_subclasses = CLASS_TO_SUBCLASSES.get('Z-DNA', [])
    
    for genome in summary:
        organism = genome['organism']
        gc_content = genome['gc_content']
        total_length_mb = genome['total_length_mb']
        phylum = genome['phylum']
        
        subclass_counts = genome.get('subclass_counts', {})
        
        row = {
            'Organism': organism,
            'Phylum': phylum,
            'GC_Content_pct': round(gc_content, 1),
            'Genome_Size_Mb': round(total_length_mb, 2),
            'Total_ZDNA': genome.get('class_counts', {}).get('Z-DNA', 0)
        }
        
        # Add each Z-DNA subclass count
        for subclass in zdna_subclasses:
            key = f'Z-DNA:{subclass}'
            count = subclass_counts.get(key, 0)
            col_name = subclass.replace(' ', '_').replace('-', '_')
            row[col_name] = count
            row[f'{col_name}_per_Mb'] = round(count / total_length_mb, 2) if total_length_mb > 0 else 0
        
        zdna_data.append(row)
    
    return pd.DataFrame(zdna_data)


def generate_curved_dna_subclass_analysis(summary):
    """Generate detailed Curved DNA subclass analysis."""
    
    curved_data = []
    
    # Define Curved DNA subclasses
    curved_subclasses = CLASS_TO_SUBCLASSES.get('Curved_DNA', [])
    
    for genome in summary:
        organism = genome['organism']
        gc_content = genome['gc_content']
        total_length_mb = genome['total_length_mb']
        phylum = genome['phylum']
        
        subclass_counts = genome.get('subclass_counts', {})
        
        row = {
            'Organism': organism,
            'Phylum': phylum,
            'GC_Content_pct': round(gc_content, 1),
            'Genome_Size_Mb': round(total_length_mb, 2),
            'Total_Curved_DNA': genome.get('class_counts', {}).get('Curved_DNA', 0)
        }
        
        # Add each Curved DNA subclass count
        for subclass in curved_subclasses:
            key = f'Curved_DNA:{subclass}'
            count = subclass_counts.get(key, 0)
            col_name = subclass.replace(' ', '_')
            row[col_name] = count
            row[f'{col_name}_per_Mb'] = round(count / total_length_mb, 2) if total_length_mb > 0 else 0
        
        curved_data.append(row)
    
    return pd.DataFrame(curved_data)


def generate_all_subclass_pivot(summary):
    """Generate a pivot table of all subclasses across all genomes."""
    
    # Collect all unique subclasses
    all_subclasses = set()
    for genome in summary:
        for key in genome.get('subclass_counts', {}).keys():
            all_subclasses.add(key)
    
    # Build pivot data
    pivot_data = []
    for genome in summary:
        organism = genome['organism']
        gc_content = genome['gc_content']
        total_length_mb = genome['total_length_mb']
        
        row = {
            'Organism': organism,
            'GC_Content_pct': round(gc_content, 1),
            'Genome_Size_Mb': round(total_length_mb, 2)
        }
        
        subclass_counts = genome.get('subclass_counts', {})
        
        for subclass_key in sorted(all_subclasses):
            count = subclass_counts.get(subclass_key, 0)
            # Clean column name
            col_name = subclass_key.replace(':', '_').replace(' ', '_').replace('/', '_').replace('-', '_')
            row[col_name] = count
        
        pivot_data.append(row)
    
    return pd.DataFrame(pivot_data)


def run_subclass_analysis(output_dir='Genomes/analysis_results'):
    """Run complete subclass-level analysis."""
    
    print(f"\n{'#'*70}")
    print("# SUBCLASS-LEVEL COMPARATIVE GENOMICS ANALYSIS")
    print(f"# Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"{'#'*70}")
    
    # Load existing analysis results
    print("\nLoading analysis results...")
    summary = load_analysis_results()
    print(f"Loaded data for {len(summary)} genomes")
    
    base_dir = os.path.dirname(os.path.abspath(__file__))
    output_path = os.path.join(base_dir, output_dir)
    os.makedirs(output_path, exist_ok=True)
    
    # Generate all subclass statistics
    print("\n" + "="*70)
    print("Generating subclass-level statistics...")
    print("="*70)
    
    # 1. Complete subclass statistics
    print("\n1. Extracting all subclass statistics...")
    subclass_df = extract_subclass_statistics(summary)
    subclass_csv = os.path.join(output_path, 'subclass_statistics.csv')
    subclass_df.to_csv(subclass_csv, index=False)
    print(f"   Saved: {subclass_csv}")
    print(f"   Total unique subclass entries: {len(subclass_df)}")
    
    # 2. Cluster analysis
    print("\n2. Generating cluster analysis...")
    cluster_df = generate_cluster_analysis(summary)
    cluster_csv = os.path.join(output_path, 'cluster_analysis.csv')
    cluster_df.to_csv(cluster_csv, index=False)
    print(f"   Saved: {cluster_csv}")
    print(f"   Total cluster entries: {len(cluster_df)}")
    
    # 3. Hybrid (overlap) analysis
    print("\n3. Generating hybrid motif analysis...")
    hybrid_df = generate_hybrid_analysis(summary)
    hybrid_csv = os.path.join(output_path, 'hybrid_analysis.csv')
    hybrid_df.to_csv(hybrid_csv, index=False)
    print(f"   Saved: {hybrid_csv}")
    print(f"   Total hybrid entries: {len(hybrid_df)}")
    
    # 4. G-Quadruplex subclass analysis
    print("\n4. Generating G-Quadruplex subclass analysis...")
    g4_df = generate_gquadruplex_subclass_analysis(summary)
    g4_csv = os.path.join(output_path, 'gquadruplex_subclass_analysis.csv')
    g4_df.to_csv(g4_csv, index=False)
    print(f"   Saved: {g4_csv}")
    
    # 5. i-Motif subclass analysis
    print("\n5. Generating i-Motif subclass analysis...")
    imotif_df = generate_imotif_subclass_analysis(summary)
    imotif_csv = os.path.join(output_path, 'imotif_subclass_analysis.csv')
    imotif_df.to_csv(imotif_csv, index=False)
    print(f"   Saved: {imotif_csv}")
    
    # 6. Z-DNA subclass analysis
    print("\n6. Generating Z-DNA subclass analysis...")
    zdna_df = generate_zdna_subclass_analysis(summary)
    zdna_csv = os.path.join(output_path, 'zdna_subclass_analysis.csv')
    zdna_df.to_csv(zdna_csv, index=False)
    print(f"   Saved: {zdna_csv}")
    
    # 7. Curved DNA subclass analysis
    print("\n7. Generating Curved DNA subclass analysis...")
    curved_df = generate_curved_dna_subclass_analysis(summary)
    curved_csv = os.path.join(output_path, 'curved_dna_subclass_analysis.csv')
    curved_df.to_csv(curved_csv, index=False)
    print(f"   Saved: {curved_csv}")
    
    # 8. Complete subclass pivot table
    print("\n8. Generating complete subclass pivot table...")
    pivot_df = generate_all_subclass_pivot(summary)
    pivot_csv = os.path.join(output_path, 'subclass_pivot_table.csv')
    pivot_df.to_csv(pivot_csv, index=False)
    print(f"   Saved: {pivot_csv}")
    
    # Print summary statistics
    print("\n" + "="*70)
    print("SUBCLASS ANALYSIS SUMMARY")
    print("="*70)
    
    # Cluster summary
    print("\n--- CLUSTER SUMMARY ---")
    if not cluster_df.empty:
        cluster_summary = cluster_df.groupby('Cluster_Type')['Count'].sum().sort_values(ascending=False)
        print(f"\nCluster types across all genomes:")
        for cluster_type, count in cluster_summary.items():
            print(f"  {cluster_type}: {count:,}")
    
    # Hybrid summary
    print("\n--- HYBRID SUMMARY (TOP 15) ---")
    if not hybrid_df.empty:
        hybrid_summary = hybrid_df.groupby('Hybrid_Type')['Count'].sum().sort_values(ascending=False).head(15)
        print(f"\nTop 15 hybrid types across all genomes:")
        for hybrid_type, count in hybrid_summary.items():
            print(f"  {hybrid_type}: {count:,}")
    
    # G4 subclass summary
    print("\n--- G-QUADRUPLEX SUBCLASS SUMMARY ---")
    if not g4_df.empty:
        g4_cols = [c for c in g4_df.columns if c not in ['Organism', 'Phylum', 'GC_Content_pct', 'Genome_Size_Mb', 'Total_G4'] and not c.endswith('_per_Mb')]
        g4_totals = g4_df[g4_cols].sum().sort_values(ascending=False)
        print(f"\nG4 subclasses across all genomes:")
        for subclass, count in g4_totals.items():
            if count > 0:
                print(f"  {subclass}: {int(count):,}")
    
    print(f"\n{'='*70}")
    print("Subclass analysis complete!")
    print(f"{'='*70}")
    
    return {
        'subclass_df': subclass_df,
        'cluster_df': cluster_df,
        'hybrid_df': hybrid_df,
        'g4_df': g4_df,
        'imotif_df': imotif_df,
        'zdna_df': zdna_df,
        'curved_df': curved_df,
        'pivot_df': pivot_df
    }


if __name__ == '__main__':
    results = run_subclass_analysis()
