#!/usr/bin/env python3
"""
Comprehensive NBST Validation Analysis
Compares NonBDNAFinder results with NBST (Non-B GFA) results on standardized test sequence.
"""

import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys

# Add current directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

import nonbscanner as nbs
from utilities import parse_fasta

# Set up plotting style
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.linewidth'] = 1.5

# Output directories
VALIDATION_DIR = Path("Genomes/NBSTVALIDATION")
OUTPUT_DIR = Path("Genomes/validation_results")
OUTPUT_DIR.mkdir(exist_ok=True)

def load_nbst_results():
    """Load NBST results from TSV files."""
    nbst_results = {}
    
    # G-Quadruplex results
    gq_file = VALIDATION_DIR / "693fc40d26a53_GQ.tsv"
    if gq_file.exists():
        df = pd.read_csv(gq_file, sep='\t')
        nbst_results['G-Quadruplex'] = df
        print(f"Loaded {len(df)} NBST G-Quadruplex motifs")
    
    # Z-DNA results
    z_file = VALIDATION_DIR / "693fc40d26a53_Z.tsv"
    if z_file.exists():
        df = pd.read_csv(z_file, sep='\t')
        nbst_results['Z-DNA'] = df
        print(f"Loaded {len(df)} NBST Z-DNA motifs")
    
    # Mirror Repeats (Triplex) results
    mr_file = VALIDATION_DIR / "693fc40d26a53_MR.tsv"
    if mr_file.exists():
        df = pd.read_csv(mr_file, sep='\t')
        nbst_results['Mirror_Repeats'] = df
        print(f"Loaded {len(df)} NBST Mirror Repeat motifs")
    
    # STR results
    str_file = VALIDATION_DIR / "693fc40d26a53_Slipped_STR.tsv"
    if str_file.exists():
        df = pd.read_csv(str_file, sep='\t')
        nbst_results['STR'] = df
        print(f"Loaded {len(df)} NBST STR motifs")
    
    # Curved DNA results
    curved_file = VALIDATION_DIR / "693fc40d26a53_curved.tsv"
    if curved_file.exists():
        df = pd.read_csv(curved_file, sep='\t')
        nbst_results['Curved'] = df
        print(f"Loaded {len(df)} NBST Curved DNA motifs")
    
    # Direct Repeats results
    dr_file = VALIDATION_DIR / "693fc40d26a53_slipped_DR.tsv"
    if dr_file.exists():
        df = pd.read_csv(dr_file, sep='\t')
        nbst_results['Direct_Repeats'] = df
        print(f"Loaded {len(df)} NBST Direct Repeat motifs")
    
    return nbst_results


def run_nonbdnafinder_analysis():
    """Run NonBDNAFinder on the validation sequence."""
    # Load sequence
    fasta_file = VALIDATION_DIR / "693fc40d26a53.fasta"
    
    with open(fasta_file, 'r') as f:
        fasta_content = f.read()
    
    # parse_fasta returns a dict by default, or generator with streaming=True
    sequences_dict = parse_fasta(fasta_content)
    seq_name = list(sequences_dict.keys())[0]
    sequence = sequences_dict[seq_name]
    
    print(f"\nAnalyzing sequence: {seq_name}")
    print(f"Sequence length: {len(sequence)} bp")
    
    # Run analysis
    motifs = nbs.analyze_sequence(sequence, seq_name)
    
    print(f"\nNonBDNAFinder detected {len(motifs)} total motifs")
    
    # Convert to DataFrame
    df = pd.DataFrame(motifs)
    
    # Print summary by class
    print("\nMotif counts by class:")
    class_counts = df['Class'].value_counts()
    for cls, count in class_counts.items():
        print(f"  {cls}: {count}")
    
    # Save detailed results
    output_file = OUTPUT_DIR / "nonbdnafinder_validation_results.json"
    with open(output_file, 'w') as f:
        json.dump(motifs, f, indent=2)
    print(f"\nSaved detailed results to {output_file}")
    
    # Save as CSV
    csv_file = OUTPUT_DIR / "nonbdnafinder_validation_results.csv"
    df.to_csv(csv_file, index=False)
    print(f"Saved CSV results to {csv_file}")
    
    return df, motifs


def compare_detections(nbst_results, nbf_df):
    """Compare NBST and NonBDNAFinder detections."""
    
    comparison_data = []
    
    # G-Quadruplex comparison
    nbst_gq = len(nbst_results.get('G-Quadruplex', []))
    nbf_gq = len(nbf_df[nbf_df['Class'] == 'G-Quadruplex'])
    comparison_data.append({
        'Motif_Class': 'G-Quadruplex',
        'NBST_Count': nbst_gq,
        'NonBDNAFinder_Count': nbf_gq,
        'Fold_Difference': nbf_gq / nbst_gq if nbst_gq > 0 else np.inf,
        'Detection_Method_NBST': 'G-island scanning',
        'Detection_Method_NBF': 'G4Hunter scoring'
    })
    
    # Z-DNA comparison
    nbst_z = len(nbst_results.get('Z-DNA', []))
    nbf_z = len(nbf_df[nbf_df['Class'] == 'Z-DNA'])
    comparison_data.append({
        'Motif_Class': 'Z-DNA',
        'NBST_Count': nbst_z,
        'NonBDNAFinder_Count': nbf_z,
        'Fold_Difference': nbf_z / nbst_z if nbst_z > 0 else np.inf,
        'Detection_Method_NBST': 'KV-score threshold',
        'Detection_Method_NBF': 'Thermodynamic modeling'
    })
    
    # Mirror Repeats/Triplex comparison
    nbst_mr = len(nbst_results.get('Mirror_Repeats', []))
    nbf_triplex = len(nbf_df[nbf_df['Class'] == 'Triplex'])
    comparison_data.append({
        'Motif_Class': 'Mirror_Repeats/Triplex',
        'NBST_Count': nbst_mr,
        'NonBDNAFinder_Count': nbf_triplex,
        'Fold_Difference': nbf_triplex / nbst_mr if nbst_mr > 0 else np.inf,
        'Detection_Method_NBST': 'Palindrome matching',
        'Detection_Method_NBF': 'Hoogsteen bond prediction'
    })
    
    # STR/Slipped DNA comparison
    nbst_str = len(nbst_results.get('STR', []))
    nbf_slipped = len(nbf_df[nbf_df['Class'] == 'Slipped_DNA'])
    comparison_data.append({
        'Motif_Class': 'STR/Slipped_DNA',
        'NBST_Count': nbst_str,
        'NonBDNAFinder_Count': nbf_slipped,
        'Fold_Difference': nbf_slipped / nbst_str if nbst_str > 0 else np.inf,
        'Detection_Method_NBST': 'Repeat unit scanning',
        'Detection_Method_NBF': 'Hairpin-forming subset'
    })
    
    # Curved DNA comparison
    nbst_curved = len(nbst_results.get('Curved', []))
    nbf_curved = len(nbf_df[nbf_df['Class'] == 'Curved_DNA'])
    comparison_data.append({
        'Motif_Class': 'Curved/A-Phased_DNA',
        'NBST_Count': nbst_curved,
        'NonBDNAFinder_Count': nbf_curved,
        'Fold_Difference': nbf_curved / nbst_curved if nbst_curved > 0 else np.inf,
        'Detection_Method_NBST': 'A-tract phasing only',
        'Detection_Method_NBF': 'Global+local curvature'
    })
    
    # Novel classes in NonBDNAFinder
    novel_classes = ['R-Loop', 'i-Motif', 'A-philic_DNA', 'Hybrid', 'Non-B_DNA_Clusters']
    for cls in novel_classes:
        nbf_count = len(nbf_df[nbf_df['Class'] == cls])
        comparison_data.append({
            'Motif_Class': cls,
            'NBST_Count': 0,
            'NonBDNAFinder_Count': nbf_count,
            'Fold_Difference': np.inf if nbf_count > 0 else 0,
            'Detection_Method_NBST': 'Not detected',
            'Detection_Method_NBF': 'Novel detection'
        })
    
    comparison_df = pd.DataFrame(comparison_data)
    
    # Save comparison table
    output_file = OUTPUT_DIR / "tool_comparison.csv"
    comparison_df.to_csv(output_file, index=False)
    print(f"\nSaved comparison table to {output_file}")
    
    return comparison_df


def analyze_position_concordance(nbst_results, nbf_df):
    """Analyze position concordance between tools."""
    
    concordance_results = []
    
    # G-Quadruplex concordance
    if 'G-Quadruplex' in nbst_results:
        nbst_gq = nbst_results['G-Quadruplex']
        nbf_gq = nbf_df[nbf_df['Class'] == 'G-Quadruplex']
        
        for _, nbst_row in nbst_gq.iterrows():
            nbst_start = nbst_row['Start']
            nbst_end = nbst_row['Stop']  # NBST uses 'Stop' not 'End'
            
            # Find overlapping NBF detections (within 50 bp tolerance)
            overlaps = nbf_gq[
                ((nbf_gq['Start'] >= nbst_start - 50) & (nbf_gq['Start'] <= nbst_end + 50)) |
                ((nbf_gq['End'] >= nbst_start - 50) & (nbf_gq['End'] <= nbst_end + 50))
            ]
            
            concordance_results.append({
                'Motif_Class': 'G-Quadruplex',
                'NBST_Start': nbst_start,
                'NBST_End': nbst_end,
                'NBF_Overlaps': len(overlaps),
                'Concordant': len(overlaps) > 0
            })
    
    # Z-DNA concordance
    if 'Z-DNA' in nbst_results:
        nbst_z = nbst_results['Z-DNA']
        nbf_z = nbf_df[nbf_df['Class'] == 'Z-DNA']
        
        for _, nbst_row in nbst_z.iterrows():
            nbst_start = nbst_row['Start']
            nbst_end = nbst_row['Stop']  # NBST uses 'Stop' not 'End'
            
            # Find overlapping NBF detections (within 20 bp tolerance for Z-DNA)
            overlaps = nbf_z[
                ((nbf_z['Start'] >= nbst_start - 20) & (nbf_z['Start'] <= nbst_end + 20)) |
                ((nbf_z['End'] >= nbst_start - 20) & (nbf_z['End'] <= nbst_end + 20))
            ]
            
            concordance_results.append({
                'Motif_Class': 'Z-DNA',
                'NBST_Start': nbst_start,
                'NBST_End': nbst_end,
                'NBF_Overlaps': len(overlaps),
                'Concordant': len(overlaps) > 0
            })
    
    concordance_df = pd.DataFrame(concordance_results)
    
    if len(concordance_df) > 0:
        # Calculate concordance statistics
        concordance_stats = concordance_df.groupby('Motif_Class').agg({
            'Concordant': ['sum', 'count', 'mean']
        }).round(3)
        
        print("\nPosition Concordance Statistics:")
        print(concordance_stats)
        
        # Save concordance data
        output_file = OUTPUT_DIR / "position_concordance.csv"
        concordance_df.to_csv(output_file, index=False)
        print(f"\nSaved concordance data to {output_file}")
    
    return concordance_df


def create_validation_visualizations(comparison_df, nbf_df):
    """Create comprehensive validation visualizations."""
    
    # 1. Comparison bar chart
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Select main motif classes for visualization
    main_classes = comparison_df[comparison_df['NBST_Count'] > 0]
    
    x = np.arange(len(main_classes))
    width = 0.35
    
    ax.bar(x - width/2, main_classes['NBST_Count'], width, label='NBST', color='#0173B2', edgecolor='black', linewidth=1.5)
    ax.bar(x + width/2, main_classes['NonBDNAFinder_Count'], width, label='NonBDNAFinder', color='#DE8F05', edgecolor='black', linewidth=1.5)
    
    ax.set_xlabel('Motif Class', fontsize=12, fontweight='bold')
    ax.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax.set_title('Motif Detection Comparison: NBST vs NonBDNAFinder', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(main_classes['Motif_Class'], rotation=45, ha='right')
    ax.legend(frameon=True, fontsize=11)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "Figure_V1_Comparison_Bar_Chart.png", dpi=300, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / "Figure_V1_Comparison_Bar_Chart.pdf", bbox_inches='tight')
    print("\nSaved Figure V1: Comparison Bar Chart")
    plt.close()
    
    # 2. Novel motif classes pie chart
    fig, ax = plt.subplots(figsize=(10, 8))
    
    novel_classes = comparison_df[comparison_df['NBST_Count'] == 0]
    novel_classes = novel_classes[novel_classes['NonBDNAFinder_Count'] > 0]
    
    colors = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#CC79A7']
    wedges, texts, autotexts = ax.pie(
        novel_classes['NonBDNAFinder_Count'],
        labels=novel_classes['Motif_Class'],
        autopct='%1.1f%%',
        colors=colors[:len(novel_classes)],
        startangle=90,
        textprops={'fontsize': 11, 'fontweight': 'bold'}
    )
    
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
    
    ax.set_title('Novel Motif Classes Detected Only by NonBDNAFinder', fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "Figure_V2_Novel_Classes.png", dpi=300, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / "Figure_V2_Novel_Classes.pdf", bbox_inches='tight')
    print("Saved Figure V2: Novel Classes Pie Chart")
    plt.close()
    
    # 3. Subclass distribution for G-Quadruplex
    fig, ax = plt.subplots(figsize=(14, 6))
    
    gq_motifs = nbf_df[nbf_df['Class'] == 'G-Quadruplex']
    subclass_counts = gq_motifs['Subclass'].value_counts()
    
    colors_subclass = sns.color_palette("husl", len(subclass_counts))
    bars = ax.bar(range(len(subclass_counts)), subclass_counts.values, color=colors_subclass, edgecolor='black', linewidth=1.5)
    
    ax.set_xlabel('G-Quadruplex Subclass', fontsize=12, fontweight='bold')
    ax.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax.set_title('G-Quadruplex Subclass Distribution (NonBDNAFinder)', fontsize=14, fontweight='bold')
    ax.set_xticks(range(len(subclass_counts)))
    ax.set_xticklabels(subclass_counts.index, rotation=45, ha='right', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    
    # Add count labels on bars
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{int(height)}',
                ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "Figure_V3_GQ_Subclass_Distribution.png", dpi=300, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / "Figure_V3_GQ_Subclass_Distribution.pdf", bbox_inches='tight')
    print("Saved Figure V3: G-Quadruplex Subclass Distribution")
    plt.close()
    
    # 4. Genome track visualization
    fig, axes = plt.subplots(6, 1, figsize=(16, 12), sharex=True)
    
    # Track 1: G-Quadruplex
    gq_motifs = nbf_df[nbf_df['Class'] == 'G-Quadruplex']
    for _, motif in gq_motifs.iterrows():
        axes[0].plot([motif['Start'], motif['End']], [0, 0], 'r-', linewidth=2, alpha=0.7)
    axes[0].set_ylabel('G4', fontweight='bold', fontsize=10)
    axes[0].set_ylim(-0.5, 0.5)
    axes[0].set_yticks([])
    axes[0].spines['left'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[0].spines['top'].set_visible(False)
    axes[0].set_title('NonBDNAFinder Motif Tracks on Validation Sequence (40,100 bp)', fontsize=14, fontweight='bold', pad=15)
    
    # Track 2: Z-DNA
    z_motifs = nbf_df[nbf_df['Class'] == 'Z-DNA']
    for _, motif in z_motifs.iterrows():
        axes[1].plot([motif['Start'], motif['End']], [0, 0], 'b-', linewidth=2, alpha=0.7)
    axes[1].set_ylabel('Z-DNA', fontweight='bold', fontsize=10)
    axes[1].set_ylim(-0.5, 0.5)
    axes[1].set_yticks([])
    axes[1].spines['left'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    axes[1].spines['top'].set_visible(False)
    
    # Track 3: Curved DNA
    curved_motifs = nbf_df[nbf_df['Class'] == 'Curved_DNA']
    for _, motif in curved_motifs.iterrows():
        axes[2].plot([motif['Start'], motif['End']], [0, 0], 'g-', linewidth=2, alpha=0.7)
    axes[2].set_ylabel('Curved', fontweight='bold', fontsize=10)
    axes[2].set_ylim(-0.5, 0.5)
    axes[2].set_yticks([])
    axes[2].spines['left'].set_visible(False)
    axes[2].spines['right'].set_visible(False)
    axes[2].spines['top'].set_visible(False)
    
    # Track 4: R-loops
    rloop_motifs = nbf_df[nbf_df['Class'] == 'R-Loop']
    for _, motif in rloop_motifs.iterrows():
        axes[3].plot([motif['Start'], motif['End']], [0, 0], 'm-', linewidth=2, alpha=0.7)
    axes[3].set_ylabel('R-Loop', fontweight='bold', fontsize=10)
    axes[3].set_ylim(-0.5, 0.5)
    axes[3].set_yticks([])
    axes[3].spines['left'].set_visible(False)
    axes[3].spines['right'].set_visible(False)
    axes[3].spines['top'].set_visible(False)
    
    # Track 5: i-Motif
    imotif_motifs = nbf_df[nbf_df['Class'] == 'i-Motif']
    for _, motif in imotif_motifs.iterrows():
        axes[4].plot([motif['Start'], motif['End']], [0, 0], 'c-', linewidth=2, alpha=0.7)
    axes[4].set_ylabel('i-Motif', fontweight='bold', fontsize=10)
    axes[4].set_ylim(-0.5, 0.5)
    axes[4].set_yticks([])
    axes[4].spines['left'].set_visible(False)
    axes[4].spines['right'].set_visible(False)
    axes[4].spines['top'].set_visible(False)
    
    # Track 6: All other classes
    other_classes = ['Triplex', 'Slipped_DNA', 'A-philic_DNA', 'Hybrid', 'Non-B_DNA_Clusters']
    colors_other = ['orange', 'purple', 'brown', 'pink', 'gray']
    for cls, color in zip(other_classes, colors_other):
        cls_motifs = nbf_df[nbf_df['Class'] == cls]
        for _, motif in cls_motifs.iterrows():
            axes[5].plot([motif['Start'], motif['End']], [0, 0], color=color, linewidth=2, alpha=0.7)
    axes[5].set_ylabel('Other', fontweight='bold', fontsize=10)
    axes[5].set_ylim(-0.5, 0.5)
    axes[5].set_yticks([])
    axes[5].spines['left'].set_visible(False)
    axes[5].spines['right'].set_visible(False)
    axes[5].spines['top'].set_visible(False)
    axes[5].set_xlabel('Genomic Position (bp)', fontweight='bold', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "Figure_V4_Genome_Tracks.png", dpi=300, bbox_inches='tight')
    plt.savefig(OUTPUT_DIR / "Figure_V4_Genome_Tracks.pdf", bbox_inches='tight')
    print("Saved Figure V4: Genome Tracks")
    plt.close()
    
    print("\nAll validation visualizations created successfully!")


def main():
    """Main validation workflow."""
    print("="*80)
    print("NBST Validation Analysis - Comprehensive Tool Comparison")
    print("="*80)
    
    # Step 1: Load NBST results
    print("\n1. Loading NBST results...")
    nbst_results = load_nbst_results()
    
    # Step 2: Run NonBDNAFinder analysis
    print("\n2. Running NonBDNAFinder analysis...")
    nbf_df, nbf_motifs = run_nonbdnafinder_analysis()
    
    # Step 3: Compare detections
    print("\n3. Comparing detections...")
    comparison_df = compare_detections(nbst_results, nbf_df)
    
    print("\n" + "="*80)
    print("COMPARISON SUMMARY")
    print("="*80)
    print(comparison_df.to_string(index=False))
    
    # Step 4: Analyze position concordance
    print("\n4. Analyzing position concordance...")
    concordance_df = analyze_position_concordance(nbst_results, nbf_df)
    
    # Step 5: Create visualizations
    print("\n5. Creating validation visualizations...")
    create_validation_visualizations(comparison_df, nbf_df)
    
    print("\n" + "="*80)
    print("VALIDATION COMPLETE!")
    print("="*80)
    print(f"\nResults saved to: {OUTPUT_DIR}")
    print("\nGenerated files:")
    print("  - nonbdnafinder_validation_results.json")
    print("  - nonbdnafinder_validation_results.csv")
    print("  - tool_comparison.csv")
    print("  - position_concordance.csv")
    print("  - Figure_V1_Comparison_Bar_Chart.png/pdf")
    print("  - Figure_V2_Novel_Classes.png/pdf")
    print("  - Figure_V3_GQ_Subclass_Distribution.png/pdf")
    print("  - Figure_V4_Genome_Tracks.png/pdf")


if __name__ == "__main__":
    main()
