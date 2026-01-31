#!/usr/bin/env python3
"""
Subclass-Level Visualization Script for Non-B DNA Comparative Genomics

This script generates publication-quality figures for subclass-level analysis
with emphasis on clusters, hybrids, and individual motif class breakdowns.

Author: NonBDNAFinder Team
Date: 2026
"""

import os
import sys
import json
import pandas as pd
import numpy as np

# Add the parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.colors import LinearSegmentedColormap
    import seaborn as sns
    PLOTTING_AVAILABLE = True
except ImportError:
    PLOTTING_AVAILABLE = False
    print("WARNING: matplotlib/seaborn not available. Install with: pip install matplotlib seaborn")


# Colorblind-friendly palette (Wong 2011)
COLORS = {
    'orange': '#E69F00',
    'sky_blue': '#56B4E9',
    'bluish_green': '#009E73',
    'yellow': '#F0E442',
    'blue': '#0072B2',
    'vermillion': '#D55E00',
    'reddish_purple': '#CC79A7',
    'black': '#000000',
    'gray': '#999999',
    'olive': '#A6761D',
    'wine': '#882255',
    'teal': '#44AA99'
}

# G4 subclass colors
G4_COLORS = {
    'Two-tetrad weak PQS': COLORS['blue'],
    'Intramolecular G-triplex': COLORS['vermillion'],
    'Extended-loop canonical': COLORS['bluish_green'],
    'Canonical intramolecular G4': COLORS['reddish_purple'],
    'Higher-order G4 array/G4-wire': COLORS['orange'],
    'Stacked G4s with linker': COLORS['sky_blue'],
    'Stacked canonical G4s': COLORS['yellow'],
    'Telomeric G4': COLORS['olive']
}

# Cluster colors
CLUSTER_COLORS = {
    'Mixed_Cluster_3_classes': COLORS['blue'],
    'Mixed_Cluster_4_classes': COLORS['vermillion'],
    'Mixed_Cluster_5_classes': COLORS['bluish_green'],
    'Mixed_Cluster_6_classes': COLORS['reddish_purple']
}


def set_publication_style():
    """Set matplotlib style for publication-quality figures."""
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 10,
        'axes.labelsize': 12,
        'axes.titlesize': 14,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 9,
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'axes.linewidth': 0.8,
        'axes.grid': True,
        'grid.alpha': 0.3,
        'grid.linestyle': '--'
    })


def load_analysis_data(results_dir='Genomes/analysis_results'):
    """Load all analysis data from CSV files."""
    base_dir = os.path.dirname(os.path.abspath(__file__))
    results_path = os.path.join(base_dir, results_dir)
    
    data = {}
    
    # Load cluster analysis
    cluster_file = os.path.join(results_path, 'cluster_analysis.csv')
    if os.path.exists(cluster_file):
        data['clusters'] = pd.read_csv(cluster_file)
    
    # Load hybrid analysis
    hybrid_file = os.path.join(results_path, 'hybrid_analysis.csv')
    if os.path.exists(hybrid_file):
        data['hybrids'] = pd.read_csv(hybrid_file)
    
    # Load G4 subclass analysis
    g4_file = os.path.join(results_path, 'gquadruplex_subclass_analysis.csv')
    if os.path.exists(g4_file):
        data['g4'] = pd.read_csv(g4_file)
    
    # Load genome statistics
    stats_file = os.path.join(results_path, 'genome_statistics.csv')
    if os.path.exists(stats_file):
        data['stats'] = pd.read_csv(stats_file)
    
    # Load curved DNA analysis
    curved_file = os.path.join(results_path, 'curved_dna_subclass_analysis.csv')
    if os.path.exists(curved_file):
        data['curved'] = pd.read_csv(curved_file)
    
    # Load Z-DNA analysis
    zdna_file = os.path.join(results_path, 'zdna_subclass_analysis.csv')
    if os.path.exists(zdna_file):
        data['zdna'] = pd.read_csv(zdna_file)
    
    return data


def plot_g4_subclass_distribution(data, output_dir):
    """Generate G4 subclass stacked bar chart."""
    if 'g4' not in data:
        print("No G4 data available")
        return
    
    g4_df = data['g4'].copy()
    
    # Sort by GC content
    g4_df = g4_df.sort_values('GC_Content_pct')
    
    # Get subclass columns (exclude metadata)
    subclass_cols = [c for c in g4_df.columns 
                     if c not in ['Organism', 'Phylum', 'GC_Content_pct', 'Genome_Size_Mb', 'Total_G4'] 
                     and not c.endswith('_per_Mb')]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # Create stacked bar chart
    organisms = g4_df['Organism'].tolist()
    x = np.arange(len(organisms))
    width = 0.7
    
    # Define colors for subclasses
    colors = list(G4_COLORS.values())[:len(subclass_cols)]
    
    bottom = np.zeros(len(organisms))
    bars = []
    
    for i, col in enumerate(subclass_cols):
        values = g4_df[col].values
        bar = ax.bar(x, values, width, bottom=bottom, label=col.replace('_', ' '), 
                    color=colors[i % len(colors)], edgecolor='white', linewidth=0.5)
        bars.append(bar)
        bottom += values
    
    # Customize plot
    ax.set_xlabel('Organism (ordered by GC content)')
    ax.set_ylabel('G-Quadruplex Count')
    ax.set_title('G-Quadruplex Subclass Distribution Across Bacterial Genomes')
    ax.set_xticks(x)
    ax.set_xticklabels([f"{org}\n({gc:.1f}% GC)" for org, gc in 
                        zip(g4_df['Organism'], g4_df['GC_Content_pct'])], 
                       rotation=45, ha='right')
    
    # Add legend
    ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1), frameon=False)
    
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.join(output_dir, 'Figure7_G4_Subclass_Distribution')
    plt.savefig(f'{output_path}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_path}.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"Saved: {output_path}.png and .pdf")


def plot_cluster_analysis(data, output_dir):
    """Generate cluster complexity analysis figure."""
    if 'clusters' not in data:
        print("No cluster data available")
        return
    
    cluster_df = data['clusters'].copy()
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Panel A: Cluster types stacked by organism
    ax1 = axes[0]
    
    # Pivot data
    pivot = cluster_df.pivot_table(
        index='Organism', 
        columns='Cluster_Type', 
        values='Count', 
        fill_value=0
    )
    
    # Sort by total clusters
    pivot['Total'] = pivot.sum(axis=1)
    pivot = pivot.sort_values('Total', ascending=True)
    pivot = pivot.drop('Total', axis=1)
    
    # Create horizontal stacked bar
    cluster_types = ['Mixed_Cluster_3_classes', 'Mixed_Cluster_4_classes', 
                     'Mixed_Cluster_5_classes', 'Mixed_Cluster_6_classes']
    cluster_types = [c for c in cluster_types if c in pivot.columns]
    
    y = np.arange(len(pivot))
    left = np.zeros(len(pivot))
    
    for i, cluster_type in enumerate(cluster_types):
        if cluster_type in pivot.columns:
            values = pivot[cluster_type].values
            ax1.barh(y, values, left=left, label=cluster_type.replace('Mixed_Cluster_', '').replace('_', ' '),
                    color=list(CLUSTER_COLORS.values())[i])
            left += values
    
    ax1.set_yticks(y)
    ax1.set_yticklabels(pivot.index)
    ax1.set_xlabel('Cluster Count')
    ax1.set_title('A. Cluster Distribution by Organism')
    ax1.legend(title='Cluster Type', loc='lower right')
    
    # Panel B: Cluster density vs GC content
    ax2 = axes[1]
    
    if 'stats' in data:
        stats_df = data['stats']
        
        for _, row in stats_df.iterrows():
            org = row['Organism']
            gc = row['GC_Content_pct']
            cluster_count = row.get('Non-B_DNA_Clusters', 0)
            genome_size = row['Genome_Size_Mb']
            density = cluster_count / genome_size if genome_size > 0 else 0
            
            phylum = row.get('Phylum', 'Unknown')
            color = COLORS['blue'] if phylum == 'Proteobacteria' else (
                COLORS['vermillion'] if phylum == 'Actinobacteria' else COLORS['bluish_green'])
            
            ax2.scatter(gc, density, s=100, c=color, edgecolor='black', linewidth=0.5)
            ax2.annotate(org.split()[0][:3], (gc, density), fontsize=8, 
                        xytext=(5, 5), textcoords='offset points')
    
    ax2.set_xlabel('GC Content (%)')
    ax2.set_ylabel('Cluster Density (per Mb)')
    ax2.set_title('B. Cluster Density vs GC Content')
    
    # Panel C: Cluster complexity pie chart
    ax3 = axes[2]
    
    cluster_totals = cluster_df.groupby('Cluster_Type')['Count'].sum().sort_values(ascending=False)
    
    colors = [CLUSTER_COLORS.get(ct, COLORS['gray']) for ct in cluster_totals.index]
    labels = [ct.replace('Mixed_Cluster_', '').replace('_', ' ') for ct in cluster_totals.index]
    
    wedges, texts, autotexts = ax3.pie(
        cluster_totals.values, 
        labels=labels,
        autopct='%1.1f%%',
        colors=colors,
        explode=[0.02] * len(cluster_totals),
        startangle=90
    )
    
    ax3.set_title('C. Overall Cluster Composition')
    
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.join(output_dir, 'Figure8_Cluster_Analysis')
    plt.savefig(f'{output_path}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_path}.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"Saved: {output_path}.png and .pdf")


def plot_hybrid_analysis(data, output_dir):
    """Generate hybrid motif analysis figure."""
    if 'hybrids' not in data:
        print("No hybrid data available")
        return
    
    hybrid_df = data['hybrids'].copy()
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    # Panel A: Top hybrid types bar chart
    ax1 = axes[0]
    
    top_hybrids = hybrid_df.groupby('Hybrid_Type')['Count'].sum().sort_values(ascending=True).tail(15)
    
    y = np.arange(len(top_hybrids))
    ax1.barh(y, top_hybrids.values, color=COLORS['reddish_purple'])
    ax1.set_yticks(y)
    ax1.set_yticklabels([h.replace('_Overlap', '').replace('_', '-') for h in top_hybrids.index], fontsize=8)
    ax1.set_xlabel('Total Count')
    ax1.set_title('A. Top 15 Hybrid Motif Types')
    
    # Panel B: Hybrid density vs GC content
    ax2 = axes[1]
    
    if 'stats' in data:
        stats_df = data['stats']
        
        for _, row in stats_df.iterrows():
            org = row['Organism']
            gc = row['GC_Content_pct']
            hybrid_count = row.get('Hybrid', 0)
            genome_size = row['Genome_Size_Mb']
            density = hybrid_count / genome_size if genome_size > 0 else 0
            
            phylum = row.get('Phylum', 'Unknown')
            color = COLORS['blue'] if phylum == 'Proteobacteria' else (
                COLORS['vermillion'] if phylum == 'Actinobacteria' else COLORS['bluish_green'])
            
            ax2.scatter(gc, density, s=100, c=color, edgecolor='black', linewidth=0.5)
            ax2.annotate(org.split()[0][:4], (gc, density), fontsize=8,
                        xytext=(5, 5), textcoords='offset points')
    
    ax2.set_xlabel('GC Content (%)')
    ax2.set_ylabel('Hybrid Density (per Mb)')
    ax2.set_title('B. Hybrid Density vs GC Content')
    
    # Panel C: Hybrid composition heatmap
    ax3 = axes[2]
    
    # Create co-occurrence matrix
    hybrid_counts = hybrid_df.groupby(['Organism', 'Hybrid_Type'])['Count'].sum().unstack(fill_value=0)
    
    # Select top 10 hybrid types for the heatmap
    top_types = hybrid_df.groupby('Hybrid_Type')['Count'].sum().sort_values(ascending=False).head(10).index
    hybrid_counts = hybrid_counts[[c for c in top_types if c in hybrid_counts.columns]]
    
    # Sort by GC content
    if 'stats' in data:
        gc_order = data['stats'].sort_values('GC_Content_pct')['Organism'].tolist()
        hybrid_counts = hybrid_counts.reindex([o for o in gc_order if o in hybrid_counts.index])
    
    # Create heatmap
    if not hybrid_counts.empty:
        sns.heatmap(hybrid_counts, ax=ax3, cmap='YlOrRd', 
                    xticklabels=[h.replace('_Overlap', '').replace('_', '-')[:20] for h in hybrid_counts.columns],
                    cbar_kws={'label': 'Count'})
        ax3.set_xticklabels(ax3.get_xticklabels(), rotation=45, ha='right', fontsize=7)
        ax3.set_yticklabels(ax3.get_yticklabels(), fontsize=8)
    
    ax3.set_title('C. Hybrid Type Heatmap')
    
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.join(output_dir, 'Figure9_Hybrid_Analysis')
    plt.savefig(f'{output_path}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_path}.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"Saved: {output_path}.png and .pdf")


def plot_curved_dna_ratio(data, output_dir):
    """Generate Curved DNA Local:Global ratio figure."""
    if 'curved' not in data:
        print("No curved DNA data available")
        return
    
    curved_df = data['curved'].copy()
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Calculate Local:Global ratio
    curved_df['Local_Global_Ratio'] = curved_df['Local_Curvature'] / curved_df['Global_Curvature'].replace(0, 1)
    
    # Sort by GC content
    curved_df = curved_df.sort_values('GC_Content_pct')
    
    # Create bar plot
    x = np.arange(len(curved_df))
    
    # Color by phylum
    colors = []
    for phylum in curved_df['Phylum']:
        if phylum == 'Proteobacteria':
            colors.append(COLORS['blue'])
        elif phylum == 'Actinobacteria':
            colors.append(COLORS['vermillion'])
        else:
            colors.append(COLORS['bluish_green'])
    
    bars = ax.bar(x, curved_df['Local_Global_Ratio'], color=colors, edgecolor='black', linewidth=0.5)
    
    # Add ratio labels on bars
    for i, (idx, row) in enumerate(curved_df.iterrows()):
        if row['Local_Global_Ratio'] > 0:
            ax.text(i, row['Local_Global_Ratio'] + 0.5, f"{row['Local_Global_Ratio']:.1f}", 
                   ha='center', fontsize=9)
    
    ax.set_xticks(x)
    ax.set_xticklabels([f"{row['Organism']}\n({row['GC_Content_pct']:.1f}% GC)" 
                        for _, row in curved_df.iterrows()], 
                       rotation=45, ha='right')
    ax.set_ylabel('Local:Global Curvature Ratio')
    ax.set_title('Curved DNA Local:Global Curvature Ratio vs GC Content')
    
    # Add legend for phyla
    legend_elements = [
        mpatches.Patch(facecolor=COLORS['blue'], edgecolor='black', label='Proteobacteria'),
        mpatches.Patch(facecolor=COLORS['vermillion'], edgecolor='black', label='Actinobacteria'),
        mpatches.Patch(facecolor=COLORS['bluish_green'], edgecolor='black', label='Firmicutes')
    ]
    ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    
    # Save figure
    output_path = os.path.join(output_dir, 'Figure10_Curved_DNA_Ratio')
    plt.savefig(f'{output_path}.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_path}.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"Saved: {output_path}.png and .pdf")


def generate_all_visualizations():
    """Generate all subclass-level visualizations."""
    
    if not PLOTTING_AVAILABLE:
        print("ERROR: matplotlib/seaborn not available. Cannot generate visualizations.")
        return
    
    set_publication_style()
    
    print("\n" + "="*70)
    print("GENERATING SUBCLASS-LEVEL VISUALIZATIONS")
    print("="*70)
    
    # Load data
    print("\nLoading analysis data...")
    data = load_analysis_data()
    print(f"Loaded data for: {list(data.keys())}")
    
    # Create output directory
    base_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(base_dir, 'Genomes')
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate figures
    print("\nGenerating figures...")
    
    print("\n1. G-Quadruplex subclass distribution...")
    plot_g4_subclass_distribution(data, output_dir)
    
    print("\n2. Cluster analysis...")
    plot_cluster_analysis(data, output_dir)
    
    print("\n3. Hybrid motif analysis...")
    plot_hybrid_analysis(data, output_dir)
    
    print("\n4. Curved DNA ratio analysis...")
    plot_curved_dna_ratio(data, output_dir)
    
    print("\n" + "="*70)
    print("Visualization generation complete!")
    print("="*70)


if __name__ == '__main__':
    generate_all_visualizations()
