# Consolidated Writeup for NonBDNAFinder Analysis

This folder contains the consolidated writeup, figures, and tables from the comprehensive analysis of disease-associated repeat expansion loci using the NonBDNAFinder computational platform.

## Contents

### Documents

1. **NonBDNA_Repeat_Expansion_Analysis.md** - Comprehensive analysis of 153 disease-associated repeat expansion loci, including:
   - Analysis of 5,721 non-B DNA structural motifs
   - Disease category classification and motif distribution
   - Top genes by motif density
   - Clinically significant findings
   - Discussion of implications for repeat expansion disorders

2. **NonBDNAFinder_Tool_Documentation.md** - Complete documentation of the NonBDNAFinder tool, including:
   - Tool architecture and design
   - Motif detection algorithms
   - Performance benchmarking
   - Scientific validation

### Figures (in `figures/` directory)

| Figure | Description | Format |
|--------|-------------|--------|
| Figure1_Motif_Distribution_Pie | Pie chart of non-B DNA motif class distribution | PNG/PDF |
| Figure2_Motif_Counts_Bar | Bar chart of motif counts by class | PNG/PDF |
| Figure3_Disease_Categories | Distribution of genes by disease category | PNG/PDF |
| Figure4_Disease_Motif_Heatmap | Heatmap of motifs across disease categories | PNG/PDF |
| Figure5_GC_vs_Motifs | Scatter plot of GC content vs total motifs | PNG/PDF |
| Figure6_Top20_Genes | Top 20 genes by non-B DNA motif density | PNG/PDF |
| Figure7_Subclass_Distribution | Distribution of top 15 motif subclasses | PNG/PDF |

### Tables (in `tables/` directory)

| Table | Description | Format |
|-------|-------------|--------|
| Table0_Summary_Statistics.csv | Overall summary statistics | CSV |
| Table1_Complete_Analysis.csv | Complete analysis for all 153 genes | CSV |
| Table2_Detailed_Motifs.csv | Detailed listing of all 5,721 motifs | CSV |
| Table3_Disease_Motif_Matrix.csv | Disease category vs motif class matrix | CSV |

## Key Findings

### Summary Statistics
- **Total Sequences Analyzed:** 153
- **Total Non-B DNA Motifs:** 5,721
- **Total Base Pairs:** 791,842
- **Average Motifs per Gene:** 37.4
- **Average Motifs per kb:** 7.2

### Motif Class Distribution
1. G-Quadruplex: 2,341 (40.9%)
2. Curved DNA: 954 (16.7%)
3. Non-B DNA Clusters: 480 (8.4%)
4. R-Loop: 447 (7.8%)
5. Triplex: 322 (5.6%)
6. Hybrid: 312 (5.5%)
7. Z-DNA: 257 (4.5%)
8. A-philic DNA: 244 (4.3%)
9. Slipped DNA: 198 (3.5%)
10. i-Motif: 161 (2.8%)
11. Cruciform: 5 (0.1%)

### Disease Categories Analyzed
- Syndromes: 43 genes
- Other: 35 genes
- Intellectual Disability: 23 genes
- Cancer: 18 genes
- Neurodegenerative: 12 genes
- Ataxias: 10 genes
- Epilepsy/Encephalopathy: 7 genes
- Cardiac Disorders: 3 genes
- Muscular Disorders: 2 genes

### Top 5 Genes by Non-B DNA Density
1. **PABPN1** (25.4/kb) - Oculopharyngeal muscular dystrophy
2. **ERF** (20.2/kb) - Craniosynostosis 4
3. **ARX** (19.0/kb) - Proud syndrome
4. **HRAS** (17.8/kb) - Thyroid carcinoma susceptibility
5. **KRT10** (17.3/kb) - Ichthyosis with epidermolytic hyperkeratosis

## Data Source

Analysis was performed on the `repeat_expansion_loci_annotated.fa` file containing 153 disease-associated repeat expansion loci from OMIM (Online Mendelian Inheritance in Man).

## Tool Information

NonBDNAFinder v2024.1
- Web server: https://NBDFinder.streamlit.app/
- GitHub: https://github.com/VRYella/NBDFinder

## Citation

If you use these results in your research, please cite:
> Yella VR et al. NonBDNAFinder: A tool for high-throughput detection and visualization of non-canonical DNA motifs. 2024.

## Contact

For questions about this analysis, please contact:
- Email: yvrajesh_bt@kluniversity.in
- Department of Biotechnology, KL University, Andhra Pradesh, India
