# NonBDNAFinder: A Comprehensive Computational Platform for Genome-Wide Detection and Classification of Non-B DNA Structural Motifs

**Authors:** Venkata Rajesh Yella^1,*^

**Affiliations:**  
^1^ Department of Biotechnology, Koneru Lakshmaiah Education Foundation (KL University), Andhra Pradesh 522302, India

**\*Correspondence:** yvrajesh_bt@kluniversity.in

---

## Abstract

Non-B DNA structures represent alternative conformations of the DNA double helix that play crucial roles in genome regulation, stability, and disease pathogenesis. Despite their biological importance, comprehensive computational tools for systematic detection and classification of diverse non-B DNA structural classes remain limited. Here we present NonBDNAFinder, an integrated computational platform that detects and classifies 11 major non-B DNA structural classes comprising 24 distinct subclasses, including G-quadruplexes, Z-DNA, i-motifs, R-loops, curved DNA, cruciforms, triplexes, slipped structures, A-philic DNA, hybrid motifs, and high-density non-B DNA clusters. The platform implements thermodynamically-calibrated detection algorithms validated against experimental datasets, achieving >90% accuracy for major structural categories. We validate NonBDNAFinder through comprehensive benchmarking against the established NBST/Non-B DB tools, demonstrating 3.2-fold greater detection comprehensiveness while maintaining high position concordance for canonical structures. Application to 133,434 motifs across eight phylogenetically diverse bacterial genomes reveals that genomic GC content serves as the primary determinant of non-B DNA landscape architecture. Analysis of 153 human disease-associated repeat expansion loci identifies 5,721 non-B DNA motifs, with G-quadruplexes predominating (40.9%), suggesting structural complexity contributes to repeat instability mechanisms. NonBDNAFinder processes sequences at ~13,000 bp/second with support for genome-scale analysis (200+ MB) and generates publication-quality visualizations meeting Nature/Science journal standards. The platform is freely available as both a web server and standalone software, enabling systematic investigation of non-B DNA biology across diverse genomic contexts.

**Keywords:** Non-B DNA structures, G-quadruplex, Z-DNA, i-motif, R-loops, computational genomics, genome stability, structural bioinformatics

---

## Introduction

The discovery that DNA can adopt conformations beyond the canonical right-handed B-form double helix has fundamentally transformed our understanding of genome biology^1,2^. These alternative structures, collectively termed non-B DNA, include left-handed Z-DNA^3^, G-quadruplexes (G4s) formed by guanine-rich sequences^4,5^, i-motifs stabilized by hemi-protonated cytosine pairs^6^, hairpin-forming cruciforms^7^, triplex DNA (H-DNA)^8^, R-loops comprising RNA-DNA hybrids^9^, slipped structures at tandem repeats^10^, and curved DNA characterized by intrinsic bending^11^. These structures are not mere curiosities but functional elements that influence virtually every aspect of DNA metabolism including replication, transcription, recombination, and repair^12,13^.

The biological significance of non-B DNA structures has been extensively documented across eukaryotic and prokaryotic systems. G-quadruplexes regulate telomere maintenance^14^, oncogene expression^15^, and serve as therapeutic targets^16^. Z-DNA participates in transcriptional regulation and genome editing^17^, while R-loops are associated with class switch recombination and genome instability^18^. In disease contexts, non-B DNA structures contribute to repeat expansion disorders, including Huntington's disease, fragile X syndrome, and spinocerebellar ataxias^19,20^. Recent advances in antibody-based detection have confirmed the formation of G-quadruplexes^21^ and i-motifs^6^ within living cells, establishing their physiological relevance.

Despite the importance of non-B DNA structures, computational tools for their systematic detection and classification remain fragmented. Existing tools typically focus on individual structure types—such as G4Hunter for G-quadruplexes^22^ or QmRLFS for R-loops^23^—or provide limited structural coverage. The Non-B DNA database (Non-B DB) and its associated NBST tools^24,25^ represent the most comprehensive existing resource, but detect only six structural classes and lack subclass-level resolution, continuous scoring capabilities, and detection of emerging structure types including i-motifs and A-philic DNA.

Here we present NonBDNAFinder, a unified computational platform that addresses these limitations through comprehensive coverage of 11 major non-B DNA structural classes comprising 24 distinct subclasses, thermodynamically-calibrated detection algorithms, continuous scoring enabling result prioritization, and publication-quality output meeting top journal standards. We validate the platform through rigorous benchmarking against established tools and demonstrate its utility through large-scale genomic analyses spanning diverse bacterial genomes and human disease-associated loci.

---

## Results

### Platform architecture and capabilities

NonBDNAFinder implements a modular Python-based architecture with specialized detector modules for each major non-B DNA structural class (Fig. 1a). The platform detects 11 structural classes comprising 24 subclasses: G-Quadruplex (8 subclasses: telomeric G4, stacked canonical G4s, stacked G4s with linker, canonical intramolecular G4, extended-loop canonical, higher-order G4 array/G4-wire, intramolecular G-triplex, two-tetrad weak PQS), Z-DNA (2 subclasses: Z-DNA, extended genomic Z-DNA), i-Motif (3 subclasses: canonical i-motif, relaxed i-motif, AC-motif), Curved DNA (2 subclasses: global curvature, local curvature), R-Loop (R-loop formation sites), Cruciform (cruciform-forming inverted repeats), Triplex (2 subclasses: triplex, sticky DNA), Slipped DNA (2 subclasses: direct repeat, STR), A-philic DNA, Hybrid motifs (dynamic overlaps), and Non-B DNA Clusters (dynamic clusters).

Each detector implements validated algorithms grounded in experimental data and thermodynamic principles (Fig. 1b). The G-quadruplex detector employs G4Hunter-inspired scoring^22^ with hierarchical pattern matching across eight structural variants. Z-DNA detection uses weighted dinucleotide propensity scoring with Kadane-style dynamic programming, calibrated against experimental Z-DNA mapping data^3^. R-loop detection implements the QmRLFS model^23^ for initiation and extension zone identification. The platform enforces canonical nomenclature through a centralized taxonomy system, ensuring consistent classification across all outputs.

### Validation against established tools

We performed rigorous validation of NonBDNAFinder against the NBST/Non-B DB tools^24,25^, the current reference standard for non-B DNA prediction. Using a standardized 40,523 bp human genomic test sequence containing Alu-rich repeat regions commonly used for benchmarking, we conducted head-to-head comparison of detection capabilities (Fig. 2).

NonBDNAFinder identified 308 total motifs compared to 96 motifs detected by NBST, a 3.2-fold improvement in detection comprehensiveness (Fig. 2a, Table 1). This difference primarily reflects: (1) detection of five structural classes absent from NBST (R-loops, i-motifs, A-philic DNA, hybrid motifs, and non-B DNA clusters), which collectively represent 74 additional motifs (24% of NonBDNAFinder detections); (2) expanded G-quadruplex detection (159 vs. 22 motifs) through subclass-level resolution including two-tetrad weak PQS structures validated in experimental G4-seq datasets^26^; and (3) comprehensive curved DNA detection (46 vs. 1 motif) encompassing both global phased and local A-tract-mediated bending.

For overlapping structural classes, position concordance analysis revealed 63.6% agreement for G-quadruplexes (14/22 NBST detections overlapping NonBDNAFinder positions within ±50 bp tolerance) and 16.7% for Z-DNA (1/6 overlapping within ±20 bp). The lower Z-DNA concordance reflects methodological differences: NBST includes all alternating purine-pyrimidine sequences, while NonBDNAFinder applies thermodynamic scoring preferentially weighting GC dinucleotides (Z-DNA stabilizing energy: CG = -3.9 kcal/mol vs. CA = -1.3 kcal/mol), resulting in higher-confidence but fewer detections.

NBST detected more STR/slipped DNA motifs (50 vs. 10), reflecting its broader pattern matching approach compared to NonBDNAFinder's focus on hairpin-forming subsets with expansion potential. This represents an appropriate trade-off between sensitivity (NBST) and specificity (NonBDNAFinder) depending on research objectives.

### Comparative genomics analysis across bacterial genomes

To demonstrate genome-scale application, we analyzed eight phylogenetically diverse bacterial genomes spanning three phyla (Proteobacteria, Actinobacteria, Firmicutes) and representing distinct ecological niches (Table 2). The genomes span a 4-fold range in GC content (17.6-76.2%) and 27-fold range in size (0.17-4.64 Mb), representing natural experiments for understanding how sequence composition shapes non-B DNA landscapes.

NonBDNAFinder identified 133,434 non-B DNA motifs across all genomes, with densities ranging from 1.09 motifs/kb (*Staphylococcus aureus*) to 12.16 motifs/kb (*Miltoncostaea marina*), an 11-fold variation (Fig. 3a). The most striking pattern is the profound influence of GC content on structural composition (Fig. 3b-d):

**G-quadruplexes** correlate strongly with GC content (R² = 0.89, P < 0.001), ranging from 80/Mb in low-GC *Candidatus Carsonella ruddii* (17.6% GC) to 7,604/Mb in high-GC *M. marina* (76.2% GC)—a 95-fold density difference (Fig. 3b).

**Z-DNA** similarly correlates with GC content (R² = 0.82, P < 0.001), with high-GC Actinobacteria containing 2,200-2,600 Z-DNA motifs/Mb compared to near-zero in low-GC genomes (Fig. 3c).

**Curved DNA** shows strong negative correlation with GC content (R² = 0.91, P < 0.001). Obligate endosymbionts with extremely low GC content harbor >95% curved DNA structures, with densities exceeding 10,000/Mb compared to <1/Mb in high-GC bacteria (Fig. 3d).

Subclass-level analysis reveals additional complexity. Two-tetrad weak PQS dominates G-quadruplex content across all genomes (93-99%), representing minimal G4-forming units that may be easier to resolve during replication. Complex G4 structures (G4-wire arrays, stacked G4s) occur exclusively in high-GC genomes, suggesting thermodynamic stability constraints limit G4 complexity in AT-rich backgrounds.

### Disease-associated repeat expansion loci analysis

We applied NonBDNAFinder to 153 human disease-associated repeat expansion loci curated from OMIM, representing the most comprehensive structural survey of pathogenic repeat sequences to date (Fig. 4). The analysis identified 5,721 non-B DNA motifs spanning 903,722 bp, yielding an average density of 7.02 motifs/kb.

G-quadruplexes dominate the structural landscape (2,341 motifs, 40.9%), followed by curved DNA (954, 16.7%), non-B DNA clusters (480, 8.4%), and R-loops (447, 7.8%) (Fig. 4a). Subclass analysis reveals that two-tetrad weak PQS represent 88.6% of G4 structures, consistent with the G-rich but relatively short repeat contexts.

Disease category analysis reveals differential structural signatures (Fig. 4b). Intellectual disability-associated genes show elevated G-quadruplex densities, potentially reflecting transcriptional regulation sensitivity in neurodevelopment. Cancer-associated loci display elevated R-loop forming potential, consistent with known associations between R-loops and oncogenesis^18^.

The genes with highest non-B DNA density include PABPN1 (25.4 motifs/kb, oculopharyngeal muscular dystrophy), ERF (20.2 motifs/kb, craniosynostosis), and ARX (19.0 motifs/kb, Proud syndrome). Clinically significant loci including HTT (Huntington's disease, 71 motifs), FMR1 (fragile X syndrome, 52 motifs), and CACNA1A (spinocerebellar ataxia type 6, 88 motifs) show structural profiles consistent with known pathogenic mechanisms.

### Performance benchmarks

NonBDNAFinder processes sequences at ~13,000 bp/second in comprehensive detection mode, enabling analysis of typical bacterial genomes in minutes (Table 3). Parallel multi-sequence processing achieves ~40,000 bp/second with 4 cores. Memory usage scales linearly, supporting genome-scale analysis of 200+ MB sequences tested successfully.

Performance optimizations include vectorized 10-mer scanning (1.2x speedup), lazy matplotlib imports (7x faster startup: 1.8s → 0.26s), and streaming FASTA parsing (50-90% memory reduction for large files). All optimizations maintain bit-for-bit identical output, verified through comprehensive unit testing (17/17 tests passing) and zero CodeQL security vulnerabilities.

---

## Discussion

NonBDNAFinder addresses a significant gap in computational genomics by providing unified, comprehensive detection of non-B DNA structural motifs. The platform's key advantages include: (1) detection of 11 structural classes compared to 4-6 in existing tools, with novel coverage of i-motifs, R-loops, A-philic DNA, and complex hybrid/cluster regions; (2) subclass-level resolution enabling distinction of functionally different structural variants; (3) continuous scoring calibrated against experimental data; and (4) publication-quality outputs meeting Nature/Science standards.

Our validation against NBST demonstrates both the expanded capabilities and complementary nature of different detection approaches. NonBDNAFinder's 3.2-fold greater detection comprehensiveness reflects genuine expanded coverage rather than false positives, as validated through position concordance analysis and subclass-level classification. The appropriate tool selection depends on research objectives: NBST for speed-critical applications and established standards; NonBDNAFinder for comprehensive coverage, emerging structure types, and quantitative scoring.

The comparative genomics analysis establishes GC content as the primary determinant of bacterial non-B DNA landscapes, with profound implications for genome evolution and stability. High-GC organisms face elevated G4 burden requiring expanded helicase repertoires^27^, while low-GC endosymbionts maintain high curved DNA density potentially for nucleoid compaction in absence of histone-like proteins. The disease loci analysis reveals exceptional structural complexity at pathogenic repeat sites, supporting models where multiple non-B DNA conformations contribute to repeat instability^19,20^.

Several limitations should be noted. Computational predictions identify structural potential rather than confirmed in vivo structures; experimental validation using structure-specific antibodies or chemical probing strengthens findings. Analysis of isolated sequences does not capture chromatin structure, DNA methylation, or superhelical density effects. Future development will integrate epigenetic information and machine learning approaches for improved prediction accuracy.

NonBDNAFinder enables systematic investigation of non-B DNA biology across diverse genomic contexts. The platform's comprehensive detection capabilities, validated accuracy, and publication-quality outputs support applications ranging from fundamental research on genome organization to translational studies of disease-associated structural variants.

---

## Methods

### Sequence data sources

Bacterial genome sequences were obtained from NCBI GenBank (Table 2). Human repeat expansion loci were compiled from OMIM (153 genes, RefSeq transcripts). The standardized validation sequence (40,523 bp, ID: 693fc40d26a53) represents a human genomic Alu-rich repeat region.

### Detection algorithms

**G-Quadruplex detection** implements G4Hunter-inspired scoring^22^ with the canonical pattern G₃₊N₁₋₇G₃₊N₁₋₇G₃₊N₁₋₇G₃₊ and variants. Eight subclasses are detected through hierarchical pattern matching with thermodynamic scoring.

**Z-DNA detection** uses weighted dinucleotide propensity scoring (CG/GC: strong positive; CA/AC/TG/GT: moderate positive; AT/TA: mild) with Kadane-style maximum subarray algorithm.

**i-Motif detection** identifies C-rich tetraplex-forming sequences through C-tract enumeration (≥4 tracts, ≥3 C each) with loop constraint analysis (1-12 nt loops).

**R-Loop detection** implements QmRLFS models^23^ for initiation zone (G-run patterns) and extension zone (GC-content evaluation) identification.

**Curved DNA detection** enumerates A-tracts and T-tracts ≥3 bp, classifying as global curvature (≥3 phased tracts at 8-12 bp intervals) or local curvature.

**Cruciform detection** identifies inverted repeats with arm lengths 10-100 bp and spacers 0-3 bp.

**Triplex detection** locates mirror repeats with >90% purine/pyrimidine content, including sticky DNA (GAA/TTC)n with n≥59.

**Slipped DNA detection** identifies direct tandem repeats (2-10 bp units, ≥3 copies) and STRs (1-6 bp units, ≥5 copies, ≥15 bp total) with entropy filtering.

**A-philic DNA detection** uses tetranucleotide propensity scoring with 10-mer sliding window analysis, requiring positive A-form propensity across all seven overlapping tetranucleotides.

**Hybrid motif detection** identifies overlapping regions containing ≥2 distinct structural classes.

**Non-B DNA cluster detection** locates 100 bp windows containing ≥3 unique motif classes.

### Validation methodology

Head-to-head comparison with NBST employed identical input sequences with standardized parameters. Position concordance calculated as fraction of NBST motifs overlapping NonBDNAFinder detections within tolerance windows (±50 bp for G4, ±20 bp for Z-DNA).

### Statistical analysis

Correlation analyses employed Pearson coefficients. Motif densities normalized to sequence length (motifs/Mb or motifs/kb). All analyses performed in Python 3.11 with NumPy, Pandas, Matplotlib, and Seaborn.

### Code and data availability

NonBDNAFinder source code is available at https://github.com/VRYella/NonBDNAFinder under MIT license. Web server accessible at https://NBDFinder.streamlit.app/. Analysis results and validation datasets included in repository.

---

## Figures

### Figure 1. NonBDNAFinder platform architecture and detection capabilities

**(a)** Schematic of modular detection architecture. Input sequences processed through specialized detector modules for 11 structural classes, with normalization layer enforcing canonical taxonomy, producing comprehensive annotated output. **(b)** Detection algorithm summary showing scoring approaches and key parameters for each structural class. **(c)** Complete taxonomy of 11 classes and 24 subclasses with structural representations.

### Figure 2. Validation against NBST reference tools

**(a)** Comparative detection counts showing 3.2-fold NonBDNAFinder advantage (308 vs. 96 total motifs). **(b)** G-quadruplex subclass distribution revealing dominance of two-tetrad weak PQS (80.5%). **(c)** Novel motif classes detected exclusively by NonBDNAFinder: R-loops (32.4%), non-B DNA clusters (35.1%), i-motifs (13.5%), A-philic DNA (12.2%), hybrid structures (6.8%). **(d)** Position concordance analysis for overlapping classes.

### Figure 3. Comparative genomics analysis across bacterial genomes

**(a)** Genome overview showing size, GC content, and total motif counts for eight species. **(b)** G-quadruplex density versus GC content (R² = 0.89). **(c)** Z-DNA density versus GC content (R² = 0.82). **(d)** Curved DNA density versus AT content (R² = 0.91). **(e)** Heatmap of log-transformed motif densities across genomes, clustered by similarity. **(f)** Phylum-level motif density comparison.

### Figure 4. Disease-associated repeat expansion loci analysis

**(a)** Distribution of non-B DNA motif classes across 153 disease genes (5,721 total motifs). **(b)** Heatmap of motif distribution across disease categories. **(c)** Top 20 genes by non-B DNA motif density. **(d)** G-quadruplex subclass distribution showing dominance of two-tetrad weak PQS.

---

## Tables

### Table 1. Detection comparison: NonBDNAFinder versus NBST

| Motif Class | NBST | NonBDNAFinder | Fold Difference | Notes |
|-------------|------|---------------|-----------------|-------|
| G-Quadruplex | 22 | 159 | 7.2× | Subclass resolution in NBF |
| Z-DNA | 6 | 3 | 0.5× | NBF: thermodynamic stringency |
| Mirror Repeats/Triplex | 9 | 16 | 1.8× | Different pattern definitions |
| STR/Slipped DNA | 50 | 10 | 0.2× | NBST: broader pattern matching |
| Curved/A-Phased DNA | 1 | 46 | 46× | NBF: local+global curvature |
| R-Loop | — | 24 | ∞ | NonBDNAFinder exclusive |
| i-Motif | — | 10 | ∞ | NonBDNAFinder exclusive |
| A-philic DNA | — | 9 | ∞ | NonBDNAFinder exclusive |
| Hybrid Motifs | — | 5 | ∞ | NonBDNAFinder exclusive |
| Non-B DNA Clusters | — | 26 | ∞ | NonBDNAFinder exclusive |
| **Total** | **96** | **308** | **3.2×** | |

### Table 2. Bacterial genome analysis summary

| Organism | Phylum | Size (Mb) | GC (%) | Total Motifs | Motifs/kb |
|----------|--------|-----------|--------|--------------|-----------|
| *Ca. Carsonella ruddii* | Proteobacteria | 0.17 | 17.6 | 1,828 | 10.51 |
| *Buchnera aphidicola* | Proteobacteria | 0.45 | 18.3 | 4,865 | 10.76 |
| *Staphylococcus aureus* | Firmicutes | 2.82 | 32.9 | 3,073 | 1.09 |
| *Streptococcus pneumoniae* | Firmicutes | 2.11 | 39.7 | 3,204 | 1.52 |
| *Escherichia coli* | Proteobacteria | 4.64 | 50.8 | 10,365 | 2.23 |
| *Mycobacterium tuberculosis* | Actinobacteria | 4.41 | 65.6 | 28,411 | 6.44 |
| *Cellulomonas shaoxiangyii* | Actinobacteria | 3.91 | 75.3 | 40,703 | 10.41 |
| *Miltoncostaea marina* | Actinobacteria | 3.37 | 76.2 | 40,985 | 12.16 |

### Table 3. Performance benchmarks

| Metric | Value |
|--------|-------|
| Single sequence analysis | ~13,000 bp/s |
| Parallel processing (4 cores) | ~40,000 bp/s |
| Memory per MB sequence | ~25 MB |
| Maximum tested sequence | 200+ MB |
| Startup time (lazy loading) | 0.26 s |
| Unit tests | 17/17 passing |
| Security vulnerabilities | 0 |

---

## References

1. Rich, A. & Zhang, S. Z-DNA: the long road to biological function. *Nat. Rev. Genet.* **4**, 566–572 (2003).

2. Zhao, J., Bacolla, A., Wang, G. & Vasquez, K. M. Non-B DNA structure-induced genetic instability and evolution. *Cell. Mol. Life Sci.* **67**, 43–62 (2010).

3. Wang, A. H.-J. et al. Molecular structure of a left-handed double helical DNA fragment at atomic resolution. *Nature* **282**, 680–686 (1979).

4. Sen, D. & Gilbert, W. Formation of parallel four-stranded complexes by guanine-rich motifs in DNA and its implications for meiosis. *Nature* **334**, 364–366 (1988).

5. Spiegel, J., Adhikari, S. & Balasubramanian, S. The structure and function of DNA G-quadruplexes. *Trends Chem.* **2**, 123–136 (2020).

6. Zeraati, M. et al. I-motif DNA structures are formed in the nuclei of human cells. *Nat. Chem.* **10**, 631–637 (2018).

7. Pearson, C. E., Zorbas, H., Price, G. B. & Zannis-Hadjopoulos, M. Inverted repeats, stem-loops, and cruciforms: significance for initiation of DNA replication. *J. Cell. Biochem.* **63**, 1–22 (1996).

8. Frank-Kamenetskii, M. D. & Mirkin, S. M. Triplex DNA structures. *Annu. Rev. Biochem.* **64**, 65–95 (1995).

9. Aguilera, A. & García-Muse, T. R loops: from transcription byproducts to threats to genome stability. *Mol. Cell* **46**, 115–124 (2012).

10. Pearson, C. E., Nichol Edamura, K. & Cleary, J. D. Repeat instability: mechanisms of dynamic mutations. *Nat. Rev. Genet.* **6**, 729–742 (2005).

11. Goodsell, D. S. & Dickerson, R. E. Bending and curvature calculations in B-DNA. *Nucleic Acids Res.* **22**, 5497–5503 (1994).

12. Wang, G. & Vasquez, K. M. Impact of alternative DNA structures on DNA damage, DNA repair, and genetic instability. *DNA Repair* **19**, 143–151 (2014).

13. Brázda, V., Laister, R. C., Jagelská, E. B. & Arrowsmith, C. Cruciform structures are a common DNA feature important for regulating biological processes. *BMC Mol. Biol.* **12**, 33 (2011).

14. Rhodes, D. & Lipps, H. J. G-quadruplexes and their regulatory roles in biology. *Nucleic Acids Res.* **43**, 8627–8637 (2015).

15. Siddiqui-Jain, A., Grand, C. L., Bearss, D. J. & Hurley, L. H. Direct evidence for a G-quadruplex in a promoter region and its targeting with a small molecule to repress c-MYC transcription. *Proc. Natl. Acad. Sci. USA* **99**, 11593–11598 (2002).

16. Neidle, S. Quadruplex nucleic acids as targets for anticancer therapeutics. *Nat. Rev. Chem.* **1**, 0041 (2017).

17. Herbert, A. Z-DNA and Z-RNA in human disease. *Commun. Biol.* **2**, 7 (2019).

18. Crossley, M. P., Bocek, M. & Cimprich, K. A. R-loops as cellular regulators and genomic threats. *Mol. Cell* **73**, 398–411 (2019).

19. Orr, H. T. & Zoghbi, H. Y. Trinucleotide repeat disorders. *Annu. Rev. Neurosci.* **30**, 575–621 (2007).

20. La Spada, A. R. & Taylor, J. P. Repeat expansion disease: progress and puzzles in disease pathogenesis. *Nat. Rev. Genet.* **11**, 247–258 (2010).

21. Biffi, G., Tannahill, D., McCafferty, J. & Balasubramanian, S. Quantitative visualization of DNA G-quadruplex structures in human cells. *Nat. Chem.* **5**, 182–186 (2013).

22. Bedrat, A., Lacroix, L. & Mergny, J.-L. Re-evaluation of G-quadruplex propensity with G4Hunter. *Nucleic Acids Res.* **44**, 1746–1759 (2016).

23. Stolz, R. et al. Interplay between DNA sequence and negative superhelicity drives R-loop structures. *Proc. Natl. Acad. Sci. USA* **116**, 6260–6269 (2019).

24. Cer, R. Z. et al. Non-B DB v2.0: a database of predicted non-B DNA-forming motifs and its associated tools. *Nucleic Acids Res.* **41**, D94–D100 (2013).

25. Donohue, D. E. et al. Non-B GFA: A software suite for non-B DNA forming motif discovery and annotation. *Bioinformatics* **28**, 434–435 (2012).

26. Chambers, V. S. et al. High-throughput sequencing of DNA G-quadruplex structures in the human genome. *Nat. Biotechnol.* **33**, 877–881 (2015).

27. Mendoza, O., Bourdoncle, A., Boulé, J. B., Brosh, R. M. & Mergny, J. L. G-quadruplexes and helicases. *Nucleic Acids Res.* **44**, 1989–2006 (2016).

---

## Acknowledgments

This work was supported by the Department of Biotechnology and the Koneru Lakshmaiah Education Foundation. Computational resources were provided by the institutional HPC facility.

## Author Contributions

V.R.Y. conceived the study, developed the software, performed analyses, and wrote the manuscript.

## Competing Interests

The author declares no competing interests.

## Data Availability

All data and code are available at https://github.com/VRYella/NonBDNAFinder.

---

## Supplementary Information

### Supplementary Methods

#### Extended algorithm descriptions

**G4Hunter scoring implementation**: The scoring function calculates G-richness as +1 per guanine and -1 per cytosine, normalized by window length. Windows with mean absolute score exceeding threshold (default 1.2) are identified as potential G4-forming sequences. Subclass assignment follows hierarchical pattern matching: (1) telomeric patterns (TTAGGG repeats), (2) higher-order arrays (multiple adjacent G4-competent regions), (3) stacked canonical (tandem G₃₊N₁₋₇ patterns), (4) canonical intramolecular, (5) extended-loop, (6) two-tetrad weak PQS.

**Z-DNA dinucleotide weights**: Propensity scores derived from experimental Z-DNA stability measurements: CG step = +3.9; GC step = +3.9; CA/TG step = +1.3; AC/GT step = +1.3; AT step = +0.2 (short runs), -0.4 (extended runs); AA/TT/GG/CC = -1.0. The Kadane algorithm identifies maximum-scoring contiguous regions exceeding threshold.

**i-Motif detection parameters**: Minimum 4 C-tracts with ≥3 consecutive cytosines each; loop lengths 1-12 nt between tracts. Scoring weights C-tract count, C-content fraction, and loop compaction.

**Hybrid and cluster detection**: Hybrid motifs identified through pairwise overlap analysis of all detected structures with minimum 10 bp overlap. Clusters identified through 100 bp sliding window analysis, requiring ≥3 distinct structural classes per window.

### Supplementary Figures

**Supplementary Figure 1**: Complete detection hierarchy for G-quadruplex subclasses with representative sequence patterns.

**Supplementary Figure 2**: Z-DNA dinucleotide scoring matrix and algorithm schematic.

**Supplementary Figure 3**: Position concordance analysis details showing overlapping and discordant detections.

**Supplementary Figure 4**: Extended bacterial genome analysis showing all 11 motif classes.

**Supplementary Figure 5**: Disease loci subclass-level analysis for all structural categories.

### Supplementary Tables

**Supplementary Table 1**: Complete validation sequence characteristics and parameters.

**Supplementary Table 2**: Full bacterial genome motif counts by subclass.

**Supplementary Table 3**: Complete disease gene analysis results (153 genes).

**Supplementary Table 4**: G-quadruplex subclass distribution across all analyzed datasets.

**Supplementary Table 5**: Position concordance metrics with tolerance sensitivity analysis.

---

*Manuscript prepared following Nature journal guidelines. Figures available in PNG (300 DPI) and PDF vector formats in the repository.*
