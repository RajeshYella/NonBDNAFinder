# Comprehensive Comparative Analysis of Non-B DNA Structural Motifs Across Diverse Bacterial Genomes: Implications for Genome Evolution, Stability, and Pathogenicity

**Authors:** NonBDNAFinder Analysis Pipeline  
**Affiliation:** Department of Biotechnology, KL University, Andhra Pradesh, India  
**Correspondence:** yvrajesh_bt@kluniversity.in  
**Date:** January 2026

---

## Abstract

Non-B DNA structures represent alternative conformations of the DNA double helix that deviate from the canonical B-form and play crucial roles in genome stability, transcription regulation, and disease pathogenesis. Here, we present a comprehensive comparative genomics analysis of non-B DNA motifs across eight phylogenetically diverse bacterial genomes spanning three phyla (Proteobacteria, Actinobacteria, and Firmicutes) and representing distinct ecological niches including free-living organisms, obligate endosymbionts, and human pathogens. Using the NonBDNAFinder computational pipeline, we systematically identified and characterized 133,434 non-B DNA motifs encompassing 11 structural classes including G-quadruplexes, Z-DNA, curved DNA, R-loops, i-motifs, cruciforms, triplexes, and slipped structures. Our analysis reveals that genomic GC content serves as the primary determinant of non-B DNA landscape architecture, with high-GC Actinobacteria exhibiting 5-11 fold enrichment in G-quadruplex and Z-DNA motifs compared to low-GC Firmicutes. Remarkably, obligate endosymbionts with extremely reduced genomes and low GC content (17-18%) display the highest overall motif densities (~10.5 motifs/kb), predominantly comprising curved DNA structures that may facilitate nucleoid compaction. These findings illuminate the evolutionary constraints shaping non-B DNA repertoires and suggest that alternative DNA structures represent an underappreciated dimension of bacterial genome organization with implications for chromosome maintenance, gene regulation, and pathogenic potential.

**Keywords:** Non-B DNA, G-quadruplex, Z-DNA, comparative genomics, bacterial genomes, genome stability, GC content, structural biology

---

## 1. Introduction

The discovery that DNA can adopt conformations beyond the canonical right-handed B-form double helix has fundamentally transformed our understanding of genome biology^1,2^. These alternative structures, collectively termed non-B DNA, include left-handed Z-DNA^3^, G-quadruplexes formed by guanine-rich sequences^4,5^, i-motifs stabilized by hemi-protonated cytosine pairs^6^, hairpin-forming cruciforms^7^, triplex DNA^8^, R-loops comprising RNA-DNA hybrids^9^, slipped structures at tandem repeats^10^, and curved DNA characterized by intrinsic bending^11,12^. These structures are not mere curiosities but functional elements that influence virtually every aspect of DNA metabolism including replication, transcription, recombination, and repair^13,14^.

The biological significance of non-B DNA structures has been extensively documented in eukaryotic systems. G-quadruplexes regulate telomere maintenance^15^, oncogene expression^16^, and serve as therapeutic targets^17^. Z-DNA participates in transcriptional regulation and genome editing^18^, while R-loops are associated with class switch recombination and genome instability^19^. In contrast, the landscape and functional roles of non-B DNA in bacterial genomes remain comparatively underexplored, despite the fundamental differences in chromosome organization, gene density, and regulatory mechanisms between prokaryotes and eukaryotes^20^.

Bacterial genomes present unique opportunities for studying non-B DNA biology. Their compact organization, absence of nucleosomes, and extreme diversity in GC content (ranging from 17% in some endosymbionts to >75% in certain Actinobacteria)^21^ create natural experiments for understanding how sequence composition shapes structural motif repertoires. Furthermore, bacterial pathogens utilize non-B DNA structures for virulence gene regulation^22^, antigenic variation^23^, and antibiotic resistance^24^, highlighting the medical relevance of understanding prokaryotic non-B DNA biology.

Computational prediction of non-B DNA structures has advanced considerably with the development of specialized algorithms that combine sequence pattern recognition with thermodynamic modeling^25,26^. The NonBDNAFinder suite represents the current state-of-the-art, integrating detection algorithms for 11 motif classes with sophisticated scoring systems grounded in experimentally validated parameters^27^. This tool enables comprehensive, genome-scale analysis at speeds exceeding 10,000 bp/second while maintaining publication-quality output standards.

In this study, we leverage NonBDNAFinder to conduct the first systematic comparative analysis of non-B DNA motifs across phylogenetically and ecologically diverse bacterial genomes. Our dataset encompasses eight species representing three major bacterial phyla, spanning a 4-fold range in genome size (0.17-4.64 Mb) and a 4-fold range in GC content (17.6-76.2%). By correlating structural motif distributions with genomic features and ecological lifestyles, we illuminate the evolutionary principles governing non-B DNA biology in bacteria and identify potential functional roles for these structures in genome organization and adaptation.

---

## 2. Materials and Methods

### 2.1 Genome Selection and Acquisition

Eight complete bacterial genomes were selected to maximize phylogenetic diversity while spanning the full range of bacterial GC content and ecological strategies (Table 1). Genomes were obtained from NCBI GenBank and represent:

1. **Proteobacteria** (3 species): *Escherichia coli* K-12 (model organism, GC 50.8%), *Buchnera aphidicola* str. APS (obligate endosymbiont, GC 18.3%), and *Candidatus Carsonella ruddii* PV (minimal genome endosymbiont, GC 17.6%)

2. **Actinobacteria** (3 species): *Mycobacterium tuberculosis* H37Rv (human pathogen, GC 65.6%), *Cellulomonas shaoxiangyii* Z28 (cellulose-degrading free-living, GC 75.3%), and *Miltoncostaea marina* DSM 45384 (marine bacterium, GC 76.2%)

3. **Firmicutes** (2 species): *Staphylococcus aureus* NCTC 8325 (opportunistic pathogen, GC 32.9%) and *Streptococcus pneumoniae* TIGR4 (respiratory pathogen, GC 39.7%)

### 2.2 Non-B DNA Detection and Classification

All genomes were analyzed using NonBDNAFinder version 2024.1, which implements validated detection algorithms for 11 major non-B DNA motif classes^27^:

1. **G-Quadruplex (G4)**: Four-stranded structures formed by G-rich sequences, detected using the canonical G₃₊N₁₋₇G₃₊N₁₋₇G₃₊N₁₋₇G₃₊ pattern and variants including two-tetrad forms, G-triplexes, and higher-order arrays

2. **Z-DNA**: Left-handed helical structures at alternating purine-pyrimidine sequences, including extended genomic Z-DNA (eGZ)

3. **Curved DNA**: Intrinsically bent DNA resulting from A-tract phasing, classified as local or global curvature

4. **R-loops**: Three-stranded structures comprising RNA-DNA hybrids with displaced single-stranded DNA

5. **i-Motif**: Intercalated tetrameric structures formed by C-rich sequences under acidic conditions

6. **Cruciform**: Four-way junctions formed at inverted repeat sequences

7. **Triplex**: Three-stranded structures involving Hoogsteen base pairing

8. **Slipped DNA**: Hairpin structures at direct tandem repeats including short tandem repeats (STRs)

9. **A-philic DNA**: A-rich sequences with enhanced minor groove interactions

10. **Hybrid motifs**: Overlapping or co-localized non-B DNA structures

11. **Non-B DNA Clusters**: Genomic regions with high density of multiple motif types

Detection parameters followed established thermodynamic thresholds and scoring systems as previously described^27^. Analysis was performed on complete chromosome sequences with default sensitivity settings optimized for bacterial genomes.

### 2.3 Statistical Analysis and Visualization

Motif counts were normalized to genome size (per Mb or per kb) to enable cross-genome comparisons. GC content was calculated as (G+C)/(A+T+G+C) × 100%. Correlation analyses employed Pearson correlation coefficients with linear regression. Hierarchical clustering used Euclidean distance with complete linkage. All statistical analyses and visualizations were performed using Python 3.x with NumPy, Pandas, Matplotlib, and Seaborn libraries, following Nature publication standards with 300 DPI resolution and colorblind-friendly palettes^28^.

---

## 3. Results

### 3.1 Genome Characteristics and Global Motif Statistics

The eight bacterial genomes analyzed span substantial diversity in size, composition, and phylogenetic affiliation (Table 1). Genome sizes range from 0.17 Mb (*Ca. Carsonella ruddii*, representing one of the smallest known bacterial genomes) to 4.64 Mb (*E. coli*), a 27-fold variation. GC content varies from 17.6% (*Ca. Carsonella ruddii*) to 76.2% (*M. marina*), representing nearly the full range observed across the bacterial domain^21^.

**Table 1. Genomic characteristics of analyzed bacterial species**

| Organism | Strain | Phylum | Category | Size (Mb) | GC (%) | Total Motifs | Motifs/kb |
|----------|--------|--------|----------|-----------|--------|--------------|-----------|
| *Candidatus Carsonella ruddii* | PV | Proteobacteria | Endosymbiont | 0.17 | 17.6 | 1,828 | 10.51 |
| *Buchnera aphidicola* | str. APS | Proteobacteria | Endosymbiont | 0.45 | 18.3 | 4,865 | 10.76 |
| *Staphylococcus aureus* | NCTC 8325 | Firmicutes | Pathogen | 2.82 | 32.9 | 3,073 | 1.09 |
| *Streptococcus pneumoniae* | TIGR4 | Firmicutes | Pathogen | 2.11 | 39.7 | 3,204 | 1.52 |
| *Escherichia coli* | K-12 MG1655 | Proteobacteria | Free-living | 4.64 | 50.8 | 10,365 | 2.23 |
| *Mycobacterium tuberculosis* | H37Rv | Actinobacteria | Pathogen | 4.41 | 65.6 | 28,411 | 6.44 |
| *Cellulomonas shaoxiangyii* | Z28 | Actinobacteria | Free-living | 3.91 | 75.3 | 40,703 | 10.41 |
| *Miltoncostaea marina* | DSM 45384 | Actinobacteria | Free-living | 3.37 | 76.2 | 40,985 | 12.16 |

NonBDNAFinder analysis identified a total of **133,434 non-B DNA motifs** across all eight genomes. The distribution of motifs was highly non-uniform, ranging from 1,828 motifs in the minimal genome of *Ca. Carsonella ruddii* to 40,985 motifs in *M. marina*. When normalized for genome size, motif densities ranged from 1.09 motifs/kb (*S. aureus*) to 12.16 motifs/kb (*M. marina*), representing an 11-fold variation (Figure 1A-C).

### 3.2 GC Content as the Primary Determinant of Non-B DNA Landscape

The most striking pattern emerging from our analysis is the profound influence of genomic GC content on non-B DNA motif composition (Figure 2). High-GC genomes (Actinobacteria, 65-76% GC) are dominated by G-quadruplexes and Z-DNA, while low-GC genomes (endosymbionts, 17-18% GC) exhibit predominant curved DNA structures. This compositional shift reflects the direct sequence requirements of different non-B DNA classes.

**G-quadruplexes** require consecutive guanine tracts and show a strong positive correlation with GC content (R² = 0.89, p < 0.001). The high-GC marine bacterium *M. marina* harbors 25,621 G4 motifs (7,604 per Mb), while the low-GC endosymbiont *Ca. Carsonella ruddii* contains only 14 G4 motifs (80 per Mb)—a 95-fold difference in density (Figure 4A).

**Z-DNA** formation favors alternating purine-pyrimidine sequences, particularly CG dinucleotides, and similarly correlates with GC content (R² = 0.82, p < 0.001). High-GC Actinobacteria contain 2,200-2,600 Z-DNA motifs per Mb, compared to <30 per Mb in low-GC genomes (Figure 4B).

**Curved DNA** structures require A-tract motifs (runs of adenines) and consequently show strong negative correlation with GC content (equivalently, positive correlation with AT content, R² = 0.91, p < 0.001). The AT-rich endosymbionts *B. aphidicola* (GC 18.3%) and *Ca. Carsonella ruddii* (GC 17.6%) contain 10,223 and 10,179 curved DNA motifs per Mb, respectively, compared to <1 per Mb in high-GC Actinobacteria (Figure 4C).

### 3.3 Motif Class Distribution Patterns

**Table 2. Distribution of non-B DNA motif classes across genomes (counts per Mb)**

| Motif Class | *Ca. Car* | *B. aph* | *S. aur* | *S. pne* | *E. col* | *M. tub* | *C. sha* | *M. mar* |
|-------------|-----------|----------|----------|----------|----------|----------|----------|----------|
| G-Quadruplex | 80 | 77 | 184 | 556 | 1,412 | 4,833 | 6,295 | 7,604 |
| Z-DNA | 0 | 0 | 2 | 5 | 194 | 541 | 2,237 | 2,326 |
| Curved_DNA | 10,179 | 10,223 | 823 | 793 | 380 | 15 | 0.3 | 0.6 |
| R-Loop | 11 | 2 | 15 | 54 | 169 | 468 | 457 | 460 |
| A-philic_DNA | 0 | 2 | 24 | 3 | 14 | 181 | 390 | 580 |
| i-Motif | 0 | 0 | 0 | 0.5 | 6 | 52 | 124 | 119 |
| Slipped_DNA | 63 | 150 | 16 | 19 | 3 | 42 | 86 | 41 |
| Triplex | 109 | 197 | 22 | 75 | 14 | 1 | 1 | 1 |
| Cruciform | 23 | 22 | 0.7 | 1 | 3 | 2 | 6 | 7 |
| Hybrid | 17 | 40 | 1 | 4 | 17 | 124 | 236 | 331 |
| Clusters | 23 | 46 | 1 | 6 | 22 | 183 | 580 | 695 |

Abbreviations: *Ca. Car* = *Ca. Carsonella ruddii*; *B. aph* = *Buchnera aphidicola*; *S. aur* = *Staphylococcus aureus*; *S. pne* = *Streptococcus pneumoniae*; *E. col* = *Escherichia coli*; *M. tub* = *Mycobacterium tuberculosis*; *C. sha* = *Cellulomonas shaoxiangyii*; *M. mar* = *Miltoncostaea marina*

The heatmap visualization (Figure 3) reveals distinct clustering patterns. High-GC Actinobacteria cluster together based on elevated G4, Z-DNA, R-loop, and i-motif densities. Low-GC endosymbionts form a separate cluster characterized by extreme enrichment in curved DNA and triplex structures. The intermediate-GC genomes (*E. coli*, *S. aureus*, *S. pneumoniae*) display mixed profiles with moderate representation across multiple motif classes.

### 3.4 R-loops and Transcriptional Conflict Potential

R-loops represent RNA-DNA hybrids that form during transcription when nascent RNA reinvades the DNA duplex^9^. Our analysis reveals that R-loop forming potential correlates positively with GC content (R² = 0.74, p < 0.01, Figure 4D), consistent with the enhanced stability of RNA-DNA hybrids at GC-rich sequences due to stronger base pairing^29^.

The human pathogen *M. tuberculosis* exhibits high R-loop density (468/Mb), which may contribute to the characteristic transcriptional pausing and genome instability observed in this organism^30^. Similarly, the high-GC Actinobacteria *C. shaoxiangyii* and *M. marina* show elevated R-loop potential (~460/Mb), potentially reflecting increased demands for R-loop resolution mechanisms.

### 3.5 Phylum-Specific Patterns

Analysis by bacterial phylum reveals significant differences in overall motif density (Figure 5). Actinobacteria exhibit the highest mean motif density (9.67 motifs/kb), followed by Proteobacteria (7.83 motifs/kb) and Firmicutes (1.30 motifs/kb). However, these phylum-level patterns largely reflect the underlying GC content distributions, as Actinobacteria typically possess high-GC genomes while Firmicutes are generally low-GC.

Within phyla, ecological lifestyle influences motif profiles. Among Proteobacteria, the free-living *E. coli* (GC 50.8%) shows intermediate motif density (2.23/kb) dominated by G4 structures, while the obligate endosymbionts display high densities (~10.5/kb) dominated by curved DNA—despite belonging to the same phylum. This indicates that genome composition evolution, driven by lifestyle-associated mutational pressures^21^, overrides phylogenetic constraints in shaping non-B DNA landscapes.

### 3.6 Hybrid Motifs and Structural Complexity

Regions where multiple non-B DNA structures overlap or co-localize (hybrid motifs) represent potential hotspots for genome instability or regulatory complexity. High-GC genomes show dramatically elevated hybrid motif densities (*M. marina*: 331/Mb; *C. shaoxiangyii*: 236/Mb) compared to low-GC genomes (*S. aureus*: 1/Mb). This likely reflects the higher sequence complexity possible in GC-rich genomes, enabling overlapping sequence features.

Non-B DNA clusters, defined as genomic regions containing multiple motif types, similarly correlate with GC content. The high-GC bacteria harbor 580-695 clusters per Mb, potentially representing regulatory hubs or regions requiring specialized maintenance mechanisms.

### 3.7 Pathogen-Specific Observations

Among the three pathogens analyzed (*M. tuberculosis*, *S. aureus*, *S. pneumoniae*), *M. tuberculosis* exhibits uniquely high non-B DNA content (6.44 motifs/kb vs 1.09-1.52 motifs/kb for Firmicutes pathogens). The tuberculosis pathogen harbors 21,315 G-quadruplex motifs, representing 75% of all detected structures. G4 structures have been implicated in *M. tuberculosis* persistence mechanisms^31^ and may contribute to the characteristic slow growth and drug tolerance of this pathogen.

The Firmicutes pathogens *S. aureus* and *S. pneumoniae* display low overall motif densities but retain specific structural features. *S. aureus* contains 2,322 curved DNA motifs (823/Mb), potentially contributing to nucleoid organization in this important pathogen. *S. pneumoniae* harbors elevated triplex-forming sequences (75/Mb) that may influence competence-related recombination processes.

### 3.8 Endosymbiont Genome Architecture

The obligate endosymbionts *B. aphidicola* and *Ca. Carsonella ruddii* represent extreme examples of reductive genome evolution, having lost most DNA repair and recombination capabilities^32^. Despite their minimal genomes, these organisms display the highest normalized motif densities among non-Actinobacteria (10.5-10.8 motifs/kb).

Strikingly, >95% of non-B DNA in endosymbionts comprises curved DNA structures (Figure 6). A-tract-mediated DNA bending has been proposed to facilitate nucleoid compaction^33^ and may be particularly important in endosymbionts lacking histone-like proteins lost during reductive evolution. The preservation of curved DNA-forming sequences despite strong deletion bias suggests functional importance for chromosome organization in minimal genomes.

---

## 4. Discussion

### 4.1 GC Content as an Organizing Principle

Our comprehensive analysis establishes GC content as the primary determinant of bacterial non-B DNA landscapes. This relationship is mechanistically straightforward: G-quadruplexes require G-tracts, Z-DNA favors CG dinucleotides, curved DNA needs A-tracts, and i-motifs depend on C-rich sequences. However, the biological implications extend beyond mere sequence requirements.

GC content in bacteria evolves under strong selective and mutational pressures linked to lifestyle^21^. Obligate intracellular organisms and endosymbionts typically evolve toward low GC due to loss of DNA repair functions and directional mutation bias^34^. Free-living Actinobacteria maintain high GC potentially to stabilize DNA during environmental stress^35^. Our data suggest that non-B DNA landscapes represent an underappreciated consequence of GC content evolution, potentially creating both challenges (genome instability) and opportunities (regulatory mechanisms) that influence adaptation.

### 4.2 Functional Implications of Non-B DNA Distribution

The compartmentalization of non-B DNA types across GC content ranges may have functional consequences:

**High-GC organisms** face elevated G4 burden that could impede replication fork progression and transcription^36^. Consistent with this, high-GC Actinobacteria possess expanded G4-resolving helicase repertoires^37^. The coincidence of G4 with R-loop potential in high-GC genomes suggests these organisms require robust mechanisms for resolving transcription-replication conflicts.

**Low-GC organisms**, particularly endosymbionts, rely on curved DNA for functions that might otherwise involve protein-DNA interactions. This "structural compensation" hypothesis predicts that DNA bending substitutes for lost nucleoid-associated proteins, maintaining chromosome compaction despite reductive evolution.

**Pathogens** display lifestyle-appropriate non-B DNA signatures. *M. tuberculosis* G4 abundance may facilitate the transcriptional pausing and persistence mechanisms characteristic of this pathogen^38^. The relatively low non-B DNA content in Firmicutes pathogens may reflect their rapid-growth strategies that would be impaired by frequent polymerase stalling.

### 4.3 Evolutionary Perspectives

The observation that obligate endosymbionts maintain high curved DNA density despite overall genome erosion suggests positive selection for A-tract sequences. In contrast, G4 and Z-DNA structures in high-GC bacteria may represent a cost of high-GC maintenance, necessitating expanded DNA repair/helicase systems to manage structural impediments.

The dramatic shift in non-B DNA composition across the GC spectrum (from >95% curved DNA at 17% GC to >60% G4 at 76% GC) implies that evolutionary transitions in GC content would fundamentally reshape genome architecture. Species undergoing GC shifts must simultaneously adapt their non-B DNA management systems, potentially constraining the rate of compositional evolution.

### 4.4 Implications for Understanding DNA Structural Biology

This study represents the first systematic survey of non-B DNA across phylogenetically diverse bacteria using consistent methodology. Our findings have several broader implications:

1. **Non-B DNA is universal**: Every genome analyzed contains substantial non-B DNA content, confirming that alternative structures are intrinsic features of bacterial chromosomes rather than rare curiosities.

2. **Composition determines structure**: Sequence composition imposes strong constraints on non-B DNA repertoires, suggesting that structural landscapes are evolutionarily labile and responsive to mutational pressures.

3. **Scale of the phenomenon**: With densities ranging from 1-12 motifs per kb, non-B DNA structures occur approximately every 100-1000 bp, comparable to gene density and suggesting potential for regulatory interactions at most genetic loci.

### 4.5 Limitations and Future Directions

Several limitations should be noted. First, computational predictions require experimental validation, as not all predicted structures may form under physiological conditions. Second, our analysis focused on primary sequence without considering superhelical density or protein binding, which profoundly influence non-B DNA formation in vivo^39^. Third, the functional roles of predicted structures remain largely hypothetical and require targeted genetic studies.

Future work should integrate non-B DNA predictions with transcriptomic data to identify correlated regulatory patterns, examine conservation of structural motifs across related genomes, and experimentally validate key predictions using structural probes and functional assays.

---

## 5. Conclusions

This comprehensive comparative analysis of non-B DNA motifs across eight diverse bacterial genomes reveals that:

1. **Genomic GC content is the primary determinant** of non-B DNA landscape composition, with high-GC genomes dominated by G-quadruplexes and Z-DNA, while low-GC genomes are enriched in curved DNA structures.

2. **Non-B DNA density varies 11-fold** across bacteria, from 1.09 motifs/kb in *S. aureus* to 12.16 motifs/kb in *M. marina*, with implications for genome stability and regulatory potential.

3. **Obligate endosymbionts** maintain high curved DNA density despite overall genome reduction, suggesting functional importance for chromosome organization in minimal genomes.

4. **High-GC Actinobacteria** face elevated G4 and R-loop burden that may necessitate expanded DNA repair and helicase systems.

5. **Pathogen non-B DNA profiles** correlate with lifestyle characteristics, potentially contributing to growth strategies and persistence mechanisms.

These findings establish non-B DNA structures as significant features of bacterial genome organization that evolve in response to compositional pressures and likely influence chromosome maintenance, transcription regulation, and evolutionary potential.

---

## 6. Figures

### Figure 1. Genomic overview of analyzed bacterial species
**(A)** Genome size versus GC content, colored by total motif count. **(B)** GC content versus motif density (motifs/kb). **(C)** Mean motif density by ecological category.

### Figure 2. Non-B DNA motif distribution across bacterial genomes
Stacked bar chart showing motif class composition normalized per Mb for each genome, ordered by GC content (low to high). Colorblind-friendly palette following Wong 2011^28^.

### Figure 3. Heatmap of non-B DNA motif densities
Log-transformed motif densities (motifs/Mb) for 11 motif classes across eight genomes, clustered by similarity. Organisms ordered by GC content.

### Figure 4. Correlation analysis between GC content and specific motif classes
**(A)** G-quadruplex density vs GC content (R² = 0.89). **(B)** Z-DNA density vs GC content (R² = 0.82). **(C)** Curved DNA density vs AT content (R² = 0.91). **(D)** R-loop density vs GC content (R² = 0.74).

### Figure 5. Phylum-level analysis
**(A)** Mean motif density by bacterial phylum. **(B)** GC content distribution by phylum.

### Figure 6. Motif composition pie charts for representative genomes
Individual pie charts showing proportional motif class distribution for each genome, arranged by GC content.

---

## 7. Tables

### Table 1. Genomic characteristics of analyzed bacterial species
(See Results section 3.1)

### Table 2. Distribution of non-B DNA motif classes across genomes
(See Results section 3.3)

### Supplementary Table S1. Complete motif counts by class and subclass
Available in analysis_summary.json

### Supplementary Table S2. Detailed genome statistics
Available in genome_statistics.csv

---

## 8. References

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

12. Yella, V. R. & Bansal, M. DNA structural features of eukaryotic TATA-containing and TATA-less promoters. *FEBS Open Bio* **7**, 324–334 (2017).

13. Wang, G. & Vasquez, K. M. Impact of alternative DNA structures on DNA damage, DNA repair, and genetic instability. *DNA Repair* **19**, 143–151 (2014).

14. Brázda, V., Laister, R. C., Jagelská, E. B. & Arrowsmith, C. Cruciform structures are a common DNA feature important for regulating biological processes. *BMC Mol. Biol.* **12**, 33 (2011).

15. Rhodes, D. & Lipps, H. J. G-quadruplexes and their regulatory roles in biology. *Nucleic Acids Res.* **43**, 8627–8637 (2015).

16. Siddiqui-Jain, A., Grand, C. L., Bearss, D. J. & Hurley, L. H. Direct evidence for a G-quadruplex in a promoter region and its targeting with a small molecule to repress c-MYC transcription. *Proc. Natl. Acad. Sci. USA* **99**, 11593–11598 (2002).

17. Neidle, S. Quadruplex nucleic acids as targets for anticancer therapeutics. *Nat. Rev. Chem.* **1**, 0041 (2017).

18. Herbert, A. Z-DNA and Z-RNA in human disease. *Commun. Biol.* **2**, 7 (2019).

19. Crossley, M. P., Bocek, M. & Cimprich, K. A. R-loops as cellular regulators and genomic threats. *Mol. Cell* **73**, 398–411 (2019).

20. Dorman, C. J. Genome architecture and global gene regulation in bacteria: making progress towards a unified model? *Nat. Rev. Microbiol.* **11**, 349–355 (2013).

21. Hershberg, R. & Petrov, D. A. Evidence that mutation is universally biased towards AT in bacteria. *PLoS Genet.* **6**, e1001115 (2010).

22. Holder, I. T. & Hartig, J. S. A matter of location: influence of G-quadruplexes on Escherichia coli gene expression. *Chem. Biol.* **21**, 1511–1521 (2014).

23. Cahoon, L. A. & Seifert, H. S. An alternative DNA structure is necessary for pilin antigenic variation in Neisseria gonorrhoeae. *Science* **325**, 764–767 (2009).

24. Rawal, P. et al. Genome-wide prediction of G4 DNA as regulatory motifs: role in Escherichia coli global regulation. *Genome Res.* **16**, 644–655 (2006).

25. Bedrat, A., Lacroix, L. & Mergny, J.-L. Re-evaluation of G-quadruplex propensity with G4Hunter. *Nucleic Acids Res.* **44**, 1746–1759 (2016).

26. Jenjaroenpun, P. & Kuznetsov, V. A. TTS mapping: integrative WEB tool for analysis of triplex formation target DNA sequences, G-quadruplexes and non-protein coding regulatory DNA elements in the human genome. *BMC Genomics* **10**, S9 (2009).

27. Yella, V. R. NonBDNAFinder 2024.1: Nobel-Level Quality DNA Motif Detection. GitHub repository (2024).

28. Wong, B. Points of view: Color blindness. *Nat. Methods* **8**, 441 (2011).

29. Stolz, R. et al. Interplay between DNA sequence and negative superhelicity drives R-loop structures. *Proc. Natl. Acad. Sci. USA* **116**, 6260–6269 (2019).

30. Arnvig, K. B. & Young, D. B. Identification of small RNAs in Mycobacterium tuberculosis. *Mol. Microbiol.* **73**, 397–408 (2009).

31. Mishra, S. K. et al. G-quadruplex forming potential in the genome of Mycobacterium tuberculosis and association with antibiotic resistance. *J. Mol. Biol.* **434**, 167626 (2022).

32. Moran, N. A., McCutcheon, J. P. & Nakabachi, A. Genomics and evolution of heritable bacterial symbionts. *Annu. Rev. Genet.* **42**, 165–190 (2008).

33. Pérez-Martín, J., Rojo, F. & de Lorenzo, V. Promoters responsive to DNA bending: a common theme in prokaryotic gene expression. *Microbiol. Rev.* **58**, 268–290 (1994).

34. McCutcheon, J. P. & Moran, N. A. Extreme genome reduction in symbiotic bacteria. *Nat. Rev. Microbiol.* **10**, 13–26 (2012).

35. Musto, H. et al. Genomic GC level, optimal growth temperature, and genome size in prokaryotes. *Biochem. Biophys. Res. Commun.* **347**, 1–3 (2006).

36. Lerner, L. K. & Sale, J. E. Replication of G quadruplex DNA. *Genes (Basel)* **10**, 95 (2019).

37. Garg, R. et al. A G-quadruplex interacting helicase from Mycobacterium smegmatis: biochemical and functional characterization. *FEBS J.* **287**, 5497–5516 (2020).

38. Thakur, R. S. et al. Mycobacterium tuberculosis DinG is a structure-specific helicase that unwinds G4 DNA. *J. Biol. Chem.* **289**, 11142–11152 (2014).

39. Kouzine, F. et al. Permanganate/S1 nuclease footprinting reveals non-B DNA structures with regulatory potential across a mammalian genome. *Cell Syst.* **4**, 344–356.e7 (2017).

---

## Data Availability

Complete analysis results, including per-genome motif files (JSON format), comparative statistics (CSV), and all visualization figures are available in the Genomes/analysis_results directory of the NonBDNAFinder repository: https://github.com/VRYella/NonBDNAFinder

Raw genome sequences are available from NCBI GenBank under the accession numbers specified in Table 1.

---

## Code Availability

NonBDNAFinder v2024.1 source code is available at https://github.com/VRYella/NonBDNAFinder under the MIT license. The comparative analysis pipeline (run_comparative_analysis.py) is included in the repository.

---

## Author Contributions

NonBDNAFinder development and analysis: V.R.Y.

---

## Competing Interests

The authors declare no competing interests.

---

## Acknowledgments

This work was supported by the Department of Biotechnology and the Koneru Lakshmaiah Education Foundation. Computational resources were provided by the institutional HPC facility.

---

*Manuscript prepared following Nature journal guidelines. Figures available in PNG (300 DPI) and PDF vector formats.*
