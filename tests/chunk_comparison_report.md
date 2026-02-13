# NonBDNAFinder Chunk Size Comparison Report

**Generated:** 2026-02-13T07:11:56.010204
**Sequence:** ecoli_like_test
**Sequence Length:** 100,000 bp

## Summary Table

| Chunk Size | Overlap | Chunks | Motifs | Time (s) | Throughput (bp/s) |
|------------|---------|--------|--------|----------|-------------------|
| 5,000 bp | 500 bp | 20 | 273 | 2.843 | 35,169 |
| 10,000 bp | 1000 bp | 10 | 274 | 2.904 | 34,430 |
| 25,000 bp | 2500 bp | 4 | 273 | 3.025 | 33,059 |
| 50,000 bp | 2500 bp | 2 | 273 | 2.903 | 34,448 |
| 100,000 bp | 2500 bp | 1 | 272 | 2.85 | 35,083 |

## Detection Consistency

Baseline chunk size: **5,000 bp**

| Chunk Size | Common | Only Baseline | Only Current | Jaccard |
|------------|--------|---------------|--------------|---------|
| 10,000 bp | 268 | 5 | 6 | 0.9606 |
| 25,000 bp | 267 | 6 | 6 | 0.9570 |
| 50,000 bp | 268 | 5 | 5 | 0.9640 |
| 100,000 bp | 267 | 6 | 5 | 0.9604 |

## Motif Distribution by Class

| Chunk Size | A-philic_DNA | Cruciform | Curved_DNA | G-Quadruplex | Hybrid | Non-B_DNA_Clusters | R-Loop | Z-DNA | i-Motif |
|------------|---------------|---------------|---------------|---------------|---------------|---------------|---------------|---------------|---------------|
| 5,000 bp | 22 | 14 | 24 | 139 | 7 | 6 | 50 | 7 | 4 |
| 10,000 bp | 22 | 14 | 24 | 139 | 9 | 5 | 50 | 7 | 4 |
| 25,000 bp | 22 | 14 | 24 | 139 | 9 | 5 | 49 | 7 | 4 |
| 50,000 bp | 22 | 14 | 24 | 139 | 9 | 5 | 49 | 7 | 4 |
| 100,000 bp | 22 | 14 | 24 | 139 | 9 | 5 | 48 | 7 | 4 |

## Analysis

### Key Observations:

1. **Best throughput:** 5,000 bp (35,169 bp/s)
2. **Most consistent:** 50,000 bp (Jaccard: 0.9640)
3. **Motif count range:** 272 - 274

### Recommendations:

Based on the analysis:
- **Recommended chunk size:** 100,000 bp
  - Provides good balance of speed and detection consistency
- **For maximum speed:** Use 5,000 bp chunks
- **For genome-scale analysis:** Consider 25,000-50,000 bp chunks with 2,500 bp overlap