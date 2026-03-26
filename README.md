# miRNA Embryonic Transitions

This repository contains an assignment investigating microRNA (miRNA) regulation during the maternal-to-zygotic transition in *Drosophila melanogaster* embryos.

The project compares miRNA expression between 2-hour and 5-hour embryonic stages using differential expression analysis, followed by target prediction and functional/pathway enrichment analysis to explore the regulatory roles of significant miRNAs during early development.

This repository includes:
- The assignment report
- R code used for downstream analysis and visualisation
- Generated figures and outputs where applicable

## Analysis Overview
The workflow includes:
- Processing differential expression results from DESeq2
- Identifying significantly up- and down-regulated miRNAs
- Generating a volcano plot of differential expression results
- Filtering predicted miRNA targets using confidence-based criteria
- Performing functional and pathway enrichment analysis
- Visualising enrichment results in R
