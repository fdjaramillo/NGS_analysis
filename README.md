# Practical NGS Processing and Analysis

This repository contains the work performed for the practical exercise of the Transcriptomics course of UVIC Omics Data Analysis Master. The practical is based on Next Generation Sequencing (NGS) data, focusing on RNA-seq, ChIP-seq, differential splicing analysis, and visualization. The goal is to integrate multi-omics datasets to explore tissue-specific gene expression and regulatory mechanisms.

This analysis focuses on a brain-liver comparison, expanding the work performed in this [hands-on session](https://public-docs.crg.es/rguigo/Data/cklein/courses/UVIC/handsOn/).

## Repository files

- **Practical NGS processing and Analysis.pdf** – Instructions and details of the practical exercise.  
- **commands_brain_liver** – Bash commands used throughout the analysis.  
- **plots/** – Directory containing all generated plots (heatmaps, bar plots, REVIGO visualizations, etc.).  
- **report.html** – Summary of results and methods used in the analysis (also available in the [wiki](https://github.com/fdjaramillo/NGS_analysis/wiki)).  

## Software used

The analysis was conducted using a Docker container. Please follow the setup instructions provided in the [hands-on session](https://public-docs.crg.es/rguigo/Data/cklein/courses/UVIC/handsOn/) material to replicate the environment.

## Results summary  

- **Differential Expression** – 225 genes identified (FDR < 0.01, |logFC| > 10), highlighting distinct expression patterns between brain and liver.  
- **GO Enrichment Analysis** – Tissue-specific enriched GO terms confirm functional differences between brain and liver.  
- **Differential Splicing** – Significant splicing events (SE, RI, MX, AF) detected using defined ΔPSI and p-value thresholds.  
- **ChIP-seq Integration** – H3K4me3 peak analysis combined with RNA-seq and ATAC-seq data provides insights into tissue-specific promoter activity.  

For detailed results and figures, refer to `report.html` or visit the [wiki](https://github.com/fdjaramillo/NGS_analysis/wiki).  

## Acknowledgments
The course instructor and the hands-on session materials provided by the CRG.
