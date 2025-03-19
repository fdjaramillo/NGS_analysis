# Practical NGS Processing and Analysis

This repository contains an extended analysis of Next Generation Sequencing (NGS) data, focusing on RNA-seq, ChIP-seq, and differential splicing between brain and liver tissues. The project integrates multi-omics datasets to explore tissue-specific gene expression and regulatory mechanisms.  

This analysis expands on the work performed in this hands-on session, now focusing on Brain-Liver comparison:  
[UVIC NGS Hands-On](https://public-docs.crg.es/rguigo/Data/cklein/courses/UVIC/handsOn/).  

## Repository Structure

- **Practical NGS processing and Analysis.pdf**  
  Instructions and details of the practical exercise.  

- **commands_brain_liver.txt**  
  Bash commands used throughout the analysis.  

- **plots/**  
  Contains all generated plots (heatmaps, barplots, REVIGO visualizations,...).  

- **report.html**  
  A summary of results and methods used in the analysis. This report is also available in the [wiki](https://github.com/fdjaramillo/NGS_analysis/wiki).  

## Prerequisites
The analysis was conducted using a Docker container, following the setup instructions provided in the hands-on materials.

## Results Summary

- **Differential Expression**  
  225 genes were identified (FDR < 0.01, |logFC| > 10), highlighting distinct expression patterns between brain and liver.  

- **GO Enrichment Analysis**  
  Tissue-specific enriched GO terms reveal functional differences between brain and liver.  

- **Differential Splicing**  
  Significant splicing events (SE, RI, MX, AF) were detected using defined ΔPSI and p-value thresholds.  

- **ChIP-seq Integration**  
  H3K4me3 peak analysis combined with RNA-seq and ATAC-seq data provides insights into tissue-specific promoter activity.  

For detailed results and figures, refer to `report.html` or visit the [wiki](https://github.com/fdjaramillo/NGS_analysis/wiki).

